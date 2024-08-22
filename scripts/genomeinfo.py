import gzip
from collections import defaultdict

import polars as pl
import pysam


class GenomeInfo:
    """Stores genome information from a large insert size analysis"""

    def __init__(
        self,
        large_insert_bam_loc: str,
        large_insert_bai_loc: str,
        coverage_loc: str,
        genome_metadata: dict = None,
        gff_loc: str = None,
        chr_list: list = None,
    ):
        self.large_insert_bam_loc = large_insert_bam_loc
        self.large_insert_bai_loc = large_insert_bai_loc
        self.coverage_loc = coverage_loc
        self.gff_loc = gff_loc
        self.genome_metadata = genome_metadata
        self.bam = None
        self.coverage = None
        self.high_coverage_regions = None
        self.gff = None
        self.chr_list = chr_list

    def load_bam(self, pareto_cutoff=0.8):
        "Load in a (filtered) BAM file + index and filter by chromosome list"
        self.bam = pysam.AlignmentFile(
            self.large_insert_bam_loc, "rb", index_filename=self.large_insert_bai_loc
        )
        # Use the chromosome list if it's there, otherwise calculate it
        if self.chr_list is None:
            self.calculate_chr_list(pareto_cutoff=pareto_cutoff)
        # self.filter_bam()

    def load_coverage(self):
        "Load in a filtered coverage file and filter by chromosome list"
        # Load in filtered coverage file
        self.coverage = pl.read_csv(self.coverage_loc, separator="\t")
        # No header, so add column names
        self.coverage.columns = ["chromosome", "position", "coverage"]
        # Filter if a chromosome list is provided
        if self.chr_list is None:
            self.calculate_chr_list()
        self.coverage = self.filter_df(self.coverage)

    def load_gff(self):
        "Load in a GFF file, keeping only the relevant lines and formatting the attributes"
        lines_to_keep = []
        # Open the GFF file, if it's gzipped use gzip.open
        is_gzipped = self.gff_loc.endswith(".gz")
        if is_gzipped:
            with gzip.open(self.gff_loc, "rt", encoding="utf-8") as f:
                for line in f:
                    if not line.startswith("#"):
                        lines_to_keep.append(line)
        else:
            with open(self.gff_loc, "r") as f:  # noqa: UP015
                for line in f:
                    if not line.startswith("#"):
                        lines_to_keep.append(line)

        linesplit = [line.split("\t") for line in lines_to_keep]
        # Turn into polars df
        gff_schema = [
            "chromosome",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ]
        self.gff = pl.DataFrame(linesplit, schema=gff_schema)
        # Filter lines that start with #
        self.gff = self.gff.filter(~pl.col("feature").str.contains("region"))
        # Filter if a chromosome list is provided
        if self.chr_list is None:
            self.calculate_chr_list()

        self.gff = self.filter_df(self.gff)

        gff_attributes = self.format_gff_attributes()

        self.gff = self.gff.drop("attribute").hstack(gff_attributes)

        # Turn start and end into integers
        self.gff = self.gff.with_columns(
            [
                (pl.col("start").str.to_integer().alias("start")),
                (pl.col("end").str.to_integer().alias("end")),
            ]
        )

    def filter_bam(self):
        """Given a BAM file, filter to only include reads from
        the chromosomes in chr_list and store back to self.bam
        """
        filtered_reads = [read for read in self.bam if read.reference_name in self.chr_list]
        # Turn into a pysam.AlignmentFile object
        filtered_reads = pysam.AlignmentFile(filtered_reads)
        self.bam = filtered_reads

    def filter_df(self, df):
        """Given a Polars dataframe, filter to include only the chromosomes in chr_list"""
        df = df.filter(pl.col("chromosome").is_in(self.chr_list))

        return df

    def parse_attribute_str(self, attribute_str: str):
        """Parses the attribute string in a GFF file into a dictionary"""
        attributes = {}
        for attribute in attribute_str.split(";"):
            key, value = attribute.split("=")
            attributes[key] = value

        return attributes

    def format_gff_attributes(self):
        """Format the GFF attributes into a dataframe

        The GFF attribute string is a semi-colon separated list of key-value pairs,
        each of which is separated by an = sign. However, the number of key-value pairs
        can vary between lines, so this creates a dictionary for each attribute string
        and then places them into a column-accurate polars dataframe
        """
        # This is SUPER janky but I guess it works well enough
        attribute_dicts = [
            self.parse_attribute_str(attribute) for attribute in self.gff["attribute"]
        ]

        unique_keys = set(key for attribute in attribute_dicts for key in attribute.keys())

        columns_dict = {key: [] for key in unique_keys}

        for attribute in attribute_dicts:
            for key in unique_keys:
                columns_dict[key].append(attribute.get(key, None))

        attribute_df = pl.DataFrame(columns_dict)

        return attribute_df

    def calculate_chr_list(self, pareto_cutoff=0.8):
        """
        Use the Pareto principle to calculate the
        chromosomes that account for 80% of the insert size data

        This is a proxy for the most important chromosomes in each genome
        for this specific analysis case, but likely doesn't capture all
        chromosome pieces very well in highly fragmented genomes
        """
        # Count the number of reads per chromosome in the bam file
        chr_counts = defaultdict(int)
        # Use pysam count to get the number of reads per chromosome
        for read in self.bam:
            chr_counts[read.reference_name] += 1

        # Sort the chromosomes by the number of reads
        chr_counts = {
            k: v for k, v in sorted(chr_counts.items(), key=lambda item: item[1], reverse=True)
        }
        # Calculate the cumulative sum of reads
        chr_counts = {k: v / sum(chr_counts.values()) for k, v in chr_counts.items()}
        chr_cumulative = 0
        chr_list = []
        for k, v in chr_counts.items():
            chr_cumulative += v
            chr_list.append(k)
            if chr_cumulative >= pareto_cutoff:
                break
        self.chr_list = chr_list

    def generate_insert_summary(self, maximum_size_threshold=100000):
        """
        Generate a summary of the positions and sizes of the
        most common inserts in the genome

        The maximum_size_threshold denoting the largest insert size
        this function will consider is derived from the parasitoid wasp
        information. For two species with different integrated polydnaviruses,
        the largest amplified regions are around 100kb
        """
        insert_summary = pl.DataFrame()
        for chromosome in self.chr_list:
            retained_read_data = []
            for read in self.bam.fetch(chromosome):
                retained_read_dict = {
                    "name": read.query_name,
                    "reference_name": read.reference_name,
                    "reference_start": read.reference_start,
                    "reference_end": read.reference_end,
                    "length": read.template_length,
                    "cigar": read.cigarstring,
                }
                retained_read_data.append(retained_read_dict)

                if len(retained_read_data) > 10000:
                    retained_read_df = pl.DataFrame(retained_read_data)
                    insert_summary = self.process_read_chunk(
                        retained_read_df,
                        insert_summary,
                        maximum_size_threshold=maximum_size_threshold,
                    )
                    retained_read_data = []

            if retained_read_data:
                retained_read_df = pl.DataFrame(retained_read_data)
                insert_summary = self.process_read_chunk(
                    retained_read_df, insert_summary, maximum_size_threshold=maximum_size_threshold
                )

        # Sort by chromosome and start position
        insert_summary = insert_summary.sort(["reference_name", "rounded_start"])

        return insert_summary

    def process_read_chunk(self, retained_reads, insert_summary, maximum_size_threshold=100000):
        "Rather than process all reads at once, process them in chunks to save memory"
        retained_reads = retained_reads.with_columns(
            [
                (pl.col("length").abs().alias("abs_length")),
                (
                    pl.col("length")
                    .abs()
                    # Round insert size to nearest 1kb
                    .map_elements(lambda x: round(x, -3), return_dtype=pl.Int64)
                    .alias("rounded_length")
                ),
                (
                    pl.col("reference_start")
                    # Round start position to nearest 1kb
                    .map_elements(lambda x: round(x, -3), return_dtype=pl.Int64)
                    .alias("rounded_start")
                ),
                (
                    pl.when(pl.col("length") > 0)
                    .then(pl.lit("+"))
                    .otherwise(pl.lit("-"))
                    .alias("strand")
                ),
            ]
        )

        # Filter out the reads that are too large
        retained_reads = retained_reads.filter(pl.col("abs_length") <= maximum_size_threshold)

        # Summarize the chunk, counting the number of inserts per
        # chromosome, strand, start and length
        chunk_summary = retained_reads.group_by(
            ["reference_name", "strand", "rounded_start", "rounded_length"]
        ).len()

        # Calculate the end position of the insert
        chunk_summary = chunk_summary.with_columns(
            pl.when(pl.col("strand") == "+")
            .then(pl.col("rounded_start") + pl.col("rounded_length"))
            .otherwise(pl.col("rounded_start") - pl.col("rounded_length"))
            .alias("rounded_end"),
        )

        if not insert_summary.is_empty():
            insert_summary = pl.concat([insert_summary, chunk_summary])
        else:
            insert_summary = chunk_summary

        return insert_summary

    def generate_insert_range(self, insert_summary_df, width_bp=20000):
        """
        Given an insert summary df, generate the genome position range
            to search coverage and annotations. This is wider than the
            insert size, because circularized dsDNA will always be smaller
            than the amplified DNA surrounding it and the insert sizes are
            probably only picking up the circularized components.
            I set the default to 20kb in both directions
        """
        insert_range = insert_summary_df.with_columns(
            [
                (
                    pl.when(pl.col("rounded_end") < pl.col("rounded_start"))
                    .then(pl.col("rounded_end") - width_bp)
                    .otherwise(pl.col("rounded_start") - width_bp)
                    .alias("start_range")
                ),
                (
                    pl.when(pl.col("rounded_end") < pl.col("rounded_start"))
                    .then(pl.col("rounded_start") + width_bp)
                    .otherwise(pl.col("rounded_end") + width_bp)
                    .alias("end_range")
                ),
            ]
        )

        # If the start or end range is negative, set it to 0
        # Otherwise the search for coverage will hang
        insert_range = insert_range.with_columns(
            [
                (
                    pl.when(pl.col("start_range") < 0)
                    .then(pl.lit(0))
                    .otherwise(pl.col("start_range"))
                    .alias("start_range")
                ),
                (
                    pl.when(pl.col("end_range") < 0)
                    .then(pl.lit(0))
                    .otherwise(pl.col("end_range"))
                    .alias("end_range")
                ),
            ]
        )

        return insert_range

    def deduplicate_insert_summary(self, insert_summary_df, z_score_threshold=3):
        """
        Given an insert summary df, deduplicate the inserts that
        overlap with each other. This is useful for finding the
        insert peaks - positions where there are many inserts, which are
        statistical outliers compared to the rest of the chromosome (which
        may have a few inserts at random positions, but not all in one place)
        """
        # Group the insert summary by chromosome, start and end range
        insert_summary_df = insert_summary_df.group_by(
            ["reference_name", "start_range", "end_range"]
        ).agg(pl.col("len").sum().alias("count"))

        # Make sure there are no duplicates
        insert_summary_df = insert_summary_df.unique(
            subset=["reference_name", "start_range", "end_range"]
        )

        # Calculate total number of inserts per chromosome
        total_inserts_by_chromosome = insert_summary_df.group_by("reference_name").agg(
            pl.col("count").sum().alias("total_count")
        )

        # Merge the total count back into the insert summary
        insert_summary_df = insert_summary_df.join(total_inserts_by_chromosome, on="reference_name")

        # Calculate the proportion of the total count per chromosome
        # This is a nice-to-have column, I find it useful for looking at
        # how many inserts I'm finding on each chromosome
        insert_summary_df = insert_summary_df.with_columns(
            [(pl.col("count") / pl.col("total_count")).alias("proportion")]
        )

        # Calculate the z-score of the number of inserts for each chromosome and insert range
        insert_summary_df = self.calculate_z_scores(insert_summary_df)

        # Filter to only include inserts that are greater than the z-score threshold
        insert_summary_df = insert_summary_df.filter(pl.col("z_score") > z_score_threshold)

        # Also filter by number, if there are fewer than 100 inserts on a chromosome don't bother
        insert_summary_df = insert_summary_df.filter(pl.col("total_count") > 100)

        return insert_summary_df

    def calculate_z_scores(self, insert_summary_df):
        """
        Given an insert summary df, calculate the z-scores of the
        insert sizes
        """
        # Calculate the z-score of the number of inserts at a given range
        grouped_df = insert_summary_df.group_by("reference_name").agg(
            [pl.col("count").mean().alias("mean"), pl.col("count").std().alias("std")]
        )

        insert_summary_df = insert_summary_df.join(grouped_df, on="reference_name")

        insert_summary_df = insert_summary_df.with_columns(
            [((pl.col("count") - pl.col("mean")) / pl.col("std")).alias("z_score")]
        )

        insert_summary_df = insert_summary_df.drop(["mean", "std"])

        return insert_summary_df

    def generate_filtering_positions(self, deduped_insert_summary):
        """
        Given an (ideally) deduplicated insert summary df, generate
        the positions to filter the coverage and annotations. Doing this by
        generating a dictionary of chromosome names and the positions within
        that chromosome
        """
        range_dict = {}
        for chromosome in deduped_insert_summary["reference_name"].unique().to_list():
            filtered_summary_df = deduped_insert_summary.filter(
                pl.col("reference_name") == chromosome
            )
            range_lists = [
                list(range(start, end + 1))
                for start, end in zip(
                    filtered_summary_df["start_range"], filtered_summary_df["end_range"], 
                    strict=False
                )
            ]
            range_list = list(set([item for sublist in range_lists for item in sublist]))
            range_dict[chromosome] = range_list

        return range_dict

    def filter_coverage(self, range_dict):
        """Filter the coverage dataframe to only include the positions
        that are within the range list per chromosome
        """
        # For each chromosome, filter the coverage dataframe
        coverage_dfs = []
        for chromosome in range_dict.keys():
            positions = range_dict[chromosome]
            filtered_coverage = self.coverage.filter(pl.col("chromosome") == chromosome)
            filtered_coverage = filtered_coverage.filter(pl.col("position").is_in(positions))
            # If dataframe is empty, skip, otherwise append
            if len(filtered_coverage) > 0:
                coverage_dfs.append(filtered_coverage)
        # If there are no dataframes, return an empty dataframe
        if len(coverage_dfs) == 0:
            return pl.DataFrame()
        else:
            return pl.concat(coverage_dfs)

    def filter_gff(self, range_dict):
        """Filter the gff dataframe to only include the positions
        that are within the range list
        """
        gff_dfs = []
        for chromosome in range_dict.keys():
            positions = range_dict[chromosome]
            filtered_gff = self.gff.filter(pl.col("chromosome") == chromosome)
            # This checks for starts of genes that are within the positions
            filtered_gff = filtered_gff.filter(pl.col("start").is_in(positions))
            if len(filtered_gff) > 0:
                gff_dfs.append(filtered_gff)

        if len(gff_dfs) == 0:
            return pl.DataFrame()
        else:
            return pl.concat(gff_dfs)
