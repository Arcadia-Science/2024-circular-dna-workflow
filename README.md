# 2024-circular-dna-workflow

This repository contains code for pulling down datasets and a workflow for automating downloading files, mapping, and formatting for getting coverage of paired-end Illumina sequencing experiments against a reference genome assembly.

This repo is part of the pub [Identifying circular DNA using short-read mapping](https://research.arcadiascience.com/pub/method-circular-DNA-ID). 


## Installation and Setup
This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.anaconda.com/free/miniconda/miniconda-install/). After installing conda and mamba, run the following command to create the pipeline run environment.
```
git clone git@github.com:Arcadia-Science/2024-parasite-capsids-workflow.git
conda create -f environment.yml -n circular_dna_finder
```
This environment includes dependencies for scripts and notebooks, as well as installing Nextflow to run the mapping/coverage parsing workflow.

To run the parasite_species_genomeinfo.ipynb notebook, you will need to set up a separate conda environment with the following command:
```
conda create -f envs/genomeinfo.yml -n genomeinfo
```

## Data
Input data:
* The input data for this workflow is a sample sheet with three columns:
    * `genome_accession`: The accession ID for the genome assembly to map against.
    * `genome_ftp_path`: The NCBI FTP path to the genome assembly file.
    * `SRA_run_accession`: The SRA run accession ID for the paired-end Illumina sequencing data to map against the genome assembly.

We've provided two example samplesheets in the `inputs` directory, one under `parasitoid_wasps/parasitoid_wasps_samplesheet.csv` and one under `parasite_species/parasite_related_species_samplesheet.csv`. The samplesheets contain the necessary information to run the workflow for the respective datasets.

If you'd like to download the Nextflow pipeline outputs for our two examples, please find them on Zenodo [here](https://zenodo.org/doi/10.5281/zenodo.13362362). Then place these outputs in `results`.

Output data:
* The Nextflow workflow outputs several files:
    *.sorted.bam: A sorted BAM file of all mapped reads
    *.large_inserts.bam: A filtered BAM file of only reads with mapped distances ≥ 1 kb
    *.coverage.txt: A tab-delimited file of coverage depth for each position in the genome
    *.average_coverage.txt: A tab-delimited file of per-scaffold average coverage depth
    *.high_coverage_regions.txt: A tab-delimited file positions with coverage depth ≥ 100× average scaffold coverage. Also available in BED format
    *.high_coverage_region_sizes.txt: A tab-delimited file of high-coverage depth regions (≥ 100× average scaffold coverage depth, with positions within 100 bp merged into a region), with columns corresponding to scaffold name, start position, end position, and total region size (bp)
    *.filtered_coverage.txt: A tab-delimited file of coverage depth only at positions where we identified mapped distances ≥ 1 kb

* If you are working with the outputs using the `GenomeInfo` class, you can also write several of the outputs to file. Examples of this can be found in `notebooks/parasite_species_genomeinfo.ipynb`. Specifically, you can write:
    * *.deduped_insert_summary.csv: A CSV file of the z-score filtered inserts per scaffold per genome
    * *.insert_filtered_coverage.csv: A CSV file of coverage depth, filtered to only positions within detected large inserts
    * *.insert_filtered_gff.csv: A CSV file of provided annotations, filtered to only start positions within detected large inserts
    * These are only example outputs - if you have unique needs for your data, you can modify the `GenomeInfo` class to write out additional files.

## Overview

### Description of the folder structure
* `inputs`: Contains genome assemblies, annotations, and sample sheets for the Nextflow pipeline
    * `parasite_species`: Contains the sample sheet and input data for the parasite/related species dataset
    * `parasitoid_wasps`: Contains the sample sheet and input data for the parasitoid wasp dataset
* `notebooks`: Contains Jupyter notebooks for analyzing the output of the Nextflow pipeline for the parasitoid wasp and parasite/related species datasets
* `results`: Should contain filtered BAM and coverage files, figures, and `GenomeInfo`-filtered example files for the parasitoid wasp and parasite/related species datasets. **NOTE**: Many of these files are too large for Github, so the samplesheets should contain the necessary information to generate these files. If you'd rather run the examples without generating the data with Nextflow, then the complete set of outputs can be found on [Zenodo](https://zenodo.org/doi/10.5281/zenodo.13362362). You can download the dataset and place the files in the appropriate directories within `results`:
    * `parasite_species`:
        * `bam_files`: Contains BAM files with mapped distances >= 1 kb
        * `coverage_files`: Contains coverage depth files filtered by the BAM files
        * `insert_filtered_results`: Contains examples of filtered file outputs from the `GenomeInfo` class. The complete set of outputs can be found on Zenodo.
        * `fig`: Figures generated for the pub.
    * `parasitoid_wasps`:
        * `bam_files`: Contains BAM files with mapped distances >= 1 kb
        * `coverage_files`: Contains coverage depth files filtered by the BAM files
        * `fig`: Figures generated for the pub.
        * `hyposoter_didymator_blastx_results`: Contains the manual BLASTx results from the _Hyposoter didymator_ search.

### Methods

#### Nextflow Workflow Documentation

The software dependencies within the workflow can be executed through either Docker containers or conda environments. The recommended way for running the workflow locally is through conda environments since this is easier to install and deploy. If you choose to go the Docker route, install Docker [according to these instructions for your operating system](https://docs.docker.com/engine/install/). The below assumes you are running via conda.


After setting up your conda environment that includes installing Nextflow, you are now ready to launch the workflow.
**IMPORTANT NOTE**: If you're using `conda` to manage your environments, please manually install `sra-tools` on your system or point your path to recognize it in your conda environment. 

The input CSV samplesheet should look like the following (also see the CSV samplesheet in `inputs/`):
```
genome_accession,genome_ftp_path,SRA_run_accession
GCF_026212275.2_iyMicDemo2.1a,https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/212/275/GCF_026212275.2_iyMicDemo2.1a,SRR3420509
GCF_026212275.2_iyMicDemo2.1a,https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/212/275/GCF_026212275.2_iyMicDemo2.1a,SRR3420507
```

The three columns `genome_accession`, `genome_ftp_path`, and `SRA_run_accession` are all required. The `genome_accession` must match the full accession name that is at the end of the FTP download path for a given assembly input. The FTP path is the full NCBI FTP download path for the entire accession, but not including the specific file to be downloaded, (such as not including GCF_026212275.2_iyMicDemo2.1a_genomic.fna.gz in the download link field) as this is done within the workflow code. The `SRA_run_accession` is just the accession and not the FTP download link.

For every pair of mapping results that you want, you must list them with each `genome_accession` and `SRA_run_accession` pairing. For example the above gives two rows of the same `M. demolitor` genome accession but different SRA run accessions. The workflow will sense duplicate genome accessions in the first column and only download it one time, but the different genome:SRA accession pairing need to be listed like the above to ensure proper mapping pairings. Therefore if you have multiple SRA run accessions to map against the same genome, you can copy columns 1 and 2 but change the SRA run accession in the third column in your samplesheet.

Now you are ready to launch the workflow:
```
nextflow run main.nf \\
    --samples samples.csv \\
    --outdir results_directory \\
    --threads 6 \\
    -profile conda
```
The only necessary inputs are your CSV samplesheet formatted such as the above, the path to you results directory, and how many threads you have available for individual processes. The workflow doesn't copy over every intermediate result file to the output directory, so if there is an intermediate file that you need or want to check out, you can look in the `work/` directory and for example with `ls */*.bam` find where all BAM files are located in the subdirectories of the `work` directory. Not all intermediate files are in the results output directory because of how Nextflow looks for existing files when resuming an incomplete run, and results are copied twice if you want them in your results output directory. This can eat up space on a machine fairly quickly when you have many samples with large FASTQ, BAM, etc. files.

If you have a run that is interrupted, or you add additional samples and don't want to rerun existing samples, you can add the `-resume` flag to your run:
```
nextflow run main.nf \\
    --samples samples.csv \\
    --outdir results_directory \\
    --threads 6 \\
    -profile conda \\
    -resume
```

This can also be used if there were SRA run accessions that failed to download due to HTTP timeout errors with sratool-kit and you need to retry downloading those samples. The workflow is configured that if they fail to download those runs will be ignored and the workflow will keep going for successful downloads. However if other steps fail for other reasons this will cause the workflow to hault.

##### Compute Specifications
The Nextflow workflow was primarily built and tested on a MacBook Pro running macOS Monterey v12.5.1 with 32 GB of RAM and using 6 threads. The workflow successfully completes for two _Microplitis demolitor_ test samples in about 10 minutes. All samples were processed running the workflow on Nextflow Tower via AWS Batch and using -profile docker.

#### `GenomeInfo` Class Documentation

The `GenomeInfo` Python class can be used to parse the outputs of the Nextflow pipeline. This class takes in several mandatory files and several optional files, including:
- `large_insert_bam_loc` (str, required): The path to the BAM file with mapped distances ≥ 1 kb
- `large_insert_bai_loc` (str, required): The path to the BAM index file
- `coverage_loc` (str, required): The path to the coverage depth file. This can either be the complete coverage data or the filtered coverage data (recommended)
- `genome_metadata` (dict, optional): A dictionary of metadata about the specific genome you're analyzing. This can be useful for keeping track of multiple genomes in a single analysis.
- `gff_loc` (str, optional): The path to the GFF file for the genome assembly, if it exists
- `chr_list` (list, optional): A list of chromosome names you're interested in analyzing in the genome assembly. If not provided, the class will calculate which scaffolds contribute to 80% of the insert size data and filter the outputs accordingly.

To instantiate the class, you can use the following code:
```python
from scripts.genomeinfo import GenomeInfo

genome_of_interest = GenomeInfo(
    large_insert_bam_loc='path/to/large_inserts.bam',
    large_insert_bai_loc='path/to/large_inserts.bam.bai',
    coverage_loc='path/to/filtered_coverage.txt',
    genome_metadata={'species': 'Microplitis demolitor'},
    gff_loc='path/to/genome.gff',
    chr_list=['scaffold_1', 'scaffold_2']
)
```

`GenomeInfo` has several methods to handle and filter your data:
- `.load_bam()`: Loads the filtered BAM file into the class. Optionally, you can provide a `pareto_cutoff` parameter (default = 0.8) to filter the scaffolds that contribute to X% of the insert size data.
- `.load_coverage()`, `.load_gff()`: Loads the coverage depth and GFF files into the class.
- `.generate_insert_summary()`: Generates a summary of the insert sizes for each scaffold in the genome assembly. Optionally, you can provide a `maximum_size_threshold` parameter (default = 100000) to only include inserts up to a certain size in base pairs.
- `.generate_insert_range()`: Generates a range surrounding inserts to use for filtering coverage and annotations. Optionally, you can provide a `width_bp` parameter (default = 20000) to only include inserts up to a certain size in base pairs.
- `.deduplicate_insert_summary()`: Finds statistical outliers in the insert summary/range data, looking for positions in the scaffold with an outsized number of large inserts. Optionally, you can provide a `z_score_threshold` parameter (default = 3) to set the threshold for z-score filtering.
- `.filter_coverage()`: Filters the coverage data to only include positions within the insert range.
- `.filter_gff()`: Filters the GFF data to only include annotations within the insert range.

Example use cases for these methods can be found in `notebooks/parasite_species_genomeinfo.ipynb`.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).