#! /usr/bin/env nextflow

// Description
// Download reference genomes and Illumina seq reads for mapping and coverage operations
// This workflow script is for running on Tower with the Docker profile


nextflow.enable.dsl=2

log.info """\

DOWNLOAD, MAP, AND PROCESS ILLUMINA MAPPING AGAINST REFERENCE GENOMES
TO RETRIEVE HIGH COVERAGE AREAS
=========================================
samples            : $params.samples
outdir             : $params.outdir
threads            : $params.threads
"""

// channels for CSV of SRA Illumina samples and genome assembly links to download, mapping queries
// this assumes reads are paired-end
sample_csv = Channel.fromPath(params.samples)
    .splitCsv(header:true)

// genome accession channel
genome_accession = sample_csv.map { [it.genome_accession, it.genome_ftp_path ]}
    .unique() // duplicates because samplesheet is for all mapping pairs

// fastq reads channel
paired_end_reads = sample_csv.map { [it.genome_accession, it.SRA_run_accession ]}


workflow {
    // download SRA files
    downloaded_reads = download_paired_SRA_accessions(paired_end_reads)

    // download genome accession
    downloaded_assembly = download_genome_accession(genome_accession)

    // QC reads
    qced_reads = fastp_filter_reads(downloaded_reads)

    // minimap2 index
    assembly_index = make_minimap2_index(downloaded_assembly)

    // map reads
    mapping_samples = qced_reads
        .combine(assembly_index, by: 0) // combines by genome accession for correct pairings

    mapped_reads = minimap2_mapping(mapping_samples)

    // samtools depth
    samtools_insert_depth_files = get_samtools_insert_depth(mapped_reads)

    // awk commands to get specific fields
    full_coverage_results = samtools_insert_depth_files.coverage_txt
    coverage_files = coverage_filtering_reformatting(full_coverage_results)

    // extract chromosome and position from large_inserts BAM file
    extracted_large_inserts = extract_chromosome_position(samtools_insert_depth_files.large_inserts_bam)

    // write matching lines between positions and coverage files
    // combine by genome accession and SRA run accession to match correct file pairings
    merged_pairings = extracted_large_inserts
        .combine(full_coverage_results, by: [0,1])

    merged_files = merge_positions_coverage(merged_pairings)
}

// download paired-end reads with SRA tools passing the SRA run accession
process download_paired_SRA_accessions {
    // download each SRA run with SRAtools
    tag "${SRA_run_accession}_download"

    errorStrategy 'ignore' // ignore failed downloads because could have NCBI FTP issues and can retry on failed ones later

    container "quay.io/staphb/sratoolkit:latest"
    // conda "envs/sratoolkit.yml"
    // If run on a Mac you will need to have a local install of sratools, if run on Linux EC2 instance run with the conda sratools package and uncomment the above

    input:
    tuple val(genome_accession), val(SRA_run_accession)

    output:
    tuple val(genome_accession), val(SRA_run_accession), path("*.fastq"), emit: fastq

    script:
    """
    fasterq-dump ${SRA_run_accession}
    """
}

// download genome assembly
process download_genome_accession {
    tag "${genome_accession}_download"
    publishDir "${params.outdir}/genome_assemblies", mode: 'copy', pattern: "*.fna"

    container "quay.io/biocontainers/wget:1.20.1"

    input:
    tuple val(genome_accession), val(genome_ftp_path)

    output:
    tuple val(genome_accession), path("*.fna"), emit: fasta

    script:
    """
    wget ${genome_ftp_path}/${genome_accession}_genomic.fna.gz

    gunzip ${genome_accession}_genomic.fna.gz
    """

}

// filter reads
process fastp_filter_reads {
    tag "${SRA_run_accession}_fastp"

    container "quay.io/biocontainers/fastp:0.23.4--hadf994f_3"
    conda 'envs/fastp.yml'

    input:
    tuple val(genome_accession), val(SRA_run_accession), path(SRA_fastq)

    output:
    tuple val(genome_accession), val(SRA_run_accession), path("*.qced.fastq.gz"), emit: filtered_fastq

    script:
    """
    fastp \\
        --in1 ${SRA_run_accession}_1.fastq \\
        --in2 ${SRA_run_accession}_2.fastq \\
        --out1 ${SRA_run_accession}_1.qced.fastq.gz \\
        --out2 ${SRA_run_accession}_2.qced.fastq.gz \\
        --json ${SRA_run_accession}.fastp.json \\
        --html ${SRA_run_accession}.fastp.html \\
        --cut_front \\
        --cut_tail \\
        --cut_mean_quality 15 \\
        -q 15 \\
        -w ${params.threads}
    """
}

// make minimap2 index of genome assembly
process make_minimap2_index {
    tag "${genome_accession}_index"

    container "quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0"
    conda 'envs/minimap.yml'

    input:
    tuple val(genome_accession), path(assembly)

    output:
    tuple val(genome_accession), path("*.mmi"), emit: index

    script:
    """
    minimap2 \\
        -t ${params.threads} \\
        -d ${genome_accession}.mmi \\
        ${assembly}

    """
}

// map reads with minimap2, sort
process minimap2_mapping {
    tag "${genome_accession}_vs_${SRA_run_accession}_mapping"

    container "quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0"
    conda 'envs/minimap.yml'

    input:
    tuple val(genome_accession), val(SRA_run_accession), path(reads), path(genome_index)

    output:
    tuple val(genome_accession), val(SRA_run_accession), path("*.bam"), path("*.bai"), emit: sorted_indexed_bam

    script:
    """
    minimap2 \\
        -ax sr \\
        ${genome_index} \\
        ${reads} | samtools sort -o ${genome_accession}_vs_${SRA_run_accession}.sorted.bam -

    samtools index ${genome_accession}_vs_${SRA_run_accession}.sorted.bam
    """
}

// samtools depth
process get_samtools_insert_depth {
    tag "${genome_accession}_vs_${SRA_run_accession}_depth"
    publishDir "${params.outdir}/bam_files", mode: 'copy', pattern: "*.bam"
    publishDir "${params.outdir}/bam_files", mode: 'copy', pattern: "*.bai"
    publishDir "${params.outdir}/coverage_results", mode: 'copy', pattern: "*.txt"

    container "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
    conda 'envs/samtools.yml'

    input:
    tuple val(genome_accession), val(SRA_run_accession), path(sorted_bam), path(indexed_bam)

    output:
    tuple val (genome_accession), val(SRA_run_accession), path("*.bam"), path("*.bai"), emit: large_inserts_bam
    tuple val(genome_accession), val(SRA_run_accession), path("*.txt"), emit: coverage_txt

    script:
    """
    samtools view -h -F 4 ${sorted_bam} | awk '(\$9 > 1000 || \$9 < -1000)' | samtools sort -o ${genome_accession}_vs_${SRA_run_accession}.large_inserts.bam

    samtools index ${genome_accession}_vs_${SRA_run_accession}.large_inserts.bam

    samtools depth -a ${sorted_bam} > ${genome_accession}_vs_${SRA_run_accession}.coverage.txt
    """
}

// coverage filtering and reformatting
process coverage_filtering_reformatting {
    tag "${genome_accession}_vs_${SRA_run_accession}_reformatting"
    publishDir "${params.outdir}/coverage_results", mode: 'copy', pattern: "*.txt"
    publishDir "${params.outdir}/coverage_results", mode: 'copy', pattern: "*.bed"

    container "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1"
    conda 'envs/bedtools.yml'

    input:
    tuple val(genome_accession), val(SRA_run_accession), path(coverage_file)

    output:
    path("*.average_coverage.txt"), emit: average_coverage_txt
    path("*.high_coverage_regions.txt"), emit: high_coverage_txt
    path("*.high_coverage_regions.bed"), emit: high_coverage_bed
    path("*.merged_high_coverage_regions.bed"), emit: merged_bed
    path("*.high_coverage_region_sizes.txt"), emit: high_coverage_sizes_txt

    script:
    """
    # get average coverage per scaffold
    awk '{cov[\$1]+=\$3; count[\$1]++} END{for (contig in cov) print contig, cov[contig]/count[contig]}' ${coverage_file} > ${genome_accession}_vs_${SRA_run_accession}.average_coverage.txt

    # find regions with high coverage
    awk 'NR==FNR {avg[\$1]=\$2; next} {if (\$3 >= 100*avg[\$1]) print \$0}' ${genome_accession}_vs_${SRA_run_accession}.average_coverage.txt ${coverage_file} > ${genome_accession}_vs_${SRA_run_accession}.high_coverage_regions.txt

    # reformat for bedtools
    awk '{print \$1 "\t" \$2-1 "\t" \$2 "\t" \$3}' ${genome_accession}_vs_${SRA_run_accession}.high_coverage_regions.txt > ${genome_accession}_vs_${SRA_run_accession}.high_coverage_regions.bed

    # merge intervals within 100bp
    bedtools merge -i ${genome_accession}_vs_${SRA_run_accession}.high_coverage_regions.bed -d 100 -c 4 -o mean > ${genome_accession}_vs_${SRA_run_accession}.merged_high_coverage_regions.bed

    # calculate region sizes
    awk '{print \$1, \$2, \$3, \$3-\$2, \$4}' ${genome_accession}_vs_${SRA_run_accession}.merged_high_coverage_regions.bed > ${genome_accession}_vs_${SRA_run_accession}.high_coverage_region_sizes.txt
    """
}

process extract_chromosome_position {
    tag "${genome_accession}_vs_${SRA_run_accession}_extract_chromosomes_pos"
    publishDir "${params.outdir}/filtered_results", mode: 'copy', pattern: "*.txt"

    container "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
    conda 'envs/samtools.yml'

    input:
    tuple val(genome_accession), val(SRA_run_accession), path(large_inserts_bam), path(indexed_inserts_bam)

    output:
    tuple val(genome_accession), val(SRA_run_accession), path("*.positions.txt"), emit: positions_file

    script:
    """
    samtools view ${large_inserts_bam} | awk '{print \$3 "\t" \$4}' > ${genome_accession}_vs_${SRA_run_accession}.positions.txt
    """
}

process merge_positions_coverage {
    tag "${genome_accession}_vs_${SRA_run_accession}_merge_tables"
    publishDir "${params.outdir}/filtered_results", mode: 'copy', pattern: "*.txt"

    container "quay.io/biocontainers/samtools:1.20--h50ea8bc_0" // process doesn't actually use samtools but Tower needs a docker container
    // consider merging this step with the above for simplicity

    input:
    tuple val(genome_accession), val(SRA_run_accession), path(extracted_large_inserts), path(coverage_txt)

    output:
    tuple val(genome_accession), val(SRA_run_accession), path("*.filtered_coverage.txt"), emit: filtered_coverage

    script:
    """
    awk 'NR==FNR {pos[\$1 "\t" \$2]; next} (\$1 "\t" \$2) in pos' ${extracted_large_inserts} ${coverage_txt} > ${genome_accession}_vs_${SRA_run_accession}.filtered_coverage.txt
    """
}
