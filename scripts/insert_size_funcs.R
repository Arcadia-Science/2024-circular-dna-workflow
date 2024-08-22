# Main functions for dealing with BAM files and visualizing insert size distribution

library(tidyverse)
library(Rsamtools)
library(data.table)

# Process a BAM file and turn it into a data table
process_bam <- function(bam_file_loc) {
    bam_file <- Rsamtools::scanBam(bam_file_loc)
    # Need to pull out dataframe before converting to dt for this to work properly
    bam_df <- data.frame(bam_file[[1]])
    bam_dt <- data.table(bam_df)
    # Convert all insert sizes to positive, since strandedness doesn't matter
    bam_dt[, isize := abs(isize)] 

    return(bam_dt)
}

# Add clearer column names to a BAM data table and
# if provided with a list of chromosome names, filter
# the data table to only include those chromosomes 
format_bam <- function(bam_dt, chromosome_names = NULL) {
    # Ensure typing of contig name and position
    bam_dt[, `:=`(chromosome = as.character(rname), position = as.numeric(pos))]

    if (!is.null(chromosome_names)) {
        bam_dt <- bam_dt[chromosome %in% chromosome_names]
    }

    return(bam_dt)
}

# Add clearer column names to a coverage depth data table
# and if provided with a list of chromosome names, filter
# the data table to only include those chromosomes
format_coverage <- function(coverage_dt, chromosome_names = NULL) {
    # Name and ensure typing of columns
    setnames(coverage_dt, c("V1", "V2", "V3"), c("chromosome", "position", "coverage"))
    coverage_dt[, c("position", "coverage") := .(as.numeric(position), as.numeric(coverage))] 

    if (!is.null(chromosome_names)) {
        coverage_dt <- coverage_dt[chromosome %in% chromosome_names]
    }

    return(coverage_dt)
}

# Generate a summary of insert sizes from a BAM data table
# This is done by breaking the insert sizes into bins and 
# counting the number of reads in each bin
# Bins are broken up by 1kb regions
# Data are grouped by chromosome, start position, and insert size
generate_insert_summary <- function(bam_dt,
                                    max_chromosome_length = 1000000000,
                                    window_size = 1000) {
    # TODO: Add a check on strandedness
    insert_summary <- bam_dt %>% # nolint
    mutate(insert_start_bin = cut(pos, breaks = seq(0,
                                                    max_chromosome_length,
                                                    by = window_size)),
           insert_size_bin = cut(isize, breaks = seq(0,
                                                     max_chromosome_length,
                                                     by = window_size))) %>%
    group_by(rname, insert_start_bin, insert_size_bin) %>%
    summarize(num_inserts = n()) %>%
    ungroup() %>%
    separate(insert_start_bin, c("min_insert_start", "max_insert_start"), sep = ",") %>%
    separate(insert_size_bin, c("min_insert_size", "max_insert_size"), sep = ",") %>%
    mutate(min_insert_start = gsub("\\(|\\)", "", min_insert_start),
           max_insert_start = gsub("\\]|\\)", "", max_insert_start),
           min_insert_start = as.numeric(min_insert_start),
           max_insert_start = as.numeric(max_insert_start)) %>%
    mutate(min_insert_size = gsub("\\(|\\)", "", min_insert_size),
           max_insert_size = gsub("\\]|\\)", "", max_insert_size),
           min_insert_size = as.numeric(min_insert_size),
           max_insert_size = as.numeric(max_insert_size))

    return(insert_summary)
}

# Find the most common insert size(s) for each chromosome given a summary table
# In theory, this should give the most common n for how_many,
# but something with the top_n function is wonky
# TODO: fix this weirdness
find_most_common_insert <- function(insert_summary, how_many = 1) {
    # TODO: Add a check on strandedness
    most_common_insert <- insert_summary %>%
    group_by(rname, min_insert_start, max_insert_start) %>%
    top_n(how_many, num_inserts) %>%
    ungroup() %>%
    group_by(rname) %>%
    top_n(1, num_inserts) %>%
    ungroup()

    return(most_common_insert)
}

# Generate a chromosome plot 
# with the most common insert as an overlaid segment
plot_chr_cov_with_insert <- function(chromosome_name, coverage_dt, most_common_insert_df) {
    # Find start position minimum, end position maximum, and calculate length
    start_pos <- filter(most_common_insert_df, rname == chromosome_name)$min_insert_start
    insert_len <- filter(most_common_insert_df, rname == chromosome_name)$max_insert_size
    end_pos <- start_pos + insert_len

    # Create the plot
    plt <- coverage_dt %>%
    filter(chromosome == chromosome_name) %>%
    ggplot(aes(x = position, y = coverage)) +
    geom_line() +
    theme_classic() +
    labs(x = "Position (bp)", y = "Coverage") +
    scale_y_continuous(labels = function(x) paste0(x, "x")) +
    annotate("segment", y = 1, x = start_pos, xend = end_pos, color = "#FFB984", size = 8) +
    ggtitle(chromosome_name) +
    theme(text = element_text(size = 18))
    
    return(plt)
}

# Iterate over list of chromosomes, plot with inserts, and store in list
plot_all_chr_cov_with_inserts <- function(chromosome_list, coverage_dt, most_common_insert_df) {
    plt_list <- list()
    
    for (chromosome in chromosome_list) {
        plt <- plot_chr_cov_with_insert(chromosome, coverage_dt, most_common_insert_df)
        plt_list[[chromosome]] <- plt
    }
    
    return(plt_list)
}
