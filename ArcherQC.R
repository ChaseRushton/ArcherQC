#!/usr/bin/Rscript
################################################################################
# Developer Notes
################################################################################
# Development Environment:
# - R version: >= 4.0.0
# - Platform: Windows/Linux
# - Required Memory: 8GB minimum
#
# Dependencies:
# - tidyverse (data manipulation)
# - RCurl (API requests)
# - jsonlite (JSON processing)
#Troubleshooting:
# - make sure folder is empty except molbar files if rerunning 
# - Ensure all input files exist in specified paths
# - Check file permissions for output directories
# - Verify historic data file format matches expected
#
# Contact:
# ArcherQCv9.R
# - Created: 05/06/2024
# - Developer: Chase Rushton
# - Email: chase.rushton@pennmedicine.upenn.edu
# - Last Updated: 07/18/2024

# Description:
# This script performs comprehensive QC analysis on Archer NGS data by:
# 1. Processing molbar API data and historical records
# 2. Applying multiple QC thresholds and validations
# 3. Generating detailed reports and visualizations
# 4. Tracking failures across multiple metrics
################################################################################

#===============================================================================
# Changelog
#===============================================================================
# Version 1.0.0 (05/06/2024)
# - Initial release
# - Core Features:
#   * Data import from molbar API and historic data files
#   * GSP2 filtering and validation
#   * Run statistics processing
#   * Historic data integration
#   * QC metric calculations and thresholds
#
# Version 1.1.0 (07/18/2024)
# - Enhanced QC Features:
#   * Added failed samples report generation
#   * Implemented QC failure graphs
#   * Added normalized value visualization
#   * Enhanced count graph generation
#   * Added below-27 threshold check
#   * Added lower limit flag        
#   * Added SD flag        
#   * Added below-SD threshold check        
#   * Added comparison with last 5 runs
#   * Added failed any QC metric flag

#===============================================================================
# Notes and To-Do List
#===============================================================================
# Future Improvements:
# 1. Compare depths with previous run only
#    - Implement direct comparison with last run's depths
#    - Add new failure column based on this comparison
#
# 2. Code Optimization
#    - Remove unused/deprecated columns
#    - Streamline data processing pipeline
#
# 3. Documentation
#    - Add detailed comments for each QC metric
#    - Document column selection rationale
#===============================================================================
# Required Libraries
#===============================================================================
# Load required libraries for data manipulation and web requests
library(tidyverse)
library(RCurl)
library(jsonlite)

#===============================================================================
# Command Line Arguments and Configuration
#===============================================================================
# Get command line arguments, with fallback defaults
args <- commandArgs(trailingOnly = TRUE)
# folder: directory containing the analysis files (default: 5154)
folder <- args[1]
# write_all_graphs: control graph output (default: "no")
write_all_graphs <- args[2]
write_all_graphs <- "no"
# Quality control parameters
number_of_SDs_flag <- 1    # Number of standard deviations for outlier detection
lower_limit_flag <- 54     # Lower limit threshold for quality checks
#===============================================================================
# Path Configuration
#===============================================================================
# Set up file paths for analysis
folder<-5162
run_name <- folder                                    # Use folder as run identifier
base_path <- ""                 # Root path for all analyses
full_path <- file.path(base_path, as.character(folder))  # Complete path to analysis folder
print(full_path)
#===============================================================================
# Validate Arguments
#===============================================================================
# Check if write_all_graphs argument is valid
if (write_all_graphs != "yes" && write_all_graphs != "no") {
  stop("Argument 4 (write out graphs) needs to be equal to yes or no, stopping.", "\n")
}
#===============================================================================
# Load and Process Historic Data
#===============================================================================
# Load historic data from file
historic_data_long2 <- read.table('', 
                                sep = "\t",  
                                header = TRUE, 
                                stringsAsFactors = FALSE, 
                                fill = TRUE, 
                                quote = "")

# Find most recent historic data file
files <- list.files(path = "/Data1/",
                   pattern = "historic_data_by.*long.*",
                   full.names = TRUE)
most_recent_file <- files[which.max(file.info(files)$mtime)]
historic_data_long2 <- read.delim(most_recent_file, sep = "\t", header = TRUE)
print(paste("Reading file:", basename(most_recent_file)))

# Filter for last 5 runs
unique_runs <- unique(historic_data_long2$run)
last_5_unique_runs <- tail(unique_runs, 5)
filtered_data <- historic_data_long2[historic_data_long2$run %in% last_5_unique_runs, ]
print(filtered_data)

historic_data_long2 <- filtered_data
data <- historic_data_long2

#===============================================================================
# Calculate Statistics per GSP2
#===============================================================================
gsp2_stats <- aggregate(counts ~ GSP2, data = data,
                       FUN = function(x) c(mean = mean(x), sd = sd(x)))

#===============================================================================
# Create Statistics Dataframe
#===============================================================================
gsp2_stats_df <- data.frame(
  GSP2 = gsp2_stats$GSP2,
  Average_Count = gsp2_stats$counts[,1],
  Standard_Deviation = gsp2_stats$counts[,2]
)
gsp2_stats_df <- gsp2_stats_df[order(-gsp2_stats_df$Average_Count), ]

#===============================================================================
# Load Historic Statistics
#===============================================================================
historic_stats <- read.table('', 
                           sep = "\t",  
                           header = TRUE, 
                           stringsAsFactors = FALSE, 
                           fill = TRUE, 
                           quote = "")

#===============================================================================
# Process Historic Statistics
#===============================================================================
historic_stats <- historic_stats %>% 
  mutate(
    counts_low_limit = mean_counts_by_primer - (mean_counts_by_primer * mean_cv_counts_by_primer * number_of_SDs_flag),
    Mean_norm_counts_low_limit = mean_normalized_counts_by_primer - (mean_normalized_counts_by_primer * mean_normalized_cv_counts_by_primer * number_of_SDs_flag)
  )

#===============================================================================
# Process Input Files
#===============================================================================
file_names <- list.files(pattern = "*counts$", full.names = TRUE)
print(file_names)

#===============================================================================
# Load and Process Input Files
#===============================================================================
listOfFiles <- lapply(file_names, function(file_names) read.table(file_names, header = FALSE))
names(listOfFiles) <- file_names
split_names <- names(listOfFiles) %>% str_split_fixed(., "_R1", n = 2)
split_names <- split_names[, 1]
split_names <- basename(split_names)
listOfFiles <- listOfFiles

#===============================================================================
# Process Files
#===============================================================================
for (z in 1:length(listOfFiles)) {
  listOfFiles[[z]] <- listOfFiles[[z]] %>% 
    dplyr::rename("GSP2" = 1, "counts" = 2) %>%
    mutate(
      mean_norm_counts = counts/mean(counts),
      library = split_names[z],
      run = run_name
    )
}

#===============================================================================
# Combine and Filter Data
#===============================================================================
job_table <- bind_rows(listOfFiles)
job_table <- job_table %>%
  filter(
    !grepl("BRAF_chr7_140477815_24_+_A1_GSP2", GSP2),
    !grepl("BRAF_chr7_140481455_22_-_A1_GSP2", GSP2),
    !grepl("CCNB3_chrX_50028185_20_-_A1_GSP2", GSP2),
    !grepl("FGFR3_chr4_1808931_19_+_A1_GSP2", GSP2),
    !grepl("KRAS_chr12_25362806_34_+_A1_GSP2", GSP2),
    !grepl("KRAS_chr12_25378592_22_-_A1_GSP2", GSP2),
    !grepl("MKL2_chr16_14345731_31_-_A1_GSP2", GSP2),
    !grepl("NRG1_chr8_32474359_32_-_A1_GSP2", GSP2),
    !grepl("PPARG_chr3_12447414_21_-_A1_GSP2", GSP2),
    !grepl("PPARG_chr3_12447547_25_+_A1_GSP2", GSP2),
    !grepl("TFE3_chrX_48891259_23_-_A1_GSP2", GSP2),
    !grepl("TMPRSS2_chr21_42870059_26_-_A1_GSP2", GSP2)
  )

#===============================================================================
# Filter Out Additional GSP2s
#===============================================================================
filterout <- read.csv("")
job_table <- job_table[!(job_table$GSP2 %in% filterout$GSP2),]

write.csv(job_table, "control.csv")


#===============================================================================
# Filter Based on Run Statistics
#===============================================================================
run_stats <- readr::read_tsv(paste0(folder, "_RunStatsFinal.txt"))
samples_to_keep <- run_stats %>% 
  filter(MeanCoverage >= 10) %>% 
  pull(SampleName)
job_table <- job_table %>% 
  filter(library %in% samples_to_keep)

#===============================================================================
# Combine with Historic Data and Calculate Statistics
#===============================================================================
historic_data_long2_job_add <- rbind(historic_data_long2, job_table)
job_table_grouped <- job_table %>% 
  group_by(GSP2) %>% 
  summarise(
    run_mean_counts = mean(counts),
    run_mean_normalized_counts = mean(mean_norm_counts)
  )

#===============================================================================
# Merge Statistics and Calculate Thresholds
#===============================================================================
stats_with_run <- merge(historic_stats, job_table_grouped)
stats_with_run_avg <- merge(stats_with_run, gsp2_stats_df)
stats_with_run_avg$Mean_Minus_SD <- stats_with_run_avg$Average_Count - stats_with_run_avg$Standard_Deviation
stats_with_run_avg$Mean_25 <- stats_with_run_avg$Average_Count * 0.25 

#===============================================================================
# Flag QC Metrics
#===============================================================================
stats_with_run_avg <- stats_with_run_avg %>% 
  mutate(
    flag_counts_below_avg = case_when(
      run_mean_counts < Mean_25 ~ "fail",
      run_mean_counts >= Mean_25 ~ "pass",
      TRUE ~ "error"
    )
  ) %>%
  mutate(
    flag_counts_below_lower_threshold = case_when(
      run_mean_counts < lower_limit_flag ~ "fail",
      run_mean_counts >= lower_limit_flag ~ "pass",
      TRUE ~ "error"
    )
  ) %>%
  mutate(
    flag_counts_below_27 = case_when(
      run_mean_counts < 27 ~ "fail",
      run_mean_counts >= 27 ~ "pass",
      TRUE ~ "error"
    )
  ) %>%
  mutate(
    failed_any_QC_metric = case_when(
      flag_counts_below_avg == "fail" | 
      flag_counts_below_lower_threshold == "fail" |
      flag_counts_below_27 == "fail" ~ "fail",
      TRUE ~ "pass"
    )
  )


#===============================================================================
# Write Output Files
#===============================================================================
historic_out_name <- paste0("historic_data_by_run_long_", run_name, ".tsv")
write.table(historic_data_long2_job_add, 
            file = historic_out_name, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)

stats_out_name <- paste0(run_name, "_run_report.tsv")
stats_write <- file.path(full_path, stats_out_name)
write.table(stats_with_run_avg, 
            file = stats_write, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)
#write.csv(stats_with_run_avg, "run_report.csv")

#===============================================================================
# Generate Failed Samples Report
#===============================================================================
# Create a summary of all failed samples with their QC metrics
failed_samples <- stats_with_run_avg %>%
  filter(failed_any_QC_metric == "fail") %>%
  select(GSP2, run_mean_counts, Mean_25, lower_limit_flag, 
         flag_counts_below_avg, flag_counts_below_lower_threshold, 
         flag_counts_below_27, failed_any_QC_metric)

write.csv(failed_samples, 
          file = paste0(run_name, "_failed_samples.csv"),
          row.names = FALSE)

#===============================================================================
# Generate QC Failure Graphs
#===============================================================================
failees <- stats_with_run_avg %>% 
  filter(failed_any_QC_metric == "fail") %>% 
  select(GSP2)

new_directory <- paste0(run_name, "_run_level_QC_fail_graphs_2")
dir_new <- file.path(full_path, new_directory)
dir.create(dir_new)

num_of_failees <- nrow(failees)
if (num_of_failees == 0 && write_all_graphs == "no") {
  stop("No primers failed QC checks, stopping.", "\n")
}

if (write_all_graphs == "yes") {
  failees <- stats_with_run %>% select(GSP2)
}

#===============================================================================
# Generate Normalized Value Graphs
#===============================================================================
for(i in 1:nrow(failees)) {
  primer_name <- as.character(failees[i, "GSP2"])
  primer_mean <- as.numeric(stats_with_run_avg %>% 
                           filter(GSP2 == primer_name) %>% 
                           select(mean_normalized_counts_by_primer))
  primer_cv <- as.numeric(stats_with_run_avg %>% 
                         filter(GSP2 == primer_name) %>% 
                         select(mean_normalized_cv_counts_by_primer))
  primer_sd <- primer_mean * primer_cv
  
  primer_data_to_graph <- historic_data_long2_job_add %>% 
    filter(GSP2 == primer_name) %>%
    mutate(run = factor(run, levels = sort(unique(run), decreasing = TRUE)))
  
  pic <- paste0(primer_name, "_mean_normalized.png")
  plot_out_name <- file.path(dir_new, pic)
  calc_bin_width <- diff(range(primer_data_to_graph$mean_norm_counts))/50
  
  ggplot(primer_data_to_graph, aes(x = run, y = mean_norm_counts)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_dotplot(binaxis = "y", 
                stackdir = "center", 
                dotsize = 0.25, 
                binwidth = calc_bin_width, 
                color = "blue3") + 
    expand_limits(y = 0) + 
    geom_hline(yintercept = primer_mean, color = "deepskyblue") + 
    geom_hline(yintercept = primer_mean - primer_sd, color = "green") + 
    geom_hline(yintercept = primer_mean + primer_sd, color = "green") + 
    geom_hline(yintercept = primer_mean - (2 * primer_sd), color = "red") + 
    geom_hline(yintercept = primer_mean + (2 * primer_sd), color = "red") + 
    labs(x = "", 
         y = "primer mean normalized counts", 
         title = primer_name) + 
    theme(axis.text.x = element_text(color = "black", 
                                   size = 4, 
                                   angle = 65, 
                                   hjust = 1))
  ggsave(filename = plot_out_name, width = 6, height = 4)
}

#===============================================================================
# Generate Count Graphs
#===============================================================================
for(i in 1:nrow(failees)) {
  primer_name <- as.character(failees[i, "GSP2"])
  primer_mean <- as.numeric(stats_with_run_avg %>% 
                           filter(GSP2 == primer_name) %>% 
                           select(mean_counts_by_primer))
  primer_cv <- as.numeric(stats_with_run_avg %>% 
                         filter(GSP2 == primer_name) %>% 
                         select(mean_cv_counts_by_primer))
  primer_sd <- primer_mean * primer_cv
  
  primer_data_to_graph <- historic_data_long2_job_add %>% 
    filter(GSP2 == primer_name) %>%
    mutate(run = factor(run, levels = sort(unique(run), decreasing = TRUE)))
  
  plot_out_name <- paste0(dir_new, run_name, primer_name, "_counts.png")
  calc_bin_width <- diff(range(primer_data_to_graph$counts))/50
  
  ggplot(primer_data_to_graph, aes(x = run, y = counts)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_dotplot(binaxis = "y", 
                stackdir = "center", 
                dotsize = 0.25, 
                binwidth = calc_bin_width, 
                color = "blue3") + 
    expand_limits(y = 0) + 
    geom_hline(yintercept = primer_mean, color = "deepskyblue") + 
    geom_hline(yintercept = primer_mean - primer_sd, color = "green") + 
    geom_hline(yintercept = primer_mean + primer_sd, color = "green") + 
    geom_hline(yintercept = primer_mean - (2 * primer_sd), color = "red") + 
    geom_hline(yintercept = primer_mean + (2 * primer_sd), color = "red") + 
    labs(x = "", 
         y = "primer count", 
         title = primer_name) + 
    theme(axis.text.x = element_text(color = "black", 
                                   size = 4, 
                                   angle = 65, 
                                   hjust = 1))
  ggsave(filename = plot_out_name, width = 6, height = 4)
}
