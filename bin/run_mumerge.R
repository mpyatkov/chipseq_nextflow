#!/usr/bin/env Rscript

# Load necessary libraries
# When this script is run, the mumerge package itself and its dependencies
# should already be installed.
suppressPackageStartupMessages(library(mumerge))
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(tidyverse)) # Ok to use here for the script
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(future))

# --- Argument Parsing ---
p <- arg_parser('Postprocessing of ChIP-seq peaks with mumerge')
p <- add_argument(p, '--files_path', help="Path to MACS2 bed/xls files.")
p <- add_argument(p, '--pattern', default = "narrow", help="File name pattern")
p <- add_argument(p, '--output_file', default="mumerge_output.bed", help="Output BED file.")
p <- add_argument(p, '--ncores', default=1, help="Number of cores for processing.")
p <- add_argument(p, '--chunk_size', default=500, help="Number of regions per chunk for parallel processing.")
argv <- parse_args(p)

# --- Main Logic ---
final_results <- mumerge::run_peak_merging(
  files_path = argv$files_path,
  pattern = argv$pattern,
  ncores = as.integer(argv$ncores),
  chunk_size = as.integer(argv$chunk_size)
)

# --- Write Output ---
vroom::vroom_write(final_results, argv$output_file, col_names = FALSE)

cat("Mumerge processing complete. Results saved to:", argv$output_file, "\n")
