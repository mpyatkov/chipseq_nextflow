#!/usr/bin/env Rscript

library(dplyr)
library(readr)

args = commandArgs(trailingOnly=TRUE)

stats_path <- args[1]
output_name <-args[2]

norm_factors <- read_csv(stats_path, col_names = F) %>%
  select(sample_id = X1, fragment_count = X2, fragment_in_peak_count = X3, fragment_in_peak_ratio = X4) %>%
  mutate(norm_factor = round(min(fragment_in_peak_count)/fragment_in_peak_count,3)) %>% 
  write_tsv(output_name, col_names = T)
