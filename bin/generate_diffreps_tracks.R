#!/usr/bin/env Rscript

library(argparser)
ParseArguments <- function() {
  p <- arg_parser('Postprocessing of diffreps outputs')
  p <- add_argument(p, '--diffreps_config', default = "Sample_Labels.csv", help = "CSV file with information about samples")
  p <- add_argument(p, '--diffreps_tracks', default = "", help = "CSV without header with sample_id and file names information")
  p <- add_argument(p, '--output_name', default = "diffreps_tracks.txt", help = "Output name with sample specific tracks")
  p <- add_argument(p, '--data_path', default = "buuser/TEST1", help = "Directory on the waxmanserver where your data will be available")
  
  return(parse_args(p))
}

argv <- ParseArguments()

DEBUG <- FALSE

if (DEBUG) {
  argv$diffreps_config <- "/projectnb/wax-dk/max/features/chipseq_nextflow_exp/work/60/4d3dac3c02a25ee319b53472ebada9/diffreps_config.csv"
  argv$diffreps_tracks <- "/projectnb/wax-dk/max/features/chipseq_nextflow_exp/work/tmp/f2/fd65859ec01cdf16eaba10d5ba80e6/collect-file.data"
  argv$data_path <- "buuser/TEST1"
}

library(dplyr)
library(stringr)
library(readr)
library(purrr)


diffreps_config <- read_csv(argv$diffreps_config, col_names = F) %>% 
  select(exp_num = 1, tr_name = 2, ctrl_name = 3, tr_samples = 4, ctrl_samples = 5, normalization = 6, window_size = 7)

diffreps_tracks <- read_csv(argv$diffreps_tracks, col_names = F) %>% 
  select(exp_num = 1, filename = 2)

combined <- left_join(diffreps_config, diffreps_tracks)

get_track_line <- function(t){
  track_path <- str_glue("http://waxmanlabvm.bu.edu/{argv$data_path}/{t$filename}")
  peakcaller <- ifelse(str_detect(t$filename, "MACS2"), "MACS2", "SICER")
  track <- str_glue("track type=bigBed name={t$exp_num}_{t$tr_name}_vs_{t$ctrl_name} ",
           "description={t$tr_name}_vs_{t$ctrl_name}_{peakcaller}_{t$normalization}_{t$window_size} visibility=dense itemRgb=on ",
           "bigDataUrl={track_path}")
  
  tibble(track = track)
}

diffreps_track_lines <- map_dfr(combined %>% mutate(r = row_number()) %>% group_by(r) %>% group_split(), get_track_line)
write_csv(diffreps_track_lines, argv$output_name, col_names = F)

