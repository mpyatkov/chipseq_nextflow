#!/usr/bin/env Rscript

library(argparser)
ParseArguments <- function() {
  p <- arg_parser('Postprocessing of diffreps outputs')
  p <- add_argument(p, '--diffreps_tracks', default = "", help = "CSV without header with sample_id and file names information")
  p <- add_argument(p, '--output_name', default = "diffreps_tracks.txt", help = "Output name with sample specific tracks")
  p <- add_argument(p, '--data_path', default = "buuser/TEST1", help = "Directory on the waxmanserver where your data will be available")
  
  return(parse_args(p))
}

argv <- ParseArguments()

DEBUG <- FALSE

if (DEBUG) {
  #argv$diffreps_config <- "/projectnb/wax-dk/max/REFACTORING/raw_configs/diffreps_config.csv"
  argv$diffreps_config <- "/projectnb/wax-dk/max/REFACTORING/work/tmp/17/f56bbd1f9278fede310d4acb7c59f3/diffreps_output.csv"
  argv$diffreps_tracks <- "/projectnb/wax-dk/max/REFACTORING/work/tmp/f4/2983095d2b562cc4fbe64374668142/collect-file.data"
  argv$data_path <- "buuser/TEST1"
}

library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(tidyr)

diffreps_tracks <- read_csv(argv$diffreps_tracks, col_names = F) %>% 
    select(group_name = 1, description = 2, filename = 3) %>%
    arrange(description)

## combined <- left_join(diffreps_config, diffreps_tracks) %>%
##   drop_na()

# uncomment code below if the order in the configuration file is random, this allows the
# tracks to be sorted so that the MANORM track for each group is followed by the
# DIFFREPS/RIPPM tracks for the same group.
# %>% 
  # mutate(gr = case_when(normalization == "MANORM" ~ str_glue("{ctrl_name}{tr_name}0{normalization}{window_size}"),
  #                       TRUE ~ str_glue("{ctrl_name}{tr_name}{normalization}{window_size}"))) %>%
  # arrange(gr) %>%
  # select(-gr)

get_track_line <- function(t){
  track_path <- str_glue("http://waxmanlabvm.bu.edu/{argv$data_path}/{t$filename}")
  track <- str_glue("track type=bigBed name='{t$description}' ",
           "description='{t$description}' visibility=dense itemRgb=on ",
           "bigDataUrl={track_path}")
  
  tibble(track = track)
}

diffreps_track_lines <- map_dfr(diffreps_tracks %>%
                                mutate(r = row_number()) %>%
                                group_by(r) %>%
                                group_split(), get_track_line)

write_csv(diffreps_track_lines, argv$output_name, col_names = F)

