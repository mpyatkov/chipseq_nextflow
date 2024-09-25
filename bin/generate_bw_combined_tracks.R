#!/usr/bin/env Rscript

library(argparser)
ParseArguments <- function() {
  p <- arg_parser('Postprocessing of diffreps outputs')
  p <- add_argument(p, '--combined_tracks', default = "", help = "CSV without header with group_name, samples string, color separated by '_' and file names information")
  p <- add_argument(p, '--output_name', default = "bw_combined_tracks.txt", help = "Output name with sample specific tracks")
  p <- add_argument(p, '--data_path', default = "buuser/TEST1", help = "Directory on the waxmanserver where your data will be available")
  return(parse_args(p))
}

argv <- ParseArguments()
print(argv)

library(dplyr)
library(stringr)
library(readr)
library(purrr)


DEBUG <- F
if(DEBUG){
  argv$combined_tracks<- "/projectnb/wax-dk/max/G157_H3K27me3_20240916/work/tmp/03/c94af41d69b79087e239d44d6d9750/collect-file.data"
  argv$data_path <- "buuser/TEST1"
}


combined_tracks <- read_csv(argv$combined_tracks, col_names = F) %>% 
  select(group_name = 1, samples_str = 2, color = 3, filename = 4) %>%
  mutate(samples_str = str_replace_all(samples_str,'__',', '), 
         color  = color %>% str_replace_all('\\"',"") %>% str_replace_all("__",","),
         description = str_glue("Combined {group_name} ({samples_str})")) %>% 
  arrange(group_name) 


get_track_line <- function(t){
  track_path <- str_glue("http://waxmanlabvm.bu.edu/{argv$data_path}/{t$filename}")
  
  track <- str_glue("track type=bigWig name={t$group_name}_COMB description='{t$description}' ",
                    "db=mm9 visibility=full color={t$color} ",
                    "autoScale=on viewLimits=0.0:100.0 yLineOnOff=off windowingFunction=mean smoothingWindow=3 maxHeightPixels=100:64:8 ",
                    "bigDataUrl={track_path}")
  
  tibble(track = track)
}

combined_track_lines <- map_dfr(combined_tracks %>% mutate(r = row_number()) %>% group_by(r) %>% group_split(), get_track_line) 
write_csv(combined_track_lines, argv$output_name, col_names = F, quote = "none", escape="none")


