#!/usr/bin/env Rscript

library(argparser)
ParseArguments <- function() {
  p <- arg_parser('Postprocessing of diffreps outputs')
  p <- add_argument(p, '--sample_labels', default = "Sample_Labels.csv", help = "CSV file with information about samples")
  p <- add_argument(p, '--sid_tracks', default = "", help = "CSV without header with sample_id and file names information")
  p <- add_argument(p, '--output_name', default = "sid_tracks.txt", help = "Output name with sample specific tracks")
  p <- add_argument(p, '--data_path', default = "buuser/TEST1", help = "Directory on the waxmanserver where your data will be available")

  return(parse_args(p))
}

argv <- ParseArguments()
print(argv)

library(dplyr)
library(stringr)
library(readr)
library(purrr)

DEBUG <- FALSE

## order of sample specific tracks in UCSC browser for each sample_id
tracks_order <- tibble(order = c(1,3,2,5,4),
                       track_type = c("bw","broad","narrow","bam","epic2"),
                       track_suffix = c("_RiPPM_norm", "_MACS2_narrow","_MACS2_broad","", "_SICER_broad"),
                       track_ucsc_type = c("bigWig", "bigBed", "bigBed","bam", "bigBed"),
                       visibility = c("full","squish","squish","hide", "hide"),
                       autoscale = c("on", "-","-","-","-"),
                       should_update_color = c(NA, "255,0,0", "0,0,0", "0,0,0","0,0,0"))


if(DEBUG){
  argv$sample_labels <- "/projectnb/wax-dk/max/features/chipseq_nextflow_exp/work/60/4d3dac3c02a25ee319b53472ebada9/sample_labels.csv"
  argv$sid_tracks <- "/projectnb/wax-dk/max/features/chipseq_nextflow_exp/work/tmp/b6/a3ce2bbbf4dc2d8c1587fc9dcabefe/collect-file.data"
  argv$data_path <- "buuser/TEST1"
}


#### MAIN
sid_tracks <- read_csv(argv$sid_tracks,col_names = F) %>% 
  select(sample_id = 1, filename = 2)

sample_labels <- read_csv(argv$sample_labels, col_names = T, col_types= as.col_spec("ccccc")) ## because last column is comma separated color

combined_df <- left_join(sid_tracks, sample_labels, by="sample_id") %>% 
  mutate(track_type = case_when(str_detect(filename, "narrow") ~ "narrow",
                                str_detect(filename, "broad") ~ "broad",
                                str_detect(filename, "bw")~ "bw",
                                str_detect(filename, "epic")~ "epic2",
                                TRUE ~ "bam")) %>% 
  left_join(., tracks_order, by = "track_type") %>% 
  mutate(color = coalesce(should_update_color, color)) %>% 
  arrange(group_id,sample_id, order) %>% 
  select(-should_update_color)

get_track_line <- function(t){
  track_path <- str_glue("http://waxmanlabvm.bu.edu/{argv$data_path}/{t$filename}")

  track <- if (t$track_type == "bw") {
    extra <- str_glue("autoScale={t$autoscale} viewLimits=0.0:100.0 yLineOnOff=off windowingFunction=mean smoothingWindow=3 maxHeightPixels=100:64:8")
    str_glue("track type={t$track_ucsc_type} name={t$sample_id}{t$track_suffix} description={t$sample_description} ",
             "db=mm9 visibility={t$visibility} color={t$color} ",
             "{extra} ",
             "bigDataUrl={track_path}")
  } else {
    str_glue("track type={t$track_ucsc_type} name={t$sample_id}{t$track_suffix} description={t$sample_description} ",
             "db=mm9 visibility={t$visibility} color={t$color} ",
             "bigDataUrl={track_path}")
  }
  
  tibble(track = track)
}

sid_track_lines <- map_dfr(combined_df %>% mutate(r = row_number()) %>% group_by(r) %>% group_split(), get_track_line) 
write_csv(sid_track_lines, argv$output_name, col_names = F, quote = "none", escape="none")
