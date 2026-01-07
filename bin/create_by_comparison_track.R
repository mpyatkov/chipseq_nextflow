#!/usr/bin/env Rscript

library(argparser)
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(readxl)
library(tidyr)

ParseArguments <- function() {
  p <- arg_parser('create individual tracks')
  p <- add_argument(p,'--report_path', default="report_name.xlsx", help="path to full report xlsx")
  p <- add_argument(p,'--track_name', default="UCSC_track", help="output name for track")
  return(parse_args(p))
}
argv <- ParseArguments()
print(argv)

DEBUG <- F
if (DEBUG) {
  setwd("/projectnb/wax-dk/max/REFACTORING/work/e5/8c8cbd97cc615441c412c3da5c6b8c")
  argv$output_prefix <- "all_groups_together"
}

df <- list.files(pattern = "*xlsx") %>%
  map_dfr(\(f){
    
    sheet_number <- ifelse(str_detect(f, "DESEQ2"),1,2)
    print(sheet_number)
    
    readxl::read_xlsx(path = f, sheet = sheet_number) %>%
      select(seqnames, start, end,log2FC, padj, delta) %>% 
      filter(!is.na(delta)) %>%
      mutate(filename = f %>% tools::file_path_sans_ext()) 
  }) 

##### ADD COLORS #####
histogram_colors <- tribble(
  ~delta, ~cols,
  "DOWN_STRONG", "red",
  "DOWN_WEAK", "pink",
  "UP_WEAK", "lightblue",
  "UP_STRONG","blue"
)

df.with_colors <- left_join(df, histogram_colors, join_by(delta))

##### CREATE BED TRACKS #####
df.with_colors %>% 
  drop_na(delta) %>% 
  select(seqnames, start, end, delta, cols) %>% 
  rowwise() %>% 
  mutate(t1 = 1000,
         t2 = ".",
         t3 = 0,
         t4 = 0,
         t5 = paste0(col2rgb(cols), collapse = ",")) %>%
  select(-cols) %>% 
  arrange(seqnames,start) %>% 
  write_tsv(file = argv$track_name, append = T, col_names = F)




