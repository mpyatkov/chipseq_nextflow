#!/usr/bin/env Rscript

library(argparser)

ParseArguments <- function() {
  p <- arg_parser('Postprocessing of diffreps outputs')
  p <- add_argument(p,'--annotated_path', default="input.annotated", help="path to diffReps annotated output")
  p <- add_argument(p, '--min_avg_count', default = 20, help = "Threshold for average number of reads per peak for CONTROL and TREATMENT")
  p <- add_argument(p, '--log2fc_cutoff', default = 1, help = "Threshold for fold change (logarithmic scale)")
  p <- add_argument(p, '--mumerge_path', help = "Path to MuMerge file which will be used instead of MACS2 union")
  p <- add_argument(p,'--output_prefix', default="NONE", help="name for group of comparisons")
  

  return(parse_args(p))
}

argv <- ParseArguments()

print(argv)

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(plyranges)

## input params
annotated_path <- argv$annotated_path

min_avg_count <- argv$min_avg_count
log2fc_cutoff <- argv$log2fc_cutoff

DEBUG <- FALSE

if (DEBUG){
  main_path <- "/projectnb/wax-dk/max/REFACTORING/work/0b/b9a21c6263fbb458b00db47119fb62/"
  
  ## make sure that you run the following command in the main_path to create the links:
  ## cat .command.run | grep ln | grep -E "mumerge|annotated" | bash
  annotated_path <- fs::dir_ls(path = main_path, glob = "*.annotated") |> as.character()
  argv$mumerge_path <- fs::dir_ls(path = main_path, glob = "*mumerge*") |> as.character()
  
  min_avg_count <- 20
  log2fc_cutoff <- 1
}



######## load diffReps output annotated or "with header" files
diffreps.df <- read_tsv(annotated_path, col_names = T, comment = "#") %>% 
  dplyr::rename(seqnames = Chrom, start = Start, end = End) %>% 
  mutate(delta = case_when(Event == "Down" & abs(log2FC) > log2fc_cutoff & Control.avg > min_avg_count & padj < 0.05  ~ "DOWN_STRONG", 
                           Event == "Down" & abs(log2FC) <= log2fc_cutoff & Control.avg > min_avg_count & padj < 0.05 ~ "DOWN_WEAK", 
                           Event == "Up" & abs(log2FC) > log2fc_cutoff & Treatment.avg > min_avg_count & padj < 0.05 ~ "UP_STRONG", 
                           Event == "Up" & abs(log2FC) <= log2fc_cutoff & Treatment.avg > min_avg_count & padj < 0.05 ~ "UP_WEAK",
                           #abs(log2FC) > log2fc_cutoff ~ "0_Low_reads_sites", 
                           .default = NA)) 

mumerge <- read_tsv(argv$mumerge_path, col_names = F) %>%
  select(seqnames = X1, start = X2, end = X3) %>%
  plyranges::as_granges() 


## only Significant and Weakest sites
diffreps_mumerge <- plyranges::join_overlap_left(mumerge, diffreps.df %>% plyranges::as_granges()) %>% 
  as_tibble() %>%
  drop_na(delta) %>% 
  as_tibble() %>% 
  filter(padj == min(padj), .by = c(seqnames,start,end)) %>% 
  filter(abs(log2FC) == max(abs(log2FC)), .by = c(seqnames,start,end)) %>% 
  filter(Treatment.avg == max(Treatment.avg) & Control.avg == max(Control.avg), .by = c(seqnames,start,end)) %>% 
  distinct()

writexl::write_xlsx(list(
  DIFFREPS = diffreps.df %>% arrange(padj),
  DIFFREPS_mumerge_centric = diffreps_mumerge), ## all quality 
  path = paste0(argv$output_prefix,".xlsx"), col_names = T)


