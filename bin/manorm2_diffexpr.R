#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

library(readxl)
library(writexl)

remotes::install_cran("MAnorm2", upgrade = "never")
library(MAnorm2)

remotes::install_cran("argparser", upgrade = "never")
library(argparser)

remotes::install_cran("plyranges", upgrade = never)
library(plyranges)


ParseArguments <- function() {
  p <- arg_parser('MAnorm2 processing')
  p <- add_argument(p,'--manorm2_profile', default="profile.xls", help="input file with profile obtained by 'profile_bins' function")
  p <- add_argument(p,'--control_samples', default="", help="list of samples separated by pipe operator. ex. 'G123_M1|G123_M2' which represent control condition")
  p <- add_argument(p,'--treatment_samples', default="", help="list of samples separated by pipe operator. ex. 'G123_M1|G123_M2' which represent treatment condition")
  p <- add_argument(p,'--control_name', default="", help="control group name")
  p <- add_argument(p,'--treatment_name', default="", help="treatment group name")
  p <- add_argument(p,'--output_name', default="", help="name for group of comparisons")
  return(parse_args(p))
}

argv <- ParseArguments()

print(argv)

DEBUG <- T
if (DEBUG) {
  argv$manorm2_profile <-  "/projectnb/wax-dk/max/G222_CHIPSEQ/G222_G156_G207/work/08/0525e292f7f92fc6205f2536a40504/Male_8wk_H3K27ac_vs_Male_3wk_H3K27ac_profile_bins.xls"
  argv$control_samples <- "G222_M16|G222_M17"
  argv$treatment_samples <- "G222_M09|G222_M10|G207_M01|G207_M02"
  argv$treatment_name <- "Male_8wk_H3K27ac"
  argv$control_name <- "Male_3wk_H3K27ac"
}

### manorm2 read and filtering input file
mx.manorm2 <- read_tsv(argv$manorm2_profile, col_names = T) %>% 
  mutate(rowsums = purrr::reduce(select(., contains("read_cnt")), `+`),
         occsums = purrr::reduce(select(., contains("occup")), `+`)) %>% 
  filter(rowsums > 20, occsums > 1) %>% 
  select(-contains("sums"))

#mx.manorm2 %>% select(all_of(contains("occupa"))) %>% summarise_all(list(sum))


## extract indexes for specific samples and marker (read_cnt,occupancy)
## "sample1|sample2", "read_cnt" -> indexes
find_indexes <- function(v, samples, marker) {
  ## find indexes for strs which contains "read_cnt"
  a <- str_detect(v, marker) %>% which(isTRUE(.))
  
  ## find indexes for strs which contains samples
  b <- str_detect(v, samples) %>% which(isTRUE(.))
  
  ## find overlap of this indexes
  intersect(a,b)
}

v <- names(mx.manorm2)
control_cnts <- find_indexes(v, argv$control_samples, "read_cnt")
control_ocuppancy <- find_indexes(v, argv$control_samples, "occupancy")
treatment_cnts <- find_indexes(v, argv$treatment_samples, "read_cnt")
treatment_ocuppancy <- find_indexes(v, argv$treatment_samples, "occupancy")


## normalization
stat5mn <- MAnorm2::normalize(mx.manorm2, count = treatment_cnts, occupancy = treatment_ocuppancy)
stat5mn <- MAnorm2::normalize(stat5mn, count = control_cnts, occupancy = control_ocuppancy)

condsmn <- list(Treatment = bioCond(stat5mn[treatment_cnts], stat5mn[treatment_ocuppancy], name = str_glue("{argv$treatment_name}.treatment")),
                Control = bioCond(stat5mn[control_cnts], stat5mn[control_ocuppancy], name = str_glue("{argv$control_name}.control")))

autosome <- !(mx.manorm2$chrom %in% c("chrX", "chrY"))
condsmn <- normBioCond(condsmn, common.peak.regions = autosome)
condsmn <- fitMeanVarCurve(condsmn, method = "parametric", occupy.only = TRUE , init.coef = c(0.1, 10))
resmn <- diffTest(condsmn[[2]],condsmn[[1]]) ## Treatment vs Control

resmn %>% as_tibble() %>% filter(padj < 0.05 & abs(Mval) > 1) %>% dplyr::summarize(n = dplyr::n())

mx.manorm2.diff <- mx.manorm2 %>% 
  mutate(rowid = as.character(row_number())) %>% 
  inner_join(.,resmn %>% tibble::rownames_to_column(var="rowid")) %>% 
  select(-rowid)

# mx.manorm2 %>% select(matches(argv$control_samples))

mx.manorm2.diff %>% write_xlsx(path = str_glue("{argv$output_name}"), col_names = T)
