#!/usr/bin/env Rscript

remotes::install_cran("writexl", upgrade = "never")
remotes::install_cran("argparser", upgrade = "never")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("plyranges", update = FALSE)
library(plyranges)

library(tidyverse)
library(argparser)

DEBUG <- F
if (DEBUG) {
  #setwd("/projectnb/wax-dk/max/G222_CHIPSEQ/G222_G156_G207/temp2/")
  setwd("/projectnb/wax-dk/max/G223_H3K27ac/work/b4/fe650147015c7651495cd00d324bfe")
}

ParseArguments <- function() {
  p <- arg_parser('Peaks processing')
  p <- add_argument(p,'--peak_prefix', default="narrow", help="MACS2 peak prefixes (broad or narrow)")
  p <- add_argument(p,'--sample_labels', default="sample_labels.csv", help="to extract addition information we can use tab separated Sample_Labels.txt")
  # p <- add_argument(p,'--output_name', default="MACS2_peaks_extradetailed_overlap_report", help="output name for MACS2 peaks overlap report")
  return(parse_args(p))
}

argv <- ParseArguments()
print(argv)

peak_prefix <- argv$peak_prefix
sample_labels_config <- argv$sample_labels
output_name <- str_glue("MACS2_{argv$peak_prefix}_peaks_extradetailed_overlap_report.xlsx")


sample_labels <- read_csv(sample_labels_config, col_names=T) %>% 
  select(peak_id = sample_id, desc = sample_description, group_id, group_description) %>% 
  mutate(ix = str_pad(row_number(), width = 2, pad = "0"),
         desc = str_glue("{ix}_{peak_id}\n{group_id}_{desc}")) %>%
  add_count(group_id,name = "n_samples") %>% 
  mutate(group=str_glue("Group_{group_id}_{group_description}\n(n={n_samples})")) %>% 
  select(-ix, -n_samples, -group_id, -group_description)

xls <- map_dfr(list.files(pattern = paste0(peak_prefix,".*\\.xls$")), function (xls_fname){
  # fname <- str_extract(xls_fname,"^G\\d+\\w+\\d+(?=_MACS)")
  tmp <- read_tsv(xls_fname, comment = "#", col_names = T) %>% 
    select(seqnames = chr, start, end, peak_id = name, length, pileup, fold_enrichment, minus_log10_qvalue)
}) 

xls_merged <- xls %>% 
  #select(seqnames,start,end,peak_id) %>% 
  plyranges::as_granges() %>% 
  GenomicRanges::reduce(., min.gapwidth = 0L) 

xls_combined <- join_overlap_left(xls_merged,xls %>% plyranges::as_granges()) %>% 
  as_tibble() %>% 
  distinct() %>% 
  mutate(peak_id = str_replace(peak_id,"(_narrow_|_broad_)MACS2_peak_\\d+","")) %>% 
  left_join(x =., y = sample_labels, join_by(peak_id)) %>% 
  select(-peak_id) %>% 
  dplyr::rename(peak_id = desc)

final <- xls_combined %>% 
  filter(fold_enrichment == max(fold_enrichment) & pileup == max(pileup) & length == max(length), 
         .by = c(seqnames, start, end, peak_id)) %>% 
  mutate(overlap = 1)
  
## column names sample_id -- values - is specific sample present in this region (0/1)
final_presence <- final %>% 
  select(seqnames,start,end,overlap, peak_id) %>% 
  add_count(seqnames,start,end, name = "number_of_overlaps") %>% 
  distinct() %>% 
  pivot_wider(names_from = peak_id, values_from = overlap, names_sort = T, values_fill = 0) %>% 
  rename_with(~paste0(.x, " (Overlap = 1)"), matches("^[[:digit:]]"))

## column names -- group names, values - how many samples from this specific group 
## present in region
final_group_presence <- final %>%
  select(seqnames,start,end, width, group) %>%
  add_count(seqnames,start,end,group, name = "samples_in_group_in_region") %>%
  distinct() %>%
  pivot_wider(names_from = group, values_from = samples_in_group_in_region, names_sort = T, values_fill = 0)


## column names -- group names, values - how many samples from this specific group 
## present in region
final_group_extra_presence <- final %>% 
  select(-strand, -overlap, -peak_id, -width) %>% 
  pivot_longer(c("length", "pileup", "fold_enrichment", "minus_log10_qvalue"), values_to = "param") %>% 
  mutate(avg = round(mean(param),2), .by = c(seqnames,start,end,group,name)) %>% 
  select(-param) %>% 
  distinct() %>% 
  pivot_wider(names_from = c(group, name), values_from = avg, names_sort = T, values_fill = NA)

#nm <- names(union_left) %>% keep(~str_detect(., "padj|log2FC"))
nm <- names(final_group_extra_presence) %>% keep(~str_detect(., "length|pileup|fold_enrichment|minus_log10_qvalue"))

# set correct order for columns: MANORM2, DIFFREPS and RIPPM
nm <- c(nm[grepl("qvalue",nm)] %>% sort,
        nm[grepl("pileup",nm)] %>% sort,
        nm[grepl("length",nm)] %>% sort,
        nm[grepl("enrichment",nm)] %>% sort)

final_group_extra_presence <- final_group_extra_presence %>% 
  relocate(all_of(nm), .after = last_col())

  
## column names -- sample_id.(length,pileup,fold_enrichment,minus_log10_qvalue) -- extended
## information about each sample  
final_extra_detailed <- final %>% 
  select(-width, -strand, -overlap, -group) %>% 
  pivot_longer(c("length", "pileup", "fold_enrichment", "minus_log10_qvalue"), values_to = "param") %>% 
  distinct() %>% 
  pivot_wider(names_from = c(peak_id, name), values_from = param, names_sort = T, values_fill = NA)

## sort columns in specific order (ordered_subcategories)
sort_sample_subcategories <- function(df_all_columns, marker_columns, ordered_subcategories) {
  ix <- str_detect(df_all_columns, marker_columns) %>% which()
  scols <- df_all_columns[ix]
  
  map_chr(ordered_subcategories,\(subcategory){
    scols[grepl(subcategory, scols)]
  }) 
}

nms <- map(sample_labels$desc, \(sample_description){
  sort_sample_subcategories(names(final_extra_detailed), sample_description, c("qvalue","pileup","length","enrichment"))  
}) %>% unlist()

final_extra_detailed <- final_extra_detailed %>% 
  relocate(all_of(nms), .after = last_col())

#### Final result
real_final <- left_join(final_group_presence, final_group_extra_presence, join_by(seqnames,start,end)) %>% 
  left_join(x = ., final_presence, join_by(seqnames,start,end)) %>%
  left_join(x = ., final_extra_detailed, join_by(seqnames,start,end)) %>% 
  arrange(desc(number_of_overlaps)) %>% 
  mutate(ucsc_coords = str_glue("{seqnames}:{start}-{end}")) %>% 
  relocate(ucsc_coords, .before = seqnames) %>% 
  relocate(number_of_overlaps, .after = width)

### excel nrows limit ~1e6, just trim table by filtering
### if the table size more then 1e6
# noverlaps <- 1
# while(nrow(real_final) > 1000000) {
#   print(noverlaps)
#   real_final <- real_final %>% 
#     filter(number_of_overlaps > noverlaps)
#   noverlaps <- noverlaps + 1
# }

real_final %>% 
  writexl::write_xlsx(path = output_name, col_names = T)
