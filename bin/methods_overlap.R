#!/usr/bin/env Rscript


library(plyranges)
library(argparser)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readxl)
library(writexl)

ParseArguments <- function() {
  p <- arg_parser('individual overlap processing')
  p <- add_argument(p,'--output_prefix', default="manorm2_vs_diffreps_overlap", help="Overlap between MAnorm2 and DIFFREPS")
  return(parse_args(p))
}
argv <- ParseArguments()
print(argv)

DEBUG <- F
if (DEBUG) {
  setwd("/projectnb/wax-dk/max/REFACTORING/work/85/a236663a8b35a1e0b4316851862d94")
  argv$output_prefix <- "all_groups_together"
}


diffreps_df <- list.files(pattern = "DIFFREPS|RIPPM") %>%
  map_dfr(\(f){
    readxl::read_xlsx(path = f, sheet = 2) %>%
      select(seqnames, start, end,log2FC, padj, delta) %>% 
      filter(!is.na(delta) & padj < 0.05) %>%
      mutate(filename = f %>% tools::file_path_sans_ext()) 
  }) 
# %>%  
#   mutate(coords = as.character(str_glue("{seqnames}:{start}-{end}")))
  
## split on upregulatation/downregulation
df_plus <- diffreps_df %>% 
  filter(log2FC >= 0) %>% 
  mutate(regulation = "positive")


df_minus <- diffreps_df %>% 
  filter(log2FC < 0)%>% 
  mutate(regulation = "negative")


## union plus
union_plus <- if (argv$output_prefix == "all_groups_together") {
  print("we are here")
  df_plus %>% 
    as_granges() %>% 
    GenomicRanges::reduce(., min.gapwidth = 0L) %>% 
    join_overlap_left(., df_plus %>% as_granges()) %>% 
    as_tibble() 
} else {
  df_plus %>% mutate(width = end - start, strand = "*")
}
  
union_plus <- union_plus %>% 
  filter(padj == min(padj), .by = c(seqnames,start,end,filename)) %>% 
  filter(abs(log2FC) == max(abs(log2FC)), .by = c(seqnames,start,end,filename)) %>% 
  distinct() %>% 
  add_count(seqnames,start,end, name = "n_any_quality_overlaps") %>% 
  ## delta starts with: 1_*,4_* - signif, 2_*,3_* - weak
  dplyr::mutate(n_signif_quality_overlaps = sum(str_detect(delta, "1_|4_")), .by = c(seqnames,start,end)) 

## union minus
union_minus <- if (argv$output_prefix == "all_groups_together") {
  df_minus %>% 
    as_granges() %>% 
    GenomicRanges::reduce(., min.gapwidth = 0L) %>% 
    join_overlap_left(., df_minus %>% as_granges()) %>% 
    as_tibble() 
} else {
  df_minus %>% mutate(width = end - start, strand = "*")
}

union_minus <- union_minus %>% 
  filter(padj == min(padj), .by = c(seqnames,start,end,filename)) %>% 
  filter(abs(log2FC) == max(abs(log2FC)), .by = c(seqnames,start,end,filename)) %>% 
  distinct() %>% 
  add_count(seqnames,start,end, name = "n_any_quality_overlaps") %>% 
  ## delta starts with: 1_*,4_* - signif, 2_*,3_* - weak, 0_* - low read regions
  dplyr::mutate(n_signif_quality_overlaps = sum(str_detect(delta, "1_|4_")), .by = c(seqnames,start,end)) 

## all together

union_left <- bind_rows(union_plus, union_minus) %>% 
  group_by(seqnames, start, end) %>%
  filter(n_distinct(regulation) == 1) %>% ## drop mixed regulation some positive, some negative
  ungroup() %>% 
  ## removing duplicated rows which contains only different coordinates
  distinct(seqnames,start,end,width,strand,regulation,n_any_quality_overlaps,
           n_signif_quality_overlaps, filename, .keep_all = T) 

## wider version of detailed data.frame
union_left_detailed <- union_left %>% 
  select(-delta) %>% 
  pivot_wider(names_from = filename, values_from = c(padj, log2FC), values_fill = NA, names_glue = "{filename}.{.value}") 

#nm <- names(union_left) %>% keep(~str_detect(., "padj|log2FC"))
## set correct order for columns: MANORM2, DIFFREPS and RIPPM
nm <- names(union_left_detailed) %>% keep(~str_detect(., "MANORM|DIFFREPS|RIPPM"))
nm <- c(nm[grepl("DIFFREPS",nm)] %>% sort,
        nm[grepl("RIPPM",nm)] %>% sort)

union_left_detailed <- union_left_detailed %>% 
  relocate(all_of(nm), .after = last_col())


## shows 0/1 overlaps for specific region for all MANORM/DIFFREPS methods
union_left_presence <- union_left %>% 
  select(seqnames,start,end,filename,regulation, delta) %>% 
  #mutate(overlap = 1) %>% 
  distinct() %>% 
  pivot_wider(names_from = filename, values_from = delta, values_fill = NA)

## change order of columns
nmp <- names(union_left_presence) %>% keep(~str_detect(., "MANORM|DIFFREPS|RIPPM"))
nmp <- c(nmp[grepl("DIFFREPS",nmp)] %>% sort,
         nmp[grepl("RIPPM",nmp)] %>% sort)
union_left_presence <- union_left_presence %>% 
  relocate(all_of(nmp),.after = last_col())

## join presence(with 0/1 integers) and detailed data.frames
union_final <- left_join(union_left_detailed, union_left_presence, join_by("seqnames","start","end","regulation")) %>% 
  relocate(all_of(nmp), .after = n_signif_quality_overlaps)


## add sorting columns to final table
final <- union_final %>% 
  mutate(padj_sort = 10^-rowMeans(-log10(select(., contains("padj"))), na.rm = TRUE),
         log2fc_sort = rowMeans(select(., contains("log2FC")), na.rm = TRUE)) %>% 
  relocate(n_signif_quality_overlaps, .after = n_any_quality_overlaps) %>% 
  arrange(desc(n_any_quality_overlaps), desc(n_signif_quality_overlaps),padj_sort, desc(abs(log2fc_sort))) %>% 
  mutate(across(contains("log2FC"), ~round(.x, digits = 2))) %>% 
  distinct()

## TOP log2FC and padj (padj -> log2fc for this padj) (ADD TO FINAL TABLE)
top_log2fc_padj <- union_left %>% 
  select(seqnames, start, end, log2FC, padj) %>% 
  mutate(best_padj = min(padj), .by = c(seqnames, start, end)) %>% 
  filter(padj == best_padj) %>% 
  mutate(log2FC = round(log2FC, digits = 2)) %>% 
  select(seqnames, start, end, `Top padj` = padj, `log2FC related to top padj` = log2FC) %>% 
  distinct()

final <- left_join(final, top_log2fc_padj, join_by(seqnames,start,end)) %>% 
  relocate(all_of(c('Top padj','log2FC related to top padj')), .after = n_signif_quality_overlaps)

top25_up <- final %>% 
  filter(log2fc_sort > 0) %>% 
  dplyr::slice_head(n = 25) %>% 
  select(chrom = seqnames, start, end, n_any_quality_overlaps, n_signif_quality_overlaps, average_padj = padj_sort, average_log2FC = log2fc_sort) 

top25_up_noXY <- final %>% 
  filter(log2fc_sort > 0 & !str_detect(seqnames,"X|Y")) %>% 
  dplyr::slice_head(n = 25) %>% 
  select(chrom = seqnames, start, end, n_any_quality_overlaps, n_signif_quality_overlaps, average_padj = padj_sort, average_log2FC = log2fc_sort) 

top25_down <- final %>% 
  filter(log2fc_sort < 0) %>% 
  dplyr::slice_head(n = 25) %>% 
  select(chrom = seqnames, start, end, n_any_quality_overlaps, n_signif_quality_overlaps, average_padj = padj_sort, average_log2FC = log2fc_sort) 

top25_down_noXY <- final %>% 
  filter(log2fc_sort < 0 & !str_detect(seqnames,"X|Y")) %>% 
  dplyr::slice_head(n = 25) %>% 
  select(chrom = seqnames, start, end, n_any_quality_overlaps, n_signif_quality_overlaps, average_padj = padj_sort, average_log2FC = log2fc_sort) 

writexl::write_xlsx(x = list(
  METHODS_OVERLAP = final %>% select(-padj_sort, -log2fc_sort, -strand),
  TOP25_UP = top25_up,
  TOP25_DOWN = top25_down,
  TOP25_UP_noXY = top25_up_noXY,
  TOP25_DOWN_noXY = top25_down_noXY
  ),
  path = str_glue("{argv$output_prefix}_overlap.xlsx"))


