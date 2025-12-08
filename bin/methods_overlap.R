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
  setwd("/projectnb/wax-dk/max/REFACTORING/work/79/80d932ed64be7a1cb91f4a978d01fe/")
  argv$output_prefix <- "all_groups_together"
}

diffreps_df <- list.files(pattern = "*xlsx") %>%
  discard(~str_detect(.x,"all_groups")) %>% 
  map_dfr(\(f){
    
    sheet_number <- ifelse(str_detect(f, "DESEQ2"),1,2)
    print(sheet_number)
    
    readxl::read_xlsx(path = f, sheet = sheet_number) %>%
      select(seqnames, start, end,log2FC, padj, delta) %>% 
      filter(!is.na(delta)) %>%
      mutate(filename = f %>% tools::file_path_sans_ext(),
             log2FC = as.numeric(log2FC)) 
  }) 

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
  dplyr::mutate(n_signif_quality_overlaps = sum(str_detect(delta, "DOWN_STRONG|UP_STRONG")), .by = c(seqnames,start,end)) 

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
  dplyr::mutate(n_signif_quality_overlaps = sum(str_detect(delta, "DOWN_STRONG|UP_STRONG")), .by = c(seqnames,start,end)) 

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



## set correct order for columns: MANORM2, DIFFREPS and RIPPM
nm <- names(union_left_detailed) %>% keep(~str_detect(., "DIFFREPS|RIPPM|DESEQ2|CSAW"))
# Extract unique comparison prefixes (e.g., "1_", "2_")
comparisons <- unique(str_extract(nm, "^\\d+_+")) %>% sort()
# For each comparison, append all methods in the desired order
nm <- unlist(lapply(comparisons, function(comp) {
  c(
    nm[grepl(paste0("^", comp, ".*DESEQ2"), nm)] %>% sort,
    nm[grepl(paste0("^", comp, ".*CSAW"), nm)] %>% sort,
    nm[grepl(paste0("^", comp, ".*DIFFREPS_200"), nm)] %>% sort,
    nm[grepl(paste0("^", comp, ".*RIPPM_200"), nm)] %>% sort,
    nm[grepl(paste0("^", comp, ".*DIFFREPS_1000"), nm)] %>% sort,
    nm[grepl(paste0("^", comp, ".*RIPPM_1000"), nm)] %>% sort
  )
}))

# nm <- names(union_left_detailed) %>% keep(~str_detect(., "DIFFREPS|RIPPM|DESEQ2|CSAW"))
# nm <- c(nm[grepl("DESEQ2",nm)] %>% sort,
#         nm[grepl("CSAW",nm)] %>% sort,
#         nm[grepl("DIFFREPS_200",nm)] %>% sort,
#         nm[grepl("RIPPM_200",nm)] %>% sort,
#         nm[grepl("DIFFREPS_1000",nm)] %>% sort,
#         nm[grepl("RIPPM_1000",nm)] %>% sort)

union_left_detailed <- union_left_detailed %>% 
  relocate(all_of(nm), .after = last_col())


## shows 0/1 overlaps for specific region for all MANORM/DIFFREPS methods
union_left_presence <- union_left %>% 
  select(seqnames,start,end,filename,regulation, delta) %>% 
  #mutate(overlap = 1) %>% 
  distinct() %>% 
  pivot_wider(names_from = filename, values_from = delta, values_fill = NA)

## change order of columns
nmp <- names(union_left_presence) %>% keep(~str_detect(., "DIFFREPS|RIPPM|DESEQ2|CSAW"))
# Extract unique comparison prefixes (e.g., "1_", "2_")
comparisons <- unique(str_extract(nm, "^\\d+_+")) %>% sort()


# For each comparison, append all methods in the desired order
nmp <- unlist(lapply(comparisons, function(comp) {
  c(
    nmp[grepl(paste0("^", comp, ".*DESEQ2"), nmp)] %>% sort,
    nmp[grepl(paste0("^", comp, ".*CSAW"), nmp)] %>% sort,
    nmp[grepl(paste0("^", comp, ".*DIFFREPS_200"), nmp)] %>% sort,
    nmp[grepl(paste0("^", comp, ".*RIPPM_200"), nmp)] %>% sort,
    nmp[grepl(paste0("^", comp, ".*DIFFREPS_1000"), nmp)] %>% sort,
    nmp[grepl(paste0("^", comp, ".*RIPPM_1000"), nmp)] %>% sort
  )
}))

# 
# nmp <- c(nmp[grepl("DESEQ2",nmp)] %>% sort,
#          nmp[grepl("CSAW",nmp)] %>% sort,
#          nmp[grepl("DIFFREPS200",nmp)] %>% sort,
#          nmp[grepl("RIPPM200",nmp)] %>% sort,
#          nmp[grepl("DIFFREPS1000",nmp)] %>% sort,
#          nmp[grepl("RIPPM1000",nmp)] %>% sort)


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


output_list <- if (argv$output_prefix == "all_groups_together") {
  tmp <- list()
  tmp[["METHODS_OVERLAP_ALL"]] <- final %>% select(-padj_sort, -log2fc_sort, -strand)
  tmp[["TOP25_UP_ALL"]]  <- top25_up
  tmp[["TOP25_DOWN_ALL"]]  <- top25_down
  tmp[["TOP25_UP_noXY_ALL"]]  <- top25_up_noXY
  tmp[["TOP25_DOWN_noXY_ALL"]]  <- top25_down_noXY
  tmp
} else {
  tmp <- list()
  tmp[["METHODS_OVERLAP"]] <- final %>% select(-padj_sort, -log2fc_sort, -strand)
  tmp[["TOP25_UP"]]  <- top25_up
  tmp[["TOP25_DOWN"]]  <- top25_down
  tmp[["TOP25_UP_noXY"]]  <- top25_up_noXY
  tmp[["TOP25_DOWN_noXY"]]  <- top25_down_noXY
  tmp
}

writexl::write_xlsx(x = output_list, path = str_glue("{argv$output_prefix}_overlap.xlsx"))

writexl::write_xlsx(x=top25_up, path = str_glue("{argv$output_prefix}_TOP25_UPREGULATED.xlsx"))
writexl::write_xlsx(x=top25_down, path = str_glue("{argv$output_prefix}_TOP25_DOWNREGULATED.xlsx"))
writexl::write_xlsx(x=top25_up_noXY, path = str_glue("{argv$output_prefix}_TOP25_NOXY_UPREGULATED.xlsx"))
writexl::write_xlsx(x=top25_down_noXY, path = str_glue("{argv$output_prefix}_TOP25_NOXY_DOWNREGULATED.xlsx"))
