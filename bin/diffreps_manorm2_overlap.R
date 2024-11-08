#!/usr/bin/env Rscript
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("plyranges", update = FALSE)
library(plyranges)

# remotes::install_cran("argparser", upgrade = "never")
library(argparser)

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readxl)
library(writexl)

# remotes::install_cran("openxlsx2", upgrade = "never")
library(openxlsx2)

ParseArguments <- function() {
  p <- arg_parser('MAnorm2 processing')
  p <- add_argument(p,'--output_prefix', default="manorm2_vs_diffreps_overlap", help="Overlap between MAnorm2 and DIFFREPS")
  return(parse_args(p))
}
argv <- ParseArguments()
print(argv)

DEBUG <- F
if (DEBUG) {
  #setwd("/projectnb2/wax-dk/max/G222_CHIPSEQ/G222_G156_G207/temp")
  #setwd("/projectnb/wax-dk/max/G223_H3K27ac/work/04/76d6113cdc09ab218d39c0a538f88b")
  
  setwd("/projectnb/wax-dk/max/G223_H3K27ac/work/87/5ef814302fb6b401f0f58fea790349")
  setwd("/projectnb/wax-dk/max/G223_H3K27ac/work/04/76d6113cdc09ab218d39c0a538f88b")
  setwd("/projectnb/wax-dk/max/G229_G207_k27me3_k9me3/work/e9/b4f0cb7a083f4e9bf2234e8fbec035")
  #setwd("/projectnb/wax-dk/max/G223_H3K27ac_COMBINED/test")
  
}

diffreps_df <- list.files(pattern = "DIFFREPS|RIPPM") %>%
  map_dfr(\(f){
    readxl::read_xlsx(path = f, skip = 5, sheet = 2) %>%
      select(seqnames, start, end, contains("avg"),log2FC, padj,peakcaller_overlap = `Overlapped with peak caller`, delta) %>%
      filter(peakcaller_overlap == 1 & !is.na(delta) & padj < 0.05) %>%
      mutate(filename = f %>% tools::file_path_sans_ext()) %>%
      select(-peakcaller_overlap)
  }) %>% 
  rename_with(~str_replace(.x, "Control\\.avg", "intensity.Control")) %>% 
  rename_with(~str_replace(.x, "Treatment\\.avg", "intensity.Treatment")) %>% 
  mutate(coords = as.character(str_glue("{seqnames}:{start}-{end}")))
  
manorm2_df <- list.files(pattern = "MANORM2") %>%
  map_dfr(\(f){
    readxl::read_xlsx(path = f, sheet = 1) %>%
      select(seqnames = chrom, start, end, contains("treatment.mean"), contains("control.mean"), log2FC = Mval, padj,delta) %>%
      rename_with(~str_replace(.x, "(.*?)\\.control\\.mean", "intensity.Control")) %>%
      rename_with(~str_replace(.x, "(.*?)\\.treatment\\.mean", "intensity.Treatment")) %>%
      #filter(padj < 0.05) %>%
      filter(!is.na(delta)) %>%  
      mutate(filename = f %>% tools::file_path_sans_ext()) %>% 
      mutate(across(contains("delta"), as.character)) %>% 
      mutate(coords = as.character(str_glue("{seqnames}:{start}-{end}")))
  }) 
#%>% 
  #mutate(across(contains("delta"), as.character)) %>% 
  # rename_with(~str_replace(.x, "(.*?)\\.control\\.mean", "intensity.Control")) %>% 
  # rename_with(~str_replace(.x, "(.*?)\\.treatment\\.mean", "intensity.Treatment")) %>% 
  #mutate(coords = as.character(str_glue("{seqnames}:{start}-{end}")))


## split on upregulatation/downregulation
diffreps_df_plus <- diffreps_df %>% 
  filter(log2FC >= 0) %>% 
  mutate(regulation = "positive")


diffreps_df_minus <- diffreps_df %>% 
  filter(log2FC < 0)%>% 
  mutate(regulation = "negative")

## if there are no any sites in manorm2
if (nrow(manorm2_df) > 0) {
  manorm2_df_plus <- manorm2_df %>% 
    filter(log2FC >=0) %>% 
    mutate(regulation = "positive")
  
  manorm2_df_minus <- manorm2_df %>% 
    filter(log2FC <0) %>% 
    mutate(regulation = "negative")
}


## if there are no any sites in manorm2
if (nrow(manorm2_df) > 0) {
  df_plus <- bind_rows(diffreps_df_plus, manorm2_df_plus) 
  df_minus <- bind_rows(diffreps_df_minus, manorm2_df_minus) 
} else {
  df_plus <- diffreps_df_plus
  df_minus <- diffreps_df_minus
}

union_plus <- df_plus %>% 
  as_granges() %>% 
  GenomicRanges::reduce(., min.gapwidth = 0L) %>% 
  join_overlap_left(., df_plus %>% as_granges()) %>% 
  as_tibble() %>% 
  filter(padj == min(padj), .by = c(seqnames,start,end,filename)) %>% 
  filter(abs(log2FC) == max(abs(log2FC)), .by = c(seqnames,start,end,filename)) %>% 
  filter(intensity.Treatment == max(intensity.Treatment) & intensity.Control == max(intensity.Control), .by = c(seqnames,start,end,filename)) %>% 
  distinct() %>% 
  add_count(seqnames,start,end, name = "n_any_quality_overlaps") %>% 
  ## delta starts with: 1_*,4_* - signif, 2_*,3_* - weak, 0_* - low read regions
  dplyr::mutate(n_signif_quality_overlaps = sum(str_detect(delta, "1_|4_")), .by = c(seqnames,start,end)) 
# %>% 
#   select(-delta) 

union_minus <- df_minus %>% 
  as_granges() %>% 
  GenomicRanges::reduce(., min.gapwidth = 0L) %>% 
  join_overlap_left(., df_minus %>% as_granges()) %>% 
  as_tibble() %>% 
  filter(padj == min(padj), .by = c(seqnames,start,end,filename)) %>% 
  filter(abs(log2FC) == max(abs(log2FC)), .by = c(seqnames,start,end,filename)) %>% 
  filter(intensity.Treatment == max(intensity.Treatment) & intensity.Control == max(intensity.Control), .by = c(seqnames,start,end,filename)) %>% 
  distinct() %>% 
  add_count(seqnames,start,end, name = "n_any_quality_overlaps") %>% 
  ## delta starts with: 1_*,4_* - signif, 2_*,3_* - weak, 0_* - low read regions
  dplyr::mutate(n_signif_quality_overlaps = sum(str_detect(delta, "1_|4_")), .by = c(seqnames,start,end)) 
# %>% 
#   select(-delta) 

union_left <- bind_rows(union_plus, union_minus) %>% 
  ## removing duplicated rows which contains only different coordinates
  distinct(seqnames,start,end,width,strand,regulation,n_any_quality_overlaps,
           n_signif_quality_overlaps, filename, .keep_all = T)

## wider version of detailed data.frame
union_left_detailed <- union_left %>% 
  select(-delta) %>% 
  pivot_wider(names_from = filename, values_from = c(coords, intensity.Control, intensity.Treatment, padj, log2FC), values_fill = NA, names_glue = "{filename}.{.value}") 

#nm <- names(union_left) %>% keep(~str_detect(., "padj|log2FC"))
## set correct order for columns: MANORM2, DIFFREPS and RIPPM
nm <- names(union_left_detailed) %>% keep(~str_detect(., "MANORM|DIFFREPS|RIPPM"))
nm <- c(nm[grepl("MANORM2",nm)] %>% sort,
        nm[grepl("DIFFREPS",nm)] %>% sort,
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
nmp <- c(nmp[grepl("MANORM2",nmp)] %>% sort,
         nmp[grepl("DIFFREPS",nmp)] %>% sort,
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
  mutate(across(contains("log2FC"), ~round(.x, digits = 2)))



## TOP log2FC and padj (padj -> log2fc for this padj) (ADD TO FINAL TABLE)
top_log2fc_padj <- union_left %>% 
  select(seqnames, start, end, log2FC, padj) %>% 
  mutate(best_padj = min(padj), .by = c(seqnames, start, end)) %>% 
  filter(padj == best_padj) %>% 
  mutate(log2FC = round(log2FC, digits = 2)) %>% 
  select(seqnames, start, end, `Top padj` = padj, `log2FC related to top padj` = log2FC)

final <- left_join(final, top_log2fc_padj, join_by(seqnames,start,end)) %>% 
  relocate(all_of(c('Top padj','log2FC related to top padj')), .after = n_signif_quality_overlaps)

## generate excel format (ex. C1:C1234, E1:E1234) ranges for "keyword" columns
## need to add scientific output format for specific columns (openxlsx2 package)
get_range <- function(df, keyword){
  col_ix <- which(str_detect(names(df),keyword))
  c(LETTERS, paste0("A",LETTERS),paste0("A",LETTERS))[col_ix] %>% 
    paste0(.,"2",":",.,nrow(final)+1)
}


## save overlap data frame to xlsx
export_final_full <- final %>% select(-padj_sort, -log2fc_sort, -strand)
wb <- wb_workbook()$add_worksheet("manorm2_diffreps_overlap", gridLines = TRUE)$add_data("manorm2_diffreps_overlap", export_final_full)
get_range(export_final_full, "padj") %>% 
  walk(\(r) wb$add_numfmt("manorm2_diffreps_overlap", dims = r, numfmt = "0.00E+00"))
wb$save(file = str_glue("{argv$output_prefix}_overlap_manorm2_vs_diffreps.xlsx"))
# writexl::write_xlsx(x = final %>% select(-padj_sort, -log2fc_sort, -strand), path = argv$output_name, col_names = T)

top25_up <- final %>% 
  filter(log2fc_sort > 0) %>% 
  dplyr::slice_head(n = 25) %>% 
  select(chrom = seqnames, start, end, n_any_quality_overlaps, n_signif_quality_overlaps, average_padj = padj_sort, average_log2FC = log2fc_sort) 

wb <- wb_workbook()$add_worksheet("top25_manorm_vs_diffreps", gridLines = TRUE)$add_data("top25_manorm_vs_diffreps", top25_up)
get_range(top25_up, "padj") %>% 
  walk(\(r) wb$add_numfmt("top25_manorm_vs_diffreps", dims = r, numfmt = "0.00E+00"))
wb$save(file = str_glue("{argv$output_prefix}_top25_UPREGULATED.xlsx"))

top25_down <- final %>% 
  filter(log2fc_sort < 0) %>% 
  dplyr::slice_head(n = 25) %>% 
  select(chrom = seqnames, start, end, n_any_quality_overlaps, n_signif_quality_overlaps, average_padj = padj_sort, average_log2FC = log2fc_sort) 

wb <- wb_workbook()$add_worksheet("top25_manorm_vs_diffreps", gridLines = TRUE)$add_data("top25_manorm_vs_diffreps", top25_down)
get_range(top25_down, "padj") %>% 
  walk(\(r) wb$add_numfmt("top25_manorm_vs_diffreps", dims = r, numfmt = "0.00E+00"))
wb$save(file = str_glue("{argv$output_prefix}_top25_DOWNREGULATED.xlsx"))

