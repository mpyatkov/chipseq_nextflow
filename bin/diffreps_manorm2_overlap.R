#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

library(readxl)
library(writexl)

remotes::install_cran("argparser", upgrade = "never")
library(argparser)

remotes::install_cran("plyranges", upgrade = never)
library(plyranges)

remotes::install_cran("openxlsx2", upgrade = never)
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
}

diffreps_df <- list.files(pattern = "DIFFREPS|RIPPM") %>%
  map_dfr(\(f){
    readxl::read_xlsx(path = f, skip = 4, sheet = 2) %>%
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
      #filter(padj < 0.05) %>%
      filter(!is.na(delta)) %>%  
      mutate(filename = f %>% tools::file_path_sans_ext())
  }) %>% 
  rename_with(~str_replace(.x, "(.*?)\\.control\\.mean", "intensity.Control")) %>% 
  rename_with(~str_replace(.x, "(.*?)\\.treatment\\.mean", "intensity.Treatment")) %>% 
  mutate(coords = as.character(str_glue("{seqnames}:{start}-{end}")))


merged_union <- bind_rows(diffreps_df, manorm2_df) %>% 
  as_granges() %>% 
  GenomicRanges::reduce(., min.gapwidth = 0L)

union_left <- join_overlap_left(merged_union, bind_rows(diffreps_df, manorm2_df) %>% as_granges()) %>% 
  as_tibble() %>% 
  filter(padj == min(padj), .by = c(seqnames,start,end,filename)) %>% 
  filter(abs(log2FC) == max(abs(log2FC)), .by = c(seqnames,start,end,filename)) %>% 
  filter(intensity.Treatment == max(intensity.Treatment) & intensity.Control == max(intensity.Control), .by = c(seqnames,start,end,filename)) %>% 
  distinct() %>% 
  add_count(seqnames,start,end, name = "n_any_quality_overlaps") %>% 
  ## delta starts with: 1_*,4_* - signif, 2_*,3_* - weak, 0_* - low read regions
  dplyr::mutate(n_signif_quality_overlaps = sum(str_detect(delta, "1_|4_")), .by = c(seqnames,start,end)) %>% 
  select(-delta) %>% 
  pivot_wider(names_from = filename, values_from = c(coords, intensity.Control, intensity.Treatment, padj, log2FC), values_fill = NA, names_glue = "{filename}.{.value}") 

#nm <- names(union_left) %>% keep(~str_detect(., "padj|log2FC"))
nm <- names(union_left) %>% keep(~str_detect(., "MANORM|DIFFREPS|RIPPM"))

# set correct order for columns: MANORM2, DIFFREPS and RIPPM
nm <- c(nm[grepl("MANORM2",nm)] %>% sort,
  nm[grepl("DIFFREPS",nm)] %>% sort,
  nm[grepl("RIPPM",nm)] %>% sort)

final <- union_left %>% 
  relocate(all_of(nm), .after = last_col()) %>% 
  mutate(padj_sort = 10^-rowMeans(-log10(select(., contains("padj"))), na.rm = TRUE),
         log2fc_sort = rowMeans(select(., contains("log2FC")), na.rm = TRUE)) %>% 
  relocate(n_signif_quality_overlaps, .after = n_any_quality_overlaps) %>% 
  arrange(desc(n_any_quality_overlaps), desc(n_signif_quality_overlaps),padj_sort, desc(abs(log2fc_sort))) %>% 
  mutate(across(contains("log2FC"), ~round(.x, digits = 2)))

## generate excel format (ex. C1:C1234, E1:E1234) ranges for "keyword" columns
## need to add scientific output format for specific columns (openxlsx2 package)
get_range <- function(df, keyword){
  col_ix <- which(str_detect(names(df),keyword))
  c(LETTERS, paste0("A",LETTERS),paste0("A",LETTERS))[col_ix] %>% 
    paste0(.,"2",":",.,nrow(final)+1)
}

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

