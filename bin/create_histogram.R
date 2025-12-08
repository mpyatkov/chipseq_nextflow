#!/usr/bin/env Rscript

library(ggplot2)
library(readxl)
library(stringr)
library(dplyr)
library(patchwork)
library(argparser)

ParseArguments <- function() {
  p <- arg_parser('Postprocessing of diffreps outputs')
  p <- add_argument(p,'--report_path', default="report_name.xlsx", help="path to full report xlsx")
  p <- add_argument(p, '--control_name', default = "CONTROL", help = "Control condition name")
  p <- add_argument(p, '--treatment_name', default = "TREATMENT", help = "Treatment condition name")
  p <- add_argument(p, '--method', default = "DIFFREPS", help = "DIFFREPS/RIPPM/DESEQ2/CSAW")
  p <- add_argument(p, '--window_size', default = "333", help = "DIFFREPS or RIPPM window size")
  p <- add_argument(p, '--treatment_samples', default = "", help = "Treatment samples (ex. G215M1,G215M2)")
  p <- add_argument(p, '--control_samples', default = "", help = "Control samples (ex. G215M3,G215M4)")
  p <- add_argument(p, '--compar_number', default = "00", help = "Comparison number, need only for output report name")
  p <- add_argument(p, '--output_name', default = "OUTPUT.pdf", help = "report name")
  return(parse_args(p))
}
argv <- ParseArguments()

print(argv)

TREATMENT_NAME <- argv$treatment_name
CONTROL_NAME <- argv$control_name
TREATMENT_SAMPLES <- argv$treatment_samples %>% str_replace_all("\\|",",")
CONTROL_SAMPLES <- argv$control_samples %>% str_replace_all("\\|",",")
METHOD <- argv$method
WINDOW_SIZE <- argv$window_size
REPORT_NAME <- argv$report_name
DF.PATH <- argv$report_path
EXPNUM <- argv$compar_number

DEBUG <- FALSE
if (DEBUG) {
  EXPNUM <- "2"
  TREATMENT_NAME <- "Male_8wk_ATAC"
  CONTROL_NAME <- "Female_8wk_ATAC"
  TREATMENT_SAMPLES <- "G241_M22|G241_M23|G228_M05|G228_M03" %>% str_replace_all("\\|",",")
  CONTROL_SAMPLES <- "G241_M45|G228_M04|G228_M06|G241_M46" %>% str_replace_all("\\|",",")
  METHOD <- "DIFFREPS"
  WINDOW_SIZE <- "1000" 
  REPORT_NAME <- "2_DIFFREPS_1000_Male_8wk_ATAC_vs_Female_8wk_ATAC" 
  DF.PATH <- "/projectnb/wax-dk/max/REFACTORING/work/c7/959e4a3d129d60bb594959ab658c87/Male_8wk_ATAC_vs_Female_8wk_ATAC_MuMerge_DIFFREPS_1000.xlsx"
}

sheet_number <- ifelse(str_detect(DF.PATH, "DEseq2"),1,2)
df <- readxl::read_xlsx(path = DF.PATH, sheet = sheet_number) %>%
  select(seqnames, start, end,log2FC, padj, delta) %>% 
  filter(!is.na(delta)) %>%
  mutate(filename = DF.PATH %>% tools::file_path_sans_ext()) 


histogram_colors <- c("red","pink", "lightblue","blue")
add_colors <- function(df_tmp, hist_colors) {
  
  df <- df_tmp %>% as_tibble()
  delta.ord <- df %>% 
    select(delta) %>% 
    distinct() %>% 
    arrange(delta) %>% 
    mutate(ord = as.numeric(str_extract(delta, "\\d+"))) 
  
  hist.ord <- tibble(cols = hist_colors) %>% 
    mutate(ord = row_number()) 
  
  delta.col <- left_join(x = delta.ord, y = hist.ord, join_by(ord)) %>% 
    distinct() %>% 
    arrange(ord) %>% 
    select(-ord)
  
  left_join(df, delta.col, join_by(delta)) 
  
}

df.with_colors <- add_colors(df, histogram_colors)


df.xy <- df.with_colors %>% 
  add_count(delta) %>% 
  mutate(condition_name = case_when(str_detect(delta,"1_") ~ str_glue("1_{CONTROL_NAME} signif. sites"),
                                    str_detect(delta,"2_") ~ str_glue("2_{CONTROL_NAME} weakest sites"),
                                    str_detect(delta,"3_") ~ str_glue("3_{TREATMENT_NAME} weakest sites"),
                                    .default = str_glue("4_{TREATMENT_NAME} signif. sites"))) %>% 
  mutate(delta_factor = as.factor(str_glue("{condition_name} ({n})")))

df.xy.arranged_colors <- df.xy %>% 
  select(delta,cols) %>% 
  distinct() %>%
  arrange(delta) %>% 
  pull(cols) 


METHOD_WINDOW_STR <- ifelse(WINDOW_SIZE == "null",METHOD,as.character(str_glue("{METHOD}_{WINDOW_SIZE}")))

title <- str_glue("Comparison: {EXPNUM}, {METHOD_WINDOW_STR}, {TREATMENT_NAME} / {CONTROL_NAME}\n",
                  "Treatment samples: {TREATMENT_SAMPLES}\n",
                  "Control samples: {CONTROL_SAMPLES}\n",
                  "Significant sites filters: |log2FC| > 1, padj < 0.05\n",
                  "Weakest sites filters: |log2FC| <= 1, padj < 0.05\n")


CHROM <- "all chromosomes"

p.xy <- ggplot(df.xy, aes(x=log2FC, fill = delta_factor))+
  geom_histogram(binwidth=.1)+
  scale_fill_manual(name = str_glue("Site_Category ({nrow(df.xy)} total sites)"), 
                    values = df.xy.arranged_colors, 
                    drop = FALSE)+
  ggtitle(str_glue("MUMERGE condition-specific sites ({CHROM})\n")) + 
  ylab("Number of Condition-specific Regions") + 
  xlab("log2(Fold Change)")+
  theme_classic()+
  theme(legend.background = element_rect(colour = "black"), 
        legend.text=element_text(size=10))


df.noxy <- df.with_colors %>% filter(!str_detect(seqnames,"chrX|chrY")) %>% 
  add_count(delta) %>% 
  arrange(delta) %>% 
  mutate(condition_name = case_when(str_detect(delta,"1_") ~ str_glue("1_{CONTROL_NAME} signif. sites"),
                                    str_detect(delta,"2_") ~ str_glue("2_{CONTROL_NAME} weakest sites"),
                                    str_detect(delta,"3_") ~ str_glue("3_{TREATMENT_NAME} weakest sites"),
                                    .default = str_glue("4_{TREATMENT_NAME} signif. sites"))) %>%
  mutate(delta_factor = as.factor(str_glue("{condition_name} ({n})")))


df.noxy.arranged_colors <- df.noxy %>% 
  select(delta,cols) %>% 
  distinct() %>%
  arrange(delta) %>% 
  pull(cols) 

CHROM <- "without XY chromosomes"
p.noxy <- ggplot(df.noxy, aes(x=log2FC, fill = delta_factor))+
  geom_histogram(binwidth=.1)+
  scale_fill_manual(name = str_glue("Site_Category ({nrow(df.noxy)} total sites)"), 
                    values = df.xy.arranged_colors, 
                    drop = FALSE)+
  ggtitle(str_glue("MUMERGE condition-specific sites ({CHROM})\n")) + 
  ylab("Number of Condition-specific Regions") + 
  xlab("log2(Fold Change)")+
  theme_classic()+
  theme(legend.background = element_rect(colour = "black"), 
        legend.text=element_text(size=10))


histograms <- p.xy + p.noxy+plot_annotation(title = title)
ggsave(argv$output_name, plot = histograms, height = 9, width = 18)
