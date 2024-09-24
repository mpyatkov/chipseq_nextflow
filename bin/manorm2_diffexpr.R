#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)

library(readxl)
library(writexl)

# remotes::install_cran("MAnorm2", upgrade = "never")
library(MAnorm2)

# remotes::install_cran("argparser", upgrade = "never")
library(argparser)

# remotes::install_cran("plyranges", upgrade = "never")
# library(plyranges)
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

BiocManager::install("plyranges")
library(plyranges)

ParseArguments <- function() {
  p <- arg_parser('MAnorm2 processing')
  p <- add_argument(p,'--manorm2_profile', default="profile.xls", help="input file with profile obtained by 'profile_bins' function")
  p <- add_argument(p,'--control_samples', default="", help="list of samples separated by pipe operator. ex. 'G123_M1|G123_M2' which represent control condition")
  p <- add_argument(p,'--treatment_samples', default="", help="list of samples separated by pipe operator. ex. 'G123_M1|G123_M2' which represent treatment condition")
  p <- add_argument(p,'--control_name', default="", help="control group name")
  p <- add_argument(p,'--treatment_name', default="", help="treatment group name")
  p <- add_argument(p,'--peakcaller', default="MACS2", help = "MACS2/EPIC2/SICER peakcaller name")
  p <- add_argument(p,'--output_prefix', default="", help="name for group of comparisons")
  p <- add_argument(p,'--exp_number', default="00", help="Experiment number")

  return(parse_args(p))
}

argv <- ParseArguments()

print(argv)

DEBUG <- F
if (DEBUG) {
  # argv$manorm2_profile <-  "/projectnb/wax-dk/max/G222_CHIPSEQ/G222_G156_G207/work/55/3314c2c968a9e297771a3f8b8cebea/Male_3wk_H3K27ac_vs_Female_3wk_H3K27ac_profile_bins.xls"
  # argv$treatment_samples <- "G222_M16|G222_M17"
  argv$manorm2_profile <- "/projectnb/wax-dk/max/G222_CHIPSEQ/G222_G156_G207/work/8e/5f595d5c7ba276620c299856439cb5/Female_8wk_H3K27ac_vs_Female_3wk_H3K27ac_profile_bins.xls"
  argv$treatment_samples <- "G207_M03|G207_M04|G222_M11|G222_M12"
  argv$control_samples <- "G222_M18|G222_M19|G222_M20"
  argv$treatment_name <- "Male_8wk_H3K27ac"
  argv$control_name <- "Male_3wk_H3K27ac"
}

min_avg_count <- 20
log2fc_cutoff <- 1
log2fc_label <- 2^log2fc_cutoff
swap_colors <- T ## swap red/blue colors on histogram and output tracks

## colors for histograms and output tracks
histogram_colors <- if (swap_colors){
  ## down - red, up - blue
  c("red","pink", "lightblue","blue")
} else {
  ## up - red, down - blue
  c("blue", "lightblue", "pink", "red")
}

### manorm2 read and filtering input file
mx.manorm2 <- read_tsv(argv$manorm2_profile, col_names = T) %>% 
  mutate(rowsums = purrr::reduce(select(., contains("read_cnt")), `+`),
         occsums = purrr::reduce(select(., contains("occup")), `+`)) %>% 
  filter(rowsums > min_avg_count, occsums > 1) %>% 
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
mn.data <- MAnorm2::normalize(mx.manorm2, count = treatment_cnts, occupancy = treatment_ocuppancy)
mn.data <- MAnorm2::normalize(mn.data, count = control_cnts, occupancy = control_ocuppancy)

condsmn <- list(Treatment = bioCond(mn.data[treatment_cnts], mn.data[treatment_ocuppancy], name = str_glue("{argv$treatment_name}.treatment")),
                Control = bioCond(mn.data[control_cnts], mn.data[control_ocuppancy], name = str_glue("{argv$control_name}.control")))

autosome <- !(mx.manorm2$chrom %in% c("chrX", "chrY"))
condsmn <- normBioCond(condsmn, common.peak.regions = autosome)
condsmn <- fitMeanVarCurve(condsmn, method = "parametric", occupy.only = TRUE , init.coef = c(0.1, 10))
resmn <- diffTest(condsmn[[2]],condsmn[[1]]) ## Treatment vs Control

#resmn %>% as_tibble() %>% filter(padj < 0.05 & abs(Mval) > 1) %>% dplyr::summarize(n = dplyr::n())

mx.manorm2.diff <- mx.manorm2 %>% 
  mutate(rowid = as.character(row_number())) %>% 
  inner_join(.,resmn %>% tibble::rownames_to_column(var="rowid")) %>% 
  select(-rowid)


mx.manorm2.diff <- mx.manorm2.diff %>% 
  mutate(delta = case_when(Mval < 0 & abs(Mval) > log2fc_cutoff & padj < 0.05 ~ str_glue("1_{argv$control_name}_Signif_sites"),
                           Mval < 0 & abs(Mval) <= log2fc_cutoff & padj < 0.05 ~ str_glue("2_{argv$control_name}_Weakest_sites"),
                           Mval > 0 & abs(Mval) > log2fc_cutoff & padj < 0.05 ~  str_glue("4_{argv$treatment_name}_Signif_sites"),
                           Mval > 0 & abs(Mval) <= log2fc_cutoff & padj < 0.05 ~  str_glue("3_{argv$treatment_name}_Weakest_sites"),
                           .default = NA))
  
# mx.manorm2.diff %>% arrange(delta) %>% writexl::write_xlsx(path = str_glue("Summary_{argv$output_prefix}.xlsx"), col_names = T)
exp_number<-str_pad(argv$exp_number, width = 2, pad = "0", side = "left")
xlsx_output_name<-str_glue("{exp_number}_{argv$output_prefix}.xlsx")
mx.manorm2.diff %>% arrange(delta) %>% writexl::write_xlsx(path = xlsx_output_name, col_names = T)
  
#### add colors
add_colors <- function(df, hist_colors) {

  delta.ord <- df %>% select(delta) %>% 
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

mx.manorm2.diff <- add_colors(mx.manorm2.diff, hist_colors = histogram_colors)

##### CREATE HISTOGRAMS #####
get_df_histogram <- function(df, filter_xy_chromosomes = F){

  df.histogram <- df %>% 
    {
      if (filter_xy_chromosomes){
        filter(., ! chrom %in% c("chrX", "chrY"))
      } else {
        .
      }
    } %>%
    filter(padj < 0.05) %>% 
    select(log2FC = Mval, delta) %>% 
    drop_na(delta) %>% 
    add_count(delta) %>% 
    arrange(delta) %>% 
    mutate(delta = as.factor(str_glue("{delta} ({n})"))) %>% 
    select(-n) #%>% 
    # mutate(ord = as.numeric(str_extract(delta, "\\d+")))

  arranged_colors <- df %>% 
    {
      if (filter_xy_chromosomes){
        filter(., ! chrom %in% c("chrX", "chrY"))
      } else {
        .
      }
    } %>%
    select(delta,cols) %>% 
    drop_na() %>%
    distinct() %>%
    arrange(delta) %>% 
    pull(cols) 

    list(df.histogram = df.histogram, arranged_colors = arranged_colors)
}

plot_histogram <- function(df.histogram, log2fc_cutoff, title_extra = "", log2fc_label = 1, arranged_colors) {

  ggplot(df.histogram, aes(x=log2FC, fill = delta))+ #factor(delta, levels = names(cols_vector))
    geom_histogram(binwidth=.1)+ ## alpha = 0.9
    scale_fill_manual(name = str_glue("Site_Category ({nrow(df.histogram)} total sites)"), 
                      values = arranged_colors, 
                      drop = FALSE)+
    ggtitle(str_glue("{title_extra}")) + 
    ylab("Count of Condition-specific Regions") + 
    xlab("log2(Fold Change)")+
    theme_classic()+
    theme(legend.background = element_rect(colour = "black"), 
          legend.text=element_text(size=10))
}

## replace separators for samples representation
str_control_samples <- str_replace_all(argv$control_samples, "\\|", ",")
str_treatment_samples <- str_replace_all(argv$treatment_samples, "\\|", ",")

hist.data_allchr <- get_df_histogram(mx.manorm2.diff, filter_xy_chromosomes = F)
title_manorm2 <-  str_glue("MANORM2 (all chromosomes) {argv$peakcaller}\n{argv$treatment_name} / {argv$control_name}\nTreatment: {str_treatment_samples}\nControl: {str_control_samples}\n",
                           "Fold Change for MAnorm2 condition-specific sites\n",
                           "Significant sites filters: |log2FC| > {log2fc_cutoff}, padj < 0.05, avg.count > {min_avg_count}\n",
                           "Weakest_sites filters: |log2FC| <= {log2fc_cutoff}, padj < 0.05, avg.count > {min_avg_count}\n")

hist.plot <- plot_histogram(hist.data_allchr$df.histogram, log2fc_cutoff = log2fc_cutoff, title_extra = title_manorm2, log2fc_label = log2fc_label, arranged_colors = hist.data_allchr$arranged_colors)
ggsave(str_glue("{argv$output_prefix}_AllChr_histogram.pdf"), plot = hist.plot, height = 9, width = 9)


hist.data_noXY <- get_df_histogram(mx.manorm2.diff, filter_xy_chromosomes = T)
title_manorm2_noxy <-  str_glue("MANORM2 (without X and Y chrom)\n {argv$peakcaller}{argv$treatment_name} / {argv$control_name}\nTreatment: {str_treatment_samples}\nControl: {str_control_samples}\n",
                           "Fold Change for MAnorm2 condition-specific sites\n",
                           "Significant sites filters: |log2FC| > {log2fc_cutoff}, padj < 0.05, avg.count > {min_avg_count}\n",
                           "Weakest_sites filters: |log2FC| <= {log2fc_cutoff}, padj < 0.05, avg.count > {min_avg_count}\n")

hist.plot <- plot_histogram(hist.data_noXY$df.histogram, log2fc_cutoff = log2fc_cutoff, title_extra = title_manorm2_noxy, log2fc_label = log2fc_label, arranged_colors = hist.data_noXY$arranged_colors)
ggsave(str_glue("{argv$output_prefix}_noXY_histogram.pdf"), plot = hist.plot, height = 9, width = 9)


##### CREATE BED TRACKS #####
ucsc_fname_filtered <- str_glue("UCSC_FILTERED_track_{argv$treatment_name}_vs_{argv$control_name}_{argv$peakcaller}_MANORM2.bed")
ucsc_header_filtered <- str_glue("track name=FILTERED_{argv$treatment_name}_vs_{argv$control_name}_{argv$peakcaller}_MANORM2 visibility=4 itemRgb=On")
write_lines(ucsc_header_filtered, ucsc_fname_filtered)
mx.manorm2.diff %>% 
  drop_na(delta) %>% 
  select(chrom, start, end, delta, cols) %>% 
  rowwise() %>% 
  mutate(t1 = 1000,
         t2 = ".",
         t3 = 0,
         t4 = 0,
         t5 = paste0(col2rgb(cols), collapse = ",")) %>%
  select(-cols) %>% 
  arrange(chrom,start) %>% 
  write_tsv(file = ucsc_fname_filtered, append = T, col_names = F)
  

