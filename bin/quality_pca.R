#!/usr/bin/env Rscript

#remotes::install_cran("writexl", upgrade = "never")
remotes::install_cran("argparser", upgrade = "never")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("plyranges", update = F)
library(plyranges)

remotes::install_cran("argparser", upgrade = "never")
library(argparser)

library(ggplot2)
remotes::install_cran("ggfortify", upgrade = "never")
library(ggfortify)
library(ggrepel)

remotes::install_cran("formattable", upgrade = "never")
library(formattable)

library(patchwork)

remotes::install_cran("ggpubr", upgrade = "never")
library(ggpubr)

library(tidyverse)

ParseArguments <- function() {
  p <- arg_parser('PCA for MACS2 narrow peaks')
  p <- add_argument(p,'--treatment_name', default="", help="treatment group name")
  p <- add_argument(p,'--control_name', default="", help="control group name")
  p <- add_argument(p,'--treatment_samples', default="", help="list of samples separated by pipe operator. ex. 'G123_M1|G123_M2' which represent treatment condition")
  p <- add_argument(p,'--control_samples', default="", help="list of samples separated by pipe operator. ex. 'G123_M1|G123_M2' which represent control condition")
  p <- add_argument(p,'--remove_chrXY', flag = T, help="remove chromosomes X and Y from analysis")
  return(parse_args(p))
}

argv <- ParseArguments()
print(argv)

DEBUG <- F
if (DEBUG) {
  setwd("/projectnb/wax-dk/max/G223_H3K27ac_COMBINED/test2/")
  # argv$treatment_samples <- "G223_M01D|G223_M02|G223_M03|G223_M04|G223_M05"
  # argv$control_samples <- "G223_M06|G223_M07|G223_M08|G223_M09|G223_M10"
  # argv$treatment_name <- "Male_2wk_H3K27ac"
  # argv$control_name <- "Female_2wk_H3K27ac"
  # 
  # argv$treatment_name <- "Male_3wk_H3K27ac"
  # argv$control_name <- "Female_3wk_H3K27ac"
  # argv$treatment_samples <- "G223_M11|G223_M12|G223_M13|G223_M14|G222_M16|G222_M17"
  # argv$control_samples <- "G223_M15|G223_M16|G223_M17|G222_M18|G222_M19|G222_M20"
  
  argv$treatment_name <- "Male_8wk_H3K27ac"
  argv$control_name <- "Female_8wk_H3K27ac"
  argv$treatment_samples <- "G223_M44|G223_M45|G223_M46|G222_M09|G222_M10"
  argv$control_samples <- "G223_M47|G223_M48|G223_M49|G222_M11|G222_M12"
  argv$remove_chrXY <- T
}

macs2_all <- map(list.files(pattern = "xls"), \(f){
  
  sample_id <- str_remove(basename(f), "_narrow_MACS2_peaks.xls")
  group <- case_when(str_detect(f,argv$treatment_samples) ~ argv$treatment_name,
                     TRUE ~ argv$control_name)

  read_tsv(f, col_names = T, comment = "#") %>% 
    select(seqnames = chr, start, end, pileup) %>% 
    mutate(group = group, sample_id = sample_id)
}) %>% list_rbind()

if (argv$remove_chrXY) {
  macs2_all <- macs2_all %>% 
    filter(!str_detect(seqnames, "chrX|chrY"))
}

## normalization coefficients (min(sum(pileup))/sum(pileup_i))
macs2_all_norm <- macs2_all %>%
  select(sample_id, pileup) %>%
  summarise(n = sum(pileup), .by = sample_id) %>% 
  mutate(n = min(n)/n) %>% 
  #print() %>%  ## show normalization factors
  left_join(macs2_all, ., join_by(sample_id)) %>% 
  mutate(pileup = pileup*n) %>% 
  select(-n)


## stats 
stats_table <- macs2_all %>%
  select(sample_id, pileup, group) %>%
  add_count(sample_id, name = "n_peaks") %>% 
  summarise(pileup_sum = sum(pileup), n_peaks = max(n_peaks), .by = c(sample_id, group)) %>% 
  mutate(norm_fct = round(min(pileup_sum)/pileup_sum, 2),
         sample_id = str_glue("{sample_id}_{group}")) %>% 
  arrange(desc(group)) %>% 
  select(-group) %>% 
  mutate(pileup_sum = formattable::comma(pileup_sum, format = 'd'))

stats_plot <- stats_table %>% 
  ggtexttable(., rows = NULL, theme = ttheme("classic", base_size = 12)) 

macs2_union <- macs2_all_norm %>% 
  plyranges::as_granges() %>% 
  GenomicRanges::reduce(., min.gapwidth = 0L)

macs2_combined <- plyranges::join_overlap_left(macs2_union, macs2_all_norm %>% plyranges::as_granges()) %>% 
  as_tibble() %>% 
  distinct() %>% 
  filter(pileup == max(pileup), .by = c(seqnames,start,end,sample_id)) %>% 
  mutate(sample_id = str_glue("{sample_id}_{group}")) %>% 
  add_count(seqnames, start,end, name = "number_of_samples_in_region") %>% 
  filter(number_of_samples_in_region > 1) 

## PCA
pca_df <- macs2_combined %>% 
  pivot_wider(names_from = sample_id, values_from = pileup, values_fill = 0) %>% 
  select(-c(1:7)) %>% 
  t() %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample_id") %>% 
  left_join(., macs2_combined %>% select(sample_id, group) %>% distinct()) %>% 
  relocate(group, .after = sample_id) 

res.pca <- prcomp(pca_df %>% select(-sample_id,-group),  scale = T, center = T)  

# comment this if you would like to have sample_id_group notation on PCA plot
pca_df <- pca_df %>% rowwise() %>% mutate(sample_id = str_remove(sample_id, paste0("_",group)))

pca_plot <- ggplot2::autoplot(res.pca, data = pca_df, label = F) +
  #ggtitle("")+
  geom_point(size=5, shape=21,aes(fill = group), color = "black" ) +
  # geom_label_repel(aes(label = sample_id, 
  #                      fill = group,  
  #                      segment.color="black"), 
  #                  size = 5, 
  #                  color = "black")+
  geom_text_repel(aes(label = sample_id, 
                      #fill = group,  
                      segment.color="black"), size = 5)+
  guides(fill = guide_legend(title = "Group", override.aes = aes(label = "")))+
  theme_bw()+
  theme(text = element_text(size=13),
        legend.title = element_text(size=12),
        legend.text=element_text(size=12),
        legend.position = "top",
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))


## correlation
order_sort <- macs2_combined %>% 
  select(sample_id, group) %>% 
  distinct() %>% 
  arrange(desc(group), sample_id) %>% 
  mutate(sort_order = row_number()) %>% 
  select(-group)

method <- "pearson"
cor_df <- macs2_combined %>%
  left_join(., order_sort, join_by(sample_id)) %>% 
  arrange(sort_order) %>% select(-sort_order) %>% 
  pivot_wider(names_from = sample_id, values_from = pileup, values_fill = 0) %>% 
  select(-c(1:7)) %>% 
  cor(., method = method)

# get only lower triangle of matrix
cor_df[upper.tri(cor_df)] <- NA

# create pairwise table from matrix [Var1, Var2, value]
td_melt <- expand.grid(dimnames(cor_df)) %>%
  cbind(., value = as.vector(cor_df)) %>%
  drop_na() 
  
mn <- round(min(cor_df, na.rm = T), 3) - 0.01
mx <- round(max(cor_df, na.rm = T), 3)

# adaptive font size for tiles
font_size <- function(x) {
  f <- floor(-0.34*x+11.36)-2
  ifelse(f<2,2,f)
}

# plot correlation
cor_plot <- ggplot(data = td_melt, aes(Var1, Var2, fill = value))+
  #ggtitle("title")+
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = str_remove(round(value, 2), "0")), size = font_size(dim(cor_df)[1]))+
  scale_fill_gradientn(colours = rainbow(4),limits=c(mn,mx), space= "Lab",
                       name=paste0(method,"\ncorrelation")) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                   size = ifelse(dim(cor_df)[1] > 35, 7,12), hjust = 1),
        axis.text.y = element_text(size=ifelse(dim(cor_df)[1] > 35, 7,12)),
        text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text=element_text(size=12))+
  coord_fixed()+
  xlab("")+
  ylab("")

chrXY_status <- ifelse(argv$remove_chrXY, "Without X and Y chromosomes.", "All chromosomes.")
title <- str_glue("{argv$treatment_name} vs {argv$control_name} (Correlation and PCA plots based on pileup of the MACS2 narrow peaks)\n",
                  "Peak union contains {formattable::comma(ncol(pca_df), format = 'd')} peaks. {chrXY_status}\n",
                  "Treatment samples: {argv$treatment_samples}\n",
                  "Control samples: {argv$control_samples}\n")

layout <- "
AAABB
AAABB
AAACC
"

final_plot <- wrap_elements(cor_plot)+wrap_elements(pca_plot)+stats_plot+
  plot_layout(design = layout)+
  plot_annotation(title = title)

chrXY_suffix <- ifelse(argv$remove_chrXY, "noXY", "")
fname <- str_glue("{argv$treatment_name}_vs_{argv$control_name}_Correlation_PCA_{chrXY_suffix}.pdf")
ggsave(fname, plot = final_plot, width = 22, height = 12)
# fname <- str_glue("{argv$treatment_name}_vs_{argv$control_name}_Correlation_PCA.png")
# ggsave(fname, device =ragg::agg_png(res = 400), plot = final_plot, width = 18, height = 10)
 

