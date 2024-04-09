#!/usr/bin/env Rscript

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("plyranges", update = F)

library(argparser)

ParseArguments <- function() {
  p <- arg_parser('Postprocessing of diffreps outputs')
  p <- add_argument(p,'--annotated_path', default="input.annotated", help="path to diffReps annotated output")
  p <- add_argument(p, '--hotspot_path', default="hotspots", help="path to diffReps hotspots output")
  p <- add_argument(p, '--blmm9_path', default="default", help="path to mm9 blacklisted regions (optional)")
  p <- add_argument(p, '--macs2_xls_dir_path', default = "./", help = "path to MACS2 xls files / SICER scoreiland files")
  p <- add_argument(p, '--sample_labels_path', default = "./Sample_Labels.txt", help = "path to Sample_Labels.txt, need for sample_id and description columns")
  p <- add_argument(p, '--min_avg_count', default = 20, help = "Threshold for average number of reads per peak for CONTROL and TREATMENT")
  p <- add_argument(p, '--log2fc_cutoff', default = 1, help = "Threshold for fold change (logarithmic scale)")
  p <- add_argument(p, '--control_name', default = "CONTROL", help = "Control condition name")
  p <- add_argument(p, '--treatment_name', default = "TREATMENT", help = "Treatment condition name")
  p <- add_argument(p, '--peak_caller', default = "MACS2", help = "MACS2 or SICER")
  p <- add_argument(p, '--normalization_caller', default = "DIFFREPS", help = "DIFFREPS or RIPPM")
  p <- add_argument(p, '--histone_mark', default = "default", help = "Dataset_id + Histone mark (ex.G215_K27ac)")
  p <- add_argument(p, '--treatment_samples', default = "", help = "Treatment samples (ex. G215M1,G215M2)")
  p <- add_argument(p, '--control_samples', default = "", help = "Control samples (ex. G215M3,G215M4)")
  
  return(parse_args(p))
}

argv <- ParseArguments()

print(argv)

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(patchwork)
library(plyranges)
library(openxlsx) #library(writexl)

## input params
annotated_path <- argv$annotated_path
hotspot_path <- argv$hotspot_path
blmm9_path <- argv$blmm9_path
macs2_xls_dir_path <- argv$macs2_xls_dir_path
sample_labels_path <- argv$sample_labels_path

min_avg_count <- argv$min_avg_count
log2fc_cutoff <- argv$log2fc_cutoff
CONTROL_NAME <- argv$control_name
TREATMENT_NAME <- argv$treatment_name

peak_caller <- argv$peak_caller
histone_mark <- argv$histone_mark
normalization_caller <- argv$normalization_caller
treatment_samples <- argv$treatment_samples %>% str_replace_all(., "\\|",",")
control_samples <- argv$control_samples %>% str_replace_all(., "\\|",",")

log2fc_label <- 2^log2fc_cutoff
swap_colors <- TRUE

DEBUG <- FALSE

if (DEBUG){
  setwd("/projectnb/wax-dk/max/exp_diffreps_genomicRanges/")
  ## input params
  # annotated_path <- "/projectnb2/wax-dk/max/G207_SG_CHIPSEQ_K27ac/Scripts/10_diffReps_1/Output_diffReps_1.diffreps_normalization.w1000/diffReps_Female_Control.vs.Male_Treatment.annotated"
  # hotspot_path <- "/projectnb2/wax-dk/max/G207_SG_CHIPSEQ_K27ac/Scripts/10_diffReps_1/Output_diffReps_1.diffreps_normalization.w1000/diffReps_Female_Control.vs.Male_Treatment.hotspot"
  # blmm9_path <- "/projectnb/wax-dk/max/G215_k4me1_k27ac_k9ac/Scripts/10_diffreps_G215K27ac_Male_vs_G215K27ac_Female_RIPPM_3/Job_Scripts/ENCODE_Blacklist/mm9-blacklist.bed.gz"
  # macs2_xls_dir_path <- "/projectnb/wax-dk/max/G207_SG_CHIPSEQ_K27ac/Scripts/10_diffReps_1/Output_diffReps_1.diffreps_normalization.w1000/XLSfiles"
  # sample_labels_path <- "/projectnb/wax-dk/max/G207_SG_CHIPSEQ_K27ac/Scripts/00_Setup_Pipeline/Sample_Labels.txt"
  
  main_path <- "/projectnb2/wax-dk/max/G222_CHIPSEQ/G222_G156_G207/work/9f/d635e5be5219b694578e9b59f0e849/"
  annotated_path <- str_glue("{main_path}/diffReps_Male_8wk_H3K27ac.vs.Male_3wk_H3K27ac_DIFFREPS_1000.annotated")
  hotspot_path <- str_glue("{main_path}/diffReps_Male_8wk_H3K27ac.vs.Male_3wk_H3K27ac_DIFFREPS_1000.hotspot")
  blmm9_path <- "default"
  macs2_xls_dir_path <- str_glue("{main_path}/XLSfiles")
  sample_labels_path <- str_glue("{main_path}/sample_labels.csv")
  
  peak_caller <- "MACS2"
  min_avg_count <- 20
  log2fc_cutoff <- 1
  CONTROL_NAME <- "Female"
  TREATMENT_NAME <- "Male"
  
  histone_mark <- "G207_K27ac"
  normalization_caller <- "DIFFREPS"
  treatment_samples <- "G207M1,G207M2"
  control_samples <- "G207M3,G207M4"
}

## "G215M1,G215M2" --> "G215M1M2"
## will not work properly for samples from different datasets "G215M1,G216M2" --> "G215_M1M2"
## in this case just remove sample names from file names, this information present inside files
collapse_sample_names <- function(sample_string) {
  dataset_number <- str_extract(sample_string,"^G\\d+")
  mnumbers <- str_extract_all(sample_string,"M\\d+", simplify = T) %>% 
    str_c(., sep = "", collapse = "")
  str_glue("{dataset_number}_{mnumbers}")
}

short_treatment_names <- collapse_sample_names(treatment_samples)
short_control_names <- collapse_sample_names(control_samples)

top_header <- str_glue("{histone_mark}, {peak_caller}, normalization: {normalization_caller}\n",
                       "Treatment: {treatment_samples} Control: {control_samples}\n")


######## Sample_labels
sample_labels <- read_csv(sample_labels_path, col_names = T) %>% 
  #select(sample_id = Sample_ID, description = Description) %>%
  select(sample_id, description = sample_description) %>% 
  mutate(description = as.character(str_glue("{sample_id}_{description}"))) %>% 
  tibble::deframe()

######## MACS2 xls for extended info
######## TODO: add similar to SICER

peakcaller_xls <- if (peak_caller == "MACS2") {
  lapply(list.files(path = macs2_xls_dir_path, full.names = T), function(fname){
    fname_prefix <- str_extract(basename(fname), "G[[:alnum:]]+_?M[[:alnum:]]+")
    
    ## select only specific columns
    tmp <- read_tsv(fname, col_names = T, comment = "#") %>%
      select(seqnames = chr, start, end, length, pileup, fold_enrichment, dplyr::any_of(c("-log10(qvalue)","minus_log10_qvalue")), peak_id = name) %>% ## -log10(qvalue)
      #rename_with(., function(x){"minus_log10_qvalue"}, starts_with("-")) %>%  ## rename -log10(qvalue) to minus_log10_qvalue if exist
      mutate(sample_id = sample_labels[fname_prefix],   ## assing sample_id_description instead of sample_id
             start = start - 1)                         ## start shift to 1, because MACS2 bed files have same shift, I just aligned to output MACS2 bed files
  }) %>% 
    purrr::reduce(bind_rows) %>% 
    as_granges()
} else {
  ## SICER
  lapply(list.files(path = macs2_xls_dir_path, full.names = T), function(fname){
    fname_prefix <- str_extract(basename(fname), "G[[:alnum:]]+_?M[[:alnum:]]+")
    
    ## select only specific columns
    tmp <- read_tsv(fname, col_names = F, comment = "#") %>%
      select(seqnames = X1, start = X2, end = X3, fold_enrichment = X4) %>% ## -log10(qvalue)
      mutate(sample_id = sample_labels[fname_prefix])
  }) %>% 
    purrr::reduce(bind_rows) %>% 
    as_granges()
}

## MACS2/SICER merging peaks
peakcaller_union  <- if (peak_caller == "MACS2") {
  peakcaller_xls %>% 
    as_tibble() %>% 
    select(seqnames,start,end,peak_id) %>% 
    as_granges() %>% 
    GenomicRanges::reduce(., min.gapwidth = 0L) %>% 
    plyranges::mutate(peakcaller_overlap = 1)
} else {
  ## SICER does not have peak_id
  peakcaller_xls %>% 
    as_tibble() %>% 
    select(seqnames,start,end) %>% 
    as_granges() %>% 
    GenomicRanges::reduce(., min.gapwidth = 0L) %>% 
    plyranges::mutate(peakcaller_overlap = 1)
  }

######## load diffReps output annotated or "with header" files
gr.annotated <- read_tsv(annotated_path, col_names = T, comment = "#") %>% 
  dplyr::rename(seqnames = Chrom, start = Start, end = End) %>% 
  # replace_na(list(strand = "*")) %>% 
  as_granges()

######## load hotspots file
gr.hotspots <- read_tsv(hotspot_path, col_names = T, comment = "#") %>% 
  dplyr::rename(seqnames = Chrom, start = Start, end = End) %>% 
  as_granges()

######## load mm9 blacklisted gz file
if (blmm9_path != "default") {
  gr.bl <- read_tsv(blmm9_path, col_names = F) %>% 
    select(seqnames = X1, start = X2, end = X3) %>% 
    mutate(strand = "*") %>% 
    as_granges()
  
  ## removing blacklisted regions from annotated file
  gr.ann.noblack <- gr.annotated %>% filter_by_non_overlaps(gr.bl)
  
  ## removing blacklisted regions from hotspots file
  gr.hotspots.noblack <- gr.hotspots %>% filter_by_non_overlaps(gr.bl)
  
} else {
  gr.ann.noblack <- gr.annotated
  gr.hotspots.noblack <- gr.hotspots
}


##### appending meta columns #####
gr.ann.noblack.extra <- gr.ann.noblack

######## create meta columns with MACS2/SICER intersection
gr.ann.noblack.extra <- join_overlap_left(gr.ann.noblack.extra, peakcaller_union) %>% 
  as_tibble() %>% 
  replace_na(list(peakcaller_overlap = 0)) %>% 
  distinct() %>% 
  as_granges()

######## create meta columns with up/down sites
gr.ann.noblack.signif.down <- gr.ann.noblack %>% 
  mutate(down_significant = 1) %>% 
  filter(Event == "Down", log2FC < -(log2fc_cutoff), Control.avg > min_avg_count, padj < 0.05) %>%   
  select(down_significant)

gr.ann.noblack.signif.up <- gr.ann.noblack %>% 
  mutate(up_significant = 1) %>% 
  filter(Event == "Up", log2FC > log2fc_cutoff, Treatment.avg > min_avg_count, padj < 0.05) %>%   
  select(up_significant)

gr.ann.noblack.extra <- join_overlap_left(gr.ann.noblack.extra, gr.ann.noblack.signif.down)
gr.ann.noblack.extra <- join_overlap_left(gr.ann.noblack.extra, gr.ann.noblack.signif.up)

gr.ann.noblack.extra <- gr.ann.noblack.extra %>% 
  as_tibble() %>% 
  replace_na(list(down_significant = 0,
                  up_significant = 0)) %>% 
  as_granges()

####### significant or not
gr.ann.noblack.extra <- gr.ann.noblack.extra %>% 
  mutate(significant = as.integer(down_significant | up_significant)) %>% 
  as_tibble() %>% 
  distinct() %>% 
  #filter(padj < 0.05) %>%
  as_granges()

###### add delta columns with information about significant, weak and low reads sites
gr.ann.noblack.extra <- gr.ann.noblack.extra %>% 
  as_tibble() %>% 
  mutate(delta = case_when(Event == "Down" & abs(log2FC) > log2fc_cutoff & Control.avg > min_avg_count & padj < 0.05  ~ str_glue("1_{CONTROL_NAME}_Signif_sites"), 
                           Event == "Down" & abs(log2FC) <= log2fc_cutoff & Control.avg > min_avg_count & padj < 0.05 ~ str_glue("2_{CONTROL_NAME}_Weakest_sites"), 
                           Event == "Up" & abs(log2FC) > log2fc_cutoff & Treatment.avg > min_avg_count & padj < 0.05 ~ str_glue("4_{TREATMENT_NAME}_Signif_sites"), 
                           Event == "Up" & abs(log2FC) <= log2fc_cutoff & Treatment.avg > min_avg_count & padj < 0.05 ~ str_glue("3_{TREATMENT_NAME}_Weakest_sites"),
                           abs(log2FC) > log2fc_cutoff ~ "0_Low_reads_sites", # low read region ## & (Treatment.avg <= min_avg_count | Control.avg <= min_avg_count)
                           .default = NA)) %>% 
  as_granges()


add_colors <- function(df_tmp, hist_colors) {
  
  df <- df_tmp %>% as_tibble()
  delta.ord <- df %>% 
    select(delta) %>% 
    distinct() %>% 
    arrange(delta) %>% 
    mutate(ord = as.numeric(str_extract(delta, "\\d+"))+1) 
  
  hist.ord <- tibble(cols = hist_colors) %>% 
    mutate(ord = row_number()) 
  
  delta.col <- left_join(x = delta.ord, y = hist.ord, join_by(ord)) %>% 
    distinct() %>% 
    arrange(ord) %>% 
    select(-ord)
  
  left_join(df, delta.col, join_by(delta)) %>% as_granges()
  
}

histogram_colors <- if (swap_colors){
  ## down - red, up - blue
  c("gray", "red","pink", "lightblue","blue")
} else {
  ## up - red, down - blue
  c("gray", "blue", "lightblue", "pink", "red")
}

gr.ann.noblack.extra <- add_colors(gr.ann.noblack.extra, histogram_colors)
# gr.ann.noblack.extra %>% as_tibble() %>% glimpse()
gr.ann.noblack.extra %>% as_tibble() %>% pull(cols) %>% table

#### PLOT DATA
#### FDR (0.05, 0.01, 0.005, 0.001)
#### input: gr.ann.noblack.extra
#### filtration: by "significant" == 1
FDRs <- c(0.05, 0.01, 0.005, 0.001)


fdrs_plot <- function(fdrs, df, peakcaller_overlap_flag = F, extra_title = "") {
  
  ## if some Up/Down is NA we have to assign them to 0
  fake <- tibble(FDR = rep(c(FDRs),2), Event = rep(c("Up","Down"), each = 4))

  df.barplot_fdr <- map_dfr(fdrs, function(fdr) {
    df %>% 
      as_tibble() %>% 
      dplyr::filter(if (!peakcaller_overlap_flag) {
        significant == 1 & padj < fdr
      } else {
        significant == 1 & padj < fdr & peakcaller_overlap == 1
      }
      ) %>% 
      #filter(significant == 1 & padj < fdr & peakcaller_overlap == 1) %>% 
      dplyr::select(Event) %>% 
      dplyr::summarise(n = dplyr::n(), .by = "Event") %>% 
      dplyr::mutate(FDR = fdr) 
  }) %>% left_join(fake,.) %>% 
    replace_na(list(n = 0))
  
  plot_by_fdr <- ggplot(df.barplot_fdr, aes(factor(FDR), n, fill=Event)) +
    geom_bar(stat="identity", position = "dodge") + 
    scale_fill_brewer(palette = "Set1", direction = -1) +
    ggtitle(extra_title) + 
    ylab("Count of Condition-specific Regions") + 
    xlab("FDR Threshold") +
    guides(fill=guide_legend(title="Event")) +
    geom_text(aes(label = n), size = 3, position=position_dodge(width=0.9), vjust=-0.25)
}

unfiltered_fdrplot_title <- str_glue("Unfiltered {TREATMENT_NAME}/{CONTROL_NAME}.\nEffect of FDR cutoff on number condition-specific sites\n",
                            "Filtering options: |log2FC| > {log2fc_cutoff}, avg.count > {min_avg_count}\n")
unfiltered_fdrplot <- fdrs_plot(FDRs, 
                                gr.ann.noblack.extra, 
                                peakcaller_overlap = F, 
                                extra_title = unfiltered_fdrplot_title)

filtered_fdrplot_title <- str_glue("{peak_caller} filtered {TREATMENT_NAME}/{CONTROL_NAME}.\nEffect of FDR cutoff on number condition-specific sites\n",
                                    "Filtering options: |log2FC| > {log2fc_cutoff}, avg.count > {min_avg_count}\n")
filtered_fdrplot <- fdrs_plot(FDRs, 
                              gr.ann.noblack.extra, 
                              peakcaller_overlap = T, 
                              extra_title = filtered_fdrplot_title)

final_fdr_barchart <- unfiltered_fdrplot+filtered_fdrplot+plot_annotation(title = top_header)

# output_name_fdrbarchart <- str_glue("FDR_Barchart_{histone_mark}_{TREATMENT_NAME}_{short_treatment_names}_vs_{CONTROL_NAME}_{short_control_names}.pdf")
output_name_fdrbarchart <- str_glue("FDR_Barchart_{TREATMENT_NAME}_vs_{CONTROL_NAME}_{peak_caller}_{normalization_caller}.pdf")
ggsave(output_name_fdrbarchart, plot = final_fdr_barchart, height = 7, width = 14)



#### Histogram
#### input: gr.ann.noblack.extra


plot_histogram_bk <- function(df, delta_column, log2fc_cutoff, min_avg_count, filter_by_peakcaller_overlap = F, title_extra = "", log2fc_label = 1, swap_colors = T) {
  
  histogram_colors <- if (swap_colors){
    ## down - red, up - blue
    c("pink", "lightblue", "grey","red", "blue")
  } else {
    ## up - red, down - blue
    c("lightblue", "pink", "grey","blue", "red")    
  }
  
  df.histogram <- df %>% 
    as_tibble() %>% 
    #filter(padj < 0.05) %>% 
    select(delta, log2FC, peakcaller_overlap) %>% 
    {
      if(filter_by_peakcaller_overlap){
        filter(., peakcaller_overlap == 1)
      } else {
        .
      }} %>%
    drop_na(delta) %>% 
    add_count(delta) %>% 
    mutate(delta = as.factor(str_glue("{delta} ({n})"))) %>% 
    select(-n)
    
  ggplot(df.histogram, aes(x=log2FC, fill = delta))+ #factor(delta, levels = names(cols_vector))
    geom_histogram(binwidth=.1)+ ## alpha = 0.9
    # scale_fill_manual(name = "Site_Category", values = as.vector(cols_vector), labels = names(cols_vector))+
    scale_fill_manual(name = str_glue("Site_Category ({nrow(df.histogram)} total sites)"), 
                      values = histogram_colors, ## c("lightblue", "pink", "grey","blue", "red")
                      drop = FALSE)+
    ggtitle(str_glue("{title_extra}")) + 
    ylab("Count of Condition-specific Regions") + 
    xlab("log2(Fold Change)")+
    theme_classic()+
    theme(legend.background = element_rect(colour = "black"), 
          legend.text=element_text(size=10))
}

plot_histogram <- function(df, filter_by_peakcaller_overlap = F, title_extra = "") { ## df should be tibble not granges
  df.histogram <- df %>% 
    select(delta, log2FC, peakcaller_overlap, cols) %>% 
    {
      if(filter_by_peakcaller_overlap){
        filter(., peakcaller_overlap == 1)
      } else {
        .
      }} %>%
    drop_na(delta) %>% 
    add_count(delta) %>% 
    arrange(delta) %>% 
    mutate(delta = as.factor(str_glue("{delta} ({n})"))) %>% 
    select(-n, -peakcaller_overlap) 
  
  arranged_colors <- df.histogram %>% 
    select(delta,cols) %>% 
    drop_na() %>%
    distinct() %>%
    arrange(delta) %>% 
    pull(cols) 
  
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


title_unfiltered = str_glue("Unfiltered {TREATMENT_NAME} / {CONTROL_NAME}.\nFold Change for diffReps condition-specific sites\n",
                            "Significant sites filters: |log2FC| > {log2fc_cutoff}, padj < 0.05, avg.count > {min_avg_count}\n",
                            "Marginal sites filters: |log2FC| > {log2fc_cutoff}, padj < 0.05, avg.count <= {min_avg_count}\n",
                            "Less_{log2fc_label}-fold filters: |log2FC| <= {log2fc_cutoff}, padj < 0.05, avg.count > {min_avg_count}\n")

hist_unfiltered <- plot_histogram(gr.ann.noblack.extra %>% as_tibble(),
                                  filter_by_peakcaller_overlap = F, 
                                  title_extra = title_unfiltered)

title_filtered = str_glue("{peak_caller} filtered {TREATMENT_NAME} / {CONTROL_NAME}.\nFold Change for diffReps condition-specific sites\n",
                          "Significant sites filters: |log2FC| > {log2fc_cutoff}, padj < 0.05, avg.count > {min_avg_count}\n",
                          "Marginal sites filters: |log2FC| > {log2fc_cutoff}, padj < 0.05, avg.count <= {min_avg_count}\n",
                          "Less_{log2fc_label}-fold filters: |log2FC| <= {log2fc_cutoff}, padj < 0.05, avg.count > {min_avg_count}\n")

hist_filtered <- plot_histogram(gr.ann.noblack.extra %>% as_tibble(),
                                  filter_by_peakcaller_overlap = T, 
                                  title_extra = title_filtered)

histograms <- hist_unfiltered + hist_filtered+plot_annotation(title = top_header)
#output_name_histograms <- str_glue("Histograms_{histone_mark}_{TREATMENT_NAME}_{short_treatment_names}_vs_{CONTROL_NAME}_{short_control_names}.pdf")
output_name_histograms <- str_glue("Histograms_{TREATMENT_NAME}_vs_{CONTROL_NAME}_{peak_caller}_{normalization_caller}.pdf")
ggsave(output_name_histograms, plot = histograms, height = 9, width = 18)

#### barplot by features
#### input: gr.ann.noblack.extra
#### filtration: by "up/down_significant" == 1
#S1_diff_data <- diffReps_output[diffReps_output$Event =="Down" & diffReps_output$log2FC < -(log2FC_cutoff) & diffReps_output$Control.avg > min_avg_count & diffReps_output$padj < 0.05,]

feature_colors <- tribble(
    ~Feature, ~color,
    "ProximalPromoter",   "dodgerblue",
    "Promoter1k",   "firebrick",
    "Promoter3k",   "forestgreen",
    "Genebody", "darkorchid4",
    "Genedesert", "dodgerblue4",
    "Pericentromere", "darkorange",
    "Subtelomere", "midnightblue",
    "OtherIntergenic", "gold")


barchart_feature <- function(colors, df, title){

  ## include all colors not only active   
  df_all_colors <- left_join(colors, df) %>% 
    replace_na(list(n = 0, pct = 0)) %>% 
    arrange(Feature)
  
  # df_all_colors %>% print
  
  color_vector <- df_all_colors %>% 
    select(Feature, color) %>% 
    tibble::deframe()
  
  # color_vector %>% print
  
  ggplot(df_all_colors, aes(x = Feature, y = pct, fill = Feature, label = factor(n))) + 
    geom_col(position = 'dodge') + 
    geom_text(position = position_dodge(width = .9),    # move to center of bars
              vjust = -0.5,    # nudge above top of bar
              size = 3) + 
    scale_y_continuous(labels = scales::percent)+
    #scale_fill_manual(values = as.vector(color_vector), labels = names(color_vector)) +
    scale_fill_manual(values = df_all_colors$color, labels = df_all_colors$Feature) +
    ylab("Percent")+
    ggtitle(title)+
    theme_light()+
    theme(axis.text.x = element_blank())
  
}

## UNFILTERED 
unf_df.up_feature <- gr.ann.noblack.extra %>% 
  as_tibble() %>% 
  filter(up_significant == 1) %>% 
  count(Feature = factor(Feature)) %>% 
  mutate(pct = prop.table(n))

title <- str_glue("Unfiltered Feature Distribution of\n{TREATMENT_NAME} Sites ({sum(unf_df.up_feature$n)})")
unfiltered_up_barchart <- barchart_feature(feature_colors, unf_df.up_feature, title = title)


unf_df.down_feature <- gr.ann.noblack.extra %>% 
  as_tibble() %>% 
  filter(down_significant == 1) %>% 
  count(Feature = factor(Feature)) %>% 
  mutate(pct = prop.table(n))

title <- str_glue("Unfiltered Feature Distribution of\n{CONTROL_NAME} Sites ({sum(unf_df.down_feature$n)})")
unfiltered_down_barchart <- barchart_feature(feature_colors, unf_df.down_feature, title = title)

## FILTERED
filt_df.up_feature <- gr.ann.noblack.extra %>% 
  as_tibble() %>% 
  filter(up_significant == 1 & peakcaller_overlap == 1) %>% 
  count(Feature = factor(Feature)) %>% 
  mutate(pct = prop.table(n))

title <- str_glue("{peak_caller} filtered Feature Distribution of\n{TREATMENT_NAME} Sites ({sum(filt_df.up_feature$n)})")
filtered_up_barchart <- barchart_feature(feature_colors, filt_df.up_feature, title = title)

filt_df.down_feature <- gr.ann.noblack.extra %>% 
  as_tibble() %>% 
  filter(down_significant == 1 & peakcaller_overlap == 1) %>% 
  count(Feature = factor(Feature)) %>% 
  mutate(pct = prop.table(n))

title <- str_glue("{peak_caller} filtered Feature Distribution of\n{CONTROL_NAME} Sites ({sum(filt_df.down_feature$n)})")
filtered_down_barchart <- barchart_feature(feature_colors, filt_df.down_feature, title = title)

#final_plot <- list(up_barchart, down_barchart) %>% keep(\(x) is.ggplot(x)) %>% purrr::reduce(`+`)
final_feature_barchart <- (unfiltered_up_barchart+unfiltered_down_barchart)/(filtered_up_barchart+filtered_down_barchart)+
  plot_layout(guides = "collect")+
  plot_annotation(title = top_header)

# output_name_barchart <- str_glue("Barchart_{histone_mark}_{TREATMENT_NAME}_{short_treatment_names}_vs_{CONTROL_NAME}_{short_control_names}.pdf")
output_name_barchart <- str_glue("Barchart_{TREATMENT_NAME}_vs_{CONTROL_NAME}_{peak_caller}_{normalization_caller}.pdf")
ggsave(output_name_barchart, plot = final_feature_barchart, width = 11, height = 8.5)



#### EXPORT DATA ####
######## export hotspots as bed file
# hotspot_fname <- str_glue("diffReps_Hotspots_{histone_mark}_{TREATMENT_NAME}_{short_treatment_names}_vs_{CONTROL_NAME}_{short_control_names}_UCSC.bed")
hotspot_fname <- str_glue("diffReps_Hotspots_UCSC_{TREATMENT_NAME}_vs_{CONTROL_NAME}_{peak_caller}_{normalization_caller}.bed")
hotspot_header <- str_glue("track name=Hotspots_{TREATMENT_NAME}_vs_{CONTROL_NAME}_{peak_caller}_{normalization_caller} visibility=4 itemRgb=On")
write_lines(hotspot_header, hotspot_fname)

gr.hotspots.noblack %>% 
  as_tibble() %>% 
  select(seqnames, start, end) %>% 
  mutate(
    t0 = "Hotspot",
    t1 = 1000,
    t2 = ".",
    t3 = 0,
    t4 = 0,
    t5 = "0,0,0") %>% 
  write_tsv(file = hotspot_fname, append = T, col_names = F)

### top200,600, bed output
# diffReps_output_sorted <- diffReps_output[order(-diffReps_output$"log2FC") & diffReps_output$padj < 0.05,]

### screenshots, diffreps output
# S1_diff_screenshot <- diffReps_output[diffReps_output$Event =="Down" & diffReps_output$log2FC < -(log2FC_cutoff) & diffReps_output$Control.avg > 100 & diffReps_output$padj < 0.05,]
# S1_diff_screenshot <- S1_diff_screenshot[order(S1_diff_screenshot$"log2FC"),]
# S2_diff_screenshot <- diffReps_output[diffReps_output$Event =="Up" & diffReps_output$log2FC > (log2FC_cutoff) & diffReps_output$Treatment.avg > 200 & diffReps_output$padj < 0.05,]
# S2_diff_screenshot <- S2_diff_screenshot[order(-S2_diff_screenshot$"log2FC"),]

### Control/Treatment bed files with different colors
ucsc_fname_unfiltered <- str_glue("UCSC_UNFILTERED_track_{TREATMENT_NAME}_vs_{CONTROL_NAME}_{peak_caller}_{normalization_caller}.bed")
ucsc_header_unfiltered <- str_glue("track name=UNFILTERED_{histone_mark}_{TREATMENT_NAME}_vs_{CONTROL_NAME}_{peak_caller}_{normalization_caller} visibility=4 itemRgb=On")
write_lines(ucsc_header_unfiltered, ucsc_fname_unfiltered)

gr.ann.noblack.extra %>% 
  as_tibble() %>% 
  filter(str_detect(delta,"1_|2_|3_|4_")) %>% 
  drop_na(delta) %>% 
  select(seqnames, start, end, delta, cols) %>% 
  rowwise() %>% 
  mutate(t1 = 1000,
         t2 = ".",
         t3 = 0,
         t4 = 0,
         t5 = paste0(col2rgb(cols), collapse = ",")) %>%
  select(-cols) %>% 
  arrange(seqnames,start) %>% 
  write_tsv(file = ucsc_fname_unfiltered, append = T, col_names = F)
  


## FILTERED, significant + peak_caller overlap
ucsc_fname_filtered <- str_glue("UCSC_FILTERED_track_{TREATMENT_NAME}_vs_{CONTROL_NAME}_{peak_caller}_{normalization_caller}.bed")
ucsc_header_filtered <- str_glue("track name=FILTERED_{histone_mark}_{TREATMENT_NAME}_vs_{CONTROL_NAME}_{peak_caller}_{normalization_caller} visibility=4 itemRgb=On")
write_lines(ucsc_header_filtered, ucsc_fname_filtered)

gr.ann.noblack.extra %>% 
  as_tibble() %>% 
  filter(peakcaller_overlap == 1 & str_detect(delta,"1_|2_|3_|4_")) %>% 
  drop_na(delta) %>% 
  select(seqnames, start, end, delta, cols) %>% 
  rowwise() %>% 
  mutate(t1 = 1000,
         t2 = ".",
         t3 = 0,
         t4 = 0,
         t5 = paste0(col2rgb(cols), collapse = ",")) %>%
  select(-cols) %>% 
  arrange(seqnames,start) %>% 
  write_tsv(file = ucsc_fname_filtered, append = T, col_names = F)

### XLSX files for filtered and unfiltered
### TODO: export in csv format
### UNFILTERED XLSX

unfiltered_xls <- gr.ann.noblack.extra %>% 
  as_tibble() %>% 
  select(-cols) %>% 
  relocate(any_of(c("peakcaller_overlap","down_significant","up_significant","significant")), .after = Control.avg) %>% 
  dplyr::rename(`Overlapped with peak caller` = peakcaller_overlap, 
                `Significant CONTROL` = down_significant, 
                `Significant TREATMENT` = up_significant, 
                `Significant by default thresholds` = significant) %>% 
  select(-Control.enr, -Treatment.enr) %>% 
  mutate(ucsc_coords = str_glue("{seqnames}:{start}-{end}")) %>% 
  relocate(ucsc_coords, .before = seqnames) %>% 
  relocate(c("Overlapped with peak caller", 
             "Significant CONTROL", 
             "Significant TREATMENT", 
             "Significant by default thresholds"), 
           .before = seqnames) %>% 
  arrange(seqnames, start, end) 

### create sheet for unfiltered data
sheet_name <- str_glue("Unfiltered_sites")
wb <- createWorkbook()
addWorksheet(wb, sheetName = sheet_name)
writeData(wb, sheet = sheet_name, str_glue("Unfiltered sites: {TREATMENT_NAME} (treatment)/{CONTROL_NAME} (control)"), startRow = 1, startCol = 1, colNames = FALSE)
writeData(wb, sheet = sheet_name, str_glue("Default thresholds: |log2FC| > {log2fc_cutoff}, padj < 0.05, avg.count > {min_avg_count}"), startRow = 2, startCol = 1, colNames = FALSE)
#writeData(wb, sheet = sheet_name, str_glue("TREATMENT samples: {short_treatment_names}, CONTROL samples: {short_control_names}"), startRow = 3, startCol = 1, colNames = FALSE)
writeData(wb, sheet = sheet_name, str_glue("TREATMENT samples: {treatment_samples}, CONTROL samples: {control_samples}"), startRow = 3, startCol = 1, colNames = FALSE)
writeData(wb, sheet = sheet_name, unfiltered_xls, startRow = 5, startCol = 1)

# write_xlsx(unfiltered_xls, path = "Unfiltered_sites_output.xlsx")

### FILTERED XLSX
### required extended information about peaks

filtered_xls <- join_overlap_left(gr.ann.noblack.extra %>% filter(peakcaller_overlap == 1), peakcaller_xls) %>%
  as_tibble() %>% 
  select(-cols) %>% 
  distinct() %>% 
  {if (peak_caller == "MACS2") {mutate(., fake_score = fold_enrichment*pileup*length)} else {mutate(.,fake_score = fold_enrichment)}} %>% 
  group_by(seqnames, start, end, sample_id) %>% 
  filter(fake_score == max(fake_score)) %>% 
  ungroup() %>% 
  select(-fake_score) %>% 
  add_count(seqnames,start,end, name = "overlap_with_n_samples") %>% 
  {
    if(peak_caller == "MACS2") {
      group_by(., seqnames,start,end) %>% 
        mutate(peak_ids = str_c(peak_id, collapse = ",")) %>% 
        ungroup() %>% 
        select(-peak_id) %>% 
        pivot_longer(c("length", "pileup", "fold_enrichment", "minus_log10_qvalue"), values_to = "param")
    } else{
      pivot_longer(., "fold_enrichment", values_to = "param")
    }
  } %>% 
  pivot_wider(names_from = c("sample_id","name"), values_from = param, names_sort = T) %>% 
  mutate(ucsc_coords = str_glue("{seqnames}:{start}-{end}")) %>% 
  relocate(ucsc_coords, .before = seqnames) %>% 
  arrange(desc(overlap_with_n_samples), seqnames, start, end) %>% 
  dplyr::rename(`Overlapped with peak caller` = peakcaller_overlap, 
                `Significant CONTROL` = down_significant, 
                `Significant TREATMENT` = up_significant, 
                `Significant by default thresholds` = significant,
                `How many samples overlap with this region` = overlap_with_n_samples) %>% 
  relocate(c("Overlapped with peak caller", 
             "Significant CONTROL", 
             "Significant TREATMENT", 
             "Significant by default thresholds", 
             "How many samples overlap with this region"), 
           .before = seqnames)

sheet_name <- str_glue("{peak_caller}_filtered_sites")
addWorksheet(wb, sheetName = sheet_name)
writeData(wb, sheet = sheet_name, str_glue("{peak_caller} filtered sites: {TREATMENT_NAME} (treatment)/{CONTROL_NAME} (control)"), startRow = 1, startCol = 1, colNames = FALSE)
writeData(wb, sheet = sheet_name, str_glue("Default thresholds: |log2FC| > {log2fc_cutoff}, padj < 0.05, avg.count > {min_avg_count}"), startRow = 2, startCol = 1, colNames = FALSE)
#writeData(wb, sheet = sheet_name, str_glue("TREATMENT samples: {short_treatment_names}, CONTROL samples: {short_control_names}"), startRow = 3, startCol = 1, colNames = FALSE)
writeData(wb, sheet = sheet_name, str_glue("TREATMENT samples: {treatment_samples}, CONTROL samples: {control_samples}"), startRow = 3, startCol = 1, colNames = FALSE)
writeData(wb, sheet = sheet_name, filtered_xls, startRow = 5, startCol = 1)

# saveWorkbook(wb, str_glue("Summary_{histone_mark}_{TREATMENT_NAME}_{short_treatment_names}_vs_{CONTROL_NAME}_{short_control_names}.xlsx"), overwrite = T)
saveWorkbook(wb, str_glue("Summary_{TREATMENT_NAME}_vs_{CONTROL_NAME}_{normalization_caller}.xlsx"), overwrite = T)




















