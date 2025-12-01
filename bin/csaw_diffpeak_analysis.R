#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(csaw)) 
suppressPackageStartupMessages(library(edgeR)) 
suppressPackageStartupMessages(library(argparser)) 
suppressPackageStartupMessages(library(plyranges)) 
suppressPackageStartupMessages(library(writexl)) 
suppressPackageStartupMessages(library(readr)) 
suppressPackageStartupMessages(library(tidyr)) 


ParseArguments <- function() {
  p <- arg_parser('CSAW processing')
  p <- add_argument(p,'--control_samples', default="ANY", help="list of samples separated by pipe operator. ex. 'G123_M1|G123_M2' which represent control condition")
  p <- add_argument(p,'--treatment_samples', default="ANY", help="list of samples separated by pipe operator. ex. 'G123_M1|G123_M2' which represent treatment condition")
  p <- add_argument(p,'--output_prefix', default="NONE", help="name for group of comparisons")
  p <- add_argument(p,'--library', default="both", help="paired/single-end information [both, none, first, second]")
  p <- add_argument(p,'--mumerge', default="./mumerge.bed", help="mumerge file for comparison")
  
  return(parse_args(p))
}

argv <- ParseArguments()

DEBUG <- FALSE
if (DEBUG) {
  argv$treatment_samples <- c("G228_M03","G228_M05","G241_M22","G241_M23") %>% paste0(collapse = "|")
  argv$control_samples <- c("G228_M04","G228_M06","G241_M45","G241_M46") %>% paste0(collapse = "|")
  argv$mumerge <- "/projectnb/wax-dk/max/G242_G241_G228_DAR_ATAC/RESULTS_G242_G241_G228_ATAC_DAR_allgroups_DESEQ2/mumerge_output/Male_8wk_ATAC_vs_Female_8wk_ATAC_mumerge.bed"
}

TREATMENT_SAMPLES <- argv$treatment_samples
CONTROL_SAMPLES <- argv$control_samples
OUTPUT_PREFIX <- argv$output_prefix
LIBRARY <- ifelse(argv$library == "paired","both","none")
MUMERGE <- argv$mumerge
  
bam.files <- list.files(pattern = "*bam$")
treatment_bam <- grep(TREATMENT_SAMPLES, bam.files, value = T)
control_bam <- grep(CONTROL_SAMPLES, bam.files, value = T)
n_treatment <- length(treatment_bam)
n_control <- length(control_bam)

## reorder bam files as [Treatment, Control]
bam.files <- c(treatment_bam, control_bam)

# Sample information
sample.info <- data.frame(
  Sample = c(paste0(rep("Treatment",n_treatment),1:n_treatment), 
             paste0(rep("Control",n_control),1:n_control)),
  Group = factor(c(rep("Treatment",n_treatment), rep("Control",n_control))),
  stringsAsFactors = FALSE
)

param <- readParam(
  pe = LIBRARY,           # Paired-end reads for ATAC-seq
  restrict = paste0("chr", c(1:22, "X", "Y"))
)

design <- model.matrix(~sample.info$Group)
colnames(design) <- levels(sample.info$Group)
rownames(design) <- sample.info$Sample
contrast <- makeContrasts(Treatment - Control, levels = design)

data <- windowCounts(
  bam.files, 
  ext = 200,      # 100 Fragment extension for ATAC-seq (typically ~200bp fragments)
  width = 150,    # Window width - smaller for ATAC-seq than ChIP-seq
  spacing = 30,   # Window spacing
  param = param
)

binned <- windowCounts(
  bam.files, 
  bin = TRUE, 
  width = 10000,    # Large bins for background estimation
  param = param
)

keep <- filterWindowsGlobal(data, binned)$filter > log2(5)
data <- data[keep,]

data <- normFactors(binned, se.out = data)

y <- asDGEList(data)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design, robust = TRUE)

results <- glmQLFTest(fit)  # Test for the cell.type coefficient

merged <- mergeResults(data, results$table, tol = 500L)

#fdr_cutoff <- 0.05
# sig.regions <- merged$regions[merged$combined$FDR <= fdr_cutoff]
# sig_idx <- which(merged$combined$FDR <= fdr_cutoff)
# rep_tests <- merged$combined$rep.test[sig_idx]
# 
# mcols(sig.regions)$num.tests <- merged$combined$num.tests[sig_idx]
# mcols(sig.regions)$num.up.logFC <- merged$combined$num.up.logFC[sig_idx]
# mcols(sig.regions)$num.down.logFC <- merged$combined$num.down.logFC[sig_idx]
# mcols(sig.regions)$logFC <- results$table$logFC[rep_tests]
# mcols(sig.regions)$logCPM <- results$table$logCPM[rep_tests]
# mcols(sig.regions)$PValue <- merged$combined$PValue[sig_idx]
# mcols(sig.regions)$FDR <- merged$combined$FDR[sig_idx]
# mcols(sig.regions)$direction <- merged$combined$direction[sig_idx]
# 
# 
# results.df <- as.data.frame(sig.regions)
# results.df <- results.df[order(results.df$FDR),]
# 
# col_order <- c("seqnames", "start", "end", "width", 
#                "num.tests", "num.up.logFC", "num.down.logFC",
#                "logFC", "logCPM", "PValue", "FDR", "direction")
# results.df <- results.df[, col_order]

# Save ALL results for exploration
all.results.df <- data.frame(
  seqnames = as.character(seqnames(merged$regions)),
  start = start(merged$regions),
  end = end(merged$regions),
  pval = merged$combined$PValue,
  padj = merged$combined$FDR,
  log2fc = results$table$logFC[merged$combined$rep.test], ## EdgeR logFC = log2FC
  #logCPM = results$table$logCPM[merged$combined$rep.test],
  #direction = merged$combined$direction,
  stringsAsFactors = FALSE
)

all.results.df <- all.results.df[order(all.results.df$pval),]
# tt <- all.results.df %>% filter(abs(log2FC) > 1 & pval < 0.05) %>% as_tibble()
# readr::write_csv(tt, "filtered_csaw.csv", col_names = T)
# readr::write_csv(all.results.df, "all_results_csaw.csv", col_names = T)

xlsx_output_name<-paste0(OUTPUT_PREFIX,".xlsx")

mu <- readr::read_tsv(MUMERGE, col_names = F) %>%
  select(seqnames = 1, start = 2, end = 3)

csaw_mu <- plyranges::join_overlap_left(mu %>% as_granges(),
                                        all.results.df %>% as_granges()) %>% 
  as_tibble() %>%
  drop_na()

writexl::write_xlsx(list(
  CSAW = all.results.df,
  CSAW_mumerge_centric = csaw_mu),
  path = xlsx_output_name, col_names = T)





