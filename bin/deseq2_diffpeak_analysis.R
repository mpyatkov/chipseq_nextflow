#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2)) 
suppressPackageStartupMessages(library(writexl)) 
suppressPackageStartupMessages(library(argparser)) 

ParseArguments <- function() {
  p <- arg_parser('DESEQ2 processing')
  p <- add_argument(p,'--coverage_file', default="coverage.csv", help="FeatureCounts output file for Treatment and Control samples")
  p <- add_argument(p,'--control_samples', default="", help="list of samples separated by pipe operator. ex. 'G123_M1|G123_M2' which represent control condition")
  p <- add_argument(p,'--treatment_samples', default="", help="list of samples separated by pipe operator. ex. 'G123_M1|G123_M2' which represent treatment condition")
  p <- add_argument(p,'--output_prefix', default="", help="name for group of comparisons")
  return(parse_args(p))
}

argv <- ParseArguments()

COVERAGE_FILE <- argv$coverage_file
TREATMENT_SAMPLES <- argv$treatment_samples
CONTROL_SAMPLES <- argv$control_samples
OUTPUT_PREFIX <- argv$output_prefix

DEBUG <- FALSE
if (DEBUG) {
  COVERAGE_FILE <- "/projectnb/wax-dk/max/src/rprojects/diffreps_manorm_deseq_overlap/Male_8wk_ATAC_vs_Female_8wk_ATAC_coverage.tsv" 
  TREATMENT_SAMPLES <- "G228_M03|G228_M05|G241_M22|G241_M23"
  CONTROL_SAMPLES <- "G228_M04|G228_M06|G241_M45|G241_M46"
  OUTPUT_PREFIX <- "OUTPUT_PREFIX"
}

count_data <- read.table(COVERAGE_FILE, 
                         header = TRUE, 
                         sep = "\t", 
                         stringsAsFactors = FALSE,
                         comment.char = "#")


# Extract count matrix (columns 7 onwards contain the count data)
count_matrix <- count_data[, 7:ncol(count_data)]
rownames(count_matrix) <- count_data$Geneid

# Define sample information
sample_names <- colnames(count_matrix)

# Remove .bam extension from sample names for cleaner output
#clean_sample_names <- gsub("_sorted_filtered.bam", "", sample_names)
#colnames(count_matrix) <- clean_sample_names

# treatment_samples <- unlist(strsplit(TREATMENT_SAMPLES, split = "\\|"))
# control_samples <- unlist(strsplit(CONTROL_SAMPLES, split = "\\|"))

## works with short names like G123_M01 without _sorted_filtered.bam tail
treatment_samples <- grep(TREATMENT_SAMPLES, sample_names, value = T)
control_samples <- grep(CONTROL_SAMPLES, sample_names, value = T)

# Create sample metadata
sample_info <- data.frame(
  sample = sample_names,
  condition = ifelse(sample_names %in% treatment_samples, "Treatment", "Control"),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_names

# Reorder sample_info to match count_matrix column order
sample_info <- sample_info[colnames(count_matrix), ]

# Convert condition to factor with Control as reference
sample_info$condition <- factor(sample_info$condition, levels = c("Control", "Treatment"))

# Filter low count genes (keep genes with at least 10 counts across all samples)
keep_genes <- rowSums(count_matrix) >= 10
count_matrix_filtered <- count_matrix[keep_genes, ]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtered,
                              colData = sample_info,
                              design = ~ condition)
dds <- DESeq(dds)
deseq2_results <- results(dds, contrast = c("condition", "Treatment", "Control"))

# Extract results with desired columns
deseq2_output <- data.frame(
  Geneid = rownames(deseq2_results),
  padj = deseq2_results$padj,
  log2fc = deseq2_results$log2FoldChange,
  stringsAsFactors = FALSE
)

deseq2_output[c("seqnames", "start", "end")] <- do.call(rbind, strsplit(as.character(deseq2_output$Geneid), "_"))
deseq2_output[c("start", "end")] <- lapply(deseq2_output[c("start", "end")], as.numeric)

# Remove rows with NA values
deseq2_output <- deseq2_output[complete.cases(deseq2_output), ]
deseq2_output <- deseq2_output[order(deseq2_output$padj),]
deseq2_output <- deseq2_output[c("seqnames","start","end","padj","log2fc")]
str(deseq2_output)
##deseq2_output <- subset(deseq2_output, padj < 0.05 & abs(log2fc) > 1, select = c(seqnames,start,end,log2fc,padj))
# write.table(deseq2_output, 
#             file = paste0(GROUP_ID,"_","DEseq2_SignifOnly.csv"),
#             sep = ",", quote = FALSE, row.names = FALSE)

xlsx_output_name<-paste0(OUTPUT_PREFIX,".xlsx")
writexl::write_xlsx(list(
  DESEQ2_mumerge_centric = deseq2_output),
  path = xlsx_output_name, col_names = T)
