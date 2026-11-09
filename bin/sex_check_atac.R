#!/usr/bin/env Rscript
#----------------------------------------------------------------------------------
# ATAC-seq sex check. Predicts each sample's sex from promoter-window counts of
# marker genes (featureCounts, SAF), normalized by RIPPM (reads in peaks per
# million), compares the prediction with the declared sex (optional "sex"
# column of the sample labels), checks that the samples of each condition
# agree, and writes sex_confirmation.xlsx with these sheets:
#  - samples:    sample_id, condition, Y_markers_RIPPM, <x_marker>_RIPPM,
#                Y_marker_reads, <x_marker>_reads, predicted_sex
#                [, sex, config_vs_prediction_sex_confirmation]
#  - conditions: condition, n_samples, n_M, n_F, n_UNCERTAIN,
#                group_sex_consistency
#  - markers:    sample_id, condition, <gene>_RIPPM per marker gene
#  - legend:     item, value, description
#
# Two votes: the summed Y markers vote M at or above --y_min_rippm (else F),
# the X marker votes F at or above --x_min_rippm (else M). Agreeing votes give
# M/F; disagreeing votes, or fewer than --min_marker_reads raw marker reads,
# give UNCERTAIN.
#
# RIPPM = raw count * norm_factor (Norm_Factors.tsv, column norm_factor).
# Default thresholds are tuned for TSS +/- 2 kb windows (y 5, x 5); for
# +/- 1 kb windows use --x_min_rippm 2.5.
#
# Every problem (NOT CONFIRMED or UNCERTAIN sample, MIXED condition) is one
# line of --warnings_file.
#----------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr)
  library(writexl)
  library(argparser)
})

p <- arg_parser("ATAC-seq sample sex prediction and confirmation workbook")
p <- add_argument(p, "--counts",           help = "featureCounts output (SAF marker windows, all samples in one file)")
p <- add_argument(p, "--norm_factors",     help = "TSV with columns sample_id and norm_factor (RIPPM factor)")
p <- add_argument(p, "--samples",          help = "CSV with sample_id, a condition column, and optional sex (M/F)")
p <- add_argument(p, "--condition_col",    help = "column of --samples used as condition", default = "group_description")
p <- add_argument(p, "--bam_suffix",       help = "regex removed from featureCounts column names to get sample ids",
                  default = "(_sorted_filtered)?\\.bam$")
p <- add_argument(p, "--y_markers",        help = "comma-separated Y-chromosome marker ids", default = "Kdm5d,Uty,Eif2s3y,Ddx3y")
p <- add_argument(p, "--x_marker",         help = "X-inactivation marker id", default = "Xist")
p <- add_argument(p, "--y_min_rippm",      help = "summed Y-marker RIPPM at or above which the Y vote is M", type = "numeric", default = 5)
p <- add_argument(p, "--x_min_rippm",      help = "X-marker RIPPM at or above which the X vote is F", type = "numeric", default = 5)
p <- add_argument(p, "--min_marker_reads", help = "minimum raw reads (Y markers + X marker) for a call", type = "numeric", default = 10)
p <- add_argument(p, "--output_file",      help = "output .xlsx path", default = "sex_confirmation.xlsx")
p <- add_argument(p, "--warnings_file",    help = "output warnings file, one line per problem", default = "sex_check_warnings.txt")
argv <- parse_args(p)
print(argv)

split_csv <- function(s) if (is.na(s) || s == "") character(0) else trimws(strsplit(s, ",", fixed = TRUE)[[1]])
y_markers <- split_csv(argv$y_markers)
x_marker  <- argv$x_marker
markers   <- c(y_markers, x_marker)

## ---- inputs ----------------------------------------------------------------
samples <- read_csv(argv$samples, col_types = cols(.default = col_character()), na = character())
if (!"sample_id" %in% names(samples)) stop("--samples has no sample_id column")
if (!argv$condition_col %in% names(samples))
  stop("--samples has no column '", argv$condition_col, "' (set --condition_col)")
samples$condition <- samples[[argv$condition_col]]

# Optional declared sex: M/F; empty means not declared for that sample.
declared <- "sex" %in% names(samples)
if (declared) {
  samples$sex <- toupper(trimws(samples$sex))
  samples$sex[samples$sex == ""] <- NA
  bad_sex <- setdiff(na.omit(unique(samples$sex)), c("M", "F"))
  if (length(bad_sex) > 0) stop("sex column must contain M or F, found: ", paste(bad_sex, collapse = ", "))
  declared <- any(!is.na(samples$sex))
}

nf <- read_tsv(argv$norm_factors, col_types = cols(.default = col_character()))
if (!all(c("sample_id", "norm_factor") %in% names(nf))) stop("--norm_factors needs columns sample_id and norm_factor")
nf <- setNames(as.numeric(nf$norm_factor), nf$sample_id)

fc <- read_tsv(argv$counts, comment = "#", col_types = cols(.default = col_character()))
count_cols <- names(fc)[7:ncol(fc)]                    # Geneid, Chr, Start, End, Strand, Length, samples...
count_ids  <- sub(argv$bam_suffix, "", basename(count_cols))
counts <- sapply(count_cols, function(col) as.numeric(fc[[col]]))
counts <- matrix(counts, nrow = nrow(fc), dimnames = list(fc$Geneid, count_ids))

absent <- setdiff(markers, rownames(counts))
if (length(absent) > 0) stop("marker(s) not found in --counts: ", paste(absent, collapse = ", "))
missing_counts <- setdiff(samples$sample_id, colnames(counts))
if (length(missing_counts) > 0) stop("no counts for sample(s): ", paste(missing_counts, collapse = ", "))
missing_nf <- setdiff(samples$sample_id, names(nf))
if (length(missing_nf) > 0) stop("no norm_factor for sample(s): ", paste(missing_nf, collapse = ", "))
extra <- setdiff(colnames(counts), samples$sample_id)
if (length(extra) > 0) message("ignored (not in --samples): ", paste(extra, collapse = ", "))

## ---- marker signals (RIPPM) -------------------------------------------------
raw   <- counts[markers, samples$sample_id, drop = FALSE]
rippm <- sweep(raw, 2, nf[samples$sample_id], `*`)
y_rippm <- colSums(rippm[y_markers, , drop = FALSE])
x_rippm <- rippm[x_marker, ]
y_reads <- colSums(raw[y_markers, , drop = FALSE])
x_reads <- raw[x_marker, ]

## ---- calls -----------------------------------------------------------------
call_sex <- function(y, x, reads) {
  if (is.na(y) || is.na(x) || reads < argv$min_marker_reads) return("UNCERTAIN")
  y_vote <- if (y >= argv$y_min_rippm) "M" else "F"
  x_vote <- if (x >= argv$x_min_rippm) "F" else "M"
  if (y_vote == x_vote) y_vote else "UNCERTAIN"
}
predicted <- mapply(call_sex, y_rippm, x_rippm, y_reads + x_reads, USE.NAMES = FALSE)

x_col <- paste0(x_marker, "_RIPPM")
x_reads_col <- paste0(x_marker, "_reads")
per_sample <- data.frame(
  sample_id       = samples$sample_id,
  condition       = samples$condition,
  Y_markers_RIPPM = round(y_rippm, 2),
  X_RIPPM         = round(x_rippm, 2),
  Y_marker_reads  = as.integer(y_reads),
  X_reads         = as.integer(x_reads),
  predicted_sex   = predicted,
  stringsAsFactors = FALSE, row.names = NULL
)
names(per_sample)[names(per_sample) == "X_RIPPM"] <- x_col
names(per_sample)[names(per_sample) == "X_reads"] <- x_reads_col

warnings_out <- character(0)

if (declared) {
  per_sample$sex <- samples$sex
  per_sample$config_vs_prediction_sex_confirmation <- ifelse(
    is.na(per_sample$sex), "NOT DECLARED",
    ifelse(per_sample$sex == per_sample$predicted_sex, "CONFIRMED", "NOT CONFIRMED"))
  bad <- per_sample[per_sample$config_vs_prediction_sex_confirmation == "NOT CONFIRMED", ]
  warnings_out <- c(warnings_out, sprintf("sex-check: sample %s NOT CONFIRMED - declared %s, predicted %s (Y_markers_RIPPM=%s, %s=%s)",
                                          bad$sample_id, bad$sex, bad$predicted_sex, bad$Y_markers_RIPPM, x_col, bad[[x_col]]))
}
# A declared sample that is UNCERTAIN was already warned as NOT CONFIRMED,
# so only undeclared UNCERTAIN samples get their own line.
unc <- per_sample[per_sample$predicted_sex == "UNCERTAIN" &
                    (if (declared) is.na(per_sample$sex) else TRUE), ]
warnings_out <- c(warnings_out, sprintf("sex-check: sample %s UNCERTAIN - Y_markers_RIPPM=%s, %s=%s, marker reads=%d",
                                        unc$sample_id, unc$Y_markers_RIPPM, x_col, unc[[x_col]],
                                        unc$Y_marker_reads + unc[[x_reads_col]]))

## ---- condition consistency ---------------------------------------------------
conditions <- unique(samples$condition)
per_condition <- do.call(rbind, lapply(conditions, function(cond) {
  pr <- per_sample$predicted_sex[per_sample$condition == cond]
  n_m <- sum(pr == "M"); n_f <- sum(pr == "F")
  data.frame(condition = cond, n_samples = length(pr), n_M = n_m, n_F = n_f,
             n_UNCERTAIN = sum(pr == "UNCERTAIN"),
             group_sex_consistency = if (n_m > 0 && n_f > 0) "MIXED" else "CONSISTENT",
             stringsAsFactors = FALSE)
}))
mixed <- per_condition[per_condition$group_sex_consistency == "MIXED", ]
warnings_out <- c(warnings_out, sprintf("sex-check: condition %s MIXED - %d predicted M, %d predicted F",
                                        mixed$condition, mixed$n_M, mixed$n_F))

## ---- per-marker sheet (values only) ------------------------------------------
marker_sheet <- data.frame(sample_id = samples$sample_id, condition = samples$condition,
                           round(t(rippm), 2), check.names = FALSE, row.names = NULL)
names(marker_sheet)[-(1:2)] <- paste0(markers, "_RIPPM")

## ---- legend ----------------------------------------------------------------
legend_rows <- list(
  c("counts file", basename(argv$counts), "featureCounts output the marker counts are read from"),
  c("Y markers", paste(y_markers, collapse = ", "), "Y-chromosome genes summed into Y_markers_RIPPM"),
  c("X marker", x_marker, paste0("X-inactivation gene reported as ", x_col)),
  c("--y_min_rippm", format(argv$y_min_rippm), "Y vote is M when Y_markers_RIPPM >= this value, otherwise F"),
  c("--x_min_rippm", format(argv$x_min_rippm), paste0("X vote is F when ", x_col, " >= this value, otherwise M")),
  c("--min_marker_reads", format(argv$min_marker_reads), "fewer raw reads on the Y markers plus the X marker than this gives UNCERTAIN"),
  c("predicted_sex", "M / F / UNCERTAIN", "M or F when the Y vote and the X vote agree; UNCERTAIN when they disagree or the marker reads are below --min_marker_reads"),
  c("Y_markers_RIPPM", "", "summed Y-marker counts x norm_factor (reads in peaks per million)"),
  c(x_col, "", "X-marker counts x norm_factor (reads in peaks per million)"),
  c("Y_marker_reads", "", "raw summed Y-marker read pairs (featureCounts, MAPQ filter as run)"),
  c(x_reads_col, "", "raw X-marker read pairs"),
  c("condition", argv$condition_col, "column of the sample labels used as condition")
)
if (declared) legend_rows <- c(legend_rows, list(
  c("sex", "M / F", "sex declared in the sample labels"),
  c("config_vs_prediction_sex_confirmation", "CONFIRMED / NOT CONFIRMED / NOT DECLARED",
    "CONFIRMED when sex equals predicted_sex; an UNCERTAIN prediction is NOT CONFIRMED; NOT DECLARED when sex is empty")
))
legend_rows <- c(legend_rows, list(
  c("group_sex_consistency", "CONSISTENT / MIXED", "conditions sheet: MIXED when the condition has both M and F predictions")
))
legend <- setNames(as.data.frame(do.call(rbind, legend_rows), stringsAsFactors = FALSE),
                   c("item", "value", "description"))

## ---- outputs ---------------------------------------------------------------
invisible(write_xlsx(list(samples = per_sample, conditions = per_condition,
                markers = marker_sheet, legend = legend),
           argv$output_file, format_headers = TRUE))
writeLines(warnings_out, argv$warnings_file)
cat(warnings_out, sep = "\n")
#----------------------------------------------------------------------------------
