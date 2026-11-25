#!/usr/bin/env Rscript

library(argparser)

ParseArguments <- function() {
    p <- arg_parser('Parsing diffReps output')
    p <- add_argument(p,'--pattern', default="_vs_", help="search pattern for diffReps output files")
    p <- add_argument(p, '--path', default="../", help="path to the dir which contains diffReps output directories")
    p <- add_argument(p, '--rippm_report', default="Norm_Factors.tsv", help="path to the file with rippm stats")
    p <- add_argument(p, '--fastq_num_reads', default="numreads.csv", help="path to the file with rippm stats")
    p <- add_argument(p, '--deseq2_factors', default="", help="dir with *_norm_factors.tsv from DESeq2 (optional)")
    p <- add_argument(p, '--csaw_factors', default="", help="dir with *_norm_factors.tsv from CSAW (optional)")
    p <- add_argument(p, '--output', default="Summary_normalization_factors.xlsx", help="output filename")
    return(parse_args(p))
}

argv <- ParseArguments()

DEBUG <- FALSE
if (DEBUG){
  setwd("/projectnb/wax-dk/max/REFACTORING/work/b5/5710b9115798de703ce30594efce5e")
  argv$path <- "./"
  #argv$fastq_num_reads <- "/projectnb/wax-dk/max/G225/work/tmp/f6/f86064f81022d98eb2cadd37bbcaff/numreads.txt"
}

library(stringr)
library(purrr)
library(readr)
library(dplyr)
library(tidyr)
library(openxlsx)

rippm_report <- read_tsv(argv$rippm_report, col_names = T) %>%
             select(-norm_factor)

## raw reads from from R1 files for each sample
fastq_num_reads <- read_csv(argv$fastq_num_reads, col_names = T)

files <- list.files(path = argv$path, pattern = argv$pattern, recursive = T, full.names = T)
#%>% 
#    keep(function(x){str_detect(x,"_vs_") && str_detect(x,"annotated|hotspot", negate = T)})

## What we need to parse (diffReps output)
# Treatment files: G220M03_fragments.bed        G220M04_fragments.bed
# Control files: G220M05_fragments.bed  G220M06_fragments.bed
# Window size: 1000
# Treatment normalization constant: 2.23        0.28
# Control normalization constant: 0.12  1.37

extract_tab_separated <- function(header, tag) {
    keep(header, function(line) {str_detect(line, tag)}) %>% 
        str_split(":", simplify = T) %>% 
        str_trim() %>% 
        str_split(pattern = "\t") %>% .[[2]]
}

parse_diffreps_file <- function(f, header) {
  normalization_caller <- ifelse(str_detect(f,"RIPPM"), "RIPPM","DIFFREPS")
  
  ## extract conditions from filename
  prefix <- str_extract(basename(f),"(\\d+_[A-Za-z]+_\\d+_)")
  fname_no_prefix <- basename(f) %>% str_replace(prefix,"")
  cond_tr <- fname_no_prefix %>% str_extract("([[:print:]]+)(?=_vs_)")
  cond_ctrl <- fname_no_prefix %>% str_extract("(?<=_vs_)([[:print:]]+)")
  
  tr_files <- extract_tab_separated(header, "Treatment files") %>% str_replace("_frag[[:print:]]+", "")
  ctrl_files <- extract_tab_separated(header, "Control files") %>% str_replace("_frag[[:print:]]+", "")
  window <- extract_tab_separated(header, "Window size")
  norm_tr <- extract_tab_separated(header, "Treatment normalization constant")
  norm_ctrl <- extract_tab_separated(header, "Control normalization constant")
  tibble(
      fname = fname_no_prefix,#basename(f), #%>% str_split(.,"_RIPPM|_DIFFREPS",simplify = T) %>% .[1],
      group = c(rep("TREATMENT", length(tr_files)), rep("CONTROL", length(ctrl_files))),
      condition = c(rep(cond_tr, length(tr_files)), rep(cond_ctrl, length(ctrl_files))),
      sample_id = c(tr_files, ctrl_files),
      normfactors = as.numeric(c(norm_tr, norm_ctrl)),
      window = as.numeric(window),
      norm_caller = normalization_caller
  )
}

tst <- map_dfr(files, function(f){
    header <- read_lines(f, n_max = 32)
    #parse_diffreps_file(f, header, argv$normalization_caller, argv$cond_tr, argv$cond_ctrl)
    parse_diffreps_file(f, header)
})

l1 <- tst %>% 
  left_join(.,rippm_report, by = c("sample_id")) %>% 
  left_join(.,fastq_num_reads, join_by(sample_id)) %>% 
  relocate(num_raw_reads, .before = fragment_count) %>% 
  mutate(window = str_glue("{norm_caller}_{window}"),
         sample_id = str_glue("{sample_id}_{group}"),
         group = str_glue("{condition}_{group}")) %>% 
  select(-norm_caller, -condition) %>% 
  group_by(fname) %>% 
  group_split()

lgroup_to_ltables <- function(df) {
  
    rippm_table <- df %>% select(sample_id, starts_with("num_raw"), starts_with("fragment")) %>% distinct()

    norm_factors_table <- df %>% 
      select(-starts_with("fragment"), -starts_with("num_raw_reads")) %>% 
      pivot_wider(names_from = window, values_from = normfactors)

    fname <- norm_factors_table %>% select(fname) %>% distinct() %>% pull(fname) %>% 
        str_replace_all(., "\\.", "_")
    
    nf <- norm_factors_table %>% 
        select(-fname, -group)
    
    ## average by factors only
    all <- nf %>% 
        summarise_at(vars(starts_with("RIPPM") | starts_with("DIFFREPS")), list(~round(mean(.), 2))) %>% 
        mutate(sample_id = "all") %>% 
        select(sample_id, everything())
    
    ## average by groups/factors
    by_groups <- norm_factors_table %>% 
        group_by(group) %>% 
        summarise_at(vars(starts_with("RIPPM") | starts_with("DIFFREPS")), list(~round(mean(.), 2))) %>% 
        mutate(group = str_glue("Group_{group}")) %>% 
        select(sample_id = group, everything())

    ratio_treat_ctrl <- by_groups %>% 
        pivot_longer(starts_with("RIPPM") | starts_with("DIFFREPS"), names_to = "factor", values_to = "values") %>% 
        mutate(sample_id = ifelse(str_detect(sample_id, "TREATMENT"), "TREATMENT", "CONTROL")) %>% 
        pivot_wider(names_from = sample_id, values_from = values) %>% 
        rowwise() %>% 
        mutate(TR_CTRL = round(TREATMENT/CONTROL,2)) %>% 
        select(-TREATMENT,-CONTROL) %>% 
        pivot_wider(names_from = factor, values_from = TR_CTRL) %>% 
        mutate(sample_id = "Treatment/Control ratio") %>% 
        select(sample_id, everything())
    
    res <- bind_rows(nf, all, by_groups, ratio_treat_ctrl)
    list(name =fname, 
         df = res, 
         rippm_df = rippm_table,
         nrow = nrow(res), ## data + header+fname
         ncol = ncol(res))
}

tt1 <- map(l1, lgroup_to_ltables)

###
### DESeq2 and CSAW factors (one *_norm_factors.tsv per comparison)
###

read_factor_tables <- function(dir) {
    if (is.na(dir) || dir == "" || !dir.exists(dir)) return(NULL)
    fs <- list.files(dir, pattern = "_norm_factors\\.tsv$", full.names = TRUE)
    if (length(fs) == 0) return(NULL)
    map_dfr(fs, function(f) read_tsv(f, col_types = cols(comparison = "c", sample_id = "c", group = "c", .default = "d"))) %>%
        mutate(cmp_n = as.integer(str_extract(comparison, "^\\d+")),
               cmp = str_remove(comparison, "^\\d+_"),
               condition = ifelse(group == "TREATMENT", str_extract(cmp, ".+(?=_vs_)"), str_extract(cmp, "(?<=_vs_).+")),
               group_label = str_glue("{condition}_{group}"),
               sample_label = str_glue("{sample_id}_{group}")) %>%
        arrange(cmp_n, desc(group), sample_id)
}

## per-comparison block: sample rows, then "all", group averages and the
## Treatment/Control ratio (same layout as the diffReps sheet)
factor_block <- function(df, value_cols, digits, side_cols) {
    samples <- df %>% select(sample_id = sample_label, all_of(value_cols))
    all <- df %>% summarise(across(all_of(value_cols), ~round(mean(.x), digits))) %>% mutate(sample_id = "all")
    grp <- df %>% group_by(group, group_label) %>%
        summarise(across(all_of(value_cols), mean), .groups = "drop") %>% arrange(group_label)
    by_groups <- grp %>% mutate(across(all_of(value_cols), ~round(.x, digits)), sample_id = str_glue("Group_{group_label}")) %>%
        select(sample_id, all_of(value_cols))
    ratio <- map_dfc(value_cols, function(v) tibble(!!v := round(grp[[v]][grp$group == "TREATMENT"] / grp[[v]][grp$group == "CONTROL"], 2))) %>%
        mutate(sample_id = "Treatment/Control ratio")
    res <- bind_rows(samples %>% mutate(across(all_of(value_cols), ~round(.x, digits))), all, by_groups, ratio) %>%
        mutate(sample_id = as.character(sample_id)) %>% select(sample_id, everything())
    list(name = df$cmp[1], df = res,
         side = df %>% select(sample_id = sample_label, all_of(side_cols)) %>% mutate(sample_id = as.character(sample_id)),
         nrow = nrow(res), ncol = ncol(res))
}

deseq2_tbl <- read_factor_tables(argv$deseq2_factors)
csaw_tbl <- read_factor_tables(argv$csaw_factors) %>%
    {if (is.null(.)) NULL else group_by(., comparison) %>%
         mutate(effective_lib_size = round(effective_lib_size),
                effective_lib_relative = effective_lib_size / mean(effective_lib_size)) %>% ungroup()}

deseq2_blocks <- if (is.null(deseq2_tbl)) NULL else
    deseq2_tbl %>% group_split(cmp_n) %>% map(factor_block, value_cols = "size_factor", digits = 3, side_cols = "reads_in_peaks")
csaw_blocks <- if (is.null(csaw_tbl)) NULL else
    csaw_tbl %>% group_split(cmp_n) %>% map(factor_block, value_cols = c("norm_factor", "effective_lib_relative"), digits = 3,
                                            side_cols = c("total_reads", "effective_lib_size"))

###
### All methods side by side: every factor divided by its average in the
### comparison, so all columns are "counts are divided by this, average 1"
###
rel <- function(x) x / mean(x)

all_methods <- list(
    tst %>%
        mutate(cmp = fname, col = str_glue("{norm_caller}_{window}")) %>%
        group_by(cmp, col) %>% mutate(value = rel(normfactors)) %>% ungroup() %>%
        select(cmp, sample_id, group, col, value),
    if (!is.null(deseq2_tbl)) deseq2_tbl %>% group_by(cmp) %>%
        mutate(value = rel(size_factor), col = "DESEQ2") %>% ungroup() %>% select(cmp, sample_id, group, col, value),
    if (!is.null(csaw_tbl)) csaw_tbl %>%
        mutate(value = effective_lib_relative, col = "CSAW") %>% select(cmp, sample_id, group, col, value)
) %>% compact() %>% bind_rows() %>%
    mutate(col = as.character(col))

## RIPPM is identical for every window size: keep one column of it
rippm_cols <- sort(unique(str_subset(all_methods$col, "^RIPPM_")))
if (length(rippm_cols) > 1) all_methods <- filter(all_methods, !col %in% rippm_cols[-1])
all_methods <- all_methods %>% mutate(col = ifelse(str_detect(col, "^RIPPM_"), "RIPPM", col))

cmp_order <- bind_rows(deseq2_tbl, csaw_tbl) %>% {if (nrow(.) == 0) tibble(cmp = character(), cmp_n = integer()) else distinct(., cmp, cmp_n)}
method_order <- c("RIPPM", sort(unique(str_subset(all_methods$col, "^DIFFREPS_"))), "DESEQ2", "CSAW")

all_blocks <- all_methods %>%
    left_join(cmp_order, by = "cmp") %>%
    mutate(condition = ifelse(group == "TREATMENT", str_extract(cmp, ".+(?=_vs_)"), str_extract(cmp, "(?<=_vs_).+")),
           group_label = str_glue("{condition}_{group}"),
           sample_label = str_glue("{sample_id}_{group}")) %>%
    pivot_wider(names_from = col, values_from = value) %>%
    arrange(cmp_n, cmp, desc(group), sample_id) %>%
    group_split(cmp_n, cmp) %>%
    map(function(df) {
        cols <- intersect(method_order, names(df))
        b <- factor_block(df, value_cols = cols, digits = 2, side_cols = character())
        ## rank agreement with RIPPM: 1 = samples in the same order
        if ("RIPPM" %in% cols) {
            agree <- map_dfc(cols, function(v) tibble(!!v := round(cor(df[[v]], df$RIPPM, method = "spearman"), 2))) %>%
                mutate(sample_id = "Same order as RIPPM (1 = identical)")
            b$df <- bind_rows(b$df, agree)
            b$nrow <- nrow(b$df)
        }
        b$side <- NULL
        b
    })

write_blocks <- function(wb, sheet, blocks) {
    addWorksheet(wb, sheetName = sheet)
    cur_row <- 1
    for (b in blocks) {
        writeData(wb, sheet = sheet, b$name, startRow = cur_row, startCol = 1, colNames = FALSE)
        writeData(wb, sheet = sheet, b$df, startRow = cur_row + 1, startCol = 1)
        if (!is.null(b$side) && ncol(b$side) > 1)
            writeData(wb, sheet = sheet, b$side, startRow = cur_row + 1, startCol = b$ncol + 2)
        cur_row <- 4 + cur_row + b$nrow
    }
    setColWidths(wb, sheet, cols = 1, widths = 38)
}

legend <- tribble(
    ~Sheet, ~Item, ~Meaning,
    "All sheets", "How factors are used", "In every sheet a sample's read counts are DIVIDED by its factor. A larger factor means the sample had more reads (or more signal) and is scaled down more.",
    "All sheets", "Block", "One block per comparison: the comparison name, a table with one row per sample, then summary rows.",
    "All sheets", "<sample>_TREATMENT / _CONTROL", "Sample ID and its group. TREATMENT is the first group in the comparison name, CONTROL the second (e.g. Male_2wk vs Female_2wk: Male = TREATMENT).",
    "All sheets", "all", "Average of the column over all samples in the comparison.",
    "All sheets", "Group_<name>", "Average of the column within one group.",
    "All sheets", "Treatment/Control ratio", "Group average of TREATMENT divided by that of CONTROL. 1 = both groups scaled the same; above 1 = TREATMENT samples are scaled down more.",
    "diffReps normalization", "DIFFREPS_<window>", "Factor diffReps computes itself from read counts in windows of this size (bp).",
    "diffReps normalization", "RIPPM_<window>", "Reads in peaks of the sample divided by the average over the samples of the comparison (average = 1). The same for every window size.",
    "diffReps normalization", "num_raw_reads", "Reads in the raw sequencing file (R1).",
    "diffReps normalization", "fragment_count", "Fragments after alignment and filtering.",
    "diffReps normalization", "fragment_in_peak_count", "Fragments inside peaks (the reads RIPPM is based on).",
    "diffReps normalization", "fragment_in_peak_ratio", "Share of fragments inside peaks (signal-to-noise; FRiP).",
    "DESeq2 normalization", "size_factor", "DESeq2 factor computed from the reads in the comparison's MuMerge peaks; a typical sample is about 1.",
    "DESeq2 normalization", "reads_in_peaks", "Reads counted in the comparison's MuMerge peaks.",
    "CSAW normalization", "total_reads", "Reads (read pairs for paired-end data) counted on the main chromosomes.",
    "CSAW normalization", "norm_factor", "Correction for signal-to-noise, computed from the sites with signal. 1 = no correction; below 1 = the sample's peaks look taller than its depth explains.",
    "CSAW normalization", "effective_lib_size", "total_reads x norm_factor: the number the counts are divided by.",
    "CSAW normalization", "effective_lib_relative", "effective_lib_size divided by its average in the comparison (average = 1); directly comparable with RIPPM.",
    "All methods", "Columns", "Each method's factor divided by its average in the comparison (average = 1), so all columns are on the same scale. If two methods treat the samples alike, their columns look alike.",
    "All methods", "Same order as RIPPM", "Rank correlation with RIPPM: 1 = the method ranks the samples in exactly the same order, 0 = unrelated."
)

wb <- createWorkbook()
addWorksheet(wb, sheetName = "Legend")
writeData(wb, sheet = "Legend", legend)
setColWidths(wb, "Legend", cols = 1:3, widths = c(24, 32, 110))
addStyle(wb, "Legend", createStyle(textDecoration = "bold"), rows = 1, cols = 1:3)
addStyle(wb, "Legend", createStyle(wrapText = TRUE, valign = "top"), rows = 2:(nrow(legend) + 1), cols = 1:3, gridExpand = TRUE)

sheet_name <- str_glue("diffReps normalization")
addWorksheet(wb, sheetName = sheet_name)

cur_row <- 1
cur_col <- 1
for(l in tt1){
    writeData(wb, sheet = sheet_name, l$name, startRow = cur_row, startCol = cur_col, colNames = FALSE)
    writeData(wb, sheet = sheet_name, l$df, startRow = cur_row+1, startCol = cur_col)
    writeData(wb, sheet = sheet_name, l$rippm_df, startRow = cur_row+1, startCol = l$ncol+2)
    cur_row <- 4+cur_row+l$nrow
}

if (!is.null(deseq2_blocks)) write_blocks(wb, "DESeq2 normalization", deseq2_blocks)
if (!is.null(csaw_blocks)) write_blocks(wb, "CSAW normalization", csaw_blocks)
if (length(all_blocks) > 0) write_blocks(wb, "All methods", all_blocks)

saveWorkbook(wb, argv$output, overwrite = T)


