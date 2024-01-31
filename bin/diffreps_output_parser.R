#!/usr/bin/env Rscript

library(argparser)

ParseArguments <- function() {
    p <- arg_parser('Parsing diffReps output')
    p <- add_argument(p,'--pattern', default="diffReps_", help="search pattern for diffReps output files")
    p <- add_argument(p, '--path', default="../", help="path to the dir which contains diffReps output directories")
    p <- add_argument(p, '--rippm_report', default="Norm_Factors.tsv", 
                      help="path to the file with rippm stats")
    p <- add_argument(p, '--output', default="Summary_normalization_factors.xlsx", help="output filename")
    return(parse_args(p))
}

argv <- ParseArguments()

library(stringr)
library(purrr)
library(readr)
library(dplyr)
library(tidyr)
library(openxlsx)

rippm_report <- read_tsv(argv$rippm_report, col_names = T) %>%
             select(-norm_factor)

files <- list.files(path = argv$path, pattern = argv$pattern, recursive = T, full.names = T) %>% 
    keep(function(x){str_detect(x,"\\.vs\\.") && str_detect(x,"annotated|hotspot", negate = T)})

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
  conditions <- basename(f) %>% str_replace("diffReps_","") %>% str_split(".vs.|_RIPPM|_DIFFREPS", simplify = T)
  cond_tr <- conditions[1]
  cond_ctrl <- conditions[2]
  tr_files <- extract_tab_separated(header, "Treatment files") %>% str_replace("_frag[[:print:]]+", "")
  ctrl_files <- extract_tab_separated(header, "Control files") %>% str_replace("_frag[[:print:]]+", "")
  window <- extract_tab_separated(header, "Window size")
  norm_tr <- extract_tab_separated(header, "Treatment normalization constant")
  norm_ctrl <- extract_tab_separated(header, "Control normalization constant")
  tibble(
      fname = basename(f) %>% str_split(.,"_RIPPM|_DIFFREPS",simplify = T) %>% .[1],
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
  mutate(window = str_glue("{norm_caller}_{window}"),
         sample_id = str_glue("{sample_id}_{group}"),
         group = str_glue("{condition}_{group}")) %>% 
  select(-norm_caller, -condition) %>% 
  group_by(fname) %>% 
  group_split()

lgroup_to_ltables <- function(df) {
  
    rippm_table <- df %>% select(sample_id, starts_with("fragment")) %>% distinct()

    norm_factors_table <- df %>% 
      select(-starts_with("fragment")) %>% 
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

sheet_name <- str_glue("diffReps normalization")
wb <- createWorkbook()
addWorksheet(wb, sheetName = sheet_name)

cur_row <- 1
cur_col <- 1
for(l in tt1){
    writeData(wb, sheet = sheet_name, l$name, startRow = cur_row, startCol = cur_col, colNames = FALSE)
    writeData(wb, sheet = sheet_name, l$df, startRow = cur_row+1, startCol = cur_col)
    writeData(wb, sheet = sheet_name, l$rippm_df, startRow = cur_row+1, startCol = l$ncol+2)
    cur_row <- 4+cur_row+l$nrow
}
saveWorkbook(wb, argv$output, overwrite = T)


