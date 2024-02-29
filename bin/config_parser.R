#!/usr/bin/env Rscript

library(readxl)
library(tidyr)
library(dplyr)
library(stringr)
library(readr)
library(argparser)

ParseArguments <- function() {
  p <- arg_parser('Parsing XLSX configuration file')
  p <- add_argument(p,'--input_xlsx', default="sample_labels.xslx", help="xlsx configuration file")
  return(parse_args(p))
}

argv <- ParseArguments()


#setwd("/projectnb/wax-dk/max/CTCF_Rad21_CHIPSEQ")

get_table <- function(path_to_file, sheet_number, marker, shift, num_cols,row_limit) {
  
  ## config to matrix for more convenient search of markers
  mx <- read_excel(path = path_to_file, sheet = sheet_number) %>% as.matrix()
  
  ## find table start, apply shift if marker column in the middle of the table
  ts <- which(mx==marker, arr.ind=TRUE) %>% 
    as_tibble(.) %>% 
    filter(row == min(row) & col == max(col)) %>% 
    mutate(row=row+1, col=col+shift) %>% 
    as.vector()
  
  ## left_upper and right_bottom corners of the table in excell letter:number notation
  corners <- paste0(LETTERS[ts$col], ts$row,":",LETTERS[ts$col+num_cols-1], row_limit)
  
  ## read {row_limit} number of rows starting by corners coords
  tmp_tab <- read_excel(path = path_to_file, sheet = sheet_number, range = corners, col_names = T) 
  
  ## detect end position
  end_pos <- which(is.na(tmp_tab), arr.ind=TRUE) %>% 
    as_tibble(.) %>% 
    filter(row == min(row) & col == max(col)) %>% 
    mutate(row=row-1) %>% as.vector()
  
  ## extract all rows until end position
  res_tab <- tmp_tab[1:end_pos$row,]
  res_tab
}

#### MAIN
### config processing
fastq <- get_table(argv$input_xlsx, sheet_number = 1, marker = "sample_id", shift = 0, num_cols = 3,row_limit = 100)
sample_labels <- get_table(argv$input_xlsx, sheet_number = 2, marker = "sample_id", shift = 0, num_cols = 3,row_limit = 50)
diffreps <- get_table(argv$input_xlsx, sheet_number = 2, marker = "exp_num", shift = 0, num_cols = 5,row_limit = 50)
groups <- get_table(argv$input_xlsx, sheet_number = 2, marker = "group_description", shift = -1, num_cols = 3,row_limit = 50)

### for each "auto" in groups config we need to assign specific color
all_colors<-c("204,0,204","204,77,0","0,127,204","146,166,196","66,162,66","194,86,2","40,150,160","255,86,0","0,153,153","130,130,210","80,80,160","20,150,88")

### remove colors which already in groups config
available_colors <- setdiff(all_colors, 
                            groups$color %>% purrr::discard(\(c) {c=="auto"})) ## drop "auto"

### assign next available color from all_colors list one by one 
final_groups <- groups %>% mutate(i=row_number(), 
                            color = case_when(color == "auto" ~ available_colors[i], TRUE ~ color)) %>% 
  select(-i)


### prepare diffreps configuration table which will be used in pipeline  
### separate (unnest) combined rows in diffreps config: A,B -> the separate rows for A and B 
t1 <- diffreps %>% 
  separate_rows(treatment_group, sep=",") %>% 
  separate_rows(control_group, sep=",") %>% 
  mutate(treatment_group = str_trim(treatment_group),
         control_group = str_trim(control_group)) %>% 
  pivot_longer(treatment_group:control_group,
               names_to = "group_type", ## control/treatment
               values_to = "group_id")  ## A,B,...

## collapse sampleids (G11, G22 sep.rows -> G11|G22)
t2 <- sample_labels %>% 
  select(sample_id, group_id) %>% 
  group_by(group_id) %>% 
  arrange(.by_group = T) %>% 
  summarise(samples = str_c(sample_id, collapse = "|")) %>% #, .by=group_id
  ungroup()

## collapse groups and make final configuration file for diffreps
final_diffreps_config <- t1 %>% left_join(t2, by = "group_id") %>% 
  left_join(groups %>% select(group_id, group_description)) %>% 
  group_by(exp_num, normalization, window, group_type) %>% 
  summarise(combined_samples = unique(samples) %>% str_c(collapse = "|"), 
            combined_group_desc = unique(group_description) %>% str_c(collapse = "_and_")) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "group_type", values_from = c("combined_samples", "combined_group_desc")) %>% 
  rename(control_samples = combined_samples_control_group,
         treatment_samples = combined_samples_treatment_group,
         control_name = combined_group_desc_control_group,
         treatment_name = combined_group_desc_treatment_group) %>% 
  select(exp_num,treatment_name,control_name,treatment_samples,control_samples,normalization, window)

  
### prepare sample_labels config
final_sample_labels <- sample_labels %>% 
  left_join(final_groups)

### limit fastq config only by samples presented in sample_labels
### even if index file is bigger we can calculate only subset of samples
fastq<-inner_join(fastq, final_sample_labels %>% select(sample_id))

write_csv(final_sample_labels, "sample_labels.csv", col_names = T)
write_csv(final_diffreps_config, "diffreps_config.csv", col_names = F)
write_csv(fastq, "fastq_config.csv", col_names = F)
