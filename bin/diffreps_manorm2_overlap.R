#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

library(readxl)
library(writexl)

remotes::install_cran("argparser", upgrade = "never")
library(argparser)

remotes::install_cran("plyranges", upgrade = never)
library(plyranges)
