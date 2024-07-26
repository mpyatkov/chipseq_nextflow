#!/usr/bin/env Rscript

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("plyranges", update = FALSE)

remotes::install_cran("argparser", upgrade = "never")
remotes::install_cran("writexl", upgrade = "never")
remotes::install_cran("openxlsx2", upgrade = "never")
remotes::install_cran("formattable", upgrade = "never")
remotes::install_cran("ggfortify", upgrade = "never")
remotes::install_cran("MAnorm2", upgrade = "never")
remotes::install_cran("openxlsx2", upgrade = "never")
