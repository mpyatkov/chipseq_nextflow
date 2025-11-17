#!/usr/bin/env Rscript

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install packages only if not installed
packages <- c("DESeq2", "plyranges")
for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
        #library(pkg, character.only = TRUE)
    }
}

# Compact version for quick scripts
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

suppressPackageStartupMessages(library(pak))
pak::pkg_install(c("argparser", "writexl", "openxlsx2", "formattable", "ggfortify", "MAnorm2"))

