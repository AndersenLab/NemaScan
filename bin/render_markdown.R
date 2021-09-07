#!/usr/bin/env Rscript

# argument information
# 1 - Nemascan report (trait) rmd file

# load arguments
args <- commandArgs(trailingOnly = TRUE)

markdown_rmd <- args[1]

rmarkdown::render(markdown_rmd)