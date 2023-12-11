library(tidyverse)
library(pbapply)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
#proc_sim_dir <- args[1]
proc_sim_dir <- "/projects/b1059/projects/Ryan/Caenorhabditis_GWAS/best_panel_subsample/20231017_CE_96.192_Alloutlier/20231204_ce96_processed/scored_sims"
#get the folder one level up
out_dir <- dirname(proc_sim_dir)

#define the date format year-mm-dd-ss
date <- format(Sys.time(), "%Y%m%d-%H%M%S")
#date <- format(Sys.Date(), "%Y%m%d")

#Function to read-in simulations
read_scored_sims <- function(file){
  data.table::fread(
    file,
    header = TRUE,
    sep = "\t",
    col.names = c(
        "QTL",
        "Simulated",
        "Detected",
        "CHROM",
        "POS",
        "RefAllele",
        "Frequency",
        "Effect",
        "Simulated.QTL.VarExp",
        "log10p",
        "aboveBF",
        "startPOS",
        "peakPOS",
        "endPOS",
        "detected.peak",
        "interval.Frequency",
        "BETA",
        "interval.log10p",
        "peak_id",
        "interval_size",
        "interval.var.exp",
        "top.hit",
        "nQTL",
        "simREP",
        "h2",
        "maf",
        "effect_distribution",
        "strain_set_id",
        "algorithm_id"
        )
    )
}


sim_files <- list.files(path = proc_sim_dir, pattern = "*.tsv", full.names = TRUE)
#print the number of files
n_simfiles <- length(sim_files)
print(
    glue::glue("Processing the following number of files: {n_simfiles}")
    )

#print head of the first file

# Read in all the files and bind them together
all_sims <- pblapply(sim_files, read_scored_sims) %>% 
  bind_rows()

print(
    glue::glue("Writing output file in {out_dir}")
    )

# Write out the file
write.table(
    all_sims,
    file = glue::glue("{out_dir}/{date}_aggregate_sims.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
