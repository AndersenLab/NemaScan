#!/usr/bin/env Rscript
library(tidyverse)

# argument information
# 1 - wormbase build (format: WSXXX, XXX has to be greater than 245)
# 2 - species: (supported options: c_elegans, c_briggsae, c_tropicalis)

# # # example args
# args <- c("WS270", "tropicalis")

# load arguments
args <- commandArgs(trailingOnly = TRUE)

wbb <- args[1]
sname <- glue::glue("c_{args[2]}")

# check wormbase build format
wbb_vec <- strsplit(wbb,split = "")[[1]]

if(wbb_vec[1] == "W" & wbb_vec[2] == "S" & as.numeric(paste(wbb_vec[3:5], collapse = "")) >= 245 & length(wbb_vec) == 5){
  print(glue::glue("WormBase build is formatted correctly"))
  
  # check if the species is supported
  if(sname %in% c("c_elegans", "c_briggsae", "c_tropicalis")){
    print(glue::glue("Input species name is supported"))
    
    # print(glue::glue("Downloading {wbb} {sname} GFF3 file from WormBase"))
    # system(glue::glue("wget ftp://ftp.wormbase.org/pub/wormbase/releases/{wbb}/species/{sname}/PRJNA13758/{sname}.PRJNA13758.{wbb}.annotations.gff3.gz"))
    
    print(glue::glue("Downloading {wbb} {sname} GTF file from WormBase"))
    if(sname == "c_elegans"){
      prj_name <- "PRJNA13758"
    } else if (sname == "c_briggsae"){
      prj_name <- "PRJNA10731"
    } else if (sname == "c_tropicalis"){
      prj_name <- "PRJNA53597"
    } else {
      print("Species name is not in the proper format, please choose one of the following: c_elegans, c_briggsae, c_tropicalis")
    }
    
    system(glue::glue("wget ftp://ftp.wormbase.org/pub/wormbase/releases/{wbb}/species/{sname}/{prj_name}/{sname}.{prj_name}.{wbb}.canonical_geneset.gtf.gz"))
    
  } else {
    print("Species name is not in the proper format, please choose one of the following: c_elegans, c_briggsae, c_tropicalis")
  }
} else {
  print("WormBase build is not in the proper format, please have in this format: WSXXX, where XXX is a number greater than or equal to 245")
}

# convert GFF3 to refFlat format
# conversion script from: (mac) http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/gtfToGenePred
# conversion script from: (linux) http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred

get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  return(tolower(os))
}

print(glue::glue("Converting the {sname} GTF file to refFlat format"))
if(get_os() == "osx"){
  system(glue::glue("{args[3]}/gtfToGenePred_mac {sname}.{prj_name}.{wbb}.canonical_geneset.gtf.gz {sname}_{wbb}_refFlat.txt"))
} else if (get_os() == "linux") {
  system(glue::glue("{args[3]}/gtfToGenePred {sname}.{prj_name}.{wbb}.canonical_geneset.gtf.gz {sname}_{wbb}_refFlat.txt"))
} else {
  print("Your operating system is not supported")
}

# append WBGeneID names to refFlat file
print(glue::glue("Extracting gene names from the {sname} GTF file"))
if(sname == "c_elegans"){
  system(glue::glue("gunzip -c {sname}.{prj_name}.{wbb}.canonical_geneset.gtf.gz | awk '$3==\"transcript\" {{print}}' | cut -f1,2 -d\";\" | cut -f2,4 -d\" \" | sed 's/\"//g' | sed 's/;//g' | tr ' ' '\t' > gene_transcripts.tsv" ))
} else if (sname == "c_briggsae"){
  system(glue::glue("gunzip -c {sname}.{prj_name}.{wbb}.canonical_geneset.gtf.gz | awk '$3==\"transcript\" {{print}}' | cut -f1,2 -d\";\" | cut -f2,4 -d\" \" | sed 's/\"//g' | sed 's/;//g' | tr ' ' '\t' > gene_transcripts.tsv" ))
} else if (sname == "c_tropicalis"){
  system(glue::glue("gunzip -c {sname}.{prj_name}.{wbb}.canonical_geneset.gtf.gz | awk '$3==\"transcript\" {{print}}' | cut -f1,3 -d\";\" | cut -f2,4 -d\" \" | sed 's/\"//g' | sed 's/;//g' | tr ' ' '\t' > gene_transcripts.tsv" ))
} else {
  print("Species name is not in the proper format, please choose one of the following: c_elegans, c_briggsae, c_tropicalis")
}

gene_transcript_conversion <- data.table::fread("gene_transcripts.tsv", header = F) %>%
  dplyr::rename(transcript = V2)

refflat_file <- data.table::fread(glue::glue("{sname}_{wbb}_refFlat.txt")) %>%
  dplyr::rename(transcript = V1)%>%
  dplyr::left_join(gene_transcript_conversion, ., by = "transcript")

system(glue::glue("rm {sname}_{wbb}_refFlat.txt"))

write.table(refflat_file, file = glue::glue("{sname}_{wbb}_refFlat.txt"), quote = F, col.names = F, row.names = F, sep = "\t")

