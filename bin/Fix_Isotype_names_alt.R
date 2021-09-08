#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

# arguments
# 1 - raw phenotype
# 2 - option to fix pheno? default true
# 3 - strain isotype lookup table

#' Resolve strain names to isotypes
resolve_isotypes <- function(...) {
    
    isotype_lookup = data.table::fread(args[3])
    strains <- unlist(list(...))
    
    purrr::map_chr(strains, function(x) {
        isotype <- isotype_lookup %>%
            dplyr::filter(
                (x == strain) |
                    (x == isotype) |
                    (x == previous_names)
            ) %>%
            dplyr::pull(isotype) %>%
            unique()
        
        if (length(isotype) == 0) {
            message(glue::glue("WARNING: {x} is not a known strain. Isotype set to NA; Please check CeNDR"))
        } else if (length(isotype) > 1) {
            message(glue::glue("WARNING: {x} resolves to multiple isotypes. Isotype set to NA; Please check CeNDR"))
        }
        if (length(isotype) != 1) {
            isotype <- NA
        }
        isotype
    })
    
}

#' Process phenotypes for mapping
process_strains <- function(df){
    
    if ( sum(grepl(colnames(df)[1], "Strain", ignore.case = T)) == 0 ) {
        message(glue::glue("WARNING: Check input data format, strain should be the first column."))
    }
    
    # ~ ~ ~ # resolve strain isotypes # ~ ~ ~ #
    # get strain isotypes
    strain_isotypes_db <- data.table::fread(args[3])
    # identify strains that were phenotyped, but are not part of an isotype
    non_isotype_strains <- dplyr::filter(df,
                                         !(strain %in% strain_isotypes_db$strain),
                                         !(strain %in% strain_isotypes_db$isotype))
    # remove any strains identified to not fall into an isotype
    if ( nrow(non_isotype_strains) > 0 ) {
        
        strains_to_remove <- unique(non_isotype_strains$strain)
        
        message(glue::glue("WARNING: Removing strain(s) {strains_to_remove} because they do not fall into a defined isotype."))
        
        df_non_isotypes_removed <- dplyr::filter( df, !( strain %in% strains_to_remove) )
    } else {
        df_non_isotypes_removed <- df
    }
    
    # ~ ~ ~ ## ~ ~ ~ ## ~ ~ ~ ## ~ ~ ~ # Resolve Isotypes # ~ ~ ~ ## ~ ~ ~ ## ~ ~ ~ ## ~ ~ ~ #
    df_isotypes_resolved <- df_non_isotypes_removed %>%
        dplyr::group_by(strain) %>%
        dplyr::mutate(isotype = resolve_isotypes(strain)) %>%
        dplyr::ungroup()
        # tidyr::gather(trait, phenotype, -strain, -isotype) %>%
        # dplyr::filter(!is.na(phenotype))
    
    # deal with multiple strains per isotype group
    test <- df_isotypes_resolved %>%
        dplyr::group_by(isotype) %>%
        dplyr::mutate(num = length(unique(strain))) %>% 
        dplyr::mutate(ref_strain = strain == isotype)
    no_issues <- test %>%
        dplyr::filter(num == 1 & ref_strain == T)
    issues <- test %>%
        dplyr::filter(num > 1 | ref_strain == F) 
    
    fixed_issues <- no_issues %>%
        dplyr::select(-num, -ref_strain)
    
    # go through each isotype issue and resolve it
    for(i in unique(issues$isotype)) {
        df <- issues %>%
            dplyr::filter(isotype == i)
        
        # if only one strain is phenotyped, just rename strain to isotype ref strain and flag
        if(length(unique(df$strain)) == 1) {
            fix <- df %>%
                dplyr::mutate(strain = isotype) %>%
                dplyr::select(-ref_strain, -num)
            message(glue::glue("WARNING: Non-isotype reference strain {df$strain[1]} renamed to isotype {i}."))
        } else {
            # remove non-isotype strains
            fix <- df %>%
                dplyr::filter(ref_strain) %>%
                dplyr::select(-ref_strain, -num)
            
            # warn the user
            if(sum(df$ref_strain) > 0) {
                message(glue::glue("WARNING: Non-isotype reference strain(s) {paste(df %>% dplyr::filter(!ref_strain) %>% dplyr::pull(strain) %>% unique(), collapse = ', ')} from isotype group {i} removed."))
            } 
            else {
                message(glue::glue("WARNING: Non-isotype reference strain(s) {paste(df %>% dplyr::filter(!ref_strain) %>% dplyr::pull(strain) %>% unique(), collapse = ', ')} from isotype group {i} removed. To include this isotype in the analysis, you can (1) phenotype {i} or (2) evaluate the similarity of these strains and choose one representative for the group."))
            }
            }
        # add to data
        fixed_issues <- dplyr::bind_rows(fixed_issues, fix)
        }

    return(fixed_issues %>% dplyr::select(strain = isotype))
}


######## process pheno

# load trait file
traits <- readr::read_tsv(args[1], col_names = F)

# fix col name
colnames(traits) <- "strain"

# print messages to file
sink("strain_issues.txt")
sink(stdout(), type = "message")
print("Strain issues: (if empty, no strain issues were found)")

# fix strain names
if(args[2] == "fix"){
    fixed_names <- process_strains(traits)
} else {
    fixed_names <- traits
}

write.table(fixed_names$strain, 
            file = glue::glue("Phenotyped_Strains.txt"),
            quote = F, col.names = F, row.names = F)

sink()

