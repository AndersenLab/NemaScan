#!/usr/bin/env Rscript
library(tidyverse)
# library(rrBLUP)
library(ggbeeswarm)
library(sommer)

# LOCAL
# setwd("~/Documents/projects/NemaScan_Performance/data/")
# setwd("~/Documents/AndersenLab/NemaScan_Performance/data/")
# argument information
# 1 - Genetoype matrix
# 2 - Phenotype data 
# 3 - Mapping data
# 4 - independent tests - eigen
# 5 - If two QTL are less than this distance from each other, combine the QTL into one, (DEFAULT = 1000)
# 6 - Number of SNVs to the left and right of the peak marker used to define the QTL confidence interval, (DEFAULT = 150)
# 7 - - String, What significance threshold to use for defining QTL, 
#       BF = Bonferroni, 
#       EIGEN = Defined by Number of independent tests from Eigen Decomposition of SNV correlation matrix, 
#       or a user defined number
# 8 - mapping label

# FOR LOCAL
# args <- c("Genotype_Matrix.tsv",
#           "pr_Emodepside_cv.EXT.tsv",
#           "temp.aggregate.mapping.tsv",
#           "total_independent_tests.txt",
#           1000, # formerly args[7], now args[5]
#           150, # formerly args[8], now args[6]
#           "EIGEN", # formerly args[11], now args[7]
#           "Emodepside_cv.EXT_AGGREGATE" #formerly args[15], now args[8]
#           )
# FOR QUEST
# load arguments
args <- commandArgs(trailingOnly = TRUE)


# load phenotpe data
phenotype_data <- data.table::fread(args[2]) %>%
  na.omit() %>%
  as.data.frame()

# load GCTA mapping data
map_df <- data.table::fread(args[3]) %>%
  dplyr::rename(marker = SNP, 
                CHROM = CHR,
                POS = POS) %>%
  dplyr::mutate(log10p = -log10(P))

# load genotype matrix
genotype_matrix <- readr::read_tsv(args[1]) %>%
  na.omit()

# define method for setting significance threshold
significance_threshold <- args[7]

if(significance_threshold == "EIGEN"){
  QTL_cutoff <- data.table::fread(args[4]) %>% dplyr::pull(V1) # Independent tests should be fed to calculation
} else if(significance_threshold == "BF"){
  QTL_cutoff <- NA
} else {
  QTL_cutoff <- as.numeric(args[7]) # Specified Threshold
}

# process mapping function
process_mapping_df <- function (mapping_df, 
                                phenotype_df, 
                                CI_size = as.numeric(args[6]),
                                snp_grouping = as.numeric(args[5]),
                                BF = NA,
                                thresh = significance_threshold,
                                geno = genotype_matrix) {
  pheno <- phenotype_df 
  
  pheno$trait <- colnames(phenotype_df)[2]
  
  colnames(pheno) <- c("strain", "value", "trait")
  
  # Determine how to make threshold
  if (is.na(QTL_cutoff)){ 
    mapping_df <- mapping_df %>% 
      dplyr::mutate(trait = colnames(phenotype_df)[2]) %>%
      dplyr::group_by(trait) %>% 
      dplyr::filter(log10p != 0) %>% 
      dplyr::mutate(BF = -log10(0.05/sum(log10p > 0, na.rm = T))) %>% 
      dplyr::mutate(aboveBF = ifelse(log10p >= BF, 1, 0))
  } 
  else if (is.numeric(QTL_cutoff) & thresh == "EIGEN"){
    mapping_df <- mapping_df %>% 
      dplyr::mutate(trait = colnames(phenotype_df)[2]) %>%
      dplyr::group_by(trait) %>% 
      dplyr::filter(log10p != 0) %>% 
      dplyr::mutate(BF = -log10(0.05/BF)) %>% 
      dplyr::mutate(aboveBF = ifelse(log10p >= BF, 1, 0))
  } 
  else {
    mapping_df <- mapping_df %>% 
      dplyr::mutate(trait = colnames(phenotype_df)[2]) %>%
      dplyr::group_by(trait) %>% 
      dplyr::filter(log10p != 0) %>% 
      dplyr::mutate(BF = BF) %>% 
      dplyr::mutate(aboveBF = ifelse(log10p >= BF, 1, 0))
  }
  
  
  
  Processed <- mapping_df %>% 
    dplyr::filter(sum(aboveBF, na.rm = T) > 0) %>% 
    dplyr::ungroup()
  
  snpsForVE <- Processed %>% 
    dplyr::filter(aboveBF == 1) %>% 
    dplyr::select(marker, trait)
  
  snpsForVE$trait <- as.character(snpsForVE$trait)
  
  if (nrow(snpsForVE) > nrow(Processed)*0.15) {
    Processed <- mapping_df %>%
      dplyr::mutate(strain = NA, value = NA, allele = NA, var.exp = NA,
                    startPOS = NA, peakPOS = NA, endPOS = NA,
                    peak_id = NA, interval_size = NA)
  } else if (nrow(snpsForVE) > 0 && nrow(snpsForVE) < nrow(Processed)*0.15) {
    
    row.names(pheno) <- gsub("-", "\\.", row.names(pheno))
    
    # pheno$trait <- gsub("-", "\\.", pheno$trait)  ## do we need this line? causes errors if trait name has hyphen...
    
    rawTr <- pheno %>%
      dplyr::left_join(., snpsForVE, by = "trait")
    
    rawTr$marker <- as.character(rawTr$marker)
    rawTr$strain <- as.character(rawTr$strain)
    
    snp_df <- geno %>% 
      dplyr::select(-REF, -ALT)
    snp_df$CHROM <- dplyr::recode(snp_df$CHROM, 
                    I = "1", 
                    II = "2",
                    III = "3",
                    IV = "4",
                    V = "5",
                    X = "6")
    
    gINFO <- snp_df %>% 
      dplyr::mutate(marker = paste(CHROM, POS, sep = ":")) %>%
      dplyr::filter(marker %in% snpsForVE$marker) %>% 
      tidyr::gather(strain, allele, -marker, -CHROM, -POS)
    
    gINFO$marker <- as.character(gINFO$marker)
    gINFO <- suppressWarnings(data.frame(gINFO) %>% 
                                dplyr::left_join(., snpsForVE, by = "marker") %>%
                                dplyr::left_join(rawTr, ., by = c("trait", "strain", "marker")))
    
    cors <- gINFO %>% dplyr::group_by(trait, marker) %>% 
      dplyr::mutate(var.exp = cor(value, allele, use = "pairwise.complete.obs", 
                                  method = "pearson")^2)
    
    CORmaps <- Processed %>%
      dplyr::mutate(CHROM = as.character(CHROM)) %>%
      dplyr::left_join(., cors, by = c("trait", "marker", "CHROM", "POS"), copy = TRUE)
    processed_mapping_df <- Processed
    correlation_df <- CORmaps
    phenotypes <- as.character(unique(processed_mapping_df$trait))
    intervals <- list()
    for (i in 1:length(phenotypes)) {
      PeakDF <- processed_mapping_df %>% 
        dplyr::filter(trait ==  phenotypes[i]) %>% 
        dplyr::group_by(CHROM, trait) %>% 
        dplyr::mutate(index = 1:n()) %>% 
        dplyr::mutate(peaks = cumsum(aboveBF)) %>% 
        dplyr::filter(aboveBF == 1) %>% 
        dplyr::group_by(CHROM, trait) %>% 
        dplyr::mutate(nBF = n()) %>% dplyr::group_by(CHROM, trait) %>% 
        dplyr::arrange(CHROM, POS)
      
      SNPindex <- processed_mapping_df %>% 
        dplyr::filter(trait == phenotypes[i]) %>% 
        dplyr::group_by(CHROM, trait) %>%
        dplyr::mutate(index = 1:n()) %>%
        dplyr::distinct(CHROM, POS, .keep_all = T) %>% 
        dplyr::select(CHROM, POS, index)
        # dplyr::filter(POS == min(POS) | POS == max(POS)) # is this necessary?
      
      findPks <- PeakDF %>% 
        dplyr::filter(trait == phenotypes[i]) %>% 
        dplyr::group_by(CHROM) %>% 
        dplyr::arrange(CHROM, POS)
      
      if (findPks$nBF == 1 & length(unique(findPks$CHROM)) ==  1) {
        findPks$pID <- 1
        findPks <- findPks %>% 
          dplyr::group_by(CHROM, pID, trait) %>%
          dplyr::mutate(start = min(index) - CI_size, end = max(index) + CI_size)
        
        for (k in 1:nrow(findPks)) {
          tSNPs <- SNPindex %>% dplyr::filter(CHROM == findPks$CHROM[k])
          if (findPks$start[k] < min(tSNPs$index)) {
            findPks$start[k] <- min(tSNPs$index)
          }
          if (findPks$end[k] > max(tSNPs$index)) {
            findPks$end[k] <- max(tSNPs$index)
          }
        }
        intervals[[i]] <- findPks %>% dplyr::ungroup()
      } else {
        findPks$pID <- 1
        for (j in 2:nrow(findPks)) {
          findPks$pID[j] <- ifelse(abs(findPks$index[j] - findPks$index[j - 1]) < snp_grouping & findPks$CHROM[j] == findPks$CHROM[j - 1],
                                   findPks$pID[j - 1],
                                   findPks$pID[j - 1] + 1)
        }
        findPks <- findPks %>% 
          dplyr::group_by(CHROM, pID, trait) %>% 
          dplyr::mutate(start = min(index) - CI_size, end = max(index) + CI_size)
        
        for (k in 1:nrow(findPks)) {
          tSNPs <- SNPindex %>% dplyr::filter(CHROM == findPks$CHROM[k])
          
          if (findPks$start[k] < min(tSNPs$index)) {
            findPks$start[k] <- min(tSNPs$index)
          }
          if (findPks$end[k] > max(tSNPs$index)) {
            findPks$end[k] <- max(tSNPs$index)
          }
        }
      }
      intervals[[i]] <- findPks %>% dplyr::ungroup()
    }
    intervalDF <- data.table::rbindlist(intervals)
    peak_df <- intervalDF
    peak_list <- intervals
    Pos_Index_Reference <- processed_mapping_df %>%
      dplyr::group_by(CHROM, trait) %>% 
      dplyr::mutate(index = 1:n()) %>% 
      dplyr::mutate(peaks = cumsum(aboveBF)) %>% 
      dplyr::select(trait, CHROM, POS, index) %>% 
      dplyr::filter(index %in% c(unique(peak_df$start), unique(peak_df$end))) %>% 
      dplyr::ungroup()
    
    Pos_Index_Reference$trait <- as.character(Pos_Index_Reference$trait)
    interval_positions <- list()
    for (i in 1:length(peak_list)) {
      print(paste(100 * signif(i/length(peak_list), 3), 
                  "%", sep = ""))
      peak_list[[i]]$trait <- as.character(peak_list[[i]]$trait)
      peak_list[[i]] <- peak_list[[i]] %>% dplyr::arrange(desc(log10p)) %>% 
        dplyr::distinct(pID, .keep_all = T)
      trait_i <- unique(peak_list[[i]]$trait)
      index_i <- c(peak_list[[i]]$start, peak_list[[i]]$end)
      CHROM_i <- peak_list[[i]]$CHROM
      PKpos <- data.frame(Pos_Index_Reference) %>% 
        dplyr::filter(trait == trait_i & index %in% index_i & CHROM %in% CHROM_i) %>% 
        dplyr::left_join(., peak_list[[i]], by = c("trait",  "CHROM")) %>% 
        dplyr::mutate(issues = ifelse(start ==  index.x | end == index.x, 1, 0)) %>% 
        dplyr::filter(issues != 0) %>% 
        dplyr::select(trait, CHROM, POS.x, POS.y, pID, log10p, index.x, index.y, start, end) %>% 
        dplyr::group_by(CHROM, pID) %>% 
        dplyr::mutate(startPOS = min(POS.x),  peakPOS = POS.y, endPOS = max(POS.x)) %>%
        dplyr::distinct(trait, CHROM, pID, peakPOS, .keep_all = T) %>% dplyr::select(trait, CHROM, POS = POS.y, startPOS, peakPOS, endPOS, peak_id = pID)
      interval_positions[[i]] <- PKpos
    }
    
    interval_pos_df <- data.frame(data.table::rbindlist(interval_positions)) %>% 
      dplyr::mutate(interval_size = endPOS - startPOS) %>%
      dplyr::mutate(CHROM = as.character(CHROM))
    
    Processed <- suppressWarnings(dplyr::left_join(correlation_df, 
                                                   interval_pos_df, by = c("trait", "CHROM", "POS"), 
                                                   copy = TRUE))
  } else {
    Processed <- mapping_df %>% 
      dplyr::mutate(strain = NA, value = NA, allele = NA, var.exp = NA, 
                    startPOS = NA, peakPOS = NA, endPOS = NA, 
                    peak_id = NA, interval_size = NA)
  }
  
  return(Processed)
}


# process mapping data, define QTL
processed_mapping <- process_mapping_df(mapping_df = map_df, 
                                        phenotype_df = phenotype_data, 
                                        CI_size = as.numeric(args[6]), 
                                        snp_grouping = as.numeric(args[5]), 
                                        BF = QTL_cutoff,
                                        thresh = significance_threshold,
                                        geno = genotype_matrix) %>%
                dplyr::mutate(CHROM = dplyr::case_when(CHROM == 1 ~ "I",
                                          CHROM == 2 ~ "II",
                                          CHROM == 3 ~ "III",
                                          CHROM == 4 ~ "IV",
                                          CHROM == 5 ~ "V",
                                          CHROM == 6 ~ "X",
                                          CHROM == 7 ~ "MtDNA",
                                          TRUE ~ "NA"),
                              marker = glue::glue("{CHROM}:{POS}"))

# save processed mapping data
readr::write_tsv(processed_mapping,
                 c(paste("processed",args[8],"mapping.tsv", sep = "_")),
                 col_names = T)

# add narrow-sense heritability point estimate
# narrow sense heritability with sommer::mmer (no bootstrap)
narrowh2 <- function(df_h){
    h2_res <- sommer::mmer(value ~ 1, random = ~sommer::vs(strain, Gu = A), data = df_h)
    h2 <- as.numeric(sommer::pin(h2_res, h2 ~ (V1) / (V1+V2))[[1]][1])
    return(h2)
}

# get trait name and rename trait column
traitname <- names(phenotype_data)[2]
names(phenotype_data) <- c("strain", "value")

# additive matrix - first filter by strain
A <- sommer::A.mat(t(genotype_matrix %>% dplyr::select(dplyr::one_of(phenotype_data$strain))))

result <- narrowh2(phenotype_data)



# extract interval information
qtl_region <- processed_mapping %>%
  na.omit() %>%
  dplyr::distinct(CHROM, marker, trait, startPOS,	peakPOS,	endPOS, peak_id) %>%
  dplyr::mutate(narrow_h2 = result)

# save processed mapping data
readr::write_tsv(qtl_region, 
                 c(paste(args[8],"qtl_region.tsv", sep = "_")),
                 col_names = T)
