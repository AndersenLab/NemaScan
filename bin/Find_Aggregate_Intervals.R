#!/usr/bin/env Rscript
library(tidyverse)
library(ggbeeswarm)

# argument information
# 1 - Genetoype matrix
# 2 - Phenotype data
# 3 - Mapping data
# 4 - independent tests - eigen
# 5 - number of qtl simulated (numeric)
# 6 - simulation replicate (numeric)
# 7 - If two QTL are less than this distance from each other, combine the QTL into one, (DEFAULT = 1000)
# 8 - Number of SNVs to the left and right of the peak marker used to define the QTL confidence interval, (DEFAULT = 150)
# 9 - simulated heritability2
# 10 - minor allele frequency (numeric)
# 11 - - String, What significance threshold to use for defining QTL,
#       BF = Bonferroni,
#       EIGEN = Defined by Number of independent tests from Eigen Decomposition of SNV correlation matrix,
#       or a user defined number
# 12 - strain set name (string)
# 13 - MAF (numeric)
# 14 - effect size range (character)
# 15 - mapping label

# load arguments
args <- commandArgs(trailingOnly = TRUE)
# args <- c("hahnel.isotypes_0.05_Genotype_Matrix.tsv",  ## TESTING
#           "1_1_0.9_0.05_0.5-5_hahnel.isotypes_sims.phen",
#           "temp.aggregate.mapping.tsv",
#           "hahnel.isotypes_0.05_total_independent_tests.txt", 1, 1, 1000, 150, 0.9, 0.05,
#           "BF", "hahnel.isotypes", 0.05, "0.5-5", "aggregate")

# define the trait name
trait_name <- glue::glue("{args[5]}_{args[6]}_{args[9]}")

# load phenotpe data
phenotype_data <- data.table::fread(args[2], col.names = c("strain", "strain_drop", trait_name)) %>%
  na.omit() %>%
  dplyr::select(-strain_drop)%>%
  as.data.frame()

# load GCTA mapping data
map_df <- data.table::fread(args[3]) %>%
  dplyr::rename(marker = SNP, CHROM = CHR) %>%
  dplyr::mutate(log10p = -log10(P))

# load genotype matrix
genotype_matrix <- readr::read_tsv(args[1]) %>%
  na.omit()

# define method for setting significance threshold
significance_threshold <- args[11]

# set significance threshold
if(significance_threshold == "EIGEN"){
  QTL_cutoff <- data.table::fread(args[4]) %>% dplyr::pull(V1) # Independent tests should be fed to calculation
} else if(significance_threshold == "BF"){
  QTL_cutoff <- NA
} else {
  QTL_cutoff <- as.numeric(args[11]) # Specified Threshold
}


# process mapping function
process_mapping_df <- function (mapping_df,
                                phenotype_df,
                                CI_size = as.numeric(args[8]),
                                snp_grouping = as.numeric(args[7]),
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
  } else if (is.numeric(QTL_cutoff) & thresh == "EIGEN"){
    mapping_df <- mapping_df %>% 
      dplyr::mutate(trait = colnames(phenotype_df)[2]) %>%
      dplyr::group_by(trait) %>% 
      dplyr::filter(log10p != 0) %>% 
      dplyr::mutate(BF = -log10(0.05/BF)) %>% 
      dplyr::mutate(aboveBF = ifelse(log10p >= BF, 1, 0))
  } else {
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

    pheno$trait <- gsub("-", "\\.", pheno$trait)

    rawTr <- pheno %>%
      dplyr::left_join(., snpsForVE, by = "trait")

    rawTr$marker <- as.character(rawTr$marker)
    rawTr$strain <- as.character(rawTr$strain)

    snp_df <- geno %>%
      dplyr::select(-REF, -ALT)

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
        dplyr::select(CHROM, POS, index) %>%
        dplyr::filter(POS == min(POS) | POS == max(POS))

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
      }
      else {
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
      dplyr::mutate(interval_size = endPOS - startPOS)


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
                                        CI_size = as.numeric(args[8]),
                                        snp_grouping = as.numeric(args[7]),
                                        BF = QTL_cutoff,
                                        thresh = significance_threshold,
                                        geno = genotype_matrix)

# save processed mapping data

readr::write_tsv(processed_mapping,
                 path = glue::glue("{trait_name}_{args[13]}_{args[14]}_{args[12]}_processed_{args[15]}_mapping.tsv"),
                 col_names = T)

# extract interval information
qtl_region <- processed_mapping %>%
  na.omit() %>%
  dplyr::distinct(CHROM, marker, trait, startPOS,	peakPOS,	endPOS, peak_id)

# save processed mapping data
readr::write_tsv(qtl_region,
                 path = glue::glue("{trait_name}_{args[13]}_{args[14]}_{args[12]}_{args[15]}_qtl_region.tsv"),
                 col_names = T)


## LD ###
# interval_pos_df_LD <- interval_pos_df %>%
#   dplyr::mutate(marker = paste(CHROM, POS, sep = ":"))
# 
# if(nrow(interval_pos_df_LD) == 1){
#   
#   marker.LD <- data.frame(interval_pos_df_LD$marker, interval_pos_df_LD$marker, NA, unique(interval_pos_df_LD$trait)) %>%
#     `colnames<-`(c("marker1","marker2","r2","trait"))
#   
# } else {
#   QTLcombos <- data.frame(t(combn(x = unique(gINFO$marker), m = 2))) %>%
#     `colnames<-`(c("marker1","marker2"))
#   LD <- list()
#   
#   for(q in 1:nrow(QTLcombos)){
#     markers.of.interest <- c(as.character(QTLcombos[q,]$marker1),
#                              as.character(QTLcombos[q,]$marker2))
#     haps <- suppressMessages(gINFO %>%
#                                dplyr::select(strain, allele, marker) %>%
#                                dplyr::filter(marker %in% markers.of.interest) %>%
#                                dplyr::mutate(allele = if_else(allele == -1, true = "REF", false = "ALT")) %>%
#                                dplyr::mutate(marker = as.factor(marker)) %>%
#                                tidyr::pivot_wider(names_from = marker, values_from = allele) %>%
#                                `colnames<-`(c("strain","A","B")) %>%
#                                tidyr::unite("hap", c(A,B), sep = "_", remove = F))
#     
#     P <- suppressMessages(haps %>%
#                             dplyr::group_by(hap) %>%
#                             dplyr::summarise(n()/nrow(haps)) %>%
#                             `colnames<-`(c("P","freq")) %>%
#                             tidyr::pivot_wider(names_from = P, values_from = freq))
#     
#     n <- suppressMessages(gINFO %>%
#                             dplyr::select(marker) %>%
#                             dplyr::mutate(marker = as.factor(marker)) %>%
#                             dplyr::group_by(marker) %>%
#                             dplyr::summarise(n()) %>%
#                             `colnames<-`(c("marker_id","total")))
#     
#     p <- suppressMessages(gINFO %>%
#                             dplyr::select(strain, allele, marker) %>%
#                             dplyr::filter(marker %in% markers.of.interest) %>%
#                             dplyr::mutate(allele = if_else(allele == -1, true = "REF", false = "ALT")) %>%
#                             dplyr::mutate(marker = as.factor(marker)) %>%
#                             dplyr::group_by(allele, marker) %>%
#                             dplyr::summarise(n()) %>%
#                             `colnames<-`(c("p","marker_id","n")) %>%
#                             dplyr::left_join(.,n) %>%
#                             dplyr::mutate(freq = n/total) %>%
#                             tidyr::unite("allele", c(p,marker_id), sep = "_") %>%
#                             dplyr::ungroup() %>%
#                             dplyr::select(allele, freq) %>%
#                             tidyr::pivot_wider(names_from = allele, values_from = freq))
#     
#     pApB <- p %>%
#       dplyr::select(contains("REF")) %>%
#       c(.) %>%
#       unlist() %>%
#       prod()
#     
#     pApBpapb <- p %>%
#       unlist() %>%
#       prod()
#     
#     D.AB <- P$REF_REF - pApB
#     r2 <- (D.AB^2)/pApBpapb
#     LD[[q]] <- gINFO %>%
#       dplyr::select(marker) %>%
#       dplyr::filter(marker %in% markers.of.interest) %>%
#       dplyr::distinct()  %>%
#       tidyr::pivot_wider(names_from = marker, values_from = marker) %>%
#       `colnames<-`(c("marker1","marker2")) %>%
#       dplyr::mutate(r2 = r2,
#                     trait = unique(gINFO$trait))
#   }
#   
#   marker.LD <- Reduce(rbind, LD)
# }


# # save mapping LD data
# readr::write_tsv(processed_mapping[[2]],
#                  file = glue::glue("{trait_name}_{args[13]}_{args[14]}_{args[12]}_{args[15]}_qtl_LD.tsv"),
#                  col_names = T)
