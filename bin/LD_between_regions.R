library(genetics) 
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
gm <- read.table(file = args[1], header = T)
processed_mapping <- read.delim(args[2], stringsAsFactors=FALSE)
TRAIT <- args[3]

snp_df <- processed_mapping %>% na.omit()

ld_snps <- dplyr::filter(gm, CHROM %in% snp_df$CHROM, POS %in% snp_df$POS)


if ( nrow(ld_snps) > 1 ) {
  
  ld_snps <- data.frame(snp_id = paste(ld_snps$CHROM, ld_snps$POS,
                                       sep = "_"), data.frame(ld_snps)[, 5:ncol(ld_snps)])
  
  sn <- list()
  
  for (i in 1:nrow(ld_snps)) {
    sn[[i]] <- genetics::genotype(as.character(gsub(1, "T/T",
                                                    gsub(-1, "A/A", ld_snps[i, 4:ncol(ld_snps)]))))
  }
  
  test <- data.frame(sn)
  colnames(test) <- (ld_snps$snp_id)
  ldcalc <- t(genetics::LD(test)[[4]])^2
  diag(ldcalc) <- 1
  
  write.table(ldcalc, paste0(TRAIT, "_LD_between_QTL_regions_", args[4], ".tsv"), quote=F, row.names = T, col.names = NA, sep="\t")
  
  ldcalc %>%
    as.data.frame() %>%
    dplyr::mutate(QTL1 = rownames(.),
                  trait = TRAIT) %>%
    tidyr::pivot_longer(cols = -c(QTL1, trait), names_to = "QTL2", values_to = "r2") %>%
    dplyr::filter(!is.na(r2)) %>%
    dplyr::select(QTL1, QTL2, everything()) %>%
    ggplot(., mapping = aes(x = QTL1, y = QTL2)) + 
    theme_classic() +
    geom_tile(aes(fill = r2),colour = "black", size = 3) + 
    geom_text(aes(label = round(r2, 4))) + 
    scale_fill_gradient(low="darkgreen", high="red", limits = c(0, 1), name = expression(r^2)) + 
    theme(axis.title = element_blank(),
          axis.text = element_text(colour = "black")) + 
    labs(title = paste0("Linkage Disequilibrium: ",TRAIT))
  
    ggsave(filename = paste0(TRAIT,"_LD.plot_", args[4], ".png"), width = 7, height = 7)
}