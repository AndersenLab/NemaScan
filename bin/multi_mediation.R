#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(tibble)
library(MultiMed)

# load arguments
args <- commandArgs(trailingOnly = TRUE)


# load genotype matrix
Genotype_Matrix <- readr::read_tsv(args[1]) #%>%
    # na.omit()

# load pheno data
trait_phenotype <- read.delim(args[3], stringsAsFactors=FALSE) 


# load eqtl data
eqtl_infor <- read.delim(args[7], stringsAsFactors=FALSE)

transcript_list <- eqtl_infor %>% 
  dplyr::select(gwtrait,e_chr,gwpeak,trait) %>% 
  dplyr::distinct()



# transcript level

texpression_pheno_raw <- data.table::fread(args[2])

texpression_pheno <- texpression_pheno_raw %>% 
  tidyr::gather(trait2,value,-strain) %>% 
  dplyr::mutate(trait=sub("(^X)(.*)","\\2",trait2)) %>% 
  dplyr::select(strain,trait,value) %>% 
  dplyr::filter(strain %in% trait_phenotype$strain) %>% 
  dplyr::filter(trait %in% transcript_list$trait) %>% 
  na.omit() 



# processed pheno data
trait_pheno_all <- trait_phenotype%>% 
  dplyr::rename(trait=tr) %>% 
  dplyr::filter(strain %in% texpression_pheno$strain)




# GWAS qtl infor
gwas_intchr = args[4]

gwas_peak = args[5] %>% as.numeric()

gwtrait = args[6]





# get the genotype at the peak marker
gwas_g_all <- Genotype_Matrix %>% 
  dplyr::filter(CHROM==gwas_intchr & POS == gwas_peak) %>% 
  dplyr::select(-(1:4)) %>% 
  tidyr::gather(strain,geno) %>% 
  dplyr::filter(strain %in% texpression_pheno$strain) %>% 
  dplyr::arrange(strain) %>%
  na.omit()




multimed_trait_list=list()

for(trss in unique(texpression_pheno$trait)) {
  
  
  
  t_lgmtpm_gwas_all <- texpression_pheno %>% 
    dplyr::filter(trait==trss,
                  strain %in% gwas_g_all$strain)
  
  
  
  t_lgmtpm_gwas <- t_lgmtpm_gwas_all %>% 
    tidyr::spread(trait,value) %>% 
    dplyr::arrange(strain)  %>% 
    dplyr::select(-strain) 
  
  
  
  
  gwas_g <- gwas_g_all %>% 
    dplyr::filter(strain %in% t_lgmtpm_gwas_all$strain)
  
  
  trait_pheno <- trait_pheno_all%>% 
    dplyr::filter(strain %in% t_lgmtpm_gwas_all$strain)
  
  
  
  if( length(unique(gwas_g$geno)) == 2 & length(unique(t_lgmtpm_gwas_all$value)) > 1 ){
    
    #pick transcripts with variation in expression
    t_lgmtpm_gwas_vari <- t_lgmtpm_gwas[vapply(t_lgmtpm_gwas, function(x) length(unique(x)) > 1, logical(1L))]
    
    
    exp_matr_transcript <- as.matrix(t_lgmtpm_gwas_vari)
    
    
    mt_multi_transcript <- medTest(gwas_g$geno, exp_matr_transcript, trait_pheno$trait, nperm = 1000)
    
    df_multi_transcript <- as.data.frame(mt_multi_transcript) 
    
    row.names(df_multi_transcript) <- (colnames(exp_matr_transcript))
    
    
    df_multi_transcript2 <- df_multi_transcript %>% 
      tibble::rownames_to_column(var="gene") 
    
    
    multimed_trait_list[[trss]] <- df_multi_transcript2
    
    
  }
  
}

df_multi_transcript2 <- dplyr::bind_rows(multimed_trait_list)  %>% 
  dplyr::arrange(p)


df_multi_med <- df_multi_transcript2 %>% 
  dplyr::left_join(eqtl_infor,by=c("gene"="trait"))


q99_S = quantile(df_multi_med$S, probs = 0.99)[[1]]


gene_qtl_genelist <- df_multi_med %>% 
  dplyr::filter(p<0.05 | S > q99_S) %>% 
  dplyr::select(gwtrait,e_chr,gwpeak,trait=gene) %>% 
  dplyr::distinct()


if(nrow(gene_qtl_genelist)>0){
  # save mapping data set
  readr::write_tsv(df_multi_med, 
                   path = glue::glue("{gwtrait}_{gwas_intchr}_{gwas_peak}_medmulti_{args[8]}.tsv"),
                   col_names = T)
  
  gene_qtl_genelist <- gene_qtl_genelist %>%
      dplyr::mutate(algorithm = args[8]) %>%
      dplyr::select(gwtrait, e_chr, gwpeak, algorithm, trait)
  # save gene list 
  readr::write_tsv(gene_qtl_genelist, 
                   path = glue::glue("{gwtrait}_{gwas_intchr}_{gwas_peak}_elist_{args[8]}.tsv"),
                   col_names = F)
  
}



