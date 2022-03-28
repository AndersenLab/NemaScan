
#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(glue)
library(purrr)
library(mediation)



# load arguments
args <- commandArgs(trailingOnly = TRUE)

# load gene
gene_expression = args[1]


# load genotype matrix
Genotype_Matrix <- readr::read_tsv(args[2]) #%>%
  # na.omit()

# load expression data
 
  expression_pheno_raw <- read.delim(args[3], stringsAsFactors=FALSE)
  
  expression_pheno <- expression_pheno_raw %>% 
    gather(trait2,value,-strain) %>% 
    dplyr::mutate(trait=sub("(^X)(.*)","\\2",trait2)) %>% 
    dplyr::select(strain,trait,value)
  
  
 


# load pheno data
trait_pheno <- read.delim(args[4], stringsAsFactors=FALSE)


# GWAS qtl infor
gwas_intchr = args[5]

gwas_peak = args[6] %>% as.numeric()

gwtrait = args[7]

# load eqtl data
eqtl_infor <- read.delim(args[8], stringsAsFactors=FALSE)

 


# mediation function 
eQTL_mediate_dQTL <- function(gwas_pchr, gwas_p, gene_exp, phenodf) {
  
  # get the genotype at the peak marker
  gwas_g <- Genotype_Matrix %>% 
    dplyr::filter(CHROM==gwas_pchr & POS == gwas_p) %>% 
    dplyr::select(-(1:4)) %>% 
    tidyr::gather(strain,geno) %>%
    na.omit()
  
  # get the expression data
  
  lgmtpm_gwas <- expression_pheno %>% 
    dplyr::filter(trait==gene_exp,
                  strain %in% gwas_g$strain) %>% 
    tidyr::spread(trait,value) %>% 
    dplyr::rename(expression=gene_exp) %>% 
    dplyr::mutate(scale_exp = (expression - mean(expression, na.rm = T)) / sd(expression, na.rm = T)) %>% 
    dplyr::select(strain,scale_exp)
  
  
  # join data
  pheno <- left_join(phenodf,lgmtpm_gwas) %>% left_join(gwas_g) %>% na.omit() %>%
    dplyr::select(strain,phenotype = tr,expression=scale_exp,geno)
  
  # mediation lm 
  model.m <- lm(expression ~ geno, data = pheno)
  model.y <- lm(phenotype ~ expression + geno, data = pheno)
  out <- mediation::mediate(model.m, model.y, sims = 1000, boot = T, treat = "geno", mediator = "expression")
  
  return(out)
}


# function to get the summary statistics from model
summarize_model <- function(model) {
  # causal mediation effect
  acme <- data.frame(var = "ACME", 
                     estimate = model$d0, 
                     ci_lower = model$d0.ci[[1]], 
                     ci_upper = model$d0.ci[[2]],
                     prob = model$d0.p)
  # direct effect
  ade <- data.frame(var = "ADE", 
                    estimate = model$z0, 
                    ci_lower = model$z0.ci[[1]], 
                    ci_upper = model$z0.ci[[2]],
                    prob = model$z0.p)
  # total effect
  total <- data.frame(var = "total", 
                      estimate = model$tau.coef, 
                      ci_lower = model$tau.ci[[1]], 
                      ci_upper = model$tau.ci[[2]],
                      prob = model$tau.p)
  # prop. mediated
  med <- data.frame(var = "MED", 
                    estimate = model$n0, 
                    ci_lower = model$n0.ci[[1]], 
                    ci_upper = model$n0.ci[[2]],
                    prob = model$n0.p)
  # make a dataframe
  df <- rbind(acme, ade, total, med)
  return(df)
}




##


model <- eQTL_mediate_dQTL(gwas_pchr=gwas_intchr, gwas_p=gwas_peak, gene_exp=gene_expression, phenodf=trait_pheno)


df <- summarize_model(model) %>%
  dplyr::mutate(var = dplyr::case_when(var == "ADE" ~ "direct",
                                       var == "MED" ~ "prop_med", 
                                       var == "ACME" ~ "med",
                                       var == "total" ~ "total",
                                       TRUE ~ "NA"))

df2 <- df %>% 
  dplyr::mutate(gene=gene_expression, gwas_qtl=paste(gwtrait,gwas_intchr,gwas_peak,sep = "_")) %>% 
  dplyr::left_join(eqtl_infor,by=c("gene"="trait"))


# save mapping data set
readr::write_tsv(df2, 
                 path = glue::glue("{gwtrait}_{gwas_intchr}_{gwas_peak}_{gene_expression}_med_{args[9]}.tsv"),
                 col_names = T)
