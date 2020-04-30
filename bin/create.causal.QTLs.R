library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
variants <- read.table(args[1]) %>%
  dplyr::filter(V1 != "MtDNA")
simulate_QTL_effects <- function(n_causal_variants){
  options(scipen = 999)
  effects <- rgamma(n = n_causal_variants, shape = 0.4, scale = 1.66)
  causal_variants <- sample(variants$V2, 
                            size = n_causal_variants, 
                            replace = F)
  direction <- sample(c(-1,1), 
                      size = n_causal_variants, 
                      replace = T)
  simulated.QTL.effects <- data.frame(causal_variants, effects*direction)
  write.table(simulated.QTL.effects, paste("causal.variants.sim",n_causal_variants,"txt",sep = "."), 
              quote = F, 
              col.names = F, 
              row.names = F)
}
n_QTLs <- c(2,10,50,100,250,500)
purrr::map(n_QTLs, simulate_QTL_effects)
