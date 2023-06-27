#!/usr/bin/env Rscript
require(tidyverse)
require(data.table)
require(GenomicRanges)

### Functions ###
var.exp <- function(data,QTL){
    #Function to calculate the variance explained by a single QTL
    aov.out <- summary(aov(data = data, trait.value ~ allele + 1))
    SST <- sum(aov.out[[1]][[2]])
    SSallele <- aov.out[[1]][[2]][[1]]
    simulated.variance.exp <- SSallele/SST
    data.frame(QTL, simulated.variance.exp)
} 




### Load in the data ###
effects_file <- "Analysis_Results-20230620/Simulations/OG0007525_OG0006567_OG0008498_OG0007923_OG0007069_c_elegans/Phenotypes/1_0.8_0.05_OG0007525_OG0006567_OG0008498_OG0007923_OG0007069_c_elegans_underground.gartersnake_sims.par"
phenos_file <- "Analysis_Results-20230620/Simulations/OG0007525_OG0006567_OG0008498_OG0007923_OG0007069_c_elegans/Phenotypes/1_0.8_0.05_OG0007525_OG0006567_OG0008498_OG0007923_OG0007069_c_elegans_underground.gartersnake_sims.phen"
gm_file <- "Analysis_Results-20230620/Genotype_Matrix/c_elegans_underground.gartersnake_0.05_Genotype_Matrix.tsv"
### Stepping through simulation metrics for a single simulation ###

# load simulated causal variants from .par files
effects <- data.table::fread(effects_file, header = T) %>% 
    tidyr::separate(QTL, c("CHROM","POS"), sep = ":", remove = F)

# load the simulated trait files from .phen files 
phenos <- data.table::fread(phenos_file, header = T) %>%
    `colnames<-`(c("strain","strain2","trait.value")) %>%
    dplyr::select(-strain2) #Remove redundant strain column

# load the genotype matrix
gm <- data.table::fread(gm_file, header = T) 

# create a complete.effects dataframe                   
complete.effects <- gm  %>% 
    tidyr::unite("QTL",c(CHROM, POS), sep = ":", remove = F) %>%
    dplyr::filter(QTL %in% effects$QTL) %>%  #Filter GM to include only simulated causal variants
    dplyr::select(-CHROM, -POS) %>% #Remove redundant columns used to create QTL column
    tidyr::pivot_longer(cols = !c(QTL, REF, ALT), #Create a long dataframe with one row per strain per QTL
                        names_to = "strain",
                        values_to = "allele")  %>% 
    dplyr::full_join(.,effects) %>%  #add the data from the simulated effects file (AF) - QTL	RefAllele	Frequency	Effect
    dplyr::full_join(.,phenos) %>%  #add the data from the simulated phenotypes file (trait.value) - trait.value
    dplyr::group_by(QTL) #group by QTL - simulated causal variants

# Check if any simulated causal variants are missing from the genotype matrix
# *Sam has some crazy function that does this by repeating the steps above and then comparing the two dataframes
# Im not sure why we would ever be missing a simulated CV from the genotype matrix
comp_effect_check <- complete.effects %>% 
    dplyr::summarise(n = n()) %>%  #count the number of strains with the simulated causal variant
    dplyr::filter(n != 1, #remove any QTL that do not have more than one strain with the causal variant
                !is.na(QTL))

if (nrow(comp_effect_check) != nrow(effects)) {
    print("Missing simulated causal variants from genotype matrix")
    print(setdiff(effects$QTL, comp_effect_check$QTL))
}

#If we pass this check then we can nest the data 
genos.effects <- nest(complete.effects)

#Calculate the simulated variance explained by each QTL
#Creates a df with one row per QTL and the phenotypic variance explained by that QTL
simQTL.variance.explained <- purrr::map2(genos.effects$data, 
                                        genos.effects$QTL, 
                                        var.exp) %>%
    Reduce(rbind,.)

# Append the simulated variance explained to the simulated causal vairants
effects <- effects %>%
    dplyr::full_join(.,simQTL.variance.explained) %>%
    dplyr::rename(Simulated.QTL.VarExp = simulated.variance.exp)

# Attempt to read in mapping .tsv files for each algo

safe.lmm.exact.inbred <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
                                                                                paste(x,"processed_LMM-EXACT-INBRED_mapping.tsv",sep = "_"), sep = "/"),header = T)},
                                        otherwise = "No successful mapping matching simulation parameters :( ")
safe.lmm.exact.inbred.pca <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
                                                                                paste(x,"processed_LMM-EXACT-INBRED_PCA_mapping.tsv",sep = "_"), sep = "/"),header = T)},
                                        otherwise = "No successful mapping matching simulation parameters :( ")   
safe.lmm.exact.loco <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
                                                                                paste(x,"processed_LMM-EXACT-LOCO_mapping.tsv",sep = "_"), sep = "/"),header = T)},
                                        otherwise = "No successful mapping matching simulation parameters :( ")
safe.lmm.exact.loco.pca <- purrr::safely(.f = function(x){data.table::fread(paste(paste(effect.range,nQTL,"Mappings", sep = "/"),
                                                                                paste(x,"processed_LMM-EXACT-LOCO_PCA_mapping.tsv",sep = "_"), sep = "/"),header = T)},
                                        otherwise = "No successful mapping matching simulation parameters :( ")

mapping.lmm.exact.inbred <- safe.lmm.exact.inbred(x)
mapping.lmm.exact.inbred.pca <- safe.lmm.exact.inbred.pca(x)
mapping.lmm.exact.loco <- safe.lmm.exact.loco(x)
mapping.lmm.exact.loco.pca <- safe.lmm.exact.loco.pca(x)

### ADD FUNCTION TO ITERATE THROUGH ALL MAPPING Algos ## - Only really need to do PCA here








