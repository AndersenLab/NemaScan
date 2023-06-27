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


args = commandArgs(trailingOnly=TRUE)


### Load in the data ###
#effects_file <- "Analysis_Results-20230627/Simulations/c_elegans/{sim_id1}/Phenotypes/2_0.8_0.05_sim_id1_c_elegans_underground.gartersnake_sims.par"
effects_file <- args[1]

#phenos_file <- "Analysis_Results-20230627/Simulations/c_elegans/{sim_id1}/Phenotypes/2_0.8_0.05_sim_id1_c_elegans_underground.gartersnake_sims.phen"
phenos_file <- args[2]

#gm_file <- "Analysis_Results-20230627/Genotype_Matrix/c_elegans_underground.gartersnake_0.05_Genotype_Matrix.tsv"
gm_file <- args[3]

proc_mapping_file <- "/projects/b1059/projects/Ryan/ortholog_sims/NemaScan/Analysis_Results-20230627/Simulations/c_elegans/sim_id1/Mappings/c_elegans_2_0.8_processed_LMM-EXACT-INBRED_PCA_mapping.tsv"
proc_mapping_file <- args[4]

sp_id <- args[5]

sim_id <- args[6]

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

#load the processed mapping file
mapping_df <- data.table::fread(proc_mapping_file, header = T)

# Identify the mapping algorithm used from the processed mapping file
file_name <- basename(proc_mapping_file)

#Check if INBRED or LOCO in the file name
if(str_detect(file_name, "INBRED"))
    mapping_algorithm <- "INBRED"
else if(str_detect(file_name, "LOCO"))
    mapping_algorithm <- "LOCO"
else
    stop("Mapping algorithm not identified")


## Start processing Files ##

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

peak.info <- mapping_df %>%
            dplyr::filter(!is.na(peak_id)) %>%
            dplyr::select(CHROM, marker, POS, AF1, BETA, log10p, startPOS, peakPOS, endPOS, peak_id, interval_size, var.exp) %>%
            dplyr::filter(!duplicated(.)) %>%
            dplyr::mutate(detected.peak = marker)

simulated.mapping.results.scores <- mapping_df %>%
            dplyr::rename(QTL = marker) %>%
            dplyr::filter(QTL %in% effects$QTL) %>%
            dplyr::select(QTL, log10p, aboveBF)

effects.scores <- effects  %>%
            dplyr::full_join(., simulated.mapping.results.scores) %>%
            dplyr::filter(!duplicated(QTL),
                            !is.na(log10p))

#Define QTL that were detected by the mapping algorithm
peaks <- GenomicRanges::GRanges(seqnames = peak.info$CHROM,
                                ranges = IRanges::IRanges(start = peak.info$startPOS,
                                                        end = peak.info$endPOS),
                                peakPOS = peak.info$peakPOS,
                                detected.peak = peak.info$detected.peak)

#Define QTL containing the simulated causal variants
real.effects <- GenomicRanges::GRanges(seqnames = effects.scores$CHROM,
                                    ranges = IRanges::IRanges(start = as.numeric(effects.scores$POS),
                                                                end = as.numeric(effects.scores$POS)),
                                    QTL = effects.scores$QTL)

#Check if the simulated causal variants were detected by the mapping algorithm ??? 
overlap <- IRanges::findOverlapPairs(real.effects, peaks) %>%
    as.data.frame() %>%
    dplyr::select(first.QTL, second.X.start, second.peakPOS, second.X.end, second.detected.peak) %>%
    `colnames<-`(c("QTL","startPOS","peakPOS","endPOS","detected.peak")) %>%
    dplyr::right_join(., peak.info) %>%
    dplyr::mutate(QTL = if_else(is.na(QTL),
                                true = detected.peak,
                                false = QTL)) %>%
    dplyr::rename(interval.log10p = log10p,
                interval.var.exp = var.exp,
                interval.Frequency = AF1) %>%
    dplyr::select(-c(CHROM, marker, POS))

#Actually score the detected QTL by detection status
all.QTL <- data.frame(c(effects.scores$QTL, overlap$QTL)) %>%
    `colnames<-`(c("QTL")) %>%
    dplyr::filter(!duplicated(QTL)) %>%
    dplyr::mutate(QTL = as.character(QTL),
                Simulated = (QTL %in% effects.scores$QTL),
                Detected = (QTL %in% overlap$QTL)) %>%
    dplyr::full_join(.,effects.scores, by = "QTL") %>%
    dplyr::full_join(.,overlap, by = "QTL") %>%
    dplyr::mutate(
                #algorithm = algorithm_id ,
                top.hit = QTL == detected.peak,
                sim = sim_id,
                species = sp_id,
                algorithm = mapping_algorithm) ### ASIGN A SIMULATION ID HERE
#Create factor columns for the outcomes

all.QTL$Simulated <- factor(all.QTL$Simulated, levels = c("TRUE","FALSE"))
all.QTL$Detected <- factor(all.QTL$Detected, levels = c("TRUE","FALSE"))

#Write the results to file
write.tsv(all.QTL, file = glue::glue("{sim_id}.{mapping_algorithm}.simulated.mapping.results.scores.tsv"), sep = "\t")