library(tidyverse)

#Define the results directoruy
results_dir <- "NemaScan/Analysis_Results-20230628" 

#For each species (c_elegans, c_briggsae, c_tropicalis)
# Find all the simulation directories (e.g. Analysis_Results-20230628/c_elegans/Simulations/*)
# 1. Read in the summary file for both algorithms (<results_dir>/<Simulations>/<species>/<sim_id>/<sim_id>.<LOCO/INBRED>.simulated.mapping,results.scores.tsv)
# 2. Aggregate the results into a single table

#species <- "c_elegans"
#sim_dirs <- list.files(file.path(results_dir, "Simulations", species), full.names = TRUE)

#For each species
all_sims <- data.frame()
for (species in c("c_elegans", "c_briggsae", "c_tropicalis")) {
  #Find all the simulation directories
  sim_dirs <- list.files(file.path(results_dir, "Simulations", species), full.names = TRUE)
  #For each simulation directory
  for (sim_dir in sim_dirs) {
    #Get the simulation id
    sim <- basename(sim_dir)
    #Read in the summary file for both algorithms
    loco_summary <- data.table::fread(file.path(sim_dir, glue::glue("{sim}.LOCO.simulated.mapping.results.scores.tsv")))
    inbred_summary <- data.table::fread(file.path(sim_dir, glue::glue("{sim}.INBRED.simulated.mapping.results.scores.tsv")))
    #Aggregate the results into a single table
    summary <- rbind(loco_summary, inbred_summary)
    #Aggregate the results into a single table
    #write_tsv(summary, file.path(sim_dir, "simulated.mapping.results.scores.tsv"))
    all_sims <- rbind(all_sims, summary)
  }
    
}

#Check all_sims df

#Write out the results

write_tsv(all_sims, file.path(results_dir, "simulated.mapping.results.scores.tsv"))


# Pull the data for all markers -- ONLY FOR FIRST REPLICATE SIMULATION
# Marker file example name c_briggsae_<sim_rep>_<h2>_process_<algorithm>_mapping.tsv
all_sim_markers <- data.frame()
for (species in c("c_elegans", "c_briggsae", "c_tropicalis")){
  sim_dirs <- list.files(file.path(results_dir, "Simulations", species), full.names = TRUE)
  for (sim_dir in sim_dirs){
    sim <- basename(sim_dir)
    # Attempt to read in the loco and inbred marker files warn if they don't exist and skip
     tryCatch({
      loco_markers <- data.table::fread(file.path(sim_dir, glue::glue("Mappings/{species}_1_0.8_processed_LMM-EXACT-LOCO_PCA_mapping.tsv")))%>%
        mutate(algorithm = "LOCO")%>%
        mutate(sim = sim)%>%
        mutate(species = species)

      inbred_markers <- data.table::fread(file.path(sim_dir, glue::glue("Mappings/{species}_1_0.8_processed_LMM-EXACT-INBRED_PCA_mapping.tsv")))%>%
        mutate(algorithm = "INBRED")%>%
        mutate(sim = sim)%>%
        mutate(species = species)%>%
        select(-N) #Remove N column unique to inbred data
      
      markers <- rbind(loco_markers, inbred_markers)
      
      all_sim_markers <- rbind(all_sim_markers, markers)
    }, error = function(e) { 
      print(glue::glue("No marker file for {sim}"))
    }) 
  }
}

write_tsv(all_sim_markers, file.path(results_dir, "simulated.mapping.results.markers.tsv"))
