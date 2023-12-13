#! usr/bin/env nextflow

if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

nextflow.preview.dsl=2
// nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd_HHmmss' )
params.bin_dir = "${workflow.projectDir}/bin" // this is different for gcp
params.out = "Analysis_Results-${date}"
sim_dir = "/projects/b1059/projects/Ryan/Caenorhabditis_GWAS/best_panel_subsample/20231102_CB_96_Alloutlier/20231102_CB_96_Alloutlier_Subsample_Simulations"

include {assess_sims_INBRED; assess_sims_LOCO} from './modules/simulations.nf'

workflow{
    //create a channel from the Simulations folder for CE - WILL ONLY WORK WITH ON NQTL PARAMETER SIMULATIONS 
    //print status messages to the console
    //process.statusCommand = "echo 'Running...'"
    effects_ch = Channel.fromPath("${sim_dir}/Simulations/gamma/5/Phenotypes/*.phen")
        // pull out the simulation ID from the file name
        .map{ file -> 
            def parts = file.getBaseName().split('_')
            return tuple(parts.take(7).join('_'), file)
        }
    
    //effects_ch.view()
    
    //channel for processed INBRED-PCA mappings - WILL ONLY WORK WITH ON NQTL PARAMETER SIMULATIONS & ONLY PULLS THE INBRED MAPPINGS
    //processed_LMM-EXACT-INBRED_PCA_mapping.tsv for CB.96 vs processed_LMM-EXACT-INBRED_mapping.tsv for CE. 96 and 192
    proc_mapping_ch = Channel.fromPath("${sim_dir}/Simulations/gamma/5/Mappings/*processed_LMM-EXACT-INBRED_PCA_mapping.tsv")
        .map{ file -> 
            def parts = file.getBaseName().split('_')
            return tuple(parts.take(7).join('_'), file)
        }

    //proc_mapping_ch.view()

    //genotype matrix input channel 
    gm_ch = Channel.fromPath(
        "${sim_dir}/Genotype_Matrix/*_Genotype_Matrix.tsv"
        )
        .map{ file -> 
            def parts = file.getBaseName().split('_')
            return tuple(parts[0], parts[1], parts[2] ,file)
        }

    //gm_ch.view()
    //causal variant input channel
    cv_ch = Channel.fromPath(
        "${sim_dir}/Simulations/gamma/5/Phenotypes/*.par"
        )
        .map{ file -> 
            def parts = file.getBaseName().split('_')
            return tuple(parts.take(7).join('_'), file)
        }

    //cv_ch.view()
    //join the effects channel to the mapping chanel using the simulation ID
    joined_ch = effects_ch
        .join(proc_mapping_ch, by: 0)
        //.join(gm_ch, by: 0)
        .join(cv_ch, by: 0)
        //.map{ tuple -> 
        //def (sim_id, effects_file, mapping_file) = tuple
        //return tuple(sim_id, effects_file, mapping_file)
       // }
    //joined_ch.view()

    //pull out simulation parameters from the simulation id
    assess_input_ch = joined_ch
        .map{ tuple -> 
            def (sim_id, effects_file, mapping_file, cv_file) = tuple
            def parts = sim_id.split('_')
                    // 1 - NQTL, 2 - SIMREP, 3 - h2,  4 - maf, 5 - effect,  6 - pop_id, 7 - strain set 
            return [parts[5], parts[6], parts[3], parts[0], parts[1], parts[2], parts[4], effects_file, mapping_file, cv_file]
        }
    //assess_input_ch.view()

    //add the genotype matrix to the channel by joining the assess_input_ch to the gm_ch using the pop_id and strain set
    gm_joined = gm_ch
        .combine(assess_input_ch, by: [0,1,2])
        // apply function to meet the input cardinality requirements of the assess_sims process
        .map{ tuple -> 
            def (pop_id, panel_id, maf, gm, nqtl, sim_rep, h2, effect_range, phenotypes, mapping_file, var_effects) = tuple
                strain_set = pop_id + '_' + panel_id
            return [strain_set, panel_id, nqtl, sim_rep, h2, maf, effect_range, var_effects, phenotypes, gm, mapping_file, "LMM-EXACT-INBRED_PCA"]
        }
    //gm_joined.view()

    gm_joined
        .combine(Channel.fromPath("${params.bin_dir}/Assess_Sim.R")) | assess_sims_INBRED | collectFile(name: "${params.out}/INBRED_PCA_all_sims.tsv")
}