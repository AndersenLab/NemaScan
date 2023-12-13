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
params.out = "test_Analysis_Results-${date}"
sim_dir = "/projects/b1059/projects/Ryan/Caenorhabditis_GWAS/best_panel_subsample/20231017_CE_96.192_Alloutlier/20231024_CE_96.192_Alloutlier"

params.group_qtl = 1000
params.sthresh = "BF"
params.ci_size = 150

include {get_gcta_intervals; assess_sims_INBRED; assess_sims_LOCO} from './modules/simulations.nf'

workflow{
    //create a channel from the Simulations folder for CE - WILL ONLY WORK WITH ON NQTL PARAMETER SIMULATIONS 
    //print status messages to the console
    //process.statusCommand = "echo 'Running...'"
    
    //pull the phenotypes file
    effects_ch = Channel.fromPath("${sim_dir}/Simulations/gamma/5/Phenotypes/*_*_*_*_*_ce.192*.phen")
        // pull out the simulation ID from the file name
        .map{ file -> 
            def parts = file.getBaseName().split('_')
            return tuple(parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6], file)
        }
    
    //effects_ch.view()
    
    //causal variant input channel

    cv_ch = Channel.fromPath(
        "${sim_dir}/Simulations/gamma/5/Phenotypes/*_*_*_*_*_ce.192*.par"
        )
        .map{ file -> 
            def parts = file.getBaseName().split('_')
            return tuple(parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6], file)
        }

    //cv_ch.view()


    //create a channel for the raw mapping files
    raw_mappings_inbred = Channel.fromPath("${sim_dir}/Simulations/gamma/5/Mappings/*_*_*_*_*_ce.192*_lmm-exact_inbred_pca.fastGWA")
        // pull out the simulation ID from the file name
        .map{ file -> 
            def parts = file.getBaseName().split('_')
                    // 1 - NQTL, 2 - SIMREP, 3 - h2,  4 - maf, 5 - effect_range,  6 - pop_id, 7 - strain set
            return tuple(parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6], file)        
        }

    raw_mappings_loco = Channel.fromPath("${sim_dir}/Simulations/gamma/5/Mappings/*_*_*_*_*_ce.192*_lmm-exact_pca.loco.mlma")
        // pull out the simulation ID from the file name
        .map{ file -> 
            def parts = file.getBaseName().split('_')
                    // 1 - NQTL, 2 - SIMREP, 3 - h2,  4 - maf, 5 - effect_range,  6 - pop_id, 7 - strain set
            return tuple(parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6], file)        
        }

    // join the raw mapping channels
    raw_mappings_joined = raw_mappings_inbred
        .join(raw_mappings_loco, by: [0,1,2,3,4,5,6])
        .join(effects_ch, by: [0,1,2,3,4,5,6])
        .join(cv_ch, by: [0,1,2,3,4,5,6])
        .map{ tuple -> 
            def (nqtl, simrep, h2, maf, effect_range, pop_id, strain_set, inbred_mapping, loco_mapping, phenotypes, causal_variants) = tuple
            return [pop_id, strain_set, maf, nqtl, simrep, h2, effect_range, inbred_mapping, loco_mapping, phenotypes, causal_variants]
        }

    //raw_mappings_joined.view()    

    //raw_mappings_inbred.view()
    // send the raw mapping file to get_gcta_intervals process 


    //channel for processed INBRED-PCA mappings - WILL ONLY WORK WITH ON NQTL PARAMETER SIMULATIONS & ONLY PULLS THE INBRED MAPPINGS
    //processed_LMM-EXACT-INBRED_PCA_mapping.tsv for CB.192 vs processed_LMM-EXACT-INBRED_mapping.tsv for CE. 192 and 192
    //proc_mapping_ch = Channel.fromPath("${sim_dir}/Simulations/gamma/5/Mappings/*processed_LMM-EXACT-INBRED_PCA_mapping.tsv")
    //    .map{ file -> 
    //        def parts = file.getBaseName().split('_')
    //        return tuple(parts.take(7).join('_'), file)
    //    }

    //proc_mapping_ch.view()

    //genotype matrix input channel 
    gm_ch = Channel.fromPath(
        "${sim_dir}/Genotype_Matrix/ce.192*_Genotype_Matrix.tsv"
        )
        .map{ file -> 
            def parts = file.getBaseName().split('_')
            return tuple(parts[0], parts[1], parts[2] ,file)
        }

    //gm_ch.view()
    
    //independent tests input chanel
    it_ch = Channel.fromPath(    
        "${sim_dir}/Genotype_Matrix/ce.192*_total_independent_tests.txt"
        )
        .map{ file -> 
            def parts = file.getBaseName().split('_')
            return tuple(parts[0], parts[1], parts[2] ,file)
        }
    
    // join the genotype matrix channel to the independent tests channel using the pop_id and strain set and MAF
    gm_joined = gm_ch
        .combine(it_ch, by: [0,1,2])
    
    //gm_joined.view()
    
    // join the genotype matrix channel to the raw mapping channel using the pop_id and strain set and MAF
    get_intervals_input = gm_joined
        .combine(raw_mappings_joined, by: [0,1,2])
        .map{ tuple -> 
            def (pop_id, panel_id, maf, gm, it, nqtl, simrep, h2, effect_range, inbred_mapping, loco_mapping, phenotypes, causal_variants) = tuple
            strain_set = pop_id + '_' + panel_id
            return [strain_set, panel_id, nqtl, simrep, h2, "fake_loci_file", gm, effect_range, it, maf, inbred_mapping, loco_mapping, causal_variants, phenotypes]
            
            
        }
        .combine(Channel.from("${params.sthresh}"))
        .combine(Channel.from("${params.group_qtl}"))
        .combine(Channel.from("${params.ci_size}")) 
        .combine(Channel.fromPath("${params.bin_dir}/Aggregate_Mappings.R"))
        .combine(Channel.fromPath("${params.bin_dir}/Find_Aggregate_Intervals.R"))
        .combine(Channel.fromPath("${params.bin_dir}/Find_GCTA_Intervals.R"))
        .combine(Channel.fromPath("${params.bin_dir}/Find_GCTA_Intervals_LOCO.R")) | get_gcta_intervals

        get_gcta_intervals.out.assess_data_inbred_pca 
            .combine(Channel.fromPath("${params.bin_dir}/Assess_Sim.R")) | assess_sims_INBRED 
    //get_intervals_input.view()
    //join the effects channel to the mapping chanel using the simulation ID
    //joined_ch = effects_ch
    //    .join(proc_mapping_ch, by: 0)
        //.join(gm_ch, by: 0)
    //    .join(cv_ch, by: 0)
        //.map{ tuple -> 
        //def (sim_id, effects_file, mapping_file) = tuple
        //return tuple(sim_id, effects_file, mapping_file)
       // }
    //joined_ch.view()

    //pull out simulation parameters from the simulation id
    // assess_input_ch = joined_ch
    //    .map{ tuple -> 
    //        def (sim_id, effects_file, mapping_file, cv_file) = tuple
    //        def parts = sim_id.split('_')
                    // 1 - NQTL, 2 - SIMREP, 3 - h2,  4 - maf, 5 - effect,  6 - pop_id, 7 - strain set 
    //        return [parts[5], parts[6], parts[3], parts[0], parts[1], parts[2], parts[4], effects_file, mapping_file, cv_file]
    //    }
    //assess_input_ch.view()

    //add the genotype matrix to the channel by joining the assess_input_ch to the gm_ch using the pop_id and strain set
    //gm_joined = gm_ch
    //    .combine(assess_input_ch, by: [0,1,2])
        // apply function to meet the input cardinality requirements of the assess_sims process
    //    .map{ tuple -> 
    //        def (pop_id, panel_id, maf, gm, nqtl, sim_rep, h2, effect_range, phenotypes, mapping_file, var_effects) = tuple
    //            strain_set = pop_id + '_' + panel_id
    //        return [strain_set, panel_id, nqtl, sim_rep, h2, maf, effect_range, var_effects, phenotypes, gm, mapping_file, "LMM-EXACT-INBRED_PCA"]
    //    }
    //gm_joined.view()

    //gm_joined
    //    .combine(Channel.fromPath("${params.bin_dir}/Assess_Sim.R")) | assess_sims_INBRED | collectFile(name: "${params.out}/INBRED_PCA_all_sims.tsv")
}