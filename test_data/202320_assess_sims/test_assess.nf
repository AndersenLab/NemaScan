#! usr/bin/env nextflow
if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

nextflow.enable.dsl=2
// nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )

params.out = "Analysis_Results-${date}"
params.bin_dir = "${workflow.projectDir}/bin" // this is different for gcp

include {prepare_simulation_files; chrom_eigen_variants_sims; collect_eigen_variants_sims; simulate_effects_loc; simulate_effects_genome; simulate_map_phenotypes; get_gcta_intervals; assess_sims} from './modules/simulations.nf'

workflow{
        assess_1 = [
            "AB1,AB2",
            "AB1,AB2",
            "10",
            "1",
            "0.8",
            "0.5",
            "gamma", 
            "/projects/b1059/projects/Ryan/NemaScan_Updates/NemaScan/20231019_asses_test_1/Simulations/gamma/5/Phenotypes/5_2_0.8_0.05_gamma_ce.96.allout15_irrepressible.grosbeak_sims.par",
            "/projects/b1059/projects/Ryan/NemaScan_Updates/NemaScan/20231019_asses_test_1/Simulations/gamma/5/Phenotypes/5_1_0.8_0.05_gamma_ce.96.allout15_irrepressible.grosbeak_sims.phen",
            "/projects/b1059/projects/Ryan/NemaScan_Updates/NemaScan/20231019_asses_test_1/Genotype_Matrix/ce.96.allout15_irrepressible.grosbeak_0.05_Genotype_Matrix.tsv",
            "/projects/b1059/projects/Ryan/NemaScan_Updates/NemaScan/20231019_asses_test_1/Simulations/gamma/5/Mappings/5_1_0.8_0.05_gamma_ce.96.allout15_irrepressible.grosbeak_processed_LMM-EXACT-INBRED_mapping.tsv",
            "${params.bin_dir}/Assess_Sim.R"
            ]
        assess_2 = [
            "AB1,AB2",
            "AB1,AB2",
            "10",
            "1",
            "0.8",
            "0.5",
            "gamma", 
            "/projects/b1059/projects/Ryan/NemaScan_Updates/NemaScan/20231019_asses_test_1/Simulations/gamma/5/Phenotypes/5_2_0.8_0.05_gamma_ce.96.allout15_irrepressible.grosbeak_sims.par",
            "/projects/b1059/projects/Ryan/NemaScan_Updates/NemaScan/20231019_asses_test_1/Simulations/gamma/5/Phenotypes/5_1_0.8_0.05_gamma_ce.96.allout15_irrepressible.grosbeak_sims.phen",
            "/projects/b1059/projects/Ryan/NemaScan_Updates/NemaScan/20231019_asses_test_1/Genotype_Matrix/ce.96.allout15_irrepressible.grosbeak_0.05_Genotype_Matrix.tsv",
            "/projects/b1059/projects/Ryan/NemaScan_Updates/NemaScan/20231019_asses_test_1/Simulations/gamma/5/Mappings/5_1_0.8_0.05_gamma_ce.96.allout15_irrepressible.grosbeak_processed_LMM-EXACT-INBRED_mapping.tsv",
            "${params.bin_dir}/Assess_Sim.R"
            ]
        
        Channel.fromList([assess_1, assess_2]) | assess_sims | collectFile(name: "${params.out}/all_sims.tsv") | view

        //both = Channel.from(assess_1, assess_2)

        //both | assess_sims

        // run assess_sims process on both channels


        
        
        
        //.collectFile(name: "${params.out}/all_sims.tsv", newLine = true)


}