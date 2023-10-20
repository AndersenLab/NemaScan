#! usr/bin/env nextflow
if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

nextflow.preview.dsl=2
// nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )

params.out = "Analysis_Results-${date}"
params.bin_dir = "${workflow.projectDir}/bin" // this is different for gcp

include {prepare_simulation_files; chrom_eigen_variants_sims; collect_eigen_variants_sims; simulate_effects_loc; simulate_effects_genome; simulate_map_phenotypes; get_gcta_intervals; assess_sims} from './modules/simulations.nf'

workflow{
        Channel.from("AB1,AB2") // strain set 
            .combine(Channel.from("AB1,AB2")) // strains
            .combine(Channel.from("10")) // NQTL
            .combine(Channel.from("1")) // SIMREP
            .combine(Channel.from("0.8")) // h2
            .combine(Channel.from("0.5")) // maf
            .combine(Channel.from("gamma")) // effect range
            .combine(Channel.fromPath("/projects/b1059/projects/Ryan/NemaScan_Updates/NemaScan/20231019_asses_test_1/Simulations/gamma/5/Phenotypes/5_2_0.8_0.05_gamma_ce.96.allout15_irrepressible.grosbeak_sims.par")) //causal variants file
            .combine(Channel.fromPath("/projects/b1059/projects/Ryan/NemaScan_Updates/NemaScan/20231019_asses_test_1/Simulations/gamma/5/Phenotypes/5_1_0.8_0.05_gamma_ce.96.allout15_irrepressible.grosbeak_sims.phen")) // 
            .combine(Channel.fromPath("/projects/b1059/projects/Ryan/NemaScan_Updates/NemaScan/20231019_asses_test_1/Genotype_Matrix/ce.96.allout15_irrepressible.grosbeak_0.05_Genotype_Matrix.tsv")) // genotype matrix
            .combine(Channel.fromPath("/projects/b1059/projects/Ryan/NemaScan_Updates/NemaScan/20231019_asses_test_1/Simulations/gamma/5/Mappings/5_1_0.8_0.05_gamma_ce.96.allout15_irrepressible.grosbeak_processed_LMM-EXACT-INBRED_mapping.tsv")) 
            .combine((Channel.fromPath("${params.bin_dir}/Assess_Sim.R")))
            | assess_sims
}