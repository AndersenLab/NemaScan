#! usr/bin/env nextflow
if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

nextflow.preview.dsl=2



params.bin_dir = "${workflow.projectDir}/bin" // this is different for gcp
params.master_snp_dir = "test_data/master_snps"
params.sparse_cut = 0.05
params.group_qtl = 1000
params.ci_size = 150
sthresh= "BF"
// params.simulate_h2 = "/projects/b1059/projects/Ryan/ortholog_sims/NemaScan/test_data/h2.csv"

include {prepare_repeated_simulation_files; chrom_eigen_variants_sims_repeated; collect_eigen_variants_sims_repeated; simulate_orthogroup_effects; simulate_map_phenotypes; get_gcta_intervals_repeated} from './modules/repeated_simulations.nf'

//ce_vcf = Channel.fromPath("/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.bcsq.vcf.gz")
//ce_vcf_index = Channel.fromPath("/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.bcsq.vcf.gz.tbi")

//cb_vcf = Channel.fromPath("/projects/b1059/data/c_briggsae/WI/variation/20210803/vcf/WI.20210803.hard-filter.isotype.bcsq.vcf.gz")
//cb_vcf_index = Channel.fromPath("/projects/b1059/data/c_briggsae/WI/variation/20210803/vcf/WI.20210803.hard-filter.isotype.bcsq.vcf.gz.tbi")

//ct_vcf = Channel.fromPath("/projects/b1059/data/c_tropicalis/WI/variation/20210901/vcf/WI.20210901.hard-filter.isotype.bcsq.vcf.gz")
//ct_vcf_index = Channel.fromPath("/projects/b1059/data/c_tropicalis/WI/variation/20210901/vcf/WI.20210901.hard-filter.isotype.bcsq.vcf.gz.tbi")
rename_key = Channel.fromPath("input_data/all_species/rename_chromosomes")
maf_file = Channel.fromPath("input_data/all_species/simulate_maf.csv").splitCsv()
workflow{

File pop_file = new File("test_data/test_orthogroup_samples.txt") ;

sp_ids = [["c_elegans", "/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.vcf.gz", "/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.vcf.gz.tbi", "/projects/b1059/projects/Ryan/ortholog_sims/NemaScan/test_data/c_elegans/underground.gartersnake/plink_files" ],
                ["c_briggsae", "/projects/b1059/data/c_briggsae/WI/variation/20210803/vcf/WI.20210803.hard-filter.isotype.vcf.gz", "/projects/b1059/data/c_briggsae/WI/variation/20210803/vcf/WI.20210803.hard-filter.isotype.vcf.gz.tbi", "/projects/b1059/projects/Ryan/ortholog_sims/NemaScan/test_data/c_briggsae/aboveground.gartersnake/plink_files" ],
                ["c_tropicalis", "/projects/b1059/data/c_tropicalis/WI/variation/20210901/vcf/WI.20210901.hard-filter.isotype.vcf.gz", "/projects/b1059/data/c_tropicalis/WI/variation/20210901/vcf/WI.20210901.hard-filter.isotype.vcf.gz.tbi", "/projects/b1059/projects/Ryan/ortholog_sims/NemaScan/test_data/c_tropicalis/inbetweenground.gartersnake/plink_files"]]

Channel.from(pop_file.collect { it.tokenize( ' ' ) })
          .map {SP, SM, STRAINS -> [SP, SM, STRAINS] }
          .join(Channel.from(sp_ids), by:[0])
          .combine(rename_key)
          .combine(maf_file)
          .map {
            tuple ->[tuple[0], //extract species ID
                    (tuple[1]), //extract sample name
                    (tuple[2]), //extract sample list
                    file(tuple[3]), // convert path to file obj for vcf
                    file(tuple[4]), // index
                    file(tuple[5]), // plink_dir FOR TESTING 
                    file(tuple[6]), // rename key 
                    tuple[7]] // MAF
        } |  prepare_repeated_simulation_files

    // eigen
    contigs = Channel.from(["1", "2", "3", "4", "5", "6"]) //Parallelize by chrom
    contigs.combine(prepare_repeated_simulation_files.out.sim_geno) // Combine with Plink files and Genotype matrix + Sim INFO
        .combine(Channel.fromPath("bin/Get_GenoMatrix_Eigen.R")) | chrom_eigen_variants_sims_repeated
    
    // Collect the eigen results
    chrom_eigen_variants_sims_repeated.out.sim_geno_eigen_join
        .groupTuple(by:[0,1,2,3]). // Collect all chromosome eigen files with the same SP, strain_set, strains, and MAF
        join(chrom_eigen_variants_sims_repeated.out.sim_geno_meta, by:[0,1,2,3]) | collect_eigen_variants_sims_repeated

    collect_eigen_variants_sims_repeated.out
        .combine(Channel.fromPath("test_data/causal_ogs.txt").splitCsv())
        .combine(Channel.from(1..2)) // number of reps per OG trait
        .combine(Channel.fromPath("${params.bin_dir}/sim_og_effects.py"))
        .combine(Channel.fromPath("${params.master_snp_dir}"))
        | simulate_orthogroup_effects
    
    //simulate_orthogroup_effects.out.view()

    sim_phen_inputs = simulate_orthogroup_effects.out.pheno_inputs 

    sim_phen_inputs
        .combine(Channel.fromPath("/projects/b1059/projects/Ryan/ortholog_sims/NemaScan/test_data/h2.csv").splitCsv()) | simulate_map_phenotypes

    simulate_map_phenotypes.out.gcta_intervals
            .combine(Channel.from("${params.sthresh}"))
            .combine(Channel.from("${params.group_qtl}"))
            .combine(Channel.from("${params.ci_size}")) 
            .combine(Channel.fromPath("${params.bin_dir}/Aggregate_Mappings.R"))
            .combine(Channel.fromPath("${params.bin_dir}/Find_Aggregate_Intervals.R"))
            .combine(Channel.fromPath("${params.bin_dir}/Find_GCTA_Intervals_Repeated.R"))
            .combine(Channel.fromPath("${params.bin_dir}/Find_GCTA_Intervals_LOCO_Repeated.R")) | get_gcta_intervals_repeated
}

