#! usr/bin/env nextflow


include {prepare_repeated_simulation_files_temp; chrom_eigen_variants_sims_repeated_temp; collect_eigen_variants_sims_repeated_temp } from './modules/repeated_simulations.nf'

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

sp_ids = [["c_elegans", "test_data/vcfs/c_elegans/WI.20220216.hard-filter.isotype.bcsq.vcf.gz", "test_data/vcfs/c_elegans/WI.20220216.hard-filter.isotype.bcsq.vcf.gz.tbi"],
                ["c_briggsae", "test_data/vcfs/c_briggsae/WI.20210803.hard-filter.isotype.bcsq.vcf.gz", "test_data/vcfs/c_briggsaeWI.20210803.hard-filter.isotype.bcsq.vcf.gz.tbi"],
                ["c_tropicalis", "test_data/vcfs/c_tropicalis/WI.20210901.hard-filter.isotype.bcsq.vcf.gz", "test_data/vcfs/c_tropicalis/WI.20210901.hard-filter.isotype.bcsq.vcf.gz.tbi"]]

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
                    file(tuple[5]), // rename key
                    tuple[6]] // MAF
        } |  prepare_repeated_simulation_files_temp

    // eigen
    contigs = Channel.from(["1", "2", "3", "4", "5", "6"]) //Parallelize by chrom
    contigs.combine(prepare_repeated_simulation_files_temp.out.sim_geno) // Combine with Plink files and Genotype matrix + Sim INFO
        .combine(Channel.fromPath("bin/Get_GenoMatrix_Eigen.R")) | chrom_eigen_variants_sims_repeated_temp
    
    // Collect the eigen results
    chrom_eigen_variants_sims_repeated_temp.out.sim_geno_eigen_join
        .groupTuple(by:[0,1,2,3]). // Collect all chromosome eigen files with the same SP, strain_set, strains, and MAF
        join(chrom_eigen_variants_sims_repeated_temp.out.sim_geno_meta, by:[0,1,2,3]) | collect_eigen_variants_sims_repeated_temp

    collect_eigen_variants_sims_repeated_temp.out
        .combine(Channel.fromPath("test_data/passing_ogs.csv").splitCsv()).view()

}

