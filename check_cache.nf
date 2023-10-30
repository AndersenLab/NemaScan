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
params.maf = 0.05

workflow{

File pop_file = new File("test_data/repeated_sim_strain_sets.txt") ;

sp_ids = [["c_elegans", "/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.vcf.gz", "/projects/b1059/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.vcf.gz.tbi"],
                ["c_briggsae", "/projects/b1059/data/c_briggsae/WI/variation/20210803/vcf/WI.20210803.hard-filter.isotype.vcf.gz", "/projects/b1059/data/c_briggsae/WI/variation/20210803/vcf/WI.20210803.hard-filter.isotype.vcf.gz.tbi"],
                ["c_tropicalis", "/projects/b1059/data/c_tropicalis/WI/variation/20210901/vcf/WI.20210901.hard-filter.isotype.vcf.gz", "/projects/b1059/data/c_tropicalis/WI/variation/20210901/vcf/WI.20210901.hard-filter.isotype.vcf.gz.tbi"]]

Channel.from(pop_file.collect { it.tokenize( ' ' ) })
          .map {SP, SM, STRAINS -> [SP, SM, STRAINS]}
          .combine()

}

process check_cache {
    input:
    tuple val(SP), val(SM), val(STRAINS), file(check_cache), path(work_flow_dir)

    output:
    stdout dir_status

    """
    python ${check_cache} ${work_flow_dir} ${SP} ${SM}
    """
}