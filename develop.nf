#! usr/bin/env nextflow

if( !nextflow.version.matches('20.0+') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )

/*
~ ~ ~ > * Parameters: common to all analyses
*/

params.help        = null
if(params.simulate) {
    params.e_mem   = "100"
} else {
    params.e_mem   = "10" // I noticed it was crashing with 100 gb for mappings... maybe too much allocation?
}
params.eigen_mem   = params.e_mem + " GB"
params.out         = "Analysis_Results-${date}"
params.debug       = null
params.species     = "elegans"
params.wbb         = "WS276"
params.data_dir    = "${workflow.projectDir}/input_data/${params.species}"
params.numeric_chrom = "${workflow.projectDir}/input_data/all_species/rename_chromosomes"
params.sparse_cut  = 0.01
params.group_qtl   = 1000
params.ci_size     = 150
params.sthresh     = "BF"
params.p3d         = "TRUE"
params.maf         = 0.05
params.genes       = "${workflow.projectDir}/bin/gene_ref_flat.Rda"
params.gcp         = null
download_vcf       = null
params.annotation  = "bcsq"
params.strains     = null
params.MAF         = 0.05

/*
~ ~ ~ > * Parameters: for burden mapping
*/

params.refflat   = "${params.data_dir}/annotations/c_${params.species}_${params.wbb}_refFlat.txt"
params.freqUpper = 0.05
params.minburden = 2


/*
~ ~ ~ > * Parameters: VCF
*/


if(params.debug) {
    println """
        *** Using debug mode ***
    """
    // debug for now with small vcf
    params.vcf = "330_TEST.vcf.gz"
    params.traitfile = "${workflow.projectDir}/input_data/elegans/phenotypes/abamectin_pheno.tsv"

    vcf_file = Channel.fromPath("${workflow.projectDir}/input_data/elegans/genotypes/330_TEST.vcf.gz")
    vcf_index = Channel.fromPath("${workflow.projectDir}/input_data/elegans/genotypes/330_TEST.vcf.gz.tbi")
        
    // debug can use same vcf for impute and normal
    impute_vcf = Channel.fromPath("${workflow.projectDir}/input_data/elegans/genotypes/330_TEST.vcf.gz")
    impute_vcf_index = Channel.fromPath("${workflow.projectDir}/input_data/elegans/genotypes/330_TEST.vcf.gz.tbi")
    
    ann_file = Channel.fromPath("${workflow.projectDir}/input_data/elegans/genotypes/WI.330_TEST.strain-annotation.bcsq.tsv")
} else if(params.gcp) { 
    // use the data directly from google on gcp
    vcf_file = Channel.fromPath("gs://caendr-data/releases/${params.vcf}/variation/WI.${params.vcf}.hard-filter.isotype.vcf.gz")
    vcf_index = Channel.fromPath("gs://caendr-data/releases/${params.vcf}/variation/WI.${params.vcf}.hard-filter.isotype.vcf.gz.tbi")

    impute_vcf = Channel.fromPath("gs://caendr-data/releases/${params.vcf}/variation/WI.${params.vcf}.impute.isotype.vcf.gz")
    impute_vcf_index = Channel.fromPath("gs://caendr-data/releases/${params.vcf}/variation/WI.${params.vcf}.impute.isotype.vcf.gz.tbi")

    ann_file = Channel.fromPath("gs://caendr-data/releases/${params.vcf}/variation/WI.${params.vcf}.strain-annotation.bcsq.tsv")
} else if(!params.vcf) {
    // if there is no VCF date provided, pull the latest vcf from cendr.
    params.vcf = "20210121"
    download_vcf = true
    
} else {
    // use the vcf data from QUEST when a cendr date is provided

    // Check that params.vcf is valid
    if("${params.vcf}" == "20210121" || "${params.vcf}" == "20200815" || "${params.vcf}" == "20180527" || "${params.vcf}" == "20170531") {
        vcf_file = Channel.fromPath("/projects/b1059/data/c_elegans/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.hard-filter.isotype.vcf.gz")
        vcf_index = Channel.fromPath("/projects/b1059/data/c_elegans/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.hard-filter.isotype.vcf.gz.tbi")

        impute_vcf = Channel.fromPath("/projects/b1059/data/c_elegans/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.impute.isotype.vcf.gz")
        impute_vcf_index = Channel.fromPath("/projects/b1059/data/c_elegans/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.impute.isotype.vcf.gz.tbi")

        // check if cendr release date is not 20210121, use snpeff annotation
        if("${params.vcf}" != "20210121") {
            println "WARNING: Using snpeff annotation. To use BCSQ annotation, please use --vcf 20210121"
            ann_file = Channel.fromPath("/projects/b1059/data/c_elegans/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.strain-annotation.snpeff.tsv")
        } else {
            ann_file = Channel.fromPath("/projects/b1059/data/c_elegans/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.strain-annotation.bcsq.tsv")
        }
    } else {
        println "Error! Cannot find ${params.vcf}.vcf. Please provide a valid CeNDR release date (20210121, 20200815, 20180527, or 20170531)."
        exit 1
    }
}

// If mapping, check that traitfile exists
if(params.map) {
    if (!file("${params.traitfile}").exists()) {
        println """
        Error: Phenotype input file (${params.traitfile}) does not exist.
        """
        System.exit(1)
    } 
}



if (params.help) {
    log.info '''

O~~~     O~~                                   O~~ ~~
O~ O~~   O~~                                 O~~    O~~
O~~ O~~  O~~   O~~    O~~~ O~~ O~~    O~~     O~~         O~~~   O~~    O~~ O~~
O~~  O~~ O~~ O~   O~~  O~~  O~  O~~ O~~  O~~    O~~     O~~    O~~  O~~  O~~  O~~
O~~   O~ O~~O~~~~~ O~~ O~~  O~  O~~O~~   O~~       O~~ O~~    O~~   O~~  O~~  O~~
O~~    O~ ~~O~         O~~  O~  O~~O~~   O~~ O~~    O~~ O~~   O~~   O~~  O~~  O~~
O~~      O~~  O~~~~   O~~~  O~  O~~  O~~ O~~~  O~~ ~~     O~~~  O~~ O~~~O~~~  O~~
    '''
    log.info "----------------------------------------------------------------"
    log.info "                         USAGE                                  "
    log.info "----------------------------------------------------------------"
    log.info ""
    log.info "nextflow develop.nf --traitfile=input_data/elegans/phenotypes/abamectin_pheno.tsv --vcf 20210121"
    log.info ""
    log.info "Profiles available:"
    log.info "-profile mappings        (Default)             Perform GWA mappings with a provided trait file"
    log.info "-profile simulations                           Perform phenotype simulations with GCTA"
    log.info "----------------------------------------------------------------"
    log.info "                         DEBUG                                  "
    log.info "----------------------------------------------------------------"
    log.info ""
    log.info "nextflow develop.nf --debug"
    log.info ""
    log.info "----------------------------------------------------------------"
    log.info "             -profile annotations USAGE"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "nextflow develop.nf --vcf 20210121 -profile annotations --species elegans"
    log.info "----------------------------------------------------------------"
    log.info "Mandatory arguments:"
    log.info "--wb_build               String                Wormbase build number, must be greater than WS270"
    log.info "--species                String                What species to download information for (elegans, briggsae, or tropicalis)"
    log.info "----------------------------------------------------------------"
    log.info "             -profile mappings USAGE"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "nextflow develop.nf --vcf 20210121 --traitfile input_data/elegans/phenotypes/PC1.tsv -profile mappings"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "Mandatory arguments:"
    log.info "--traitfile              String                Name of file that contains phenotypes. File should be tab-delimited with the columns: strain trait1 trait2 ..."
    log.info "--vcf                    String                CeNDR release date of VCF to extract variants from. If none is provided, the pipeline will download the latest VCF from CeNDR"
    log.info "Optional arguments:"
    log.info "--maf                    String                Minimum minor allele frequency to use for single-marker mapping (Default: 0.05)"
    log.info "--lmm                    String                Perform GCTA mapping with --fastGWA-lmm algorithm (Default: RUN, option to not run is null)"
    log.info "--lmm-exact              String                Perform GCTA mapping with --fastGWA-lmm-exact algorithm (Default: RUN, option to not run is null)"
    log.info "--sparse_cut             String                Any off-diagonal value in the genetic relatedness matrix greater than this is set to 0 (Default: 0.05)"
    log.info "----------------------------------------------------------------"
    log.info "             -profile simulations USAGE"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "nextflow develop.nf --vcf input_data/elegans/genotypes/WI.20180527.impute.vcf.gz -profile simulations"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "Mandatory arguments:"
    log.info "--simulate_nqtl          File.                 A CSV file with the number of QTL to simulate per phenotype, one value per line (Default is located: input_data/all_species/simulate_nqtl.csv)"
    log.info "--simulate_h2            File                  A CSV file with phenotype heritability, one value per line (Default is located: input_data/all_species/simulate_h2.csv)"
    log.info "Optional arguments:"
    log.info "--simulate_reps          String                The number of replicates to simulate per number of QTL and heritability (Default: 2)"
    log.info "--simulate_maf           File                  A CSV file where each line is a minor allele frequency threshold to test for simulations (Default: input_data/all_species/simulate_maf.csv)"
    log.info "--simulate_eff           File                  A CSV file where each line is an effect size to test for simulations (Default: input_data/all_species/simulate_effect_sizes.csv)"
    log.info "--simulate_strains       File                  A TSV file with two columns: the first is a name for the strain set and the second is a comma-separated strain list without spaces (Default: input_data/all_species/simulate_strains.csv)"
    log.info "--simulate_qtlloc        File                  A BED file with three columns: chromosome name (numeric 1-6), start postion, end postion. The genomic range specified is where markers will be pulled from to simulate QTL (Default: null [which defaults to using the whole genome to randomly simulate a QTL])"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "Optional arguments (General):"
    log.info "--out                    String                Name of folder that will contain the results"
    log.info "--e_mem                  String                Value that corresponds to the amount of memory to allocate for eigen decomposition of chromosomes (DEFAULT = 100)"
    log.info "Optional arguments (Marker):"
    log.info "--sthresh                String                Significance threshold for QTL - Options: BF - for bonferroni correction, EIGEN - for SNV eigen value correction, or another number e.g. 4"
    log.info "--group_qtl              Integer               If two QTL are less than this distance from each other, combine the QTL into one, (DEFAULT = 1000)"
    log.info "--ci_size                Integer               Number of SNVs to the left and right of the peak marker used to define the QTL confidence interval, (DEFAULT = 150)"
    log.info "Optional arguments (Burden):"
    log.info "--freqUpper              Float                 Maximum allele frequency for a variant to be considered for burden mapping, (DEFAULT = 0.05)"
    log.info "--minburden              Interger              Minimum number of strains to have a variant for the variant to be considered for burden mapping, (DEFAULT = 2)"
    log.info "--genes                  String                refFlat file format that contains start and stop genomic coordinates for genes of interest, (DEFAULT = bin/gene_ref_flat.Rda)"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Information describing the stucture of the input files can be located in input_files/README.txt"
    log.info ""
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "Required software packages to be in users path:"
    log.info "BCFtools               v1.9"
    log.info "plink                  v1.9"
    log.info "R-cegwas2              Found on GitHub"
    log.info "R-tidyverse            v1.2.1"
    log.info "R-correlateR           Found on GitHub"
    log.info "R-rrBLUP               v4.6"
    log.info "R-sommer               v3.5"
    log.info "R-RSpectra             v0.13-1"
    log.info "R-ggbeeswarm           v0.6.0"
    log.info "                                               "
    log.info "Options to run without downloading individual software:"
    log.info "(1) On QUEST, a local conda environment and corresponding R packages are already supplied."
    log.info "(2) On QUEST or outside of QUEST, you can also supply a docker/singularity container with -profile mappings_docker. Remember to load singularity first."
    log.info "--------------------------------------------------------"
    exit 1
} else {
    log.info '''
O~~~     O~~                                   O~~ ~~
O~ O~~   O~~                                 O~~    O~~
O~~ O~~  O~~   O~~    O~~~ O~~ O~~    O~~     O~~         O~~~   O~~    O~~ O~~
O~~  O~~ O~~ O~   O~~  O~~  O~  O~~ O~~  O~~    O~~     O~~    O~~  O~~  O~~  O~~
O~~   O~ O~~O~~~~~ O~~ O~~  O~  O~~O~~   O~~       O~~ O~~    O~~   O~~  O~~  O~~
O~~    O~ ~~O~         O~~  O~  O~~O~~   O~~ O~~    O~~ O~~   O~~   O~~  O~~  O~~
O~~      O~~  O~~~~   O~~~  O~  O~~  O~~ O~~~  O~~ ~~     O~~~  O~~ O~~~O~~~  O~~
    '''
log.info ""
log.info "Trait File                              = ${params.traitfile}"
log.info "VCF                                     = ${params.vcf}"
log.info "Significance Threshold                  = ${params.sthresh}"
log.info "Result Directory                        = ${params.out}"
log.info "-----------------------------------------------------------"
log.info "Run mapping profile?                    = ${params.maps}"
log.info "Run simulation profile?                 = ${params.simulate}"
log.info "Run annotation profile?                 = ${params.annotate}"
log.info "-----------------------------------------------------------"
log.info "Git info:                               = $workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""
}



/*
~ ~ ~ > * WORKFLOW
*/

workflow {

    // if no VCF is provided, download the latest version from CeNDR
    if(download_vcf) {
        pull_vcf()

        vcf_file = pull_vcf.out.hard_vcf
        vcf_index = pull_vcf.out.hard_vcf_index
        impute_vcf = pull_vcf.out.impute_vcf
        impute_vcf_index = pull_vcf.out.impute_vcf_index
        ann_file = pull_vcf.out.ann_vcf
    }

    // for mapping
    if(params.maps) {

        // Fix strain names
        Channel.fromPath("${params.traitfile}")
            .combine(Channel.fromPath("${params.data_dir}/isotypes/strain_isotype_lookup.tsv")) | fix_strain_names_bulk        

        traits_to_map = fix_strain_names_bulk.out.fixed_strain_phenotypes
                .flatten()
                .map { file -> tuple(file.baseName.replaceAll(/pr_/,""), file) }

        // Genotype matrix
        pheno_strains = fix_strain_names_bulk.out.phenotyped_strains_to_analyze

        vcf_file.spread(vcf_index)
                .combine(pheno_strains) | vcf_to_geno_matrix

        // EIGEN
        contigs = Channel.from(["I", "II", "III", "IV", "V", "X"])
        contigs.combine(vcf_to_geno_matrix.out) | chrom_eigen_variants
        chrom_eigen_variants.out.collect() | collect_eigen_variants

        // GWAS mapping
        pheno_strains
            .spread(traits_to_map)
            .spread(vcf_file.spread(vcf_index))
            .spread(Channel.fromPath("${params.numeric_chrom}")) | prepare_gcta_files | gcta_grm | gcta_lmm_exact_mapping

        // process GWAS mapping
        traits_to_map
            .spread(collect_eigen_variants.out)
            .spread(vcf_to_geno_matrix.out)
            .combine(Channel.from("${params.p3d}"))
            .combine(Channel.from("${params.sthresh}"))
            .combine(Channel.from("${params.group_qtl}"))
            .combine(Channel.from("${params.ci_size}"))
            .join(gcta_lmm_exact_mapping.out) | gcta_intervals_maps

        // plot
        gcta_intervals_maps.out.to_plots | generate_plots 

        // LD b/w regions
        gcta_intervals_maps.out.to_plots | LD_between_regions


        // summarize all peaks
        peaks = gcta_intervals_maps.out.qtl_peaks
            .collectFile(keepHeader: true, name: "QTL_peaks.tsv", storeDir: "${params.out}/Mapping/Processed")

        // prep LD files
        peaks
            .splitCsv(sep: '\t', skip: 1)
            .join(generate_plots.out.maps_from_plot, by: 2)
            .spread(impute_vcf.spread(impute_vcf_index))
            .spread(pheno_strains)
            .spread(Channel.fromPath("${params.numeric_chrom}")) | prep_ld_files

        //fine mapping
        prep_ld_files.out.finemap_preps
            .combine(ann_file) | gcta_fine_maps

        // divergent regions and haplotypes
        peaks
            .combine(Channel.fromPath("${params.data_dir}/isotypes/divergent_bins.bed"))
            .combine(Channel.fromPath("${params.data_dir}/isotypes/divergent_df_isotype.bed"))
            .combine(Channel.fromPath("${params.data_dir}/isotypes/haplotype_df_isotype.bed"))
            .combine(Channel.fromPath("${params.data_dir}/isotypes/div_isotype_list.txt")) | divergent_and_haplotype

        // generate main html report
        peaks
            .spread(traits_to_map)
            .combine(divergent_and_haplotype.out.div_done)
            //.combine(gcta_fine_maps.out.finemap_done) | html_report_main
            .join(gcta_fine_maps.out.finemap_done, by: 1, remainder: true) | html_report_main

    } else if(params.annotate) {

        // what does annotate do?? just this one process?
        save_dir = "${params.input_data}/${params.species}/annotations"

        Channel.fromPath("${params.script_loc}")
            .combine(save_dir) | update_annotations

    } else if(params.matrix) {

        // only run geno matrix step - and fix isotype names if needed
        Channel.fromPath("${params.strains}")
            .combine(Channel.fromPath("${params.data_dir}/isotypes/strain_isotype_lookup.tsv")) | fix_strain_names_alt
        
        pheno_strains = fix_strain_names_alt.out.phenotyped_strains_to_analyze

        vcf_file.spread(vcf_index)
                .combine(pheno_strains) | vcf_to_geno_matrix

    } else if(params.simulate) {

        // vcf = Channel.fromPath("input_data/elegans/genotypes/WI.20180527.impute.vcf.gz")
        // vcf_index = Channel.fromPath("input_data/elegans/genotypes/WI.20180527.impute.vcf.gz.tbi")

        // for simulations
        File pop_file = new File(params.simulate_strains);

        Channel.from(pop_file.collect { it.tokenize( ' ' ) })
            .map { SM, STRAINS -> [SM, STRAINS] }
            .spread(vcf_file.spread(vcf_index))
            .spread(Channel.fromPath("${params.numeric_chrom}"))
            .spread(Channel.fromPath("${params.simulate_maf}").splitCsv()) | prepare_simulation_files

        // eigen
        contigs = Channel.from(["1", "2", "3", "4", "5", "6"])
        contigs.combine(prepare_simulation_files.out.sim_geno) | chrom_eigen_variants_sims

        chrom_eigen_variants_sims.out.sim_geno_eigen_join
            .groupTuple(by:[0,1,2]).
            join(chrom_eigen_variants_sims.out.sim_geno_meta, by:[0,1,2]) | collect_eigen_variants_sims

        // simulate qtl or genome
        if(params.simulate_qtlloc){

            collect_eigen_variants_sims.out
                .spread(Channel.fromPath("${params.simulate_nqtl}").splitCsv())
                .spread(Channel.fromPath("${params.simulate_qtlloc}"))
                .spread(Channel.fromPath("${params.simulate_eff}").splitCsv())
                .combine(Channel.from(1..params.simulate_reps)) | simulate_effects_loc

            sim_phen_inputs = simulate_effects_loc.out

        } else {

            collect_eigen_variants_sims.out
                .spread(Channel.fromPath("${params.simulate_nqtl}").splitCsv())
                .spread(Channel.fromPath("${params.simulate_eff}").splitCsv())
                .combine(Channel.from(1..params.simulate_reps)) | simulate_effects_genome

            sim_phen_inputs = simulate_effects_genome.out

        }

        sim_phen_inputs
            .spread(Channel.fromPath("${params.simulate_h2}").splitCsv()) | simulate_map_phenotypes

        // simulation mappings
        simulate_map_phenotypes.out.gcta_intervals
            .spread(Channel.from("${params.sthresh}"))
            .spread(Channel.from("${params.group_qtl}"))
            .spread(Channel.from("${params.ci_size}")) | get_gcta_intervals

        // LD from simulations
        get_gcta_intervals.out.processed_gcta | LD_simulated_maps

    }

}

/*
==============================================
~ > *                                    * < ~
~ ~ > *                                * < ~ ~
~ ~ ~ > *  DOWNLOAD VCF FROM CENDR   * < ~ ~ ~
~ ~ > *                                * < ~ ~
~ > *                                    * < ~
==============================================
*/

process pull_vcf {

    tag {"PULLING VCF FROM CeNDR"}

    output:
        path "*hard-filter.isotype.vcf.gz", emit: hard_vcf 
        path "*hard-filter.isotype.vcf.gz.tbi", emit: hard_vcf_index 
        path "*impute.isotype.vcf.gz", emit: impute_vcf 
        path "*impute.isotype.vcf.gz.tbi", emit: impute_vcf_index 
        path "*.strain-annotation*.tsv", emit: ann_vcf

    """
        wget https://storage.googleapis.com/elegansvariation.org/releases/${params.vcf}/variation/WI.${params.vcf}.small.hard-filter.isotype.vcf.gz
        tabix -p vcf WI.${params.vcf}.small.hard-filter.isotype.vcf.gz

        wget https://storage.googleapis.com/elegansvariation.org/releases/${params.vcf}/variation/WI.${params.vcf}.impute.isotype.vcf.gz
        tabix -p vcf WI.${params.vcf}.impute.isotype.vcf.gz

        wget https://storage.googleapis.com/elegansvariation.org/releases/${params.vcf}/variation/WI.${params.vcf}.strain-annotation.bcsq.tsv

    """
}


/*
==============================================
~ > *                                    * < ~
~ ~ > *                                * < ~ ~
~ ~ ~ > *  UPDATE ANNOTATION INPUTS  * < ~ ~ ~
~ ~ > *                                * < ~ ~
~ > *                                    * < ~
==============================================
*/

process update_annotations {

    publishDir "${save_dir}", mode: 'copy'

    input:
        tuple val(gtf_to_refflat), val(save_dir)

    output:
        tuple file("*canonical_geneset.gtf.gz"), file("c_${params.species}_${params.wb_build}_refFlat.txt")

    """
        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/update_annotations.R > update_annotations.R
        Rscript --vanilla update_annotations.R ${params.wb_build} ${params.species} ${gtf_to_refflat}
    """

}   

/*
==============================================================
~ > *                                                    * < ~
~ ~ > *                                                * < ~ ~
~ ~ ~ > *  FIX STRAIN NAMES TO MATCH THOSE ON CENDR  * < ~ ~ ~
~ ~ > *                                                * < ~ ~
~ > *                                                    * < ~
==============================================================
*/

/*
THIS WILL NEED TO BE UPDATED TO HANDLE OTHER SPECIES
*/


process fix_strain_names_bulk {

    tag {"BULK TRAIT"}

    publishDir "${params.out}/Phenotypes", mode: 'copy', pattern: "*pr_*.tsv"
    publishDir "${params.out}/Phenotypes", mode: 'copy', pattern: "strain_issues.txt"

    input:
        tuple file(phenotypes), file(isotype_lookup)

    output:
        path "pr_*.tsv", emit: fixed_strain_phenotypes 
        path "Phenotyped_Strains.txt", emit: phenotyped_strains_to_analyze 
        file("strain_issues.txt")

    """
        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Fix_Isotype_names_bulk.R > Fix_Isotype_names_bulk.R 
        Rscript --vanilla Fix_Isotype_names_bulk.R ${phenotypes} fix $isotype_lookup

        # check to make sure there are more than 40 strains for a mapping.
        if [[ \$(wc -l <Phenotyped_Strains.txt) -le 40 ]]
        then
            echo "ERROR: Please provide at least 40 strains for a GWAS mapping."
            exit 1
        fi

    """

}

process fix_strain_names_alt {

    publishDir "${params.out}/Phenotypes", mode: 'copy', pattern: "*.txt"

    input:
        tuple file(phenotypes), file(isotype_lookup)

    output:
        path "Phenotyped_Strains.txt", emit: phenotyped_strains_to_analyze 
        file("strain_issues.txt")

    """
        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Fix_Isotype_names_alt.R > Fix_Isotype_names_alt.R 
        Rscript --vanilla Fix_Isotype_names_alt.R ${phenotypes} fix $isotype_lookup

    """

}


/*
===========================================================
~ > *                                                 * < ~
~ ~ > *                                             * < ~ ~
~ ~ ~ > *  CONVERT THE VCF TO A GENOTYPE MATRIX   * < ~ ~ ~
~ ~ > *                                             * < ~ ~
~ > *                                                 * < ~
===========================================================
*/

process vcf_to_geno_matrix {

    publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

    cpus 1

    input:
        tuple file(vcf), file(index), file(strains)

    output:
        file("Genotype_Matrix.tsv") 

    """
        bcftools view -S ${strains} ${vcf} |\\
        bcftools filter -i N_MISSING=0 -Oz -o Phenotyped_Strain_VCF.vcf.gz
        tabix -p vcf Phenotyped_Strain_VCF.vcf.gz
        plink --vcf Phenotyped_Strain_VCF.vcf.gz \\
            --snps-only \\
            --biallelic-only \\
            --maf ${params.MAF} \\
            --set-missing-var-ids @:# \\
            --indep-pairwise 50 10 0.8 \\
            --geno \\
            --allow-extra-chr
        awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
        sort -k1,1d -k2,2n > markers.txt
        bcftools query -l Phenotyped_Strain_VCF.vcf.gz |\\
        sort > sorted_samples.txt 
        bcftools view -v snps \\
        -S sorted_samples.txt \\
        -R markers.txt \\
        Phenotyped_Strain_VCF.vcf.gz |\\
        bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' |\\
        sed 's/[[# 0-9]*]//g' |\\
        sed 's/:GT//g' |\\
        sed 's/0|0/-1/g' |\\
        sed 's/1|1/1/g' |\\
        sed 's/0|1/NA/g' |\\
        sed 's/1|0/NA/g' |\\
        sed 's/.|./NA/g'  |\\
        sed 's/0\\/0/-1/g' |\\
        sed 's/1\\/1/1/g'  |\\
        sed 's/0\\/1/NA/g' |\\
        sed 's/1\\/0/NA/g' |\\
        sed 's/.\\/./NA/g' > Genotype_Matrix.tsv
    """

}



/*
============================================================
~ > *                                                  * < ~
~ ~ > *                                              * < ~ ~
~ ~ ~ > *  EIGEN DECOMPOSITION OF GENOTYPE MATRIX  * < ~ ~ ~
~ ~ > *                                              * < ~ ~
~ > *                                                  * < ~
============================================================
*/


/*
------------ Decomposition per chromosome
*/

process chrom_eigen_variants {

    tag { CHROM }

    cpus 6
    memory params.eigen_mem

    input:
        tuple val(CHROM), file(genotypes)


    output:
        file("${CHROM}_independent_snvs.csv")


    """
        cat Genotype_Matrix.tsv |\\
        awk -v chrom="${CHROM}" '{if(\$1 == "CHROM" || \$1 == chrom) print}' > ${CHROM}_gm.tsv
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Get_GenoMatrix_Eigen.R > Get_GenoMatrix_Eigen.R
        Rscript --vanilla Get_GenoMatrix_Eigen.R ${CHROM}_gm.tsv ${CHROM}
    """

}

/*
------------ Sum independent tests for all chromosomes
*/

process collect_eigen_variants {

    publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

    cpus 1

    input:
        file(chrom_tests) 

    output:
        file("total_independent_tests.txt")

    """
        cat *independent_snvs.csv |\\
        grep -v inde |\\
        awk '{s+=\$1}END{print s}' > total_independent_tests.txt
    """

}


/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  RUN GWAS MAPPING  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/

/*
------------ GCTA
*/

process prepare_gcta_files {

    cpus 4

    input:
        tuple file(strains), val(TRAIT), file(traits), file(vcf), file(index), file(num_chroms)

    output:
        tuple val(TRAIT), file("plink_formated_trats.tsv"), file("${TRAIT}.bed"), file("${TRAIT}.bim"), file("${TRAIT}.fam"), file("${TRAIT}.map"), file("${TRAIT}.nosex"), file("${TRAIT}.ped"), file("${TRAIT}.log")

    """
    bcftools annotate --rename-chrs ${num_chroms} ${vcf} |\\
    bcftools view -S ${strains} |\\
    bcftools filter -i N_MISSING=0 -Oz -o renamed_chroms.vcf.gz
    tabix -p vcf renamed_chroms.vcf.gz
    plink --vcf renamed_chroms.vcf.gz \\
    --snps-only \\
    --biallelic-only \\
    --maf ${params.MAF} \\
    --set-missing-var-ids @:# \\
    --indep-pairwise 50 10 0.8 \\
    --geno \\
    --not-chr MtDNA \\
    --allow-extra-chr
    tail -n +2 ${traits} | awk 'BEGIN {OFS="\\t"}; {print \$1, \$1, \$2}' > plink_formated_trats.tsv
    plink --vcf renamed_chroms.vcf.gz \\
    --make-bed \\
    --snps-only \\
    --biallelic-only \\
    --maf ${params.MAF} \\
    --set-missing-var-ids @:# \\
    --extract plink.prune.in \\
    --geno \\
    --recode \\
    --out ${TRAIT} \\
    --allow-extra-chr \\
    --pheno plink_formated_trats.tsv
    """
}

process gcta_grm {

    cpus 4

    input:
        tuple val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log)

    output:
        tuple val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file("${TRAIT}_gcta_grm.grm.bin"), file("${TRAIT}_gcta_grm.grm.id"), file("${TRAIT}_gcta_grm.grm.N.bin"), file("${TRAIT}_heritability.hsq"), file("${TRAIT}_heritability.log"), file("${TRAIT}_gcta_grm_inbred.grm.bin"), file("${TRAIT}_gcta_grm_inbred.grm.id"), file("${TRAIT}_gcta_grm_inbred.grm.N.bin"), file("${TRAIT}_heritability_inbred.hsq"), file("${TRAIT}_heritability_inbred.log")

    when:
        params.maps

    """
    gcta64 --bfile ${TRAIT} --autosome --maf ${params.MAF} --make-grm --out ${TRAIT}_gcta_grm --thread-num 10
    gcta64 --bfile ${TRAIT} --autosome --maf ${params.MAF} --make-grm-inbred --out ${TRAIT}_gcta_grm_inbred --thread-num 10
    gcta64 --grm ${TRAIT}_gcta_grm --pheno plink_formated_trats.tsv --reml --out ${TRAIT}_heritability --thread-num 10
    gcta64 --grm ${TRAIT}_gcta_grm_inbred --pheno plink_formated_trats.tsv --reml --out ${TRAIT}_heritability_inbred --thread-num 10
    """
}


process gcta_lmm_exact_mapping {

    cpus 4

    publishDir "${params.out}/Mapping/Raw", pattern: "*fastGWA", overwrite: true
    publishDir "${params.out}/Mapping/Raw", pattern: "*loco.mlma", overwrite: true

    errorStrategy 'ignore'

    input:
    tuple val(TRAIT), file(traits), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(grm_bin), file(grm_id), file(grm_nbin), file(h2), file(h2log), file(grm_bin_inbred), file(grm_id_inbred), file(grm_nbin_inbred), file(h2_inbred), file(h2log_inbred)

    output:
    tuple val(TRAIT), file("${TRAIT}_lmm-exact_inbred.fastGWA"), file("${TRAIT}_lmm-exact.loco.mlma")


    """
    gcta64 --grm ${TRAIT}_gcta_grm --make-bK-sparse ${params.sparse_cut} --out ${TRAIT}_sparse_grm
    gcta64 --mlma-loco \\
        --grm ${TRAIT}_sparse_grm \\
        --bfile ${TRAIT} \\
        --out ${TRAIT}_lmm-exact \\
        --pheno ${traits} \\
        --maf ${params.maf}
    gcta64 --grm ${TRAIT}_gcta_grm_inbred --make-bK-sparse ${params.sparse_cut} --out ${TRAIT}_sparse_grm_inbred
    gcta64 --fastGWA-lmm-exact \\
        --grm-sparse ${TRAIT}_sparse_grm \\
        --bfile ${TRAIT} \\
        --out ${TRAIT}_lmm-exact_inbred \\
        --pheno ${traits} \\
        --maf ${params.maf}
    """
}



process gcta_intervals_maps {

    publishDir "${params.out}/Mapping/Processed", mode: 'copy', pattern: "*AGGREGATE_mapping.tsv"
    publishDir "${params.out}/Mapping/Processed", mode: 'copy', pattern: "*AGGREGATE_qtl_region.tsv" //would be nice to put all these files per trait into one file

    memory '48 GB'

    input:
        tuple val(TRAIT), file(pheno), file(tests), file(geno), val(P3D), val(sig_thresh), val(qtl_grouping_size), val(qtl_ci_size), file(lmmexact_inbred), file(lmmexact_loco)

    output:
        tuple file(geno), file(pheno), file(tests), val(TRAIT), file("*AGGREGATE_mapping.tsv"), emit: to_plots
        path "*AGGREGATE_qtl_region.tsv", emit: qtl_peaks

    """
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Aggregate_Mappings.R > Aggregate_Mappings.R
    Rscript --vanilla Aggregate_Mappings.R ${lmmexact_loco} ${lmmexact_inbred}

    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Find_Aggregate_Intervals_Maps.R > Find_Aggregate_Intervals_Maps.R
    Rscript --vanilla Find_Aggregate_Intervals_Maps.R ${geno} ${pheno} temp.aggregate.mapping.tsv ${tests} ${qtl_grouping_size} ${qtl_ci_size} ${sig_thresh} ${TRAIT}_AGGREGATE
    """
}


process generate_plots {

    publishDir "${params.out}/Plots/EffectPlots", mode: 'copy', pattern: "*_effect.plot.png"
    publishDir "${params.out}/Plots/ManhattanPlots", mode: 'copy', pattern: "*_manhattan.plot.png"

    input:
        tuple file(geno), file(pheno), file(tests), val(TRAIT), file(aggregate_mapping)

    output:
        tuple file(geno), file(pheno), val(TRAIT), file(aggregate_mapping), emit: maps_from_plot
        file("*.png")

    """
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/pipeline.plotting.mod.R > pipeline.plotting.mod.R
    Rscript --vanilla pipeline.plotting.mod.R ${aggregate_mapping} ${tests}
    """
}


process LD_between_regions {

  publishDir "${params.out}/Mapping/Processed", mode: 'copy', pattern: "*LD_between_QTL_regions.tsv"
  publishDir "${params.out}/Plots/LDPlots", mode: 'copy', pattern: "*_LD.plot.png"

  input:
        tuple file(geno), file(pheno), file(tests), val(TRAIT), file(aggregate_mapping)

  output:
        tuple val(TRAIT), path("*LD_between_QTL_regions.tsv"), file("*.png") optional true
        val TRAIT, emit: linkage_done

  """
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/LD_between_regions.R > LD_between_regions.R 
    Rscript --vanilla LD_between_regions.R ${geno} ${aggregate_mapping} ${TRAIT}
  """
}


/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  RUN FINE MAPPING  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/

/*
------------ Extract QTL interval genotype matrix
*/



process prep_ld_files {

    tag {TRAIT}

    publishDir "${params.out}/Fine_Mappings/Data", mode: 'copy', pattern: "*ROI_Genotype_Matrix.tsv"
    publishDir "${params.out}/Fine_Mappings/Data", mode: 'copy', pattern: "*LD.tsv"

    input:
        tuple val(TRAIT), val(CHROM), val(marker), val(start_pos), val(peak_pos), val(end_pos), val(peak_id), val(h2), file(geno), \
        file(pheno), file(aggregate_mapping), file(imputed_vcf), file(imputed_index), file(phenotype), file(num_chroms)

    output:
        tuple val(TRAIT), file(pheno), file("*ROI_Genotype_Matrix.tsv"), file("*LD.tsv"), file("*.bim"), file("*.bed"), file("*.fam"), emit: finemap_preps

    """
        echo "HELLO"
        cat ${aggregate_mapping} |\\
        awk '\$0 !~ "\\tNA\\t" {print}' |\\
        awk '!seen[\$1,\$12,\$19,\$20,\$21]++' |\\
        awk 'NR>1{print \$1, \$12, \$19, \$20, \$21}' OFS="\\t" > ${TRAIT}_QTL_peaks.tsv
        filename='${TRAIT}_QTL_peaks.tsv'
        echo Start
        while read p; do 
            chromosome=`echo \$p | cut -f1 -d' '`
            trait=`echo \$p | cut -f2 -d' '`
            start_pos=`echo \$p | cut -f3 -d' '`
            peak_pos=`echo \$p | cut -f4 -d' '`
            end_pos=`echo \$p | cut -f5 -d' '`
        
        cat ${phenotype} | awk '\$0 !~ "strain" {print}' | cut -f1 > phenotyped_samples.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S phenotyped_samples.txt |\\
        bcftools filter -i N_MISSING=0 |\\
        bcftools annotate --rename-chrs rename_chromosomes |\\
        awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt
        bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${imputed_vcf} \
        -S phenotyped_samples.txt |\\
        bcftools filter -i N_MISSING=0 -Oz |\\
        bcftools annotate --rename-chrs ${num_chroms} -o finemap.vcf.gz
        plink --vcf finemap.vcf.gz \\
            --snps-only \\
            --maf ${params.MAF} \\
            --biallelic-only \\
            --allow-extra-chr \\
            --set-missing-var-ids @:# \\
            --geno \\
            --make-bed \\
            --recode vcf-iid bgz \\
            --extract \$trait.\$chromosome.\$start_pos.\$end_pos.txt \\
            --out \$trait.\$chromosome.\$start_pos.\$end_pos
        nsnps=`wc -l \$trait.\$chromosome.\$start_pos.\$end_pos.txt | cut -f1 -d' '`
        chrom_num=`cat rename_chromosomes | grep -w \$chromosome | cut -f2 -d' '`
        plink --r2 with-freqs \\
            --allow-extra-chr \\
            --snps-only \\
            --ld-window-r2 0 \\
            --ld-snp \$chrom_num:\$peak_pos \\
            --ld-window \$nsnps \\
            --ld-window-kb 6000 \\
            --chr \$chrom_num \\
            --out \$trait.\$chromosome:\$start_pos-\$end_pos.QTL \\
            --set-missing-var-ids @:# \\
            --vcf \$trait.\$chromosome.\$start_pos.\$end_pos.vcf.gz
        cut \$trait.\$chromosome:\$start_pos-\$end_pos.QTL.ld -f2-10 > \$trait.\$chromosome.\$start_pos.\$end_pos.LD.tsv
        bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' finemap.vcf.gz |\\
            sed 's/[[# 0-9]*]//g' |\\
            sed 's/:GT//g' |\\
            sed 's/0|0/-1/g' |\\
            sed 's/1|1/1/g' |\\
            sed 's/0|1/NA/g' |\\
            sed 's/1|0/NA/g' |\\
            sed 's/.|./NA/g'  |\\
            sed 's/0\\/0/-1/g' |\\
            sed 's/1\\/1/1/g'  |\\
            sed 's/0\\/1/NA/g' |\\
            sed 's/1\\/0/NA/g' |\\
            sed 's/.\\/./NA/g' |\\
            sed 's/^23/X/g' > \$trait.\$chromosome:\$start_pos-\$end_pos.ROI_Genotype_Matrix.tsv
        done < \$filename
    """
}


process gcta_fine_maps {

    publishDir "${params.out}/Fine_Mappings/Data", mode: 'copy', pattern: "*.fastGWA"
    publishDir "${params.out}/Fine_Mappings/Data", mode: 'copy', pattern: "*_genes.tsv"
    publishDir "${params.out}/Fine_Mappings/Plots", mode: 'copy', pattern: "*.pdf"

    memory '48 GB'
    
    //errorStrategy 'ignore'

    input:
        tuple val(TRAIT), file(pheno), file(ROI_geno), file(ROI_LD), file(bim), file(bed), file(fam), file(annotation)

    output:
        tuple file("*.fastGWA"), val(TRAIT), file("*.prLD_df.tsv"), file("*.pdf"), file("*_genes.tsv")
        //val true, emit: finemap_done
        tuple file("*_genes.tsv"), val(TRAIT), emit: finemap_done

    """
    tail -n +2 ${pheno} | awk 'BEGIN {OFS="\\t"}; {print \$1, \$1, \$2}' > plink_finemap_traits.tsv
    for i in *ROI_Genotype_Matrix.tsv;
        do
        chr=`echo \$i | cut -f2 -d "." | cut -f1 -d ":"`
        start=`echo \$i | cut -f2 -d "." | cut -f2 -d ":" | cut -f1 -d "-"`
        stop=`echo \$i | cut -f2 -d "." | cut -f2 -d ":" | cut -f2 -d "-"`
        gcta64 --bfile ${TRAIT}.\$chr.\$start.\$stop --autosome --maf ${params.MAF} --make-grm-inbred --out ${TRAIT}.\$chr.\$start.\$stop.FM_grm_inbred --thread-num 10
        gcta64 --grm ${TRAIT}.\$chr.\$start.\$stop.FM_grm_inbred --make-bK-sparse ${params.sparse_cut} --out ${TRAIT}.\$chr.\$start.\$stop.sparse_FM_grm_inbred
        gcta64 --fastGWA-lmm-exact \\
        --grm-sparse ${TRAIT}.\$chr.\$start.\$stop.sparse_FM_grm_inbred \\
        --bfile ${TRAIT}.\$chr.\$start.\$stop  \\
        --out ${TRAIT}.\$chr.\$start.\$stop.finemap_inbred \\
        --pheno plink_finemap_traits.tsv \\
        --maf ${params.maf}

        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Finemap_QTL_Intervals.R  > Finemap_QTL_Intervals.R 
        Rscript --vanilla Finemap_QTL_Intervals.R  ${TRAIT}.\$chr.\$start.\$stop.finemap_inbred.fastGWA \$i ${TRAIT}.\$chr.\$start.\$stop.LD.tsv
        
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/plot_genes.R  > plot_genes.R 
        Rscript --vanilla plot_genes.R  ${TRAIT}.\$chr.\$start.\$stop.prLD_df.tsv ${pheno} ${params.genes} ${annotation}
        done
    """
}



/*
------ Slice out the QTL region for plotting divergent region and haplotype data.
*/


process divergent_and_haplotype {

  publishDir "${params.out}/Divergent_and_haplotype", mode: 'copy'


  input:
    tuple file("QTL_peaks.tsv"), file("divergent_bins"), file(divergent_df_isotype), file(haplotype_df_isotype), file(div_isotype_list)

  output:
    tuple file("all_QTL_bins.bed"), file("all_QTL_div.bed"), file("haplotype_in_QTL_region.txt"), file("div_isotype_list2.txt") //, emit: div_hap_table
    val true, emit: div_done


  """
  awk NR\\>1 QTL_peaks.tsv | awk -v OFS='\t' '{print \$1,\$4,\$6}' > QTL_region.bed
  bedtools intersect -wa -a ${divergent_bins} -b QTL_region.bed | sort -k1,1 -k2,2n | uniq > all_QTL_bins.bed
  bedtools intersect -a ${divergent_df_isotype} -b QTL_region.bed | sort -k1,1 -k2,2n | uniq > all_QTL_div.bed
  bedtools intersect -a ${haplotype_df_isotype} -b QTL_region.bed -wo | sort -k1,1 -k2,2n | uniq > haplotype_in_QTL_region.txt
  cp ${div_isotype_list} ./div_isotype_list2.txt
  """

}

// generate trait-specific html reports
process html_report_main {

  errorStrategy 'ignore'

  tag {TRAIT}
  memory '16 GB'
  

  publishDir "${params.out}/Reports", mode: 'copy'


  input:
    tuple val(TRAIT), file("QTL_peaks.tsv"), file(pheno), val(div_done), file("genes.tsv")

  output:
    tuple file("NemaScan_Report_*.Rmd"), file("NemaScan_Report_*.html")


  """
    cat "${workflow.projectDir}/bin/NemaScan_Report_main.Rmd" | sed "s/TRAIT_NAME_HOLDER/${TRAIT}/g" > NemaScan_Report_${TRAIT}_main.Rmd 
    cat "${workflow.projectDir}/bin/NemaScan_Report_region_template.Rmd" > NemaScan_Report_region_template.Rmd 
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" > .Rprofile
    Rscript -e "rmarkdown::render('NemaScan_Report_${TRAIT}_main.Rmd', knit_root_dir='${workflow.launchDir}/${params.out}')"
  """
}



/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  GWAS SIMULATIONS  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/


process prepare_simulation_files {

    cpus 4

    input:
        tuple val(strain_set), val(strains), file(vcf), file(index), file(num_chroms), val(MAF)

    output:
        tuple val(strain_set), val(strains), file("TO_SIMS.bed"), file("TO_SIMS.bim"), file("TO_SIMS.fam"), file("TO_SIMS.map"), file("TO_SIMS.nosex"), file("TO_SIMS.ped"), file("TO_SIMS.log"), file("${strain_set}_${MAF}_Genotype_Matrix.tsv"), val(MAF), emit: sim_geno
        tuple val(strain_set), val(strains), val(MAF), file("renamed_chroms.vcf.gz"), file("renamed_chroms.vcf.gz.tbi"), emit: renamed_chrom_vcf_to_ld


    """
    bcftools annotate --rename-chrs rename_chromosomes ${vcf} |\\
    bcftools view -s `echo ${strains} | tr -d '\\n'` |\\
    bcftools filter -i N_MISSING=0 -Oz -o renamed_chroms.vcf.gz
    tabix -p vcf renamed_chroms.vcf.gz
    plink --vcf renamed_chroms.vcf.gz \\
    --snps-only \\
    --biallelic-only \\
    --maf ${MAF} \\
    --set-missing-var-ids @:# \\
    --indep-pairwise 50 10 0.8 \\
    --geno \\
    --not-chr MtDNA \\
    --allow-extra-chr
    plink --vcf renamed_chroms.vcf.gz \\
    --make-bed \\
    --snps-only \\
    --biallelic-only \\
    --maf ${MAF} \\
    --set-missing-var-ids @:# \\
    --extract plink.prune.in \\
    --geno \\
    --recode \\
    --out TO_SIMS \\
    --allow-extra-chr
    awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
    sort -k1,1d -k2,2n > markers.txt
    bcftools query -l renamed_chroms.vcf.gz |\\
    sort > sorted_samples.txt
    bcftools view -v snps \\
    -S sorted_samples.txt \\
    -R markers.txt \\
    renamed_chroms.vcf.gz |\\
    bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' |\\
    sed 's/[[# 0-9]*]//g' |\\
    sed 's/:GT//g' |\\
    sed 's/0|0/-1/g' |\\
    sed 's/1|1/1/g' |\\
    sed 's/0|1/NA/g' |\\
    sed 's/1|0/NA/g' |\\
    sed 's/.|./NA/g'  |\\
    sed 's/0\\/0/-1/g' |\\
    sed 's/1\\/1/1/g'  |\\
    sed 's/0\\/1/NA/g' |\\
    sed 's/1\\/0/NA/g' |\\
    sed 's/.\\/./NA/g' > ${strain_set}_${MAF}_Genotype_Matrix.tsv
    """
}


/*
------------ Decomposition per chromosome
*/

process chrom_eigen_variants_sims {

    tag { CHROM }

    cpus 6
    memory params.eigen_mem

    input:
        tuple val(CHROM), val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), val(MAF)

    output:
        tuple val(strain_set), val(strains), val(MAF), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), emit: sim_geno_meta
        tuple val(strain_set), val(strains), val(MAF), file("${CHROM}_${strain_set}_${MAF}_independent_snvs.csv"), emit: sim_geno_eigen_join


    """
        cat ${geno} |\\
        awk -v chrom="${CHROM}" '{if(\$1 == "CHROM" || \$1 == chrom) print}' > ${CHROM}_gm.tsv
        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Get_GenoMatrix_Eigen.R > Get_GenoMatrix_Eigen.R
        Rscript --vanilla Get_GenoMatrix_Eigen.R ${CHROM}_gm.tsv ${CHROM}
        mv ${CHROM}_independent_snvs.csv ${CHROM}_${strain_set}_${MAF}_independent_snvs.csv
    """

}

/*
------------ Sum independent tests for all chromosomes
*/

process collect_eigen_variants_sims {

    //executor 'local'

    publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

    cpus 1

    input:
        tuple val(strain_set), val(strains), val(MAF), file(tests), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno)


    output:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(sex), file(ped), file(log), file(geno), val(MAF), file("${strain_set}_${MAF}_total_independent_tests.txt")

    """
        cat *independent_snvs.csv |\\
        grep -v inde |\\
        awk '{s+=\$1}END{print s}' > ${strain_set}_${MAF}_total_independent_tests.txt
    """

}


process simulate_effects_loc {

    tag {NQTL}

    cpus 4

    input:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), file(qtl_loc_bed), val(effect_range), val(SIMREP)

    output:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(SIMREP), val(effect_range), file("causal.variants.sim.${NQTL}.${SIMREP}.txt")

    """
        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/create_causal_QTLs.R > create_causal_QTLs.R
        Rscript --vanilla create_causal_QTLs.R ${bim} ${NQTL} ${effect_range} ${qtl_loc_bed}
        mv causal.variants.sim.${NQTL}.txt causal.variants.sim.${NQTL}.${SIMREP}.txt
    """
}


process simulate_effects_genome {

    tag {NQTL}

    cpus 4

    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(effect_range), val(SIMREP)

    output:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(SIMREP), val(effect_range), file("causal.variants.sim.${NQTL}.${SIMREP}.txt")


    """
        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/create_causal_QTLs.R > create_causal_QTLs.R
        Rscript --vanilla create_causal_QTLs.R ${bim} ${NQTL} ${effect_range}
        mv causal.variants.sim.${NQTL}.txt causal.variants.sim.${NQTL}.${SIMREP}.txt
    """
}


process simulate_map_phenotypes {

    tag {"${NQTL} - ${SIMREP} - ${H2} - ${MAF}"}

    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", pattern: "*fastGWA", overwrite: true
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", pattern: "*loco.mlma", overwrite: true
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Phenotypes", pattern: "*.phen", overwrite: true
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Phenotypes", pattern: "*.par", overwrite: true

    cpus 4

    errorStrategy 'ignore'

    input:
        tuple val(strain_set), val(strains), file(bed), file(bim), file(fam), file(map), file(nosex), file(ped), file(log), file(gm), val(MAF), file(n_indep_tests), val(NQTL), val(SIMREP), val(effect_range), file(loci), val(H2)

    output:
        tuple file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.bed"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.bim"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.fam"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.map"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.nosex"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.ped"), file("TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set}.log"), val(NQTL), val(SIMREP), file(loci), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.par"), emit: sim_phen_output
        tuple file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.fastGWA"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.log"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.log"), emit: sim_GCTA_mapping_results
        path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.fastGWA", emit: lmm_exact_inbred_analyze_sims
        path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma", emit: lmm_exact_loco_analyze_sims
        path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen", emit: simphen_analyze_sims
        path "${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.par", emit: simgen_analyze_sims
        tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), val(MAF), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred.fastGWA"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact.loco.mlma"), file("${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen"), emit: gcta_intervals

    """
    gcta64 --bfile TO_SIMS \\
         --simu-qt \\
         --simu-causal-loci ${loci} \\
         --simu-hsq ${H2} \\
         --simu-rep 1 \\
         --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims
    plink --bfile TO_SIMS \\
        --make-bed \\
        --snps-only \\
        --biallelic-only \\
        --maf ${MAF} \\
        --set-missing-var-ids @:# \\
        --geno \\
        --recode \\
        --out TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
        --allow-extra-chr \\
        --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen
    gcta64 --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} --autosome --maf ${MAF} --make-grm --out TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm --thread-num 10
    gcta64 --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} --autosome --maf ${MAF} --make-grm-inbred --out TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm_inbred --thread-num 10
    gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm_inbred --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen --reml --out check_vp --thread-num 10
    vp=`grep Vp check_vp.hsq | head -1 | cut -f2`
    if (( \$(echo "0.00001 > \$vp" |bc -l) ));
    then
    awk '{print \$1, \$2, \$3*1000}' ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen > temp.phen;
    rm ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen
    mv temp.phen ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen
    fi
    gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm --make-bK-sparse ${params.sparse_cut} --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm
    gcta64 --mlma-loco \\
        --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
        --grm ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm \\
        --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact \\
        --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen \\
        --maf ${MAF}
    gcta64 --grm TO_SIMS_${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_gcta_grm_inbred --make-bK-sparse ${params.sparse_cut} --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm_inbred
    gcta64 --fastGWA-lmm-exact \\
        --grm-sparse ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sparse_grm_inbred \\
        --bfile TO_SIMS_${NQTL}_${SIMREP}_${MAF}_${effect_range}_${strain_set} \\
        --out ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_lmm-exact_inbred \\
        --pheno ${NQTL}_${SIMREP}_${H2}_${MAF}_${effect_range}_${strain_set}_sims.phen \\
        --maf ${MAF}
    """
}


process get_gcta_intervals {

    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*processed_aggregate_mapping.tsv"
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*processed_LMM-EXACT-INBRED_mapping.tsv"
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*processed_LMM-EXACT-LOCO_mapping.tsv"
    publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*qtl_region.tsv"

    memory '48 GB'

    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), val(MAF), file(lmmexact_inbred), file(lmmexact_loco), file(phenotypes), val(THRESHOLD), val(QTL_GROUP_SIZE), val(QTL_CI_SIZE)

    output:
    tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), file(phenotypes), val(THRESHOLD), file("*processed_aggregate_mapping.tsv"), file("*processed_LMM-EXACT-INBRED_mapping.tsv"), file("*processed_LMM-EXACT-LOCO_mapping.tsv"), emit: processed_gcta
    tuple val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file("*aggregate_qtl_region.tsv"), emit: gcta_qtl_to_ld
    tuple val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file(loci), file(phenotypes), emit: simulated_phenotypes

    """
        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Aggregate_Mappings.R > Aggregate_Mappings.R
        Rscript --vanilla Aggregate_Mappings.R ${lmmexact_loco} ${lmmexact_inbred}
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Find_Aggregate_Intervals.R > Find_Aggregate_Intervals.R
        Rscript --vanilla Find_Aggregate_Intervals.R ${gm} ${phenotypes} temp.aggregate.mapping.tsv ${n_indep_tests} ${NQTL} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} aggregate
        
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Find_GCTA_Intervals.R > Find_GCTA_Intervals.R
        Rscript --vanilla Find_GCTA_Intervals.R ${gm} ${phenotypes} ${lmmexact_inbred} ${n_indep_tests} ${NQTL} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} LMM-EXACT-INBRED
        
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Find_GCTA_Intervals_LOCO.R > Find_GCTA_Intervals_LOCO.R
        Rscript --vanilla Find_GCTA_Intervals_LOCO.R ${gm} ${phenotypes} ${lmmexact_loco} ${n_indep_tests} ${NQTL} ${SIMREP} ${QTL_GROUP_SIZE} ${QTL_CI_SIZE} ${H2} ${params.maf} ${THRESHOLD} ${strain_set} ${MAF} ${effect_range} LMM-EXACT-LOCO
    """
}

process LD_simulated_maps {

  publishDir "${params.out}/Simulations/${effect_range}/${NQTL}/Mappings", mode: 'copy', pattern: "*LD_between_QTL_regions.tsv"

  errorStrategy 'ignore'
  
  input:
        tuple val(strain_set), val(strains), val(NQTL), val(SIMREP), val(H2), file(loci), file(gm), val(effect_range), file(n_indep_tests), file(phenotypes), val(THRESHOLD), file(aggregate_mapping), file(inbred_mapping), file(loco_mapping)

  output:
        path("*LD_between_QTL_regions.tsv") optional true

  """
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/LD_between_regions_sims.R > LD_between_regions_sims.R 
    Rscript --vanilla LD_between_regions_sims.R ${gm} ${aggregate_mapping} ${NQTL} ${SIMREP} ${H2} ${params.maf} ${effect_range} ${strain_set}
  """
}

/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *         LD         * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
    gcta_qtl_to_ld
        .combine(emmma_qtl_to_ld, by: [0,1,2,3,4,5,6])
        .combine(simulated_phenotypes, by: [0,1,2,3,4,5,6])
        .combine(renamed_chrom_vcf_to_ld, by: [0,1,2])
        .into{peaks_to_ld}
    process extract_qtl_ld {
    cpus 1
    input:
    set val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file(lmmexact_peaks), file(lmmexact_inbred_peaks), file(emma_peaks), file(loci), file(phenotypes), file(vcf), file(vcfindex) from peaks_to_ld
    output:
    set val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file(lmmexact_peaks), file(lmmexact_inbred_peaks), file(emma_peaks), file(loci), file(phenotypes), file(vcf), file(vcfindex), file("*ROI_Genotype_Matrix.tsv"), file("*LD.tsv") into LD_files_to_finemap
    """
    cat *_qtl_region.tsv |\\
    awk '\$0 !~ "CHROM" {print}' > QTL_peaks.tsv
    filename='QTL_peaks.tsv'
    echo Start
    while read p; do
        chromosome=`echo \$p | cut -f1 -d' '`
        trait=`echo \$p | cut -f3 -d' '`
        start_pos=`echo \$p | cut -f4 -d' '`
        peak_pos=`echo \$p | cut -f5 -d' '`
        end_pos=`echo \$p | cut -f6 -d' '`
        map_algo=`echo \$p | cut -f8 -d' '`
    bcftools view --regions \$chromosome:\$start_pos-\$end_pos -s ${strains} ${vcf} |\\
    bcftools filter -i N_MISSING=0 |\\
    awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt
    bcftools view --regions \$chromosome:\$start_pos-\$end_pos -s ${strains} ${vcf} |\\
    bcftools filter -i N_MISSING=0 -Oz -o finemap.vcf.gz
    plink --vcf finemap.vcf.gz \\
        --snps-only \\
        --maf ${MAF} \\
        --biallelic-only \\
        --allow-extra-chr \\
        --set-missing-var-ids @:# \\
        --geno \\
        --make-bed \\
        --recode vcf-iid bgz \\
        --extract \$trait.\$chromosome.\$start_pos.\$end_pos.txt \\
        --out \$trait.\$chromosome.\$start_pos.\$end_pos
    nsnps=`wc -l \$trait.\$chromosome.\$start_pos.\$end_pos.txt | cut -f1 -d' '`
    plink --r2 with-freqs \\
        --allow-extra-chr \\
        --snps-only \\
        --ld-window-r2 0 \\
        --ld-snp \$chromosome:\$peak_pos \\
        --ld-window \$nsnps \\
        --ld-window-kb 6000 \\
        --chr \$chromosome \\
        --out \$trait.\$chromosome:\$start_pos-\$end_pos.QTL \\
        --set-missing-var-ids @:# \\
        --vcf \$trait.\$chromosome.\$start_pos.\$end_pos.vcf.gz
    sed 's/  \\t/g' \$trait.\$chromosome:\$start_pos-\$end_pos.QTL.ld |\\
    cut -f2-10 |\\
    sed 's/^23/X/g' | sed 's/\\t23\\t/\\tX\\t/g' > \$trait.\$chromosome.\$start_pos.\$end_pos.\$map_algo.LD.tsv
    bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' finemap.vcf.gz |\\
        sed 's/[[# 0-9]*]//g' |\\
        sed 's/:GT//g' |\\
        sed 's/0|0/-1/g' |\\
        sed 's/1|1/1/g' |\\
        sed 's/0|1/NA/g' |\\
        sed 's/1|0/NA/g' |\\
        sed 's/.|./NA/g'  |\\
        sed 's/0\\/0/-1/g' |\\
        sed 's/1\\/1/1/g'  |\\
        sed 's/0\\/1/NA/g' |\\
        sed 's/1\\/0/NA/g' |\\
        sed 's/.\\/./NA/g' |\\
        sed 's/^23/X/g' > \$trait.\$chromosome:\$start_pos-\$end_pos.\$map_algo.ROI_Genotype_Matrix.tsv
    done < \$filename
    """
}
*/

/*
------------ Run fine mapping
process sim_fine_maps {
    echo 'true'
    input:
        set val(strain_set), val(strains), val(MAF), val(NQTL), val(SIMREP), val(H2), val(effect_range), file(lmmexact_peaks), file(lmmexact_inbred_peaks), file(emma_peaks), file(loci), file(phenotypes), file(vcf), file(vcfindex), file(roi_geno_matrix), file(roi_ld) from LD_files_to_finemap
    output:
    """
        for i in *ROI_Genotype_Matrix.tsv;
        do
            echo \$i
        done
    """
}
*/

/*
=====================================
~ > *                           * < ~
~ ~ > *                       * < ~ ~
~ ~ ~ > *  GENERATE REPORT  * < ~ ~ ~
~ ~ > *                       * < ~ ~
~ > *                           * < ~
=====================================
*/

workflow.onComplete {

    summary = """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
    { Parameters }
    ---------------------------
    Phenotype File                          = ${params.traitfile}
    VCF                                     = ${params.vcf}
    Significance Threshold                  = ${params.sthresh}
    P3D                                     = ${params.p3d}
    Max AF for Burden Mapping               = ${params.freqUpper}
    Min Strains with Variant for Burden     = ${params.minburden}
    Threshold for grouping QTL              = ${params.group_qtl}
    Number of SNVs to define CI             = ${params.ci_size}
    Eigen Memory allocation                 = ${params.eigen_mem}
    Path to R libraries.                    = ${params.R_libpath}
    Mapping                                 = ${params.maps}
    Simulation                              = ${params.simulate}
    Simulate QTL effects                    = ${params.simulate_qtlloc}
    Annotation                              = ${params.annotate}
    Result Directory                        = ${params.out}
    """

    println summary

    def outlog = new File("${params.out}/log.txt")
    outlog.newWriter().withWriter {
        outlog << param_summary
        outlog << summary
    }
}

    // mail summary
    //if (params.email) {
    //    ['mail', '-s', 'cegwas2-nf', params.email].execute() << summary
    //}




