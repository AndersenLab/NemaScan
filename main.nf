#! usr/bin/env nextflow

if( !nextflow.version.matches('>23.0') ) {
    println "This workflow requires Nextflow version 23.0 or greater -- You are running version $nextflow.version"
    if ( !params.matches("Local") ) {
        println "On ${params.platform}, you can use `module load python/${params.anaconda}; source activate ${params.softwareDir}/conda_envs/nf23_env`"
    } else {
        println "Locally, you can create and activate a conda environment with 'nextflow>=23.0'"
    }
    exit 1
}

nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )

/*
~ ~ ~ > * Parameters setup - GENERAL
*/

params.bin_dir = "${workflow.projectDir}/bin" // this is different for gcp
params.data_dir = "${workflow.projectDir}/input_data" // this is different for gcp
params.out = "Analysis_Results-${date}"
// params.algorithm = 'inbred' //options: inbred, loco - now run both


/*
~ ~ ~ > * Parameters setup - MAPPING
*/
params.genes = "${params.data_dir}/${params.species}/annotations/${params.species}.gff"

// mediation only with c_elegans
if(params.species == "c_briggsae" || params.species == "c_tropicalis") {
    med = false
} else {
    med = params.mediation
}

// VCF parameters
if(params.debug) {
    println """
        *** Using debug mode ***
    """
    // debug for now with small vcf
    params.vcf = "${params.species}.test.vcf.gz"
    params.traitfile = "${params.data_dir}/${params.species}/phenotypes/test_pheno.tsv"
    
    vcf_file = Channel.fromPath("${params.data_dir}/${params.species}/genotypes/${params.vcf}")
    vcf_index = Channel.fromPath("${params.data_dir}/${params.species}/genotypes/${params.vcf}.tbi")
    
    // debug can use same vcf for impute and normal
    impute_file = "${params.species}.test.vcf.gz" // just to print out for reference
    impute_vcf = Channel.fromPath("${params.data_dir}/${params.species}/genotypes/${params.vcf}")
    impute_vcf_index = Channel.fromPath("${params.data_dir}/${params.species}/genotypes/${params.vcf}.tbi")
    
    ann_file = Channel.fromPath("${params.data_dir}/${params.species}/genotypes/WI.330_TEST.strain-annotation.tsv")

    // for genomatrix profile
    params.strains = "${params.data_dir}/${params.species}/phenotypes/strain_file.tsv"
    download_vcf = false
} else if(params.gcp) { 
    // use the data directly from google on gcp - switch to elegansvariation.org for now?
    // vcf_file = Channel.fromPath("gs://cendr-site-public-bucket/dataset_release/${params.species}/${params.vcf}/variation/WI.${params.vcf}.small.hard-filter.isotype.vcf.gz")
    // vcf_index = Channel.fromPath("gs://cendr-site-public-bucket/dataset_release/${params.species}/${params.vcf}/variation/WI.${params.vcf}.small.hard-filter.isotype.vcf.gz.tbi")

    vcf_file = Channel.fromPath("gs://caendr-site-public-bucket/dataset_release/${params.species}/${params.vcf}/variation/WI.${params.vcf}.small.hard-filter.isotype.vcf.gz")
    vcf_index = Channel.fromPath("gs://caendr-site-public-bucket/dataset_release/${params.species}/${params.vcf}/variation/WI.${params.vcf}.small.hard-filter.isotype.vcf.gz.tbi")

    impute_file = "WI.${params.vcf}.impute.isotype.vcf.gz" // just to print out for reference
    // impute_vcf = Channel.fromPath("gs://cendr-site-public-bucket/dataset_release/${params.species}/${params.vcf}/variation/WI.${params.vcf}.impute.isotype.vcf.gz")
    // impute_vcf_index = Channel.fromPath("gs://cendr-site-public-bucket/dataset_release/${params.species}/${params.vcf}/variation/WI.${params.vcf}.impute.isotype.vcf.gz.tbi")

    impute_vcf = Channel.fromPath("gs://caendr-site-public-bucket/dataset_release/${params.species}/${params.vcf}/variation/WI.${params.vcf}.impute.isotype.vcf.gz")
    impute_vcf_index = Channel.fromPath("gs://caendr-site-public-bucket/dataset_release/${params.species}/${params.vcf}/variation/WI.${params.vcf}.impute.isotype.vcf.gz.tbi")

    // ann_file = Channel.fromPath("gs://cendr-site-public-bucket/dataset_release/${params.species}/${params.vcf}/variation/WI.${params.vcf}.strain-annotation.tsv")
    ann_file = Channel.fromPath("gs://caendr-site-public-bucket/dataset_release/${params.species}/${params.vcf}/variation/WI.${params.vcf}.strain-annotation.tsv")

    params.strains = "input_data/${params.species}/phenotypes/strain_file.tsv"
    download_vcf = false
} else if(!params.vcf) {
    // if there is no VCF date provided, pull the latest vcf from caendr.
    params.vcf = "20231213"
    vcf_file = "20231213 - CaeNDR"
    vcf_index = "20231213 - CaeNDR"
    impute_file = "20231213 - CaeNDR"
    download_vcf = true
} else {
    download_vcf = false
    // Check that params.vcf is valid
    if("${params.vcf}" == "20231213" || "${params.vcf}" == "20220216" || "${params.vcf}" == "20210121" || "${params.vcf}" == "20200815" || "${params.vcf}" == "20180527" || "${params.vcf}" == "20170531" || "${params.vcf}" == "20210901" || "${params.vcf}" == "20210803") {
        // if("${params.vcf}" in ["20210121", "20200815", "20180527", "20170531", "20210901"]) {
        // check to make sure 20210901 is tropicalis
        if("${params.vcf}" == "20210901") {
            if("${params.species}" == "c_elegans" || "${params.species}" == "c_briggsae") {
                println """
                Error: VCF file (${params.vcf}) does not match species ${params.species} (should be c_tropicalis). Please enter a new vcf date or a new species to continue.
                """
                System.exit(1)
            }
        }
        // check to make sure vcf matches species for briggsae
        if("${params.vcf}" == "20210803") {
            if("${params.species}" == "c_elegans" || "${params.species}" == "c_tropicalis") {
                println """
                Error: VCF file (${params.vcf}) does not match species ${params.species} (should be c_briggsae). Please enter a new vcf date or a new species to continue.
                """
                System.exit(1)
            }
        }
        // check to make sure vcf matches species for elegans
        if("${params.vcf}" == "20231213" || "${params.vcf}" == "20220216" || "${params.vcf}" == "20210121" || "${params.vcf}" == "20200815" || "${params.vcf}" == "20180527" || "${params.vcf}" == "20170531") {
            if("${params.species}" == "c_briggsae" || "${params.species}" == "c_tropicalis") {
                println """
                Error: VCF file (${params.vcf}) does not match species ${params.species} (should be c_elegans). Please enter a new vcf date or a new species to continue.
                """
                System.exit(1)
            }
        }
        // use the vcf data from QUEST when a caendr date is provided
        vcf_file = Channel.fromPath("${params.dataDir}/${params.species}/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.small.hard-filter.isotype.vcf.gz")
        vcf_index = Channel.fromPath("${params.dataDir}/${params.species}/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.small.hard-filter.isotype.vcf.gz.tbi")


        impute_file = "WI.${params.vcf}.impute.isotype.vcf.gz" // just to print out for reference
        impute_vcf = Channel.fromPath("${params.dataDir}/${params.species}/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.impute.isotype.vcf.gz")
        impute_vcf_index = Channel.fromPath("${params.dataDir}/${params.species}/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.impute.isotype.vcf.gz.tbi")

        // check if caendr release date is before 20210121, use snpeff annotation
        if("${params.vcf}" == "20200815" || "${params.vcf}" == "20180527" || "${params.vcf}" == "20170531") {
            println "WARNING: Using snpeff annotation. To use BCSQ annotation, please use a newer vcf (2021 or later)"
            ann_file = Channel.fromPath("${params.dataDir}/${params.species}/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.strain-annotation.snpeff.tsv")
        } else {
            ann_file = Channel.fromPath("${parms.dataDir}/${params.species}/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.strain-annotation.tsv")
        }
    } else {
        // check that vcf file exists, if it does, use it. If not, throw error
        if (!file("${params.vcf}").exists()) {
            println """
            Error: VCF file (${params.vcf}) does not exist. Please provide a valid filepath or a valid CaeNDR release date (i.e. 20210121)
            """
            System.exit(1)
        } else {
            // if it DOES exist
            println """
            WARNING: Using a non-CaeNDR VCF for analysis. Same VCF will be used for both GWA and fine mapping. 
            """
            vcf_file = Channel.fromPath("${params.vcf}")
            vcf_index = Channel.fromPath("${params.vcf}.tbi")

            impute_file = "${params.vcf}" // just to print out for reference
            impute_vcf = Channel.fromPath("${params.vcf}")
            impute_vcf_index = Channel.fromPath("${params.vcf}.tbi")

            //choose default caendr date based on species for ann_file
            if(params.species == "c_elegans") {
                default_date = "20231213"
            } else if(params.species == "c_briggsae") {
                default_date = "20231213"
            } else {
                default_date = "20231213"
            }

            // this does not work for another species...
            ann_file = Channel.fromPath("${params.dataDir}/${params.species}/WI/variation/${default_date}/vcf/WI.${default_date}.strain-annotation.tsv")
        }
    }
}

if (params.matrix || params.mapping){
    simulation = false
} else {
    simulation = true
}


if (params.help) {
    log.info '''
O~~~     O~~                                      O~ O~~
O~ O~~   O~~                                    O~~    O~~
O~~ O~~  O~~    O~~   O~~~ O~~ O~~     O~~       O~~          O~~~    O~~     O~~ O~~
O~~  O~~ O~~  O~   O~  O~~  O~  O~~  O~~  O~~      O~~      O~~     O~~  O~~   O~~  O~~
O~~   O~ O~~ O~~~~~ O~ O~~  O~  O~~ O~~   O~~         O~~  O~~     O~~   O~~   O~~  O~~
O~~    O~O~~ O~        O~~  O~  O~~ O~~   O~~   O~~    O~~  O~~    O~~   O~~   O~~  O~~
O~~      O~~   O~~~~  O~~~  O~  O~~   O~~ O~~~    O~ O~~      O~~~   O~~ O~~~ O~~~  O~~
    '''
    log.info "----------------------------------------------------------------"
    log.info "                      USAGE                                     "
    log.info "----------------------------------------------------------------"
    log.info ""
    log.info "nextflow main.nf --debug"
    log.info "nextflow main.nf --traitfile input_data/${params.species}/phenotypes/PC1.tsv --vcf 20231213"
    log.info ""
    log.info "Profiles available:"
    log.info "standard              Profile                Perform selected analysis on Rockfish (default simulation)"
    log.info "rockfish              Profile                Perform selected analysis on Rockfish (default simulation)"
    log.info "quest                 Profile                Perform selected analysis on QUEST (default simulation)"
    log.info "gcp                   Profile                Perform selected analysis on GCP (default GWA mappings)"
    log.info "local                 Profile                Perform selected analysis using docker on local machine"
    log.info "----------------------------------------------------------------"
    log.info "Optional arguments (General):"
    log.info "--out                    String                Name of folder that will contain the results"
    log.info "Optional arguments (Marker):"
    log.info "--sthresh                String                Significance threshold for QTL - Options: BF - for bonferroni correction, EIGEN - for SNV eigen value correction, or another number e.g. 4"
    log.info "--group_qtl              Integer               If two QTL are less than this distance from each other, combine the QTL into one, (DEFAULT = 1000)"
    log.info "--ci_size                Integer               Number of SNVs to the left and right of the peak marker used to define the QTL confidence interval, (DEFAULT = 150)"
    log.info ""    
    log.info "Information describing the stucture of the input files can be located in input_files/README.txt"
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "             for simulation (default)"
    log.info "----------------------------------------------------------------"
    log.info "nextflow main.nf --vcf 20231213"
    log.info "----------------------------------------------------------------"
    log.info "Mandatory arguments:"
    log.info "--simulate_nqtl          File.                 A CSV file with the number of QTL to simulate per phenotype, one value per line (Default is located: input_data/all_species/simulate_nqtl.csv)"
    log.info "--simulate_h2            File                  A CSV file with phenotype heritability, one value per line (Default is located: input_data/all_species/simulate_h2.csv)"
    log.info "Optional arguments:"
    log.info "--simulate_reps          String                The number of replicates to simulate per number of QTL and heritability (Default: 2)"
    log.info "--simulate_maf           File                  A CSV file where each line is a minor allele frequency threshold to test for simulations (Default: input_data/all_species/simulate_maf.csv)"
    log.info "--simulate_eff           File                  A CSV file where each line is an effect size range (e.g. 0.2-0.3) to test for simulations (Default: input_data/all_species/simulate_effect_sizes.csv)"
    log.info "--simulate_strains       File                  A TSV file with two columns: the first is a name for the strain set and the second is a comma-separated strain list without spaces (Default: input_data/all_species/simulate_strains.csv)"
    log.info "--simulate_qtlloc        File                  A BED file with three columns: chromosome name (numeric 1-6), start postion, end postion. The genomic range specified is where markers will be pulled from to simulate QTL (Default: null [which defaults to using the whole genome to randomly simulate a QTL])"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "             for GWAS mappings (--mapping)"
    log.info "----------------------------------------------------------------"
    log.info "nextflow main.nf --vcf 20231213 --traitfile input_data/${params.species}/phenotypes/PC1.tsv --mapping"
    log.info "----------------------------------------------------------------"
    log.info "Mandatory arguments:"
    log.info "--traitfile              String                Name of file that contains phenotypes. File should be tab-delimited with the columns: strain trait1 trait2 ..."
    log.info "--vcf                    String                Generally a CaeNDR release date (i.e. 20231213). Can also provide a user-specified VCF with index in same folder."
    log.info "Optional arguments:"
    log.info "--MAF, --maf             String                Minimum minor allele frequency to use for single-marker mapping (Default: 0.05)"
    log.info "--lmm                    String                Perform GCTA mapping with --fastGWA-lmm algorithm (Default: RUN, option to not run is null)"
    log.info "--lmm-exact              String                Perform GCTA mapping with --fastGWA-lmm-exact algorithm (Default: RUN, option to not run is null)"
    log.info "--sparse_cut             String                Any off-diagonal value in the genetic relatedness matrix greater than this is set to 0 (Default: 0.05)"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "             for vcf to geno matrix (--matrix)"
    log.info "----------------------------------------------------------------"
    log.info "nextflow main.nf --vcf 20210121 --matrix"
    log.info "----------------------------------------------------------------"
    exit 1
} else {
    log.info '''
O~~~     O~~                                      O~ O~~
O~ O~~   O~~                                    O~~    O~~
O~~ O~~  O~~    O~~   O~~~ O~~ O~~     O~~       O~~          O~~~    O~~     O~~ O~~
O~~  O~~ O~~  O~   O~  O~~  O~  O~~  O~~  O~~      O~~      O~~     O~~  O~~   O~~  O~~
O~~   O~ O~~ O~~~~~ O~ O~~  O~  O~~ O~~   O~~         O~~  O~~     O~~   O~~   O~~  O~~
O~~    O~O~~ O~        O~~  O~  O~~ O~~   O~~   O~~    O~~  O~~    O~~   O~~   O~~  O~~
O~~      O~~   O~~~~  O~~~  O~  O~~   O~~ O~~~    O~ O~~      O~~~   O~~ O~~~ O~~~  O~~
'''
log.info ""
log.info "Trait File                              = ${params.traitfile}"
log.info "Strain File                             = ${params.strains}"
log.info "Species                                 = ${params.species}"
log.info ""
log.info "VCF                                     = ${params.vcf}"
log.info "Impute VCF                              = ${impute_file}"
log.info ""
log.info "Significance Threshold                  = ${params.sthresh}"
log.info "Result Directory                        = ${params.out}"
log.info "Minor allele frequency                  = ${params.maf}"
log.info "Mediation run?                          = ${med}"
log.info ""
}

// Includes
include {pull_vcf; fix_strain_names_bulk; fix_strain_names_alt; vcf_to_geno_matrix; chrom_eigen_variants; collect_eigen_variants} from './modules/setup.nf'
include {prepare_gcta_files; gcta_grm; gcta_lmm_exact_mapping; gcta_lmm_exact_mapping_nopca; gcta_intervals_maps} from './modules/mapping.nf'
include {summarize_mapping; generate_plots; LD_between_regions; divergent_and_haplotype} from './modules/post-mapping.nf'
include {gcta_fine_maps; html_report_main} from './modules/post-mapping.nf'
include {prep_ld_files as prep_ld_files_inbred; prep_ld_files as prep_ld_files_loco} from './modules/post-mapping.nf'
include {mediation_data; multi_mediation; simple_mediation; summary_mediation} from './modules/mediation.nf'
include {prepare_simulation_files; chrom_eigen_variants_sims; collect_eigen_variants_sims; simulate_effects_loc; simulate_effects_genome; simulate_map_phenotypes; get_gcta_intervals; assess_sims_INBRED; assess_sims_LOCO} from './modules/simulations.nf'

/*
~ ~ ~ > * WORKFLOW
*/
workflow {

    // if no VCF is provided, download the latest version from CaeNDR
    if(download_vcf) {
        pull_vcf()

        vcf_file = pull_vcf.out.hard_vcf
        vcf_index = pull_vcf.out.hard_vcf_index
        impute_vcf = pull_vcf.out.impute_vcf
        impute_vcf_index = pull_vcf.out.impute_vcf_index
        ann_file = pull_vcf.out.ann_vcf
    }

    // for mapping
    if(params.mapping) {

        // Fix strain names
         Channel.fromPath("${params.traitfile}")
                .combine(Channel.fromPath("${params.data_dir}/${params.species}/isotypes/strain_isotype_lookup.tsv"))
                .combine(Channel.fromPath("${params.bin_dir}/Fix_Isotype_names_bulk.R"))
                .combine(Channel.from("${params.fix}")) | fix_strain_names_bulk
        traits_to_map = fix_strain_names_bulk.out.fixed_strain_phenotypes
                .flatten()
                .map { file -> tuple(file.baseName.replaceAll(/pr_/,""), file) }

        // Genotype matrix
        pheno_strains = fix_strain_names_bulk.out.phenotyped_strains_to_analyze

        vcf_file.combine(vcf_index)
                .combine(pheno_strains) | vcf_to_geno_matrix

        // EIGEN
        contigs = Channel.from(["I", "II", "III", "IV", "V", "X"])
        contigs.combine(vcf_to_geno_matrix.out)
                .combine(Channel.fromPath("${params.bin_dir}/Get_GenoMatrix_Eigen.R")) | chrom_eigen_variants
        chrom_eigen_variants.out.collect() | collect_eigen_variants

        // GWAS mapping
        if(params.pca) {
            mapping_output = pheno_strains
                .combine(traits_to_map)
                .combine(vcf_file.combine(vcf_index))
                .combine(Channel.fromPath("${params.data_dir}/all_species/rename_chromosomes")) | prepare_gcta_files | gcta_grm | gcta_lmm_exact_mapping
        } else {
            mapping_output = pheno_strains
                .combine(traits_to_map)
                .combine(vcf_file.combine(vcf_index))
                .combine(Channel.fromPath("${params.data_dir}/all_species/rename_chromosomes")) | prepare_gcta_files | gcta_grm | gcta_lmm_exact_mapping_nopca
        }
        
        // process GWAS mapping
        traits_to_map
            .combine(collect_eigen_variants.out)
            .combine(vcf_to_geno_matrix.out)
            .combine(Channel.from("${params.p3d}"))
            .combine(Channel.from("${params.sthresh}"))
            .combine(Channel.from("${params.group_qtl}"))
            .combine(Channel.from("${params.ci_size}"))
            .join(mapping_output)
            .combine(Channel.fromPath("${params.bin_dir}/Find_Aggregate_Intervals_Maps.R")) | gcta_intervals_maps

        // plot
        gcta_intervals_maps.out.maps_to_plot
            .combine(Channel.fromPath("${params.bin_dir}/pipeline.plotting.mod.R")) | generate_plots 

        // LD b/w regions
        gcta_intervals_maps.out.maps_to_plot
            .combine(Channel.fromPath("${params.bin_dir}/LD_between_regions.R")) | LD_between_regions

        // summarize all peaks
        peaks_inbred = gcta_intervals_maps.out.qtl_peaks_inbred
            .collectFile(keepHeader: true, name: "QTL_peaks_inbred.tsv", storeDir: "${params.out}/INBRED/Mapping/Processed")

        peaks_loco = gcta_intervals_maps.out.qtl_peaks_loco
            .collectFile(keepHeader: true, name: "QTL_peaks_loco.tsv", storeDir: "${params.out}/LOCO/Mapping/Processed")

        peaks_inbred
            .combine(peaks_loco)
            .combine(Channel.fromPath("${params.data_dir}/${params.species}/genotypes/${params.species}_chr_lengths.tsv"))
            .combine(Channel.fromPath("${params.bin_dir}/summarize_mappings.R")) | summarize_mapping

        // // run mediation with gaotian's eqtl
        if(med) {

            File transcripteqtl_all = new File("${params.data_dir}/${params.species}/phenotypes/expression/eQTL6545forMed.tsv")
            transcript_eqtl = transcripteqtl_all.getAbsolutePath()

            traits_to_mediate = fix_strain_names_bulk.out.fixed_strain_phenotypes
                .flatten()
                .map { file -> tuple(file.baseName.replaceAll(/pr_/,""), file) }

            // combine inbred and loco qtl peaks for mediation
            peaks_inbred
                .splitCsv(sep: '\t', skip: 1)
                .map { tch,marker,logPvalue,TRAIT,tstart,tpeak,tend,peak_id,h2 -> [TRAIT,tch,tstart,tpeak,tend,logPvalue,peak_id,h2,marker] }
                .combine(Channel.of("inbred"))
                .mix(peaks_loco
                    .splitCsv(sep: '\t', skip: 1)
                    .map { tch,marker,logPvalue,TRAIT,tstart,tpeak,tend,peak_id,h2  -> [TRAIT,tch,tstart,tpeak,tend,logPvalue,peak_id,h2,marker] }
                    .combine(Channel.of("loco")))
                .combine(traits_to_mediate, by: 0)
                .combine(Channel.of(transcript_eqtl))
                .combine(Channel.fromPath("${params.bin_dir}/mediaton_input.R")) | mediation_data

            mediation_data.out
                .combine(vcf_to_geno_matrix.out)
                .combine(Channel.fromPath("${params.data_dir}/${params.species}/phenotypes/expression/tx5291exp_st207.tsv"))
                .combine(Channel.fromPath("${params.bin_dir}/multi_mediation.R")) | multi_mediation

            multi_mediation.out.eQTL_gene
                 .splitCsv(sep: '\t')
                 .combine(mediation_data.out, by: [0,1,2,3])
                 .combine(vcf_to_geno_matrix.out) 
                 .combine(Channel.fromPath("${params.data_dir}/${params.species}/phenotypes/expression/tx5291exp_st207.tsv"))
                 .combine(Channel.fromPath("${params.bin_dir}/simple_mediation.R")) | simple_mediation

            multi_mediation.out.result_multi_mediate
                .groupTuple(by: [0,1])
                .join(simple_mediation.out.groupTuple(by: [0,1]), by: [0,1], remainder: true)
                .combine(Channel.fromPath("${params.bin_dir}/summary_mediation.R")) | summary_mediation
        }


        // easiest case: don't run finemap, divergent, or html if params.finemap = false
        if(params.finemap) {
            // prep LD files
            peaks_inbred
                .splitCsv(sep: '\t', skip: 1)
                .combine(Channel.of("inbred"))
                .join(generate_plots.out.maps_from_plot_inbred, by: 3)
                .combine(impute_vcf.combine(impute_vcf_index))
                .combine(pheno_strains)
                .combine(Channel.fromPath("${params.data_dir}/all_species/rename_chromosomes")) | prep_ld_files_inbred
            peaks_loco
                .splitCsv(sep: '\t', skip: 1)
                .combine(Channel.of("loco"))
                .join(generate_plots.out.maps_from_plot_loco, by: 3)
                .combine(impute_vcf.combine(impute_vcf_index))
                .combine(pheno_strains)
                .combine(Channel.fromPath("${params.data_dir}/all_species/rename_chromosomes")) | prep_ld_files_loco

            //fine mapping
            prep_inbred = prep_ld_files_inbred.out.finemap_preps
                .join(gcta_grm.out)
            prep_loco = prep_ld_files_loco.out.finemap_preps
                .join(gcta_grm.out)

            prep_inbred.mix(prep_loco)
                .combine(ann_file)
                .combine(Channel.fromPath("${params.genes}"))
                .combine(Channel.fromPath("${params.bin_dir}/Finemap_QTL_Intervals.R"))
                .combine(Channel.fromPath("${params.bin_dir}/plot_genes.R")) | gcta_fine_maps


            // divergent regions and haplotypes
            // only for elegans right now
            if(params.species == "c_elegans") {
                peaks_inbred
                    .combine(Channel.of("inbred"))
                    .mix(
                        peaks_loco
                        .combine(Channel.of("loco")))
                    .combine(Channel.fromPath("${params.data_dir}/${params.species}/isotypes/divergent_bins.bed"))
                    .combine(Channel.fromPath("${params.data_dir}/${params.species}/isotypes/divergent_df_isotype.bed"))
                    .combine(Channel.fromPath("${params.data_dir}/${params.species}/isotypes/haplotype_df_isotype.bed"))
                    .combine(Channel.fromPath("${params.data_dir}/${params.species}/isotypes/div_isotype_list.txt")) | divergent_and_haplotype

                if(med) {
                    // generate main html report

                    traits_to_map                           // [[t1, t1.tsv], [t2, t2.tsv], ..., [tN, tN.tsv]]
                        .combine(peaks_inbred)                              // QTL_peaks_inbred.tsv
                        .combine(peaks_loco)                                // QTL_peaks_inbred.tsv
                        .combine(fix_strain_names_bulk.out.strain_issues)   // strain_issues.txt
                        .combine(collect_eigen_variants.out)                // total_independent_tests.txt
                        .combine(vcf_to_geno_matrix.out)                    // Genotype_Matrix.tsv
                        .combine(Channel.fromPath("${params.bin_dir}/NemaScan_Report_main.Rmd"))
                        .combine(Channel.fromPath("${params.bin_dir}/NemaScan_Report_region_template.Rmd"))
                        .combine(Channel.fromPath("${params.bin_dir}/render_markdown.R"))
                        .combine(Channel.fromPath("${params.bin_dir}/NemaScan_Report_algorithm_template.Rmd"))
                        .combine(Channel.of(med))                           // val(true|false) - plot mediaton?
                        .combine(Channel.of("${params.species}"))           // val(species)
                        .combine(divergent_and_haplotype.out.div_hap_table_inbred)  // 4 files
                        .combine(divergent_and_haplotype.out.div_hap_table_loco)    // 4 files
                        .join(gcta_intervals_maps.out.for_html)             // *AGGREGATE_mapping_inbred.tsv *AGGREGATE_mapping_loco.tsv
                        .join(gcta_fine_maps.out.finemap_GWA_inbred, remainder: true)   // fine mapping data
                        .join(gcta_fine_maps.out.finemap_df_inbred, remainder: true)    // fine mapping data
                        .join(gcta_fine_maps.out.finemap_genes_inbred, remainder: true) // fine mapping data
                        .join(gcta_fine_maps.out.finemap_GWA_loco, remainder: true)     // fine mapping data
                        .join(gcta_fine_maps.out.finemap_df_loco, remainder: true)      // fine mapping data
                        .join(gcta_fine_maps.out.finemap_genes_loco, remainder: true)   // fine mapping data
                        .join(prep_ld_files_inbred.out.finemap_LD, remainder: true)     // LD files
                        .join(prep_ld_files_inbred.out.finemap_ROI, remainder: true)    // LD files
                        .join(prep_ld_files_loco.out.finemap_LD, remainder: true)       // LD files
                        .join(prep_ld_files_loco.out.finemap_ROI, remainder: true)      // LD files
                        .join(summary_mediation.out.final_mediation_inbred, remainder: true)
                        .join(summary_mediation.out.final_mediation_loco, remainder: true) | html_report_main
                } else {
                    // generate main html report

                    report_input = traits_to_map                           // [[t1, t1.tsv], [t2, t2.tsv], ..., [tN, tN.tsv]]
                        .combine(peaks_inbred)                              // QTL_peaks_inbred.tsv
                        .combine(peaks_loco)                                // QTL_peaks_inbred.tsv
                        .combine(fix_strain_names_bulk.out.strain_issues)   // strain_issues.txt
                        .combine(collect_eigen_variants.out)                // total_independent_tests.txt
                        .combine(vcf_to_geno_matrix.out)                    // Genotype_Matrix.tsv
                        .combine(Channel.fromPath("${params.bin_dir}/NemaScan_Report_main.Rmd"))
                        .combine(Channel.fromPath("${params.bin_dir}/NemaScan_Report_region_template.Rmd"))
                        .combine(Channel.fromPath("${params.bin_dir}/render_markdown.R"))
                        .combine(Channel.fromPath("${params.bin_dir}/NemaScan_Report_algorithm_template.Rmd"))
                        .combine(Channel.of(med))                           // val(true|false) - plot mediaton?
                        .combine(Channel.of("${params.species}"))           // val(species)
                        .combine(divergent_and_haplotype.out.div_hap_table_inbred)  // 4 files
                        .combine(divergent_and_haplotype.out.div_hap_table_loco)    // 4 files
                        .join(gcta_intervals_maps.out.for_html)             // *AGGREGATE_mapping_inbred.tsv *AGGREGATE_mapping_loco.tsv
                        .join(gcta_fine_maps.out.finemap_GWA_inbred, remainder: true)   // fine mapping data
                        .join(gcta_fine_maps.out.finemap_df_inbred, remainder: true)    // fine mapping data
                        .join(gcta_fine_maps.out.finemap_genes_inbred, remainder: true) // fine mapping data
                        .join(gcta_fine_maps.out.finemap_GWA_loco, remainder: true)     // fine mapping data
                        .join(gcta_fine_maps.out.finemap_df_loco, remainder: true)      // fine mapping data
                        .join(gcta_fine_maps.out.finemap_genes_loco, remainder: true)   // fine mapping data
                        .join(prep_ld_files_inbred.out.finemap_LD, remainder: true)     // LD files
                        .join(prep_ld_files_inbred.out.finemap_ROI, remainder: true)    // LD files
                        .join(prep_ld_files_loco.out.finemap_LD, remainder: true)       // LD files
                        .join(prep_ld_files_loco.out.finemap_ROI, remainder: true)      // LD files
                        
                    report_input.combine(Channel.of(null)).combine(Channel.of(null)) | html_report_main
                }  
            } else {
                // generate main html report

                report_input = traits_to_map                           // [[t1, t1.tsv], [t2, t2.tsv], ..., [tN, tN.tsv]]
                    .combine(peaks_inbred)                              // QTL_peaks_inbred.tsv
                    .combine(peaks_loco)                                // QTL_peaks_inbred.tsv
                    .combine(fix_strain_names_bulk.out.strain_issues)   // strain_issues.txt
                    .combine(collect_eigen_variants.out)                // total_independent_tests.txt
                    .combine(vcf_to_geno_matrix.out)                    // Genotype_Matrix.tsv
                    .combine(Channel.fromPath("${params.bin_dir}/NemaScan_Report_main.Rmd"))
                    .combine(Channel.fromPath("${params.bin_dir}/NemaScan_Report_region_template.Rmd"))
                    .combine(Channel.fromPath("${params.bin_dir}/render_markdown.R"))
                    .combine(Channel.fromPath("${params.bin_dir}/NemaScan_Report_algorithm_template.Rmd"))
                    .combine(Channel.of(med))                           // val(true|false) - plot mediaton?
                    .combine(Channel.of("${params.species}"))           // val(species)
                    .combine(Channel.of(null)).combine(Channel.of(null))  // dummy files for divergent_and_haplotype
                    .combine(Channel.of(null)).combine(Channel.of(null))  // dummy files for divergent_and_haplotype
                    .combine(Channel.of(null)).combine(Channel.of(null))  // dummy files for divergent_and_haplotype
                    .combine(Channel.of(null)).combine(Channel.of(null))  // dummy files for divergent_and_haplotype
                    .join(gcta_intervals_maps.out.for_html)             // *AGGREGATE_mapping_inbred.tsv *AGGREGATE_mapping_loco.tsv
                    .join(gcta_fine_maps.out.finemap_GWA_inbred, remainder: true)   // fine mapping data
                    .join(gcta_fine_maps.out.finemap_df_inbred, remainder: true)    // fine mapping data
                    .join(gcta_fine_maps.out.finemap_genes_inbred, remainder: true) // fine mapping data
                    .join(gcta_fine_maps.out.finemap_GWA_loco, remainder: true)     // fine mapping data
                    .join(gcta_fine_maps.out.finemap_df_loco, remainder: true)      // fine mapping data
                    .join(gcta_fine_maps.out.finemap_genes_loco, remainder: true)   // fine mapping data
                    .join(prep_ld_files_inbred.out.finemap_LD, remainder: true)     // LD files
                    .join(prep_ld_files_inbred.out.finemap_ROI, remainder: true)    // LD files
                    .join(prep_ld_files_loco.out.finemap_LD, remainder: true)       // LD files
                    .join(prep_ld_files_loco.out.finemap_ROI, remainder: true)      // LD files
                    
                report_input.combine(Channel.of(null)).combine(Channel.of(null)) | html_report_main
            }
        }
    
    }
    if(params.matrix) {

        // only run geno matrix step - and fix isotype names if needed
        Channel.fromPath("${params.strains}")
            .combine(Channel.fromPath("${params.data_dir}/${params.species}/isotypes/strain_isotype_lookup.tsv"))
            .combine(Channel.fromPath("${params.bin_dir}/Fix_Isotype_names_alt.R"))
            .combine(Channel.of("${params.fix}")) | fix_strain_names_alt
        
        pheno_strains = fix_strain_names_alt.out.phenotyped_strains_to_analyze

        vcf_file.combine(vcf_index)
                .combine(pheno_strains) | vcf_to_geno_matrix

    }
    if(simulation) {

        // for simulations
        Channel.fromPath("${params.data_dir}/${params.simulate_strains}")
            .splitCsv(sep:" ")
            .map { SM, STRAINS -> [SM, STRAINS] }
            .combine(vcf_file.combine(vcf_index))
            .combine(Channel.fromPath("${params.data_dir}/all_species/rename_chromosomes"))
            .combine(Channel.fromPath("${params.data_dir}/${params.simulate_maf}").splitCsv()) | prepare_simulation_files

        // eigen
        contigs = Channel.of(1, 2, 3, 4, 5, 6)
        contigs.combine(prepare_simulation_files.out.sim_geno)
            .combine(Channel.fromPath("${params.bin_dir}/Get_GenoMatrix_Eigen.R")) | chrom_eigen_variants_sims

        chrom_eigen_variants_sims.out.sim_geno_eigen_join
            .groupTuple(by:[0,1,2]).
            join(chrom_eigen_variants_sims.out.sim_geno_meta, by:[0,1,2]) | collect_eigen_variants_sims

        // simulate qtl or genome
        if(params.simulate_qtlloc){

            collect_eigen_variants_sims.out
                .combine(Channel.fromPath("${params.data_dir}/${params.simulate_nqtl}").splitCsv())
                .combine(Channel.fromPath("${params.data_dir}/${params.simulate_qtlloc}"))
                .combine(Channel.fromPath("${params.data_dir}/${params.simulate_eff}").splitCsv())
                .combine(Channel.of(1..params.simulate_reps))
                .combine(Channel.fromPath("${params.bin_dir}/create_causal_QTLs.R")) | simulate_effects_loc

            sim_phen_inputs = simulate_effects_loc.out

        } else {

            collect_eigen_variants_sims.out
                .combine(Channel.fromPath("${params.data_dir}/${params.simulate_nqtl}").splitCsv())
                .combine(Channel.fromPath("${params.data_dir}/${params.simulate_eff}").splitCsv())
                .combine(Channel.of(1..params.simulate_reps))
                .combine(Channel.fromPath("${params.bin_dir}/create_causal_QTLs.R")) | simulate_effects_genome

            sim_phen_inputs = simulate_effects_genome.out

        }

        sim_phen_inputs
            .combine(Channel.fromPath("${params.data_dir}/${params.simulate_h2}").splitCsv()) 
            .combine(Channel.fromPath("${params.bin_dir}/check_vp.py")) | simulate_map_phenotypes

        // simulation mappings
        simulate_map_phenotypes.out.gcta_intervals
            .combine(Channel.of("${params.sthresh}"))
            .combine(Channel.of("${params.group_qtl}"))
            .combine(Channel.of("${params.ci_size}")) 
            .combine(Channel.fromPath("${params.bin_dir}/Aggregate_Mappings.R"))
            .combine(Channel.fromPath("${params.bin_dir}/Find_Aggregate_Intervals.R"))
            .combine(Channel.fromPath("${params.bin_dir}/Find_GCTA_Intervals.R"))
            .combine(Channel.fromPath("${params.bin_dir}/Find_GCTA_Intervals_LOCO.R")) | get_gcta_intervals
     
        get_gcta_intervals.out.assess_data_inbred_pca 
            .combine(Channel.fromPath("${params.bin_dir}/Assess_Sim.R")) | assess_sims_INBRED | collectFile(name: "${params.out}/INBRED_PCA_all_sims.tsv") | view
        get_gcta_intervals.out.assess_data_loco_pca 
            .combine(Channel.fromPath("${params.bin_dir}/Assess_Sim.R")) | assess_sims_LOCO | collectFile(name: "${params.out}/LOCO_PCA_all_sims.tsv") | view

    
    }
}



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
    Threshold for grouping QTL              = ${params.group_qtl}
    Number of SNVs to define CI             = ${params.ci_size}
    Mapping                                 = ${params.mapping}
    Simulation                              = ${params.simulate}
    Simulate QTL effects                    = ${params.simulate_qtlloc}
    Result Directory                        = ${params.out}
    """

    // println summary

}
