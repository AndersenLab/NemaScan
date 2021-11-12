/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *      MEDIATION      * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/


process mediation_data {
 
    executor 'local'
    tag {TRAIT}
    label "mediation"

    input:
        tuple val(TRAIT),val(tch),val(tstart),val(tpeak),val(tend),val(logPvalue), val(peak_id),val(h2), val(marker), file(t_file), \
        val(transcript_eqtl), file(mediation_input)

    output:
        tuple val(TRAIT),val(tch),val(tpeak),val(tstart),val(tend), file("${TRAIT}_scaled_mapping.tsv"),file("${TRAIT}_${tch}_${tpeak}_eqtl.tsv")

    """
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${mediation_input} > mediation_input 
    Rscript --vanilla mediation_input ${TRAIT} ${t_file} ${tch} ${tstart} ${tend} ${tpeak} ${transcript_eqtl}

    """
}


process multi_mediation {

    cpus 1
    memory '2 GB'
    label "mediation"

    tag {"${TRAIT}_${tch}_${tpeak}"}

    input:
        tuple val(TRAIT),val(tch),val(tpeak), val(tstart),val(tend), file(pheno), file(tr_eqtl), file(geno), file(texpression), file("multi_mediation")


    output:
        path "${TRAIT}_${tch}_${tpeak}_medmulti.tsv", emit: result_multi_mediate optional true
        path "${TRAIT}_${tch}_${tpeak}_elist.tsv", emit: eQTL_gene optional true


    """
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${multi_mediation} > multi_mediation_file 
    Rscript --vanilla multi_mediation_file ${geno} ${texpression} ${pheno} ${tch} ${tpeak} ${TRAIT} ${tr_eqtl}
    
    """
}


process simple_mediation {
 
    cpus 1
    memory '2 GB'
    tag {"${TRAIT}_${gene}"}
    label "mediation"

    input:
        tuple val(TRAIT),val(tch),val(tpeak),val(gene), val(tstart),val(tend), file(pheno), file(tr_eqtl), file(geno), file(expression), file(simple_mediation)

    output:
        file("${TRAIT}_${tch}_${tpeak}_${gene}_med.tsv") 

    """
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${simple_mediation} > simple_mediation_file 
    Rscript --vanilla simple_mediation_file ${gene} ${geno} ${expression} ${pheno} ${tch} ${tpeak} ${TRAIT} ${tr_eqtl}

    """
}


process summary_mediation {

    cpus 2
    memory '32 GB'
    label "mediation"

    publishDir "${params.out}/Mediation/file_summary", mode: 'copy', pattern: "*mediation.tsv"
    publishDir "${params.out}/Mediation/plot_summary", mode: 'copy', pattern: "*plot.png"

    input:
     tuple val(TRAIT),val(tch),val(tstart),val(tpeak),val(tend),val(logPvalue), val(peak_id),val(h2), val(marker), \
     file(summary_mediation), file("*"), file("*")//file("*_medmulti.tsv"), file("*_med.tsv")


    output:
        tuple val(TRAIT), file("${TRAIT}_mediation.tsv"), emit: final_mediation
        file("*plot.png") optional true


    """
    cat ${TRAIT}_*medmulti.tsv > ${TRAIT}_multi_mediation_analysis.tsv
    cat ${TRAIT}_*med.tsv  > ${TRAIT}_indiv_mediation_analysis.tsv

    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${summary_mediation} > summary_mediation_file 
    Rscript --vanilla summary_mediation_file ${TRAIT}_multi_mediation_analysis.tsv ${TRAIT}_indiv_mediation_analysis.tsv ${TRAIT}

    """
}