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
 
    tag {TRAIT}

    label "mediation"
    label "xs"

    input:
        tuple val(TRAIT),val(tch),val(tstart),val(tpeak),val(tend),val(logPvalue), val(peak_id),val(h2), val(marker), val(algorithm), file(t_file), \
        file(transcript_eqtl), file(mediation_input)

    output:
        tuple val(TRAIT),val(tch),val(tpeak), val(algorithm),val(tstart),val(tend), file("${TRAIT}_scaled_mapping_${algorithm}.tsv"),file("${TRAIT}_${tch}_${tpeak}_eqtl_${algorithm}.tsv")

    """
    Rscript --vanilla ${mediation_input} ${TRAIT} ${t_file} ${tch} ${tstart} ${tend} ${tpeak} ${transcript_eqtl} ${algorithm}
    """
}


process multi_mediation {

    label "xs"
    label "mediation"

    tag {"${TRAIT}_${tch}_${tpeak}"}

    input:
        tuple val(TRAIT),val(tch),val(tpeak), val(algorithm), val(tstart),val(tend), file(pheno), file(tr_eqtl), file(geno), file(texpression), file("multi_mediation")


    output:
        tuple val(TRAIT), val(algorithm), path("${TRAIT}_${tch}_${tpeak}_medmulti_${algorithm}.tsv"), emit: result_multi_mediate, optional: true
        path "${TRAIT}_${tch}_${tpeak}_elist_${algorithm}.tsv", emit: eQTL_gene, optional: true
        // val 'algorithm', emit: med_alg


    """
    Rscript --vanilla ${multi_mediation} ${geno} ${texpression} ${pheno} ${tch} ${tpeak} ${TRAIT} ${tr_eqtl} ${algorithm}
    
    """
}


process simple_mediation {
 
    tag {"${TRAIT}_${gene}"}

    label "mediation"
    label 'xs'

    input:
        tuple val(TRAIT),val(tch),val(tpeak), val(algorithm), val(gene), val(tstart),val(tend), file(pheno), file(tr_eqtl), file(geno), file(expression), file(simple_mediation)

    output:
        tuple val(TRAIT), val(algorithm), file("${TRAIT}_${tch}_${tpeak}_${gene}_med_${algorithm}.tsv") 

    """
    Rscript --vanilla ${simple_mediation} ${gene} ${geno} ${expression} ${pheno} ${tch} ${tpeak} ${TRAIT} ${tr_eqtl} ${algorithm}

    """
}


process summary_mediation {

    label "ml"
    label "mediation"

    publishDir "${params.out}/INBRED/Mediation", mode: 'copy', pattern: "*mediation_inbred.tsv"
    publishDir "${params.out}/LOCO/Mediation", mode: 'copy', pattern: "*mediation_loco.tsv"
    publishDir "${params.out}/INBRED/Mediation", mode: 'copy', pattern: "*plot_inbred.png"
    publishDir "${params.out}/LOCO/Mediation", mode: 'copy', pattern: "*plot_loco.png"

    input:
        tuple val(TRAIT), val(algorithm), file("*"), file("*"), file(summary_mediation)
    //  tuple val(TRAIT),val(tch),val(tstart),val(tpeak),val(tend),val(logPvalue), val(peak_id),val(h2), val(marker), val(algorithm), \
    //  file(summary_mediation), file("*"), file("*")//file("*_medmulti.tsv"), file("*_med.tsv")


    output:
        tuple val(TRAIT), file("${TRAIT}_mediation_inbred.tsv"), emit: final_mediation_inbred, optional: true
        tuple val(TRAIT), file("${TRAIT}_mediation_loco.tsv"), emit: final_mediation_loco, optional: true
        file("*plot_*.png") optional true


    """
    cat ${TRAIT}_*medmulti_${algorithm}.tsv > ${TRAIT}_multi_mediation_analysis_${algorithm}.tsv
    cat ${TRAIT}_*med_${algorithm}.tsv  > ${TRAIT}_indiv_mediation_analysis_${algorithm}.tsv

    Rscript --vanilla ${summary_mediation} ${TRAIT}_multi_mediation_analysis_${algorithm}.tsv ${TRAIT}_indiv_mediation_analysis_${algorithm}.tsv ${TRAIT} ${algorithm}
    
    """
}
