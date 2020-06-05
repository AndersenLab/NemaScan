#!/bin/bash
export BASE=/projects/b1059/projects/Sam/NemaScan
export BIN=${BASE}/bin
# CURRENTLY HARD-CODED VCF - in the future maybe this can be fed through NemaScan parameters
# OR - in the future eliminate completely and use NemaScan genotype matrix to pick variants, i.e.
export VCF=${BASE}/cegwas2-nf/bin/WI.HARD-FILTERED.isotype.vcf.gz
export VCF_NAME=$(ls ${VCF} | cut -f9 -d "/" | cut -f1-3 -d "." | cat)
export GCTABIN=/projects/b1059/software/gcta_1.93.0beta

# Set Parameters For Simulations
export HSQ=$1
export MAF=$2
export LD=$3

mkdir HSQ.${HSQ}.MAF.${MAF}.LD.${LD}.simulations
cd HSQ.${HSQ}.MAF.${MAF}.LD.${LD}.simulations

# LD and MAF Pruning for Causal Variants
plink --vcf ${VCF} \
--snps-only \
--biallelic-only \
--maf ${MAF} \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 ${LD} \
--geno \
--allow-extra-chr \
--make-bed

plink --bfile plink \
--extract plink.prune.in \
--allow-extra-chr \
--make-bed \
--out ${VCF_NAME}.pruned

export BIMNAME=${BASE}/HSQ.${HSQ}.MAF.${MAF}.LD.${LD}.simulations/${VCF_NAME}.pruned
Rscript --vanilla ${BIN}/recode.bim.R ${BIMNAME}.bim
mv ${BIMNAME}.bed recoded.bed
mv ${BIMNAME}.fam recoded.fam
export RECODED=${BASE}/HSQ.${HSQ}.MAF.${MAF}.LD.${LD}.simulations/recoded

for i in {1..1000}; \
  do sbatch ${BIN}/gcta.simulation.prunedvariants.sh $i ${MAF} ${HSQ}; \
done
