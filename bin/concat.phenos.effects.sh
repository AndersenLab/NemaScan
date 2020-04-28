#!/bin/bash
#SBATCH -A b1042                        # Allocation
#SBATCH -p genomicsguestA               # Queue
#SBATCH -t 12:00:00                     # Walltime/duration of the job
#SBATCH -N 1                            # Number of Nodes
#SBATCH --job-name="GCTA.cat"    # Name of job

BASE=/projects/b1059/projects/Sam/NemaScan
NQTL=$1
HSQ=$2
MAF=$3
LD=$4

cd ${BASE}/HSQ.${HSQ}.MAF.${MAF}.LD.${LD}.simulations
mkdir ${NQTL}.QTLs
mkdir logs

cd ${NQTL}.QTLs
mkdir causal.QTLs
mkdir simulated.phenotypes

DESTINATION_DIR=${BASE}/HSQ.${HSQ}.MAF.${MAF}.LD.${LD}.simulations/${NQTL}.QTLs
echo "Heritability: ${HSQ}"
echo "Minor Allele Frequency: ${MAF}"
echo "LD Cutoff: ${LD}"
echo "Architecture: ${NQTL} QTL"

for i in {1..1000}; \
  do cd ${BASE}/HSQ.${HSQ}.MAF.${MAF}.LD.${LD}.simulations/simulation.rep.$i; \

  # Causal QTLs and Effects
  mv variants.sim.${NQTL}.out.par causal.QTLs.rep.$i; \
  mv causal.QTLs.rep.$i ${DESTINATION_DIR}/causal.QTLs/; \

  # Simulated Phenotypes
  mv variants.sim.${NQTL}.out.phen sim.phenos.QTLs.rep.$i; \
  mv sim.phenos.QTLs.rep.$i ${DESTINATION_DIR}/simulated.phenotypes/; \

done
