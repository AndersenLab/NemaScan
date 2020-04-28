#!/bin/bash
#SBATCH -A b1042                        # Allocation
#SBATCH -p genomicsguestA               # Queue
#SBATCH -t 12:00:00                     # Walltime/duration of the job
#SBATCH -N 1                            # Number of Nodes
#SBATCH --job-name="GCTA.SJW"    # Name of job

export REP=$1
# Example:
#

conda activate sam_env
mkdir simulation.rep.${REP}
cd simulation.rep.${REP}


echo "Generating Causal QTLs"
Rscript --vanilla ${SCRIPTS}/create.causal.QTLs.R ${RECODED}.bim
echo "REPLICATE ${REP} TEST STATEMENT"

# Simulated Traits, Supplied Heritability
export TEXTS=$(ls causal.variants.sim.*.txt | cat)
export SIMQTLS=$(ls ${TEXTS} | cut -f2-4 -d "." | cat)
COUNT=1; \
for i in ${SIMQTLS}; \
  do ${GCTABIN}/./gcta64 \
  --bfile ${RECODED} \
  --simu-qt \
  --simu-causal-loci causal.$i.txt \
  --simu-hsq ${HSQ} \
  --simu-rep 1 \
  --out $i.out; \
  let COUNT=$COUNT+1; \
done \
