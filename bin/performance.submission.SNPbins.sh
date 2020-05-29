#!/bin/bash
#SBATCH -A b1042                              # Allocation
#SBATCH -p genomicsguestA                     # Queue
#SBATCH -t 8:00:00                            # Walltime/duration of the job
#SBATCH -N 1                                  # Number of Nodes
#SBATCH -n 16
#SBATCH --mem=16G
#SBATCH --job-name="GWA.perf.bp.bins"

cd /projects/b1059/projects/Sam/NemaScan/
module load R/3.6.3
Rscript bin/performance.assessment.BPbins.R $1 $2
