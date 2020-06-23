#!/bin/bash
#SBATCH -A b1042                              # Allocation
#SBATCH -p genomicsguestA                     # Queue
#SBATCH -t 24:00:00                            # Walltime/duration of the job
#SBATCH -N 1                                  # Number of Nodes
#SBATCH -n 8
#SBATCH --mem=32G
#SBATCH --job-name="GWA.perf"

cd /projects/b1059/projects/Sam/NemaScan/
module load R/3.6.3
Rscript bin/performance.assessment.20200609.R $1 $2

# usage: sbatch bin/performance.submission.20200609.sh [directory] [bin size]
