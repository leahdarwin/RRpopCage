#!/bin/bash

# ============================================================================
# Script: Haplotype frequency job submission
#
# Description:
#   Slurm job array script that runs haplo_frqs.R once per population
#   (f1-f60, one job per array task) to estimate founder haplotype
#   frequencies in sliding windows from the variant frequency table.
#
# Author: Leah Darwin
# ============================================================================

#SBATCH -p batch
#SBATCH --mem=5G
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=1-60
#SBATCH -o %j.out
#SBATCH -e %j.err

##cluster specific module loads
module load r/4.5.1-iikl

varfrq="/users/drand/data/RR_popcage_poolseq/aligned_reads_6.32/var_frq.tsv"

##alternatively use data/foundergt_extended.names for calling haplotypes for additional founders
founder="data/foundergt.names"

output="data/haplo_frq/"

##each population is submitted as part of a job array 
pop="f${SLURM_ARRAY_TASK_ID}"
echo "Running job for population = $pop"

Rscript haplo_frqs.R "$varfrq" "$founder" "${output}${pop}.tsv" "$pop"

