#!/bin/bash

# ============================================================================
# Script: GLM job submission for treatment-time-replicate allele frequency model
#
# Description:
#   Slurm batch script that runs poolFreqDiff on the joined sync file to fit
#   the treatment/time/replicate GLM, then converts the generated R script
#   output into a GLM results table.
#
# Author: Leah Darwin
# ============================================================================

#SBATCH -p batch
#SBATCH --mem=10G
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o %j.out
#SBATCH -e %j.err

code_dir="/tools/poolFreqDiff/"

mts="_RR"

model="treatment_time_repl${mts}"
script="poolFreqDiff_${model}.py"

glm_dir="data/"
sync_dir="aligned_reads_6.32/"

##this file is too large to upload to github and will need to be regenerated using raw reads and popoolation2 for the sync file  
input_file="${sync_dir}joined.sync"
output_file="${glm_dir}${model}.glm"

nsamps=60

##compute cluster specific module loads
module load r/4.5.1-iikl
module load miniforge3/25.3.0-3
source ${MAMBA_ROOT_PREFIX}/etc/profile.d/conda.sh
##yaml file for this conda environment is given on the github site
conda activate py27

python "${code_dir}${script}" -filename "$input_file" -npops $nsamps -nlevels 1 -n 200 -mincnt 10 -minc 30 -maxc 200 -rescale nr -zeroes 1 > "${output_file}.rin"

Rscript "${output_file}.rin" > "$output_file"
