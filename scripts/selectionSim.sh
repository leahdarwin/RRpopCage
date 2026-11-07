#!/bin/bash

# ============================================================================
# Script: Selection simulation job array
#
# Description:
#   Slurm job array. Each task takes one selection coefficient from S_VALUES,
#   runs the Wright-Fisher + pool-seq simulation (selectionSim.R) N_SIMS times,
#   applies the treatment x time GLM each time, and writes the p-values to
#   output/selectionSim/pvals_s<s>.csv.
#
# Usage:
#   sbatch scripts/selectionSim.sh
#
# Author: Leah Darwin
# ============================================================================

#SBATCH -p batch
#SBATCH --mem=4G
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=0-9
#SBATCH -o selectionSim_%A_%a.out
#SBATCH -e selectionSim_%A_%a.err

S_VALUES=(0 0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.5)
N_SIMS=10000

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}/scripts"
OUT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}/output/selectionSim"
mkdir -p "$OUT_DIR"

S=${S_VALUES[$SLURM_ARRAY_TASK_ID]}

# Distinct seed per task so runs are independent but reproducible
Rscript "$SCRIPT_DIR/run_selectionSim.R" "$S" "$N_SIMS" \
  "$OUT_DIR/pvals_s${S}.csv" $((1000 + SLURM_ARRAY_TASK_ID))
