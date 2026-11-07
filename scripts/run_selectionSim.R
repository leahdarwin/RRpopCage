# ============================================================================
# Script: Run selection simulations for one selection coefficient
#
# Description:
#   Runs the simulator in selectionSim.R n_sims times for a single selection
#   coefficient and writes the treatment x time GLM p-values (columns sim, s,
#   p_value) to a csv. Called once per selection coefficient by
#   selectionSim.sh.
#
# Usage: Rscript run_selectionSim.R <s> <n_sims> <out_csv> [seed]
#
# Author: Leah Darwin
# ============================================================================

args   <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript run_selectionSim.R <s> <n_sims> <out_csv> [seed]")
s      <- as.numeric(args[1])
n_sims <- as.integer(args[2])
out    <- args[3]
seed   <- if (length(args) >= 4) as.integer(args[4]) else 1L

script_dir <- dirname(normalizePath(sub("^--file=", "",
                grep("^--file=", commandArgs(), value = TRUE)[1])))
source(file.path(script_dir, "selectionSim.R"))

set.seed(seed)
p_values <- vapply(seq_len(n_sims), function(i) {
  simulate_wf_poolseq(
    N           = 2150,
    s           = s,
    h           = 0.5,
    n_gens      = 25,
    n_reps      = 5,
    p0          = 0.5,
    n_ind       = 200,
    depth       = 50,
    sample_gens = c(2, 10, 15, 25)
  )$p_value
}, numeric(1))

write.csv(data.frame(sim = seq_len(n_sims), s = s, p_value = p_values),
          out, row.names = FALSE)
cat("s =", s, ": wrote", n_sims, "p-values to", out,
    "| fraction p < 0.05:", mean(p_values < 0.05, na.rm = TRUE), "\n")
