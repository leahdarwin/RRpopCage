# ============================================================================
# Script: Selection simulation (Wright-Fisher + pool-seq sampling + GLM)
#
# Description:
#   Defines the simulator used for the selection power simulations: a diploid
#   Wright-Fisher model with selection followed by two-stage pool-seq sampling,
#   tested with the treatment x time GLM from
#   poolFreqDiff_treatment_time_repl_RR.py. Sourced by run_selectionSim.R; can
#   also be run directly for a single example simulation.
#
# Author: Leah Darwin
# ============================================================================
#
# Model details
# -------------
# One biallelic locus, one selection coefficient `s` per run. Three populations
# (treatments, as in the BEI/YAK/ZIM mito design): population 1 evolves with
# selection coefficient s, populations 2 and 3 are identical and neutral
# (s = 0). Each population has `n_reps` replicate cages. The locus starts at p0
# in every cage. The simulation returns the p-value of the treatment:time term
# from the quasibinomial GLM used in the real analysis.
#
# Per generation:
#   1. Selection (deterministic):
#        p' = [p^2 (1+s) + p q (1+hs)] / wbar
#        wbar = p^2 (1+s) + 2pq (1+hs) + q^2
#   2. Drift: p_{t+1} = Binomial(2N, p') / 2N
#
# At each sampled generation, two-step binomial sampling (Spitzer et al. 2020,
# Fig. 1b / Appendix A3):
#   1. Individual sampling: n_ind diploids -> g = 2*n_ind chromosomes
#        K1 ~ Binomial(g, p)
#   2. Read sampling: coverage D is fixed; K2 ~ Binomial(D, K1 / g)
#        p_obs = K2 / D
# Var(K2) = D p (1-p) (1 + (D-1)/g)  (law of total variance, eq. A21)
# =============================================================================


# ---- Selection step (diploid, vectorised) -----------------------------------
select_diploid <- function(p, s, h) {
  q    <- 1 - p
  wbar <- p^2 * (1 + s) + 2 * p * q * (1 + h * s) + q^2
  (p^2 * (1 + s) + p * q * (1 + h * s)) / wbar
}


# ---- Two-step binomial sampling ---------------------------------------------
# p:     vector of true population allele frequencies
# n_ind: number of diploid individuals sampled (2*n_ind chromosomes)
# depth: read depth (coverage), fixed
poolseq_sample <- function(p, n_ind, depth) {
  m <- length(p)
  n_chrom <- 2 * n_ind

  # Step 1: sample chromosomes from the population
  k      <- rbinom(m, n_chrom, p)
  p_pool <- k / n_chrom

  # Step 2: sample reads from the sampled chromosomes
  r     <- rbinom(m, depth, p_pool)
  p_obs <- r / depth

  list(p_pool = p_pool, depth = rep(depth, m), alt_reads = r, p_obs = p_obs)
}


# ---- GLM test (same model as poolFreqDiff_treatment_time_repl_RR.py) ---------
# obs: data.frame with columns tr_l, rep, gen, alt_reads, depth
# Mirrors get_dat() in tools/poolFreqDiff_AF_RR.R: if any allele count is 0,
# 1 is added to every allele count (zeroes = 1). Returns the LRT p-value of the
# tr_l:time term (treatment-by-time interaction).
glm_test <- function(obs, zeroes = TRUE) {
  A   <- obs$alt_reads
  Ref <- obs$depth - obs$alt_reads
  if (zeroes && (any(A == 0) || any(Ref == 0))) {
    A   <- A + 1
    Ref <- Ref + 1
  }
  dat <- data.frame(A_Cnt   = A,
                    Tot_Cnt = A + Ref,
                    tr_l    = factor(obs$tr_l),
                    rep     = factor(obs$rep),
                    time    = obs$gen)
  res <- glm(cbind(A_Cnt, Tot_Cnt - A_Cnt) ~ tr_l * time + rep + tr_l:rep,
             family = "quasibinomial", data = dat)
  mod <- anova(res, test = "LRT")
  mod[["Pr(>Chi)"]][rownames(mod) == "tr_l:time"]
}


# ---- Main simulator ---------------------------------------------------------
# Returns a list with:
#   $p_value : p-value of the treatment:time test
#   $obs     : data.frame (tr_l, rep, s, gen, p_true, p_pool, depth,
#              alt_reads, p_obs) for the sampled generations
#   $true    : data.frame (tr_l, rep, s, gen, p_true) for every generation
#   $params
simulate_wf_poolseq <- function(N           = 1000,
                                s            = 0.05,
                                h            = 0.5,
                                n_gens       = 25,
                                n_reps       = 5,
                                p0           = 0.5,
                                n_ind        = 50,
                                depth        = 100,
                                sample_gens  = seq(0, n_gens, by = 5),
                                seed         = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample_gens <- sort(unique(sample_gens[sample_gens >= 0 & sample_gens <= n_gens]))

  # Cages: treatment 1 = selected, treatments 2 and 3 = neutral
  design <- expand.grid(rep_in_tr = seq_len(n_reps), tr_l = 1:3)
  design$rep <- seq_len(nrow(design))            # unique cage id (1..3*n_reps)
  design$s   <- ifelse(design$tr_l == 1, s, 0)
  n_cages <- nrow(design)

  p    <- rep(p0, n_cages)
  traj <- matrix(NA_real_, nrow = n_gens + 1, ncol = n_cages)
  traj[1, ] <- p

  for (t in seq_len(n_gens)) {
    p_sel <- select_diploid(p, design$s, h)
    p <- rbinom(n_cages, 2 * N, p_sel) / (2 * N)
    traj[t + 1, ] <- p
  }

  true_df <- data.frame(
    tr_l   = rep(design$tr_l, each = n_gens + 1),
    rep    = rep(design$rep,  each = n_gens + 1),
    s      = rep(design$s,    each = n_gens + 1),
    gen    = rep(0:n_gens, times = n_cages),
    p_true = as.vector(traj)
  )

  obs_df <- true_df[true_df$gen %in% sample_gens, ]
  ps <- poolseq_sample(obs_df$p_true, n_ind, depth)
  obs_df$p_pool    <- ps$p_pool
  obs_df$depth     <- ps$depth
  obs_df$alt_reads <- ps$alt_reads
  obs_df$p_obs     <- ps$p_obs
  rownames(obs_df) <- NULL

  list(p_value = glm_test(obs_df), obs = obs_df, true = true_df,
       params = list(N = N, s = s, h = h, n_gens = n_gens, n_reps = n_reps,
                     p0 = p0, n_ind = n_ind, depth = depth,
                     sample_gens = sample_gens))
}


# ---- Example ----------------------------------------------------------------
# Only runs when this file is executed directly, not when source()d
if (sys.nframe() == 0) {
  sim <- simulate_wf_poolseq(
    N           = 5000,
    s           = 0.1,
    h           = 0.5,
    n_gens      = 25,
    n_reps      = 5,
    p0          = 0.5,
    n_ind       = 200,
    depth       = 50,
    sample_gens = c(2, 10, 15, 25),
    seed        = 99
  )
  head(sim$obs)
  cat("treatment:time p-value:", sim$p_value, "\n")

  # Power at a given s: fraction of runs with p < 0.05
  # mean(replicate(200, simulate_wf_poolseq(s = 0.02)$p_value) < 0.05)
}
