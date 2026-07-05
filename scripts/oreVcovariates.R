# ============================================================================
# Script: Founder haplotype frequency vs genomic covariates
#
# Description:
#   Examines relationships between estimated founder haplotype frequencies
#   (from haplo_frqs.R output) and two genomic covariates, recombination
#   rate (Comeron 2012 map) and CDS density, at a chosen generation.
#   Produces scatter/quantile plots by chromosome and genome-wide, plus
#   Spearman and partial Spearman correlation tables.
#
# Author: Leah Darwin
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)
library(ppcor)
library(knitr)
library(kableExtra)

# ============================================================
# Config
# ============================================================
founders        <- read.csv("data/foundergt.names", header = FALSE, col.names = "name")
recomb          <- read.csv("data/Comeron_100kb_genome.csv")
gtf_path        <- "data/dmel-all-r6.57.gtf" ##gtf file can be obtained from flybase

##files containing estimates of haplotype frequencies for each population generate with haplo_frqs.sh
data_dir        <- "data/haplo_frq/"
chr_order       <- c("2L", "2R", "3L", "3R", "4", "X")
chosen_founders <- founders$name[1:2]
chosen_gen      <- 25
outdir          <- "output"

q_labels_5 <- paste0(seq(0, 0.8, 0.2), "-", seq(0.2, 1.0, 0.2))

# ============================================================
# Population labels
# ============================================================
tr_l = c(rep("Bei", 20), rep("Yak", 20), rep("Zim", 20))
tr_l = as.character(tr_l)
rep  = c(rep(rep(1:5, each = 4), 3))
rep  = as.character(rep)
rep  = paste(tr_l, rep, sep = "")
time = c(rep(c(10, 15, 25, 2), times = 15))

poplabs            <- data.frame(tr_l, rep, time)
poplabs$population <- paste("f", rownames(poplabs), sep = "")
poplabs            <- filter(poplabs, rep != "Bei3")

# ============================================================
# Load haplotype frequency data
# ============================================================
pop_ids <- if (!is.null(chosen_gen)) {
  filter(poplabs, time == chosen_gen)$population
} else {
  poplabs$population
}

pop <- pop_ids %>%
  map_dfr(~ {
    f <- file.path(data_dir, paste0(.x, ".tsv"))
    read.csv(f, sep = "\t") %>% mutate(population = .x)
  }) %>%
  unique() %>%
  filter(NSNPs > 0) 

long_all <- pop %>%
  separate_rows(frequencies, sep = ";") %>%
  mutate(frequencies = as.numeric(frequencies)) %>%
  group_by(chr, pos, population) %>%
  mutate(founder = founders$name[seq_len(n())]) %>%
  ungroup() %>%
  left_join(poplabs, by = join_by(population), relationship = "many-to-one") %>%
  filter(founder %in% chosen_founders, chr %in% chr_order) %>%
  mutate(
    chr  = factor(chr, levels = chr_order),
    tr_l = factor(tr_l, levels = c("Bei", "Zim", "Yak"))
  )

# Filter to Comeron 2012 normal-recombination intervals (dm6 liftover).
# Positions outside each chromosome's interval, and all of chr4 (no row in
# the file), are dropped entirely from this analysis.
chrmask <- read.csv("chrmask_comeron2012_dm6liftover.csv")

long_all <- long_all %>%
  left_join(chrmask, by = "chr") %>%
  filter(!is.na(start), pos >= start, pos <= end) %>%
  dplyr::select(-start, -end)

# ============================================================
# Reference windows: non-overlapping 100kb (from recomb map)
# ============================================================
recomb_wins <- recomb %>%
  mutate(win_start = midpoint - 50000, win_end = midpoint + 50000 - 1)

# ============================================================
# CDS counts per 100kb window (from GTF) — computed once
# ============================================================
genes <- read.table(gtf_path, sep = "\t", quote = "", comment.char = "#",
                    col.names = c("chr","source","feature","start","end",
                                  "score","strand","frame","attributes")) %>%
  filter(feature == "CDS") %>%
  mutate(chr = sub("^chr", "", chr)) %>%
  dplyr::select(chr, gene_start = start, gene_end = end)

cds_counts <- recomb_wins %>%
  left_join(genes, by = join_by(chr, win_start <= gene_end, win_end >= gene_start)) %>%
  group_by(chr, midpoint) %>%
  summarise(n_cds = sum(!is.na(gene_start)), .groups = "drop")

# ============================================================
# Join windows once for all founders
# ============================================================
joined_wins_all <- long_all %>%
  mutate(chr = as.character(chr)) %>%
  left_join(recomb_wins, by = join_by(chr, pos >= win_start, pos <= win_end)) %>%
  filter(!is.na(midpoint))

win_pop_counts <- joined_wins_all %>%
  filter(founder == chosen_founders[1]) %>%
  group_by(chr, midpoint) %>%
  summarise(n_pops = n_distinct(population), .groups = "drop")

n_total   <- nrow(win_pop_counts)
n_full    <- sum(win_pop_counts$n_pops == length(pop_ids))
n_dropped <- n_total - n_full
cat(sprintf(
  "Recomb windows: %d total | %d with all %d populations | %d dropped (%.1f%%)\n",
  n_total, n_full, length(pop_ids), n_dropped,
  100 * n_dropped / n_total
))

# ============================================================
# Compute per-founder window frequencies
# ============================================================
all_founder_data <- lapply(chosen_founders, function(fnd) {
  jw <- filter(joined_wins_all, founder == fnd)

  # Per-population window frequencies; filter to windows with all populations
  win_freq_pop <- jw %>%
    group_by(chr, midpoint, c, population) %>%
    summarise(pop_freq = mean(frequencies, na.rm = TRUE), .groups = "drop") %>%
    group_by(chr, midpoint) %>%
    filter(n_distinct(population) == length(pop_ids)) %>%
    ungroup() %>%
    left_join(cds_counts, by = c("chr", "midpoint")) %>%
    filter(!is.na(n_cds))

  # Average across populations for visualization
  win_freq_fnd <- win_freq_pop %>%
    group_by(chr, midpoint, c, n_cds) %>%
    summarise(win_freq = mean(pop_freq, na.rm = TRUE), .groups = "drop")

  gw_mean_fnd <- mean(win_freq_fnd$win_freq, na.rm = TRUE)
  gw_sd_fnd   <- sd(win_freq_fnd$win_freq,   na.rm = TRUE)

  cat(sprintf("Founder %s: genome-wide mean = %.4f, SD = %.4f\n", fnd, gw_mean_fnd, gw_sd_fnd))

  list(
    founder      = fnd,
    data         = win_freq_fnd,
    data_noX     = filter(win_freq_fnd, chr != "X"),
    pop_data     = win_freq_pop,
    pop_data_noX = filter(win_freq_pop, chr != "X"),
    gw_mean      = gw_mean_fnd,
    gw_sd        = gw_sd_fnd
  )
})
names(all_founder_data) <- chosen_founders

# 95% CI for Spearman rho via Fisher Z transformation
fisher_ci <- function(rho, n, gp = 0, level = 0.95) {
  se   <- 1 / sqrt(n - gp - 3)
  crit <- qnorm((1 + level) / 2)
  tanh(atanh(rho) + c(-1, 1) * crit * se)
}


# ============================================================
# Shared plot functions
# ============================================================
make_scatter <- function(dat, x_col, x_lab, chr_name = NULL,
                         y_col = "win_freq", y_lab_arg = "OreR frequency") {
  base_sz <- if (is.null(chr_name)) 11 else 10
  ggplot(dat, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(alpha = 0.2, size = 1.2, colour = "slateblue") +
    geom_smooth(method = "loess", span = 0.75, se = FALSE,
                colour = "slateblue", linewidth = 0.8) +
    labs(x = x_lab, y = y_lab_arg, title = chr_name) +
    theme_bw(base_size = base_sz) +
    theme(panel.grid = element_blank())
}

make_quantile <- function(dat, x_col, x_lab, chr_name = NULL,
                          y_col = "win_freq", y_lab_arg = "OreR frequency") {
  dat_q <- dat %>% mutate(q_bin = ntile(.data[[x_col]], 5))
  qs <- dat_q %>%
    group_by(q_bin) %>%
    summarise(mean_y = mean(.data[[y_col]], na.rm = TRUE),
              se_y   = sd(.data[[y_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[y_col]]))),
              .groups = "drop") %>%
    arrange(q_bin)
  base_sz <- if (is.null(chr_name)) 11 else 10
  ggplot(dat_q, aes(x = q_bin, y = .data[[y_col]], group = q_bin)) +
    geom_jitter(colour = "lightgray", width = 0.15,
                alpha = 0.4, size = 1, show.legend = FALSE) +
    geom_point(data = qs, aes(x = q_bin, y = mean_y),
               size = 2.5, colour = "slateblue", inherit.aes = FALSE) +
    geom_errorbar(data = qs,
                  aes(x = q_bin, ymin = mean_y - se_y, ymax = mean_y + se_y),
                  width = 0.1, colour = "slateblue", linewidth = 0.7,
                  inherit.aes = FALSE) +
    scale_x_continuous(breaks = 1:5, labels = q_labels_5) +
    labs(x = x_lab, y = y_lab_arg, title = chr_name) +
    theme_bw(base_size = base_sz) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
}

save_chr_patchworks <- function(dat, x_col, x_lab, stem,
                                y_col = "win_freq", y_lab_arg = "OreR frequency") {
  chr_plots <- chr_order %>%
    keep(~ nrow(filter(dat, chr == .x)) >= 10) %>%
    map(~ {
      d_chr   <- filter(dat, chr == .x)
      chr_lbl <- sprintf("%s: %.0f%% avg.", .x, 100 * mean(d_chr[[y_col]], na.rm = TRUE))
      list(
        scatter  = make_scatter( d_chr, x_col, x_lab, chr_name = chr_lbl,
                                 y_col = y_col, y_lab_arg = y_lab_arg),
        quantile = make_quantile(d_chr, x_col, x_lab, chr_name = chr_lbl,
                                 y_col = y_col, y_lab_arg = y_lab_arg)
      )
    })
  n_cols      <- ceiling(length(chr_plots) / 2)
  pw_scatter  <- wrap_plots(map(chr_plots, "scatter"),  ncol = n_cols, nrow = 2) +
    plot_layout(axis_titles = "collect")
  pw_quantile <- wrap_plots(map(chr_plots, "quantile"), ncol = n_cols, nrow = 2) +
    plot_layout(axis_titles = "collect")

  f_s <- file.path(outdir, paste0(stem, "_bychr_scatter.pdf"))
  f_q <- file.path(outdir, paste0(stem, "_bychr_quantile.pdf"))
  ggsave(f_s, pw_scatter,  width = 1.7 * n_cols, height = 4)
  ggsave(f_q, pw_quantile, width = 1.7 * n_cols, height = 4)
  message("Saved: ", f_s)
  message("Saved: ", f_q)

  r_s <- file.path(outdir, paste0(stem, "_bychr_scatter.rds"))
  r_q <- file.path(outdir, paste0(stem, "_bychr_quantile.rds"))
  saveRDS(pw_scatter,  r_s)
  saveRDS(pw_quantile, r_q)
  message("Saved: ", r_s)
  message("Saved: ", r_q)
}

n_founders <- length(chosen_founders)

# ============================================================
# Section 1: Recombination rate vs founder frequency
# ============================================================
recomb_scatter_list  <- list()
recomb_quantile_list <- list()

for (fd in all_founder_data) {
  fnd      <- fd$founder
  dat_nX   <- fd$data_noX
  y_lab    <- paste0(fnd, " frequency")

  recomb_scatter_list[[fnd]]  <- make_scatter( dat_nX, "c", "Recombination rate",
                                               chr_name = fnd, y_col = "win_freq",
                                               y_lab_arg = y_lab)
  recomb_quantile_list[[fnd]] <- make_quantile(dat_nX, "c", "Recombination rate quantile",
                                               chr_name = fnd, y_col = "win_freq",
                                               y_lab_arg = y_lab)
}

p_recomb_scatter  <- wrap_plots(recomb_scatter_list,  ncol = n_founders) +
  plot_annotation(title = "Founder frequency vs recombination rate")
p_recomb_quantile <- wrap_plots(recomb_quantile_list, ncol = n_founders) +
  plot_annotation(title = "Founder frequency vs recombination rate quantile")

ggsave(file.path(outdir, "recombVfreq_scatter.pdf"),  p_recomb_scatter,  width = 5 * n_founders, height = 4)
ggsave(file.path(outdir, "recombVfreq_quantile.pdf"), p_recomb_quantile, width = 5 * n_founders, height = 4)
saveRDS(p_recomb_quantile, file.path(outdir, "recombVfreq_quantile.rds"))
message("Saved: recombination rate plots")

for (fd in all_founder_data) {
  is_orer <- fd$founder == "Ore.Ore"
  y_lab   <- if (is_orer) "Estimated Frequency (OreR)" else paste0(fd$founder, " frequency")
  x_lab   <- if (is_orer) "Recombination Rate Quantile" else "Recombination rate"
  save_chr_patchworks(fd$data, "c", x_lab,
                      paste0("recombVfreq_", fd$founder),
                      y_col = "win_freq", y_lab_arg = y_lab)
}

# ============================================================
# Genome-wide Manhattan plots per founder
# ============================================================
chr_order_noX    <- chr_order[chr_order != "X"]
manhattan_list   <- list()

for (fd in all_founder_data) {
  fnd    <- fd$founder
  dat_nX <- fd$data_noX
  y_lab  <- paste0(fnd, " frequency")

  chr_offsets <- dat_nX %>%
    group_by(chr) %>%
    summarise(chr_len = max(midpoint, na.rm = TRUE), .groups = "drop") %>%
    arrange(match(chr, chr_order_noX)) %>%
    mutate(offset = lag(cumsum(chr_len), default = 0))

  dat_cum <- dat_nX %>%
    left_join(dplyr::select(chr_offsets, chr, offset), by = "chr") %>%
    mutate(cum_pos = midpoint + offset)

  chr_label_pos <- dat_cum %>%
    group_by(chr) %>%
    summarise(label_pos = mean(cum_pos, na.rm = TRUE), .groups = "drop")

  manhattan_list[[fnd]] <- ggplot(dat_cum, aes(x = cum_pos, y = win_freq, colour = chr)) +
    geom_point(alpha = 0.3, size = 0.8) +
    scale_colour_manual(values = rep(c("slateblue", "gray50"), length(chr_order_noX)), guide = "none") +
    scale_x_continuous(breaks = chr_label_pos$label_pos, labels = chr_label_pos$chr) +
    labs(x = "Chromosome", y = y_lab, title = fnd) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank())
}

ggsave(file.path(outdir, "freqVgenome.pdf"),
       wrap_plots(manhattan_list, ncol = 1), width = 10, height = 4 * n_founders)
message("Saved: genome-wide Manhattan plots")

# ============================================================
# Section 2: CDS density vs founder frequency
# ============================================================
cds_scatter_list  <- list()
cds_quantile_list <- list()

for (fd in all_founder_data) {
  fnd    <- fd$founder
  dat_nX <- fd$data_noX
  y_lab  <- paste0(fnd, " frequency")

  cds_scatter_list[[fnd]]  <- make_scatter( dat_nX, "n_cds", "CDS count per 100kb window",
                                            chr_name = fnd, y_col = "win_freq",
                                            y_lab_arg = y_lab)
  cds_quantile_list[[fnd]] <- make_quantile(dat_nX, "n_cds", "CDS quantile",
                                            chr_name = fnd, y_col = "win_freq",
                                            y_lab_arg = y_lab)
}

p_cds_scatter  <- wrap_plots(cds_scatter_list,  ncol = n_founders) +
  plot_annotation(title = "Founder frequency vs CDS density")
p_cds_quantile <- wrap_plots(cds_quantile_list, ncol = n_founders) +
  plot_annotation(title = "Founder frequency vs CDS density quantile")

ggsave(file.path(outdir, "cdsVfreq_scatter.pdf"),  p_cds_scatter,  width = 5 * n_founders, height = 4)
ggsave(file.path(outdir, "cdsVfreq_quantile.pdf"), p_cds_quantile, width = 5 * n_founders, height = 4)
saveRDS(p_cds_quantile, file.path(outdir, "cdsVfreq_quantile.rds"))
message("Saved: CDS density plots")

for (fd in all_founder_data) {
  y_lab <- if (fd$founder == "Ore.Ore") "Estimated Frequency (OreR)" else paste0(fd$founder, " frequency")
  save_chr_patchworks(fd$data, "n_cds", "CDS count per 100kb window",
                      paste0("cdsVfreq_", fd$founder),
                      y_col = "win_freq", y_lab_arg = y_lab)
}

# ============================================================
# Partial correlations and correlation table (per founder)
# ============================================================
fmt_p <- function(p) {
  ifelse(p < 0.001,
         formatC(p, format = "e", digits = 2),
         formatC(p, format = "f", digits = 3))
}

for (fd in all_founder_data) {
  fnd    <- fd$founder
  dat_nX <- fd$data_noX

  res_pc1 <- tryCatch(
    pcor.test(dat_nX$win_freq, dat_nX$n_cds, dat_nX$c, method = "spearman"),
    error = function(e) list(estimate = NA_real_, p.value = NA_real_))
  res_pc2 <- tryCatch(
    pcor.test(dat_nX$win_freq, dat_nX$c, dat_nX$n_cds, method = "spearman"),
    error = function(e) list(estimate = NA_real_, p.value = NA_real_))

  cat(sprintf("\n=== Founder: %s ===\n", fnd))
  cat("--- Partial Spearman (averaged replicates): frequency ~ CDS density | recomb rate ---\n")
  cat(sprintf("rho = %.4f,  p = %.3e\n", res_pc1$estimate, res_pc1$p.value))
  cat("--- Partial Spearman (averaged replicates): frequency ~ recomb rate | CDS density ---\n")
  cat(sprintf("rho = %.4f,  p = %.3e\n", res_pc2$estimate, res_pc2$p.value))

  chrs_for_table <- chr_order %>% keep(~ nrow(filter(fd$data, chr == .x)) >= 10)

  cor_stats <- map_dfr(c("All", chrs_for_table), function(ch) {
    d_avg  <- if (ch == "All") fd$data_noX else filter(fd$data, chr == ch)
    n_wins <- nrow(d_avg)

    rho_c   <- cor(d_avg$c,     d_avg$win_freq, method = "spearman", use = "complete.obs")
    rho_cds <- cor(d_avg$n_cds, d_avg$win_freq, method = "spearman", use = "complete.obs")
    p_c     <- cor.test(d_avg$c,     d_avg$win_freq, method = "spearman")$p.value
    p_cds   <- cor.test(d_avg$n_cds, d_avg$win_freq, method = "spearman")$p.value
    ci_c    <- fisher_ci(rho_c,   n_wins)
    ci_cds  <- fisher_ci(rho_cds, n_wins)

    pc_c   <- tryCatch(
      pcor.test(d_avg$win_freq, d_avg$c,     d_avg$n_cds, method = "spearman"),
      error = function(e) list(estimate = NA_real_, p.value = NA_real_))
    pc_cds <- tryCatch(
      pcor.test(d_avg$win_freq, d_avg$n_cds, d_avg$c,     method = "spearman"),
      error = function(e) list(estimate = NA_real_, p.value = NA_real_))
    ci_pc_c   <- fisher_ci(pc_c$estimate,   n_wins, gp = 1)
    ci_pc_cds <- fisher_ci(pc_cds$estimate, n_wins, gp = 1)

    data.frame(
      chr           = ch,
      c_rho         = rho_c,             c_p           = p_c,
      c_rho_lo      = ci_c[1],           c_rho_hi      = ci_c[2],
      c_pc_rho      = pc_c$estimate,     c_pc_p        = pc_c$p.value,
      c_pc_rho_lo   = ci_pc_c[1],        c_pc_rho_hi   = ci_pc_c[2],
      cds_rho       = rho_cds,           cds_p         = p_cds,
      cds_rho_lo    = ci_cds[1],         cds_rho_hi    = ci_cds[2],
      cds_pc_rho    = pc_cds$estimate,   cds_pc_p      = pc_cds$p.value,
      cds_pc_rho_lo = ci_pc_cds[1],      cds_pc_rho_hi = ci_pc_cds[2],
      stringsAsFactors = FALSE
    )
  })

  n_chr <- length(c("All", chrs_for_table))

  tbl <- bind_rows(
    cor_stats %>%
      transmute(chr,
                rho    = sprintf("%.3f [%.3f, %.3f]", c_rho,     c_rho_lo,     c_rho_hi),
                p      = fmt_p(c_p),
                pc_rho = sprintf("%.3f [%.3f, %.3f]", c_pc_rho,  c_pc_rho_lo,  c_pc_rho_hi),
                pc_p   = fmt_p(c_pc_p)),
    cor_stats %>%
      transmute(chr,
                rho    = sprintf("%.3f [%.3f, %.3f]", cds_rho,    cds_rho_lo,    cds_rho_hi),
                p      = fmt_p(cds_p),
                pc_rho = sprintf("%.3f [%.3f, %.3f]", cds_pc_rho, cds_pc_rho_lo, cds_pc_rho_hi),
                pc_p   = fmt_p(cds_pc_p))
  )

  kt <- kable(tbl, format = "latex", booktabs = TRUE, escape = FALSE,
              col.names = c("Chr", "$\\rho$ [95\\% CI]", "$p$",
                            "$\\rho_{\\text{partial}}$ [95\\% CI]", "$p_{\\text{partial}}$"),
              caption = paste0(
                "Spearman rank correlations between ", fnd, " frequency and recombination rate or CDS density. ",
                "Partial correlations control for the other covariate. ",
                "Gen ", chosen_gen, ", ", fnd, ".")) %>%
    pack_rows("Recombination rate", 1,         n_chr, bold = FALSE, italic = TRUE) %>%
    pack_rows("CDS density",        n_chr + 1, 2 * n_chr, bold = FALSE, italic = TRUE)

  out_tex <- file.path(outdir, paste0("freqVcovariates_correlations_", fnd, ".tex"))
  writeLines(kt, out_tex)
  message("Saved: ", out_tex)
}
