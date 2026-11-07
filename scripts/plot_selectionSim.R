# ============================================================================
# Script: Selection simulation p-value boxplots
#
# Description:
#   Reads the per-selection-coefficient p-value files written by selectionSim.sh
#   (output/selectionSim/pvals_s<s>.csv; columns sim, s, p_value) and draws one
#   boxplot of -log10(p) for each value of s. A horizontal dotted line marks the
#   same Bonferroni threshold used in plot_glm.R (alpha / number of tests), where
#   the number of tests is the number of lines in the observed GLM results file.
#
# Usage (from the project root):
#   Rscript scripts/plot_selectionSim.R
#
# Author: Leah Darwin
# ============================================================================

library(ggplot2)
library(dplyr)

pval_dir    <- "output/selectionSim"
output_file <- "output/selectionSim_boxplot.pdf"

# ---- Significance threshold (same definition as plot_glm.R) ----------------
n_tests <- 5902239 
alpha <- 0.05
bonferonni <- alpha / n_tests

# ---- Read p-value files -----------------------------------------------------
# Non-recursive listing of the single output directory
files <- list.files(pval_dir, pattern = "^pvals_s.*\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No pvals_s*.csv files found in ", pval_dir)

# s is taken from the second column (it matches the file name)
df <- bind_rows(lapply(files, read.csv))
colnames(df)[1:3] <- c("sim", "s", "p_value")

# Drop p-values that cannot be shown on a -log10 scale (NA or exactly 0)
n_bad <- sum(!is.finite(-log10(df$p_value)))
if (n_bad > 0) message(n_bad, " non-finite -log10(p) values dropped")
df <- df %>% filter(is.finite(-log10(p_value)))

# Order s numerically so the boxes run from weakest to strongest selection
df$s <- factor(df$s, levels = sort(unique(df$s)))

# ---- Plot -------------------------------------------------------------------
p <- ggplot(df, aes(x = s, y = -log10(p_value))) +
  geom_boxplot(fill = "grey90", outlier.size = 0.5, outlier.alpha = 0.4) +
  geom_hline(yintercept = -log10(bonferonni), linetype = "dotted",
             color = "black", linewidth = 0.8) +
  labs(
    x = "Selection coefficient (s)",
    y = expression(-log[10](italic(p)))
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10))

ggsave(output_file, plot = p, width = 5, height = 3)
