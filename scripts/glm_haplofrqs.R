# ============================================================================
# Script: GLM Manhattan plot and haplotype frequency trace
#
# Description:
#   Combines treatment GLM p-values with estimated haplotype frequencies for
#   a chosen founder and generation, producing a Manhattan plot aligned to a
#   haplotype frequency trace across the genome for each mitochondrial
#   background (Bei, Zim, Yak). Recombination-masked regions are shown with
#   reduced line opacity.
#
# Author: Leah Darwin
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(ggrastr)
library(patchwork)

# ============================================================
# Config
# ============================================================

##this file is too large to host on github and will need to be regenerated using the script glm.sh
glm_file  <- "data/treatment_time_repl_RR_F25.glm"

##scripts containing haplotype frequency estimates for each population
data_dir  <- "data/haplo_frq/"
outdir    <- "output"
outfile   <- file.path(outdir, "glm_haplofrqs_RR_F25.pdf")

founders       <- read.csv("foundergt.names", header = FALSE, col.names = "name")
chosen_founder <- founders$name[1]
chosen_gen     <- 25

# Chromosome set shared by both data sources; GLM output has no chr4 calls,
# so it is dropped here too to keep the cumulative x-axis aligned.
chr_order <- c("2L", "2R", "3L", "3R", "X")

color_palette <- c(
  "Bei" = "#1C448E",
  "Yak" = "#FCAB10",
  "Zim" = "#52C2BA"
)

eps <- 1e-6

# ============================================================
# Population labels
# ============================================================
tr_l <- c(rep("Bei", 20), rep("Yak", 20), rep("Zim", 20))
tr_l <- as.character(tr_l)
rep  <- c(rep(rep(1:5, each = 4), 3))
rep  <- as.character(rep)
rep  <- paste(tr_l, rep, sep = "")
time <- c(rep(c(10, 15, 25, 2), times = 15))

poplabs            <- data.frame(tr_l, rep, time)
poplabs$population <- paste("f", rownames(poplabs), sep = "")

pop_ids <- if (!is.null(chosen_gen)) {
  filter(poplabs, time == chosen_gen)$population
} else {
  poplabs$population
}

# ============================================================
# GLM data (Manhattan plot)
# ============================================================
df_treatment <- read.csv(glm_file, header = FALSE, sep = "\t")
colnames(df_treatment) <- c("CHR", "BP", "REF", "tr_l", "P", "rep", "tr_l:rep")

df_treatment <- df_treatment %>%
  filter(if_all(everything(), ~ !is.infinite(.) & !is.na(.))) %>%
  filter(CHR %in% chr_order) %>%
  na.omit() %>%
  mutate(CHR = factor(CHR, levels = chr_order))

alpha      <- 0.05
bonferonni <- alpha / nrow(df_treatment)

# ============================================================
# Haplotype frequency data
# ============================================================
pop <- pop_ids %>%
  map_dfr(~ {
    f <- file.path(data_dir, paste0(.x, ".tsv"))
    read.csv(f, sep = "\t") %>% mutate(population = .x)
  }) %>%
  unique() %>%
  filter(NSNPs > 0)

long <- pop %>%
  separate_rows(frequencies, sep = ";") %>%
  mutate(frequencies = as.numeric(frequencies) + eps) %>%
  group_by(chr, pos, population) %>%
  mutate(founder = founders$name[seq_len(n())]) %>%
  ungroup() %>%
  left_join(poplabs, by = join_by(population), relationship = "many-to-one") %>%
  filter(founder == chosen_founder, chr %in% chr_order) %>%
  mutate(
    chr  = factor(chr, levels = chr_order),
    tr_l = factor(tr_l, levels = c("Bei", "Zim", "Yak"))
  )

# Comeron 2012 recombination map (dm6 liftover): one normal-recombination
# interval per chromosome. Positions outside their chromosome's interval are
# masked by lowering line alpha, not by removing them.
chrmask <- read.csv("data/chrmask_comeron2012_dm6liftover.csv")

long <- long %>%
  left_join(chrmask, by = "chr") %>%
  mutate(masked = is.na(start) | pos < start | pos > end) %>%
  select(-start, -end)

# ============================================================
# Unified cumulative genome offsets, shared by both panels so the
# Manhattan plot and haplotype frequency traces line up on the same x-axis.
# ============================================================
chr_lengths_glm   <- df_treatment %>%
  group_by(chr = CHR) %>% summarise(m = max(BP), .groups = "drop")
chr_lengths_haplo <- long %>%
  group_by(chr) %>% summarise(m = max(pos), .groups = "drop")

chr_sizes <- bind_rows(chr_lengths_glm, chr_lengths_haplo) %>%
  group_by(chr) %>%
  summarise(chr_len = max(m), .groups = "drop") %>%
  arrange(match(chr, chr_order)) %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))

chr_bounds <- chr_sizes %>%
  mutate(end = offset + chr_len, mid = offset + chr_len / 2)

df_treatment <- df_treatment %>%
  left_join(chr_sizes %>% select(chr, offset), by = c("CHR" = "chr")) %>%
  mutate(cum_pos = (BP + offset) / 1e6)

long <- long %>%
  left_join(chr_sizes %>% select(chr, offset), by = "chr") %>%
  mutate(
    cum_pos   = (pos + offset) / 1e6,
    color_grp = if_else(masked, "masked", as.character(tr_l))
  )

x_breaks <- chr_bounds$mid / 1e6
x_labels <- chr_bounds$chr
x_vlines <- chr_bounds$end[-nrow(chr_bounds)] / 1e6

# ============================================================
# Panel 1: Manhattan plot
# ============================================================
df_treatment <- df_treatment %>%
  mutate(
    sig = P < bonferonni,
    chr_idx = as.integer(factor(CHR, levels = chr_order)),
    color_group = if_else(chr_idx %% 2 == 1, "dark", "light")
  )

n_sig <- sum(df_treatment$sig)
cat(sprintf("Significant SNPs (p < %.3e): %d\n", bonferonni, n_sig))

manhattan_plot <- ggplot(df_treatment, aes(x = cum_pos, y = -log10(P), color = color_group)) +
  geom_point_rast(size = 0.05, alpha = 0.6, raster.dpi = 700) +
  geom_hline(yintercept = -log10(bonferonni), linetype = "dotted", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = x_vlines, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  scale_color_manual(
    values = c("dark" = "#454545", "light" = "#818181"),
    guide = "none"
  ) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
  labs(x = NULL, y = expression(-log[10](italic(p)))) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# ============================================================
# Panel 2: Haplotype frequencies
# ============================================================
p_freq <- ggplot(long, aes(x = cum_pos, y = frequencies,
                           color = color_grp, group = interaction(rep, chr))) +
  rasterise(geom_line(aes(alpha = masked)), dpi = 700) +
  scale_alpha_manual(values = c(`FALSE` = 0.6, `TRUE` = 0.15), guide = "none") +
  facet_wrap(~tr_l, ncol = 1, strip.position = "left") +
  geom_vline(xintercept = x_vlines, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01, 0)) +
  theme_bw() +
  labs(
    x = "Chromosome",
    y = paste0("Estimated Frequency (", chosen_founder, ")")
  ) +
  scale_color_manual(values = c(color_palette, masked = "grey70")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black")
  ) +
  ylim(0, 1) + guides(color = "none")

# ============================================================
# Combine, aligned on cumulative genome position
# ============================================================
p_combined <- manhattan_plot/wrap_plots(p_freq) +
  plot_layout(axis_titles = "collect", ncol=1, heights=c(1.5,3)) + 
  plot_annotation(tag_level="a")

saveRDS(p_combined, paste0(outfile, ".rds"))
ggsave(p_combined, filename = outfile, width = 6, height = 4)
message("Saved: ", outfile)
