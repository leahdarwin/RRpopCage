# ============================================================================
# Script: Allele frequency trajectories and change at significant 2L SNPs
#
# Description:
#   Subsets the pooled allele frequency table to the GLM-significant SNPs
#   (plot_glm.R output) on chromosome 2L, then visualizes per-SNP frequency
#   trajectories across generations and the distribution of frequency change
#   (F25-F2) by mitochondrial background.
#
# Author: Leah Darwin
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# Significant SNPs (Bonferroni-filtered GLM results, see plot_glm.R)
sig = read.csv("data/treatment_time_repl_RR_F25_sig.csv")

# Restrict pooled-sequencing allele frequencies to those significant SNPs on 2L
frq = read.csv("/users/drand/data/RR_popcage_poolseq/aligned_reads_6.32/joined.sync.MAF01.frq", sep="\t") %>%
    right_join(sig, by= join_by("chrom"=="CHR", "pos"=="BP")) %>%
    filter(chrom == "2L")

##################################
## hard code treatment ordering
## {BEI,YAK,ZIM} ~ {1,2,3}
tr_l = c(rep("Bei", 20), rep("Yak", 20), rep("Zim", 20))
tr_l = as.character(tr_l)

## hard code replicate ordering
rep = c(rep(rep(1:5, each = 4), 3))
rep = as.character(rep)
rep = paste(tr_l, rep, sep = "")

## hard code time point ordering
## {F2,F10} ~ {2,10}
time = c(rep(c(10, 15, 25, 2), times = 15))
##################################

color_palette <- c(
  "Bei" = "#1C448E",
  "Yak" = "#FCAB10",
  "Zim" = "#52C2BA"
)


poplabs = data.frame(tr_l, rep, time)
poplabs$population = paste("f", rownames(poplabs), sep = "")

# Reshape to one row per SNP/population/generation, with population metadata attached
frq_long <- frq %>%
  pivot_longer(
    cols = starts_with("f"),
    names_to = "population",
    values_to = "frequency"
  ) %>%
  left_join(poplabs, by = "population")

## Repolarize SNPs that are on average decreasing in Zim across replicates
## so that Zim-specific changes are uniformly increasing
snps_to_flip <- frq_long %>%
  filter(tr_l.y == "Zim", time.y %in% c(2, 25)) %>%
  pivot_wider(
    id_cols    = c(chrom, pos, rep.y),
    names_from = time.y,
    values_from = frequency,
    names_prefix = "t"
  ) %>%
  group_by(chrom, pos) %>%
  summarise(mean_zim_delta = mean(t25 - t2, na.rm = TRUE), .groups = "drop") %>%
  filter(mean_zim_delta < 0) %>%
  select(chrom, pos)

frq_long <- frq_long %>%
  left_join(snps_to_flip %>% mutate(flip = TRUE), by = c("chrom", "pos")) %>%
  mutate(frequency = if_else(!is.na(flip), 1 - frequency, frequency)) %>%
  select(-flip)

# Per-SNP, per-replicate frequency trajectory across generations
snp_traj <- ggplot(frq_long %>% filter(!is.na(frequency)), aes(
    x     = time.y,
    y     = frequency,
    group = interaction(chrom, pos, rep.y),
    color = tr_l.y
  )) +
  geom_line(alpha = 0.3, linewidth = 0.4) +
  scale_color_manual(values = color_palette) +
  scale_x_continuous(breaks = c(2, 10, 15, 25)) +
  labs(x = "Generation", y = "Allele Frequency", color = "Treatment") +
  theme_bw()

ggsave("snp_trajectories.pdf", snp_traj, width = 8, height = 6)

# Frequency change per SNP/replicate between the first (F2) and last (F25) generation
delta_frq <- frq_long %>%
  filter(time.y %in% c(2, 25)) %>%
  pivot_wider(
    id_cols = c(chrom, pos, rep.y, tr_l.y),
    names_from  = time.y,
    values_from = frequency,
    names_prefix = "t"
  ) %>%
  mutate(delta = t25 - t2) %>%
  filter(is.finite(delta))

delta_frq$tr_l.y <- factor(delta_frq$tr_l.y, levels = c("Bei", "Zim", "Yak"))

# Distribution of |F25-F2| frequency change by mitochondrial background
combined <- ggplot(delta_frq, aes(x = abs(delta), fill = tr_l.y)) +
  geom_histogram(aes(color = after_scale(alpha(fill, 0.9))), bins = 20, alpha = 0.4, position = "identity") +
  scale_fill_manual(values = color_palette) +
  labs(x = "|F25-F2| Frequencies", y = "Count", fill = "Population mtDNA") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "top",
        axis.text = element_text(size = 13), axis.title = element_text(size = 14))

ggsave("output/delta_histograms.pdf", combined, width = 4.5, height = 4.5)
