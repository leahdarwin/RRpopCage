# ============================================================================
# Script: Extended founder haplotype frequency plots (contamination check)
#
# Description:
#   Plots estimated haplotype frequencies per population using the extended
#   founder set (data/foundergt_extended.names), which includes founders
#   beyond OreR and DGRP-375 (w1118, Bei11, ZW142). Produces per-population
#   frequency traces and a stacked bar plot of mean founder frequency in a
#   defined region, used to check for contamination from unexpected founder
#   sources rather than just the two founders (OreR, 375) used in the main
#   nuclear-background analysis.
#
# Usage:
#   Rscript plot_haplofrqs_ext.R <data_dir> <outfile> [chosen_gen]
#
# Author: Leah Darwin
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)

# Extended founder set: OreR, DGRP-375, w1118, Bei11, ZW142
founders = read.csv("data/foundergt_extended.names", header = FALSE, col.names = "name")

##################################
## hard code treatment ordering
## {BEI,YAK,ZIM} ~ {1,2,3}
## mito = treatment = tr_l
tr_l = c(rep("Bei", 20), rep("Yak", 20), rep("Zim", 20))
tr_l = as.character(tr_l)

## hard code replicate ordering
## {B1,...,B5,Y1,...,Y5,Z1,...Z5} ~ {1,...5,6,...,10,11,...15}
rep = c(rep(rep(1:5, each = 4), 3))
rep = as.character(rep)
rep = paste(tr_l, rep, sep = "")

## hard code time point ordering
## {F2,F10} ~ {2,10}
time = c(rep(c(10, 15, 25, 2), times = 15))
##################################

poplabs = data.frame(tr_l, rep, time)
poplabs$population = paste("f", rownames(poplabs), sep = "")

print(poplabs)

color_palette <- c(
  "Ore.Ore"  = "#993d98",
  "DGRP.375" = "#ee4237",
  "w1118"    = "#E8DBC5",
  "Bei11"    = "#1C448E",
  "ZW142"    = "#52C2BA"
)

eps = 1e-6  # small offset added to frequencies to avoid exact 0s on the log/line plots

chr_order <- c("2L", "2R", "3L", "3R", "4", "X")

# Set region_chr to zoom into a single chromosome/region; leave NULL for whole genome
#region_chr   = NULL  # e.g. "2L"
#region_start = NULL  # start position in bp, e.g. 1e6
#region_end   = NULL  # end position in bp, e.g. 5e6

# Currently zoomed to a ~2Mb window padded around a region of interest on 2L
region_chr   = "2L"  # e.g. "2L"
region_start = 7164140-1000000  # start position in bp, e.g. 1e6
region_end   = 8974571+1000000  # end position in bp, e.g. 5e6


args = commandArgs(trailingOnly = TRUE)

data_dir   = args[1]
outfile    = args[2]
chosen_gen = if (length(args) >= 3) as.numeric(args[3]) else NULL

pop_ids = if (!is.null(chosen_gen)) {
  filter(poplabs, time == chosen_gen)$population
} else {
  poplabs$population
}

# Load per-population extended haplotype frequency estimates (haplo_frqs.R output, extended founder set)
pop = pop_ids %>%
  map_dfr(~ {
    f <- file.path(data_dir, paste0(.x, "_ext.tsv"))
    read.csv(f, sep = "\t") %>% mutate(population = .x)
  }) %>%
  unique() %>%
  filter(NSNPs > 0)

# Unpack the semicolon-delimited per-window frequency vector into one row per founder
long = pop %>%
  separate_rows(frequencies, sep = ";") %>%
  mutate(frequencies = as.numeric(frequencies) + eps) %>%
  group_by(chr, pos, population) %>%
  mutate(founder = founders$name[seq_len(n())]) %>%
  ungroup() %>%
  left_join(poplabs, by = join_by(population), relationship = "many-to-one") %>%
  filter(chr %in% chr_order) %>%
  mutate(
    chr  = factor(chr, levels = chr_order),
    tr_l = factor(tr_l, levels = c("Bei", "Zim", "Yak"))
  )

if (!is.null(region_chr)) {
  long <- long %>% filter(chr == region_chr)
  if (!is.null(region_start)) long <- long %>% filter(pos >= region_start)
  if (!is.null(region_end))   long <- long %>% filter(pos <= region_end)
}

print(unique(long$population))
print(unique(long$founder))

# Per-chromosome offsets for cumulative x-axis (whole-genome view only; unused when zoomed to a region)
chr_sizes <- long %>%
  group_by(chr) %>%
  summarise(chr_len = max(pos), .groups = "drop") %>%
  arrange(match(chr, chr_order)) %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))

chr_bounds <- chr_sizes %>%
  mutate(
    end = offset + chr_len,
    mid = offset + chr_len / 2
  )

chr_vlines <- chr_bounds$end[-nrow(chr_bounds)] / 1e6

x_scale <- if (!is.null(region_chr)) {
  scale_x_continuous(labels = function(x) paste0(x, " Mb"), expand = c(0.01, 0))
} else {
  scale_x_continuous(breaks = chr_bounds$mid / 1e6, labels = chr_bounds$chr, expand = c(0.01, 0))
}

long <- long %>%
  left_join(chr_sizes %>% select(chr, offset), by = "chr") %>%
  mutate(cum_pos = (pos + offset) / 1e6)

# Per-population founder frequency trace across the chosen chromosome/region
plot_frqs = function(pop, df) {
  df %>% filter(rep == pop) %>%
  ggplot(aes(x = cum_pos, y = frequencies, color = founder,
             group = interaction(founder, chr))) +
  geom_vline(xintercept = chr_vlines, linetype = "dashed",
             color = "gray50", linewidth = 0.4) +
  geom_line(alpha = 0.8) +
  theme_bw() +
  labs(
    color = "Founder",
    x = "",
    y = "",
    title = pop
  ) +
  scale_color_manual(values = color_palette) +
  x_scale +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ylim(0, 1)
}

uniqpops = unique(long$rep)
print(uniqpops)

plots = lapply(uniqpops, plot_frqs, df = long)

p = wrap_plots(plots, ncol = 5) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Stacked bar plot: average founder frequencies over the defined region, per
# population. A non-trivial contribution from founders other than OreR/375
# (w1118, Bei11, ZW142) here would indicate contamination.
avg_frq <- long %>%
  group_by(founder, rep, tr_l) %>%
  summarise(mean_freq = mean(frequencies, na.rm = TRUE), .groups = "drop") %>%
  mutate(founder = factor(founder, levels = c("DGRP.375","Ore.Ore","w1118","ZW142","Bei11")))



bar_plot <- ggplot(avg_frq, aes(x = rep, y = mean_freq, fill = founder)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ tr_l, scales = "free_x") +
  scale_fill_manual(values = color_palette,
                    labels = c("Ore.Ore" = "OreR", "DGRP.375" = "Ral375",
                               "Bei11" = "Bei", "ZW142" = "Zim", "w1118" = "w1118 (Yak)")) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text       = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle = 90, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(x = "", y = "Mean Estimated Frequency", fill = "Founder", title=paste0("Generation F",chosen_gen))

ggsave(bar_plot, filename = paste0("output/stacked_haplo_F", chosen_gen, ".pdf"), width = 5, height = 3)
saveRDS(bar_plot, paste0("output/stacked_haplo_F", chosen_gen, ".rds"))
