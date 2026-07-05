# ============================================================================
# Script: GLM p-value Manhattan plot
#
# Description:
#   Plots a genome-wide Manhattan plot of p-values from the treatment/time/
#   replicate GLM fit to pooled allele frequencies (glm.sh output), with a
#   Bonferroni significance threshold. Also writes the significant SNPs
#   (sorted by p-value) to a CSV.
#
# Author: Leah Darwin
# ============================================================================

library(ggplot2)
library(dplyr)
library(ggrastr)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Store args as a variable
input_file = "data/treatment_time_repl_RR.glm"
output_file = "output/glmplot.pdf"
sorted_file = "output/glm_sorted.csv"

chrs = c("2L", "2R", "3L", "3R", "X")

df_treatment = read.csv(input_file, header=FALSE, sep="\t")

##replace the term you want to plot with "P", in this case I am plotting the interaction effect
#col_names = c("CHR", "BP", "REF", "tr_l", "time", "rep", "tr_l:rep")
col_names = c("CHR", "BP", "REF", "tr_l", "time", "rep", "P")

colnames(df_treatment) = col_names

# Drop rows with non-finite p-values and restrict to the major chromosome arms
df_treatment = df_treatment %>%
	filter(if_all(everything(), ~ !is.infinite(.) & !is.na(.))) %>%
	filter(CHR %in% chrs)

##calculate sig threshold
alpha=0.05
bonferonni=alpha/nrow(df_treatment)

##sort pvalues and save to csv
df_sorted = df_treatment[order(df_treatment$P),]
df_sorted = df_sorted %>% filter(P<bonferonni)
write.csv(df_sorted,sorted_file,row.names=FALSE)

df_treatment = df_treatment %>% na.omit()

# Build cumulative base-pair positions per chromosome
df_treatment$CHR <- factor(df_treatment$CHR, levels = chrs)

# Per-chromosome offset so BP can be laid out along one continuous x-axis
chr_lengths <- df_treatment %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP), .groups = "drop") %>%
  arrange(match(CHR, chrs)) %>%
  mutate(cum_offset = lag(cumsum(as.numeric(max_bp)), default = 0))

df_treatment <- df_treatment %>%
  left_join(chr_lengths[, c("CHR", "cum_offset")], by = "CHR") %>%
  mutate(cum_bp = BP + cum_offset)

# Axis tick positions at midpoint of each chromosome
axis_df <- df_treatment %>%
  group_by(CHR) %>%
  summarise(center = (min(cum_bp) + max(cum_bp)) / 2, .groups = "drop")

# Significance flag and alternating grey palette (odd/even chromosomes get
# different greys so adjacent arms are visually distinguishable; significant
# SNPs are highlighted regardless of chromosome)
df_treatment <- df_treatment %>%
  mutate(
    sig = P < bonferonni,
    chr_idx = as.integer(factor(CHR, levels = chrs)),
    color_group = case_when(
      #sig & CHR == "2L" ~ "sig", ##optionally only highlight sig snps on chromosome 2L
      sig ~ "sig",
      chr_idx %% 2 == 1 ~ "dark",
      TRUE ~ "light"
    )
  )

# Rasterize points (geom_point_rast) since this is typically a genome-wide
# SNP-level plot with far too many points for a vector PDF to stay lightweight
manhattan_plot <- ggplot(df_treatment, aes(x = cum_bp, y = -log10(P), color = color_group)) +
  geom_point_rast(size = 0.8, alpha = 0.6, raster.dpi = 700) +
  geom_hline(yintercept = -log10(bonferonni), linetype = "dotted", color = "black", linewidth = 0.8) +
  scale_color_manual(
    values = c("sig" = "slateblue", "dark" = "#454545", "light" = "#818181"),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR,
    expand = c(0.01, 0)
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(p)))
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(output_file, plot = manhattan_plot, width = 5, height = 3)
