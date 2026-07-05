# ============================================================================
# Script: PCA of pooled-sequencing allele frequencies
#
# Description:
#   Runs a PCA on genome-wide allele frequencies (from the joined,
#   MAF-filtered sync file) across all populations and generations, then
#   plots PC1-3 against generation, colored by mitochondrial background,
#   with IQR-based outlier populations labeled by replicate.
#
# Author: Leah Darwin
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
})

##this file is too large to upload to github and can be reproduced from short reads and popoolation2 and calfreq script (https://github.com/Yiguan/popoolation2helper)
file_name   <- "data/joined.sync.MAF01.frq"
output_name <- "output/pca_F25_PC1-3.pdf"
labels_file <- "data/sync_labs.csv"

# Color palette for mitotypes (B = Bei, Z = Zim, Y = Yak)
color_palette <- c(
  "B"   = "#1C448E",
  "Z"  = "#52C2BA",
  "Y"  = "#FCAB10"
)

cat("File name:", file_name, "\n")

# Load per-site allele frequencies; drop the mitochondrial "chromosome" (chrom == "Y")
df <- read.table(file_name, sep = "\t", header = TRUE, check.names = FALSE) %>%
  filter(chrom != "Y")

cat("Loaded data.\n")
print(head(df))

# Drop character/position columns, leaving one column of frequencies per sample
df <- df[, !(names(df) %in% c("chrom", "pos", "ref"))]
print(head(df))

# Sample metadata: mitotype, generation, and replicate for each column of df
labels <- read.csv(labels_file)
print(head(labels))

rows_before <- nrow(df)
cat("Number of rows:", rows_before, "\n")

# Drop sites with missing frequency in any sample so all samples share the same site set
df <- df[complete.cases(df), ]
rows_after <- nrow(df)

dropped_rows <- rows_before - rows_after
cat("Number of dropped rows:", dropped_rows, "\n")

# Transpose so samples are rows and sites are columns (required for prcomp)
df_t <- as.data.frame(t(df))

# Center the data (scaling left off since all variables are already frequencies)
centered_data <- scale(df_t, center = TRUE, scale = FALSE)

# Only the first 6 PCs are needed downstream, so cap the decomposition there
n_comp <- 6
pca_result <- prcomp(centered_data, center = FALSE, scale. = FALSE, rank. = n_comp)

# Percent variance explained by each of the first 3 PCs, for axis labels
var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
pc1_pct <- round(var_explained[1] * 100, 2)
pc2_pct <- round(var_explained[2] * 100, 2)
pc3_pct <- round(var_explained[3] * 100, 2)

# Combine PC scores with sample metadata for plotting
scores <- as.data.frame(pca_result$x)
plot_df <- data.frame(
  PC1        = scores[, 1],
  PC2        = scores[, 2],
  PC3        = scores[, 3],
  Mito       = labels$Mito,
  Generation = labels$Generation,
  Rep        = labels$Rep
)

# Flags populations whose PC score falls outside 1.5*IQR within a generation,
# so they can be labeled by replicate on the plot. Generations/mitotypes to
# check are chosen per-PC based on where outliers are visually apparent.
iqr_outliers <- function(df, pc_col, gens, exclude_mito = NULL) {
  df %>%
    filter(Generation %in% gens) %>%
    { if (!is.null(exclude_mito)) filter(., !Mito %in% exclude_mito) else . } %>%
    group_by(Generation) %>%
    mutate(
      q1  = quantile(.data[[pc_col]], 0.25),
      q3  = quantile(.data[[pc_col]], 0.75),
      iqr = q3 - q1
    ) %>%
    filter(.data[[pc_col]] < q1 - 1.5 * iqr | .data[[pc_col]] > q3 + 1.5 * iqr) %>%
    ungroup()
}

outliers_pc1 <- iqr_outliers(plot_df, "PC1", gens = c(25))
outliers_pc2 <- iqr_outliers(plot_df, "PC2", gens = c(2, 10, 25))
outliers_pc3 <- iqr_outliers(plot_df, "PC3", gens = c(2, 10, 15, 25), exclude_mito = "Z")

# PC1 vs generation, colored by mitotype, with outlier populations labeled
p1 <- ggplot(plot_df, aes(x = as.numeric(Generation), y = PC1, color = Mito)) +
  geom_point(size = 4, position = position_jitter(width = 1, seed = 99), alpha = 0.6) +
  geom_text_repel(
    data = outliers_pc1,
    aes(x = as.numeric(Generation), y = PC1, label = Rep, color = Mito),
    size = 3, show.legend = FALSE, seed = 99
  ) +
  labs(
    x     = "Generation",
    y     = paste0("PC1 (", pc1_pct, "%)"),
    color = "Population Type (Mito)"
  ) +
  scale_color_manual(values = color_palette) +
  theme_classic(base_size = 14) +
  scale_x_continuous(breaks = c(2, 10, 15, 25))

# PC2 vs generation, colored by mitotype, with outlier populations labeled
p2 <- ggplot(plot_df, aes(x = as.numeric(Generation), y = PC2, color = Mito)) +
  geom_point(size = 4, position = position_jitter(width = 1, seed = 99), alpha = 0.6) +
  geom_text_repel(
    data = outliers_pc2,
    aes(x = as.numeric(Generation), y = PC2, label = Rep, color = Mito),
    size = 3, show.legend = FALSE, seed = 99
  ) +
  labs(
    x     = "Generation",
    y     = paste0("PC2 (", pc2_pct, "%)"),
    color = "Population Type (Mito)"
  ) +
  scale_color_manual(values = color_palette) +
  theme_classic(base_size = 14) +
  scale_x_continuous(breaks = c(2, 10, 15, 25))

# PC3 vs generation, colored by mitotype, with outlier populations labeled
p3 <- ggplot(plot_df, aes(x = as.numeric(Generation), y = PC3, color = Mito)) +
  geom_point(size = 4, position = position_jitter(width = 1, seed = 99), alpha = 0.6) +
  geom_text_repel(
    data = outliers_pc3,
    aes(x = as.numeric(Generation), y = PC3, label = Rep, color = Mito),
    size = 3, show.legend = FALSE, seed = 99
  ) +
  labs(
    x     = "Generation",
    y     = paste0("PC3 (", pc3_pct, "%)"),
    color = "Population Type (Mito)"
  ) +
  scale_color_manual(values = color_palette) +
  theme_classic(base_size = 14) +
  scale_x_continuous(breaks = c(2, 10, 15, 25))

# Combine the three PC panels into one figure with a shared, bottom-anchored legend
final <- p1 + p2 + p3 + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "bottom")
ggsave(output_name, plot = final, width = 12, height = 4)
