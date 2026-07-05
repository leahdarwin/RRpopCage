# ============================================================================
# Script: Population size estimation and harmonic oscillation analysis
#
# Description:
#   Analyzes population dynamics across generations by fitting harmonic
#   (sinusoidal) models to estimated population size data derived from
#   fly weight measurements. Uses leave-one-out cross-validation to identify
#   optimal oscillation periods. Generates plots with confidence intervals
#   and per-cage harmonic means.
#
# Author: Leah Darwin
# ============================================================================

library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(broom)

# Load population size data (derived from weight measurements)
df <- read.csv("data/popSize_weight.csv") %>%
  select(c(Mito, Cage, Gen, Weight..g.)) %>%
  mutate(
    Gen = as.integer(sub("^F0*", "", Gen)),
    N   = Weight..g. * 1000
  )

# Define grid for predictions and period candidates
gen_seq     <- seq(min(df$Gen), max(df$Gen), length.out = 300)
period_grid <- seq(2, diff(range(df$Gen)), by = 1)

# LOO-CV R² via hat matrix — avoids refitting n times per period candidate
cv_r2 <- function(d, p) {
  d$sin_term <- sin(2 * pi * d$Gen / p)
  d$cos_term <- cos(2 * pi * d$Gen / p)
  fit    <- lm(N ~ sin_term + cos_term, data = d)
  h      <- pmin(hatvalues(fit), 0.999)
  e      <- residuals(fit)
  ss_cv  <- sum((e / (1 - h))^2)
  ss_tot <- sum((d$N - mean(d$N))^2)
  1 - ss_cv / ss_tot
}

# Find period with highest LOO-CV R²
best_period <- function(d) {
  r2s <- map_dbl(period_grid, ~ cv_r2(d, .x))
  period_grid[which.max(r2s)]
}

# Fit harmonic model with given period
fit_with_period <- function(d, p) {
  d$sin_term <- sin(2 * pi * d$Gen / p)
  d$cos_term <- cos(2 * pi * d$Gen / p)
  lm(N ~ sin_term + cos_term, data = d)
}

# Generate predictions with confidence intervals for harmonic model
pred_with_period <- function(model, p) {
  newdata <- data.frame(
    sin_term = sin(2 * pi * gen_seq / p),
    cos_term = cos(2 * pi * gen_seq / p)
  )
  pred <- predict(model, newdata = newdata, interval = "confidence")
  data.frame(Gen = gen_seq, fit = pred[, "fit"], lwr = pred[, "lwr"], upr = pred[, "upr"])
}

# Fit harmonic models for each mitotype-cage combination
models <- df %>%
  group_by(Mito, Cage) %>%
  nest() %>%
  mutate(
    period = map_dbl(data, best_period),
    model  = map2(data, period, fit_with_period)
  )

# Extract model fit metrics
fit_metrics <- models %>%
  mutate(metrics = map(model, glance)) %>%
  select(Mito, Cage, period, metrics) %>%
  unnest(metrics) %>%
  select(Mito, Cage, period, r.squared, adj.r.squared, sigma, p.value)

# Generate fitted curves with confidence intervals
fits <- models %>%
  mutate(fitted = map2(model, period, pred_with_period)) %>%
  select(Mito, Cage, fitted) %>%
  unnest(fitted)

# Order panels by mitotype then cage name for consistent display
mito_order   <- c("Bei", "Zim", "Yak")
panel_levels <- df %>%
  distinct(Mito, Cage) %>%
  mutate(Mito = factor(Mito, levels = mito_order)) %>%
  arrange(Mito, Cage) %>%
  mutate(panel = paste(Mito, Cage)) %>%
  pull(panel)

df   <- df   %>% mutate(panel = factor(paste(Mito, Cage), levels = panel_levels))
fits <- fits %>% mutate(panel = factor(paste(Mito, Cage), levels = panel_levels))

# Calculate mean of fitted harmonic curve per cage
hmeans <- fits %>%
  group_by(Mito, Cage) %>%
  summarise(hmean_N = mean(fit), .groups = "drop")

# Prepare per-panel annotations with model statistics
annot <- fit_metrics %>%
  left_join(hmeans, by = c("Mito", "Cage")) %>%
  mutate(
    panel = factor(paste(Mito, Cage), levels = panel_levels),
    label = paste0(
      "R²=", round(r.squared, 2),
      "  p=", signif(p.value, 2),
      "\nHmean=", round(hmean_N)
    )
  )

# Color palette for mitotypes
color_palette <- c(
  "Bei" = "#1C448E",
  "Zim" = "#52C2BA",
  "Yak" = "#FCAB10"
)

# Create main plot with observed data and fitted harmonic curves
p = ggplot() +
  geom_point(
    data = df,
    aes(x = Gen, y = N, color = Mito),
    alpha = 0.4, size = 1.5
  ) +
  geom_ribbon(
    data = fits,
    aes(x = Gen, ymin = lwr, ymax = upr, fill = Mito),
    alpha = 0.2
  ) +
  geom_line(
    data = fits,
    aes(x = Gen, y = fit, color = Mito)
  ) +
  geom_text(
    data  = annot,
    aes(label = label),
    x     = -Inf, y = -Inf,
    hjust = -0.05, vjust = -0.3,
    size  = 2.8, color = "black"
  ) +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~ panel, nrow=3) +
  labs(x = "Generation", y = "Estimated N") +
  theme_linedraw() +
  theme(legend.position = "none",
        strip.text = element_text(color = "black"),
        strip.background = element_blank(),
        panel.grid = element_line(color = "gray90", linewidth = 0.5))

# Save plot and harmonic means
ggsave("output/popSize.pdf", p, width=7, height=5)
write.csv(hmeans, "output/hmeans.csv", row.names = F, quote = F)
