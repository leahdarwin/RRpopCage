# Scripts

This directory contains R scripts for analyzing and visualizing phenotypic traits and population dynamics in the RRpopCage experiment.

## Scripts

### `plot_phenos.R`

**Description:** Analyzes and visualizes phenotypic traits (development time and fecundity) across different mitochondrial backgrounds (Bei, Zim, Yak) and population cages.

**Method:** Performs mixed-model ANOVA with mitotype as a fixed effect and cage/replicate as random effects. Calculates estimated marginal means for each mitotype.

**Output Files:**
- `../output/pheno_plot.pdf` — Combined boxplot figure showing development time and fecundity distributions by mitotype and cage
- `../output/emmeans_table.tex` — LaTeX table of estimated marginal means with 95% confidence intervals for both traits

### `plot_popSize.R`

**Description:** Analyzes population dynamics by fitting harmonic (sinusoidal) models to population size estimates derived from weight measurements across generations.

**Method:** Uses leave-one-out cross-validation (LOO-CV) to identify the optimal oscillation period for each cage, then fits a harmonic model with sine and cosine terms. Calculates R², p-value, and harmonic mean for each population.

**Output Files:**
- `../output/popSize.pdf` — Multi-panel figure showing observed population sizes (points) with fitted harmonic curves and 95% confidence intervals, annotated with R², p-value, and harmonic mean for each cage
- `../output/hmeans.csv` — Table of harmonic means by mitotype and cage

## Running the Scripts

To run these scripts from the project root, open the `RRpopCage.Rproj` file in RStudio to set the working directory correctly, then:

```R
source("scripts/plot_phenos.R")
source("scripts/plot_popSize.R")
```

Ensure that the required data files are present in the `data/` directory and that the `output/` directory exists before running the scripts.

## Dependencies

Both scripts require the following R packages:
- `ggplot2` — plotting
- `dplyr` — data manipulation
- `tidyr` — data tidying
- `lme4` — mixed-effects models
- `lmerTest` — ANOVA for mixed models
- `emmeans` — estimated marginal means
- `purrr` — functional programming (popSize only)
- `broom` — extracting model summaries (popSize only)
- `patchwork` — combining plots (phenos only)
- `knitr` and `kableExtra` — table formatting (phenos only)
