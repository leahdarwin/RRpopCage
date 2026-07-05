# RRpopCage

Scripts and data to reproduce results and figures for the paper "Do divergent mitochondrial backgrounds alter the fate of nuclear alleles during hybridization and admixture?"

## Directory Structure

```
RRpopCage/
├── data/                          # Input data files
│   ├── PupalDev_2.0.csv          # Development time measurements
│   ├── Feucndity2.0.csv          # Female fecundity counts
│   └── popSize_weight.csv        # Weight measurements for population size estimation
├── scripts/                       # Analysis scripts
│   ├── plot_phenos.R             # Phenotypic trait analysis and visualization
│   ├── plot_popSize.R            # Population size and harmonic oscillation analysis
│   └── README.md                 # Detailed description of scripts and outputs
├── output/                        # Generated figures and tables
├── RRpopCage.Rproj               # RStudio project file
└── README.md                      # This file
```

## Getting Started

### Prerequisites

- R (version 3.6 or later)
- RStudio (recommended)
- Required R packages (automatically installed if missing):
  - Data manipulation: `dplyr`, `tidyr`, `purrr`
  - Statistical modeling: `lme4`, `lmerTest`, `emmeans`
  - Visualization: `ggplot2`, `patchwork`
  - Table formatting: `knitr`, `kableExtra`
  - Model summaries: `broom`

### Reproducing the Figures

1. **Open the RProject:** Open the `RRpopCage.Rproj` file in RStudio. This automatically sets your working directory to the project root.

2. **Run the analysis scripts:**
   ```R
   source("scripts/plot_phenos.R")
   source("scripts/plot_popSize.R")
   ```

3. **Output files** will be generated in the `output/` directory:
   - `pheno_plot.pdf` — Figure of development time and fecundity by mitotype and cage
   - `emmeans_table.tex` — LaTeX table of estimated marginal means
   - `popSize.pdf` — Figure of population size dynamics with fitted harmonic models
   - `hmeans.csv` — Table of harmonic means by population

## Data Description

### `PupalDev_2.0.csv`
Development time (hours) for individuals grouped by mitotype, cage, replicate, and count.

### `Feucndity2.0.csv`
Female fecundity counts (number of offspring) by mitotype and cage.

### `popSize_weight.csv`
Population weight (grams) measurements by generation (Gen), mitotype, and cage. Population size (N) is estimated as weight × 1000.

## Script Details

For detailed information about each script's functionality, outputs, and methods, see [scripts/README.md](scripts/README.md).

## Author

Leah Darwin
