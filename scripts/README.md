# Scripts

This directory contains R and shell scripts for analyzing and visualizing phenotypic traits, population dynamics, and allele frequency GLM results in the RRpopCage experiment.

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

### `oreVcovariates.R`

**Description:** Examines relationships between estimated founder haplotype frequencies (from `haplo_frqs.R` output) and two genomic covariates, recombination rate (Comeron 2012 map) and CDS density, at a chosen generation.

**Method:** Joins per-population haplotype frequencies to 100kb recombination-map windows and CDS counts (from a Flybase GTF), then produces scatter and quantile plots of founder frequency vs. each covariate (by chromosome and genome-wide), plus Spearman and partial Spearman rank correlations controlling for the other covariate.

**Output Files:**
- `output/recombVfreq_scatter.pdf`, `output/recombVfreq_quantile.pdf` — Founder frequency vs. recombination rate, all founders
- `output/recombVfreq_<founder>_bychr_scatter.pdf`, `output/recombVfreq_<founder>_bychr_quantile.pdf` — Per-chromosome panels for a given founder
- `output/freqVgenome.pdf` — Genome-wide Manhattan-style plot of founder frequency per chromosome
- `output/cdsVfreq_scatter.pdf`, `output/cdsVfreq_quantile.pdf` — Founder frequency vs. CDS density, all founders
- `output/cdsVfreq_<founder>_bychr_scatter.pdf`, `output/cdsVfreq_<founder>_bychr_quantile.pdf` — Per-chromosome panels for a given founder
- `output/freqVcovariates_correlations_<founder>.tex` — LaTeX table of Spearman and partial Spearman correlations by chromosome
- Corresponding `.rds` files for saved plot objects

### `pca.R`

**Description:** Runs a PCA on genome-wide allele frequencies across all populations and generations (from the joined, MAF-filtered pooled-sequencing sync file) and plots the first three PCs against generation, colored by mitochondrial background.

**Method:** Transposes the site-by-sample frequency table so samples are rows, mean-centers it, and computes the first 6 principal components with `prcomp`. Populations whose PC1-3 scores fall outside 1.5×IQR within a generation are labeled by replicate on the corresponding panel.

**Output Files:**
- `output/pca_F25_PC1-3.pdf` — Three-panel figure of PC1, PC2, and PC3 vs. generation, colored by mitotype

### `join_vcfSync.sh`

**Description:** Slurm batch script that builds the sample-level variant table used by `haplo_frqs.R` and `poolFreqDiff`. Combines the former `vcf2tab.sh` and `join_var_sync.sh` scripts into one job.

**Method:** Subsets `filtered_snps.pass.vcf` to the samples listed in `vcfsubsetNames_ext.txt` with `bcftools`, converts each sample's genotype to a numeric call (1 = homozygous REF, 0 = homozygous ALT, NA = missing/heterozygous), then joins the resulting variant table with the pooled-sequencing `joined.sync.MAF01.frq` file on CHROM/POS, keeping only sites present in both.

**Output Files:**
- `filtered_var_ext.table` — Intermediate CHROM/POS/REF/ALT table with recoded genotypes per sample
- `/users/drand/data/RR_popcage_poolseq/aligned_reads_6.32/var_frq_ext.tsv` — Joined variant/frequency table

### `glm.sh`

**Description:** Slurm batch script that fits the treatment/time/replicate GLM to pooled allele frequencies using `poolFreqDiff`.

**Method:** Runs the `poolFreqDiff_treatment_time_repl_RR.py` script (py27 conda environment) on the joined MAF-filtered sync file to generate an R script, then executes that script with `Rscript` to produce the GLM results table.

**Output Files:**
- `output/treatment_time_repl_RR.glm` — GLM results table (CHR, BP, REF, treatment, p-value, replicate, interaction terms)

### `plot_glm.R`

**Description:** Plots a genome-wide Manhattan plot of p-values from the treatment/time/replicate GLM (`glm.sh` output) for a chosen model term, with a Bonferroni significance threshold.

**Method:** Lays out SNPs from the major chromosome arms (2L, 2R, 3L, 3R, X) along one cumulative x-axis, computes a Bonferroni threshold from the number of tested SNPs, and rasterizes points for a lightweight PDF. Also writes SNPs below the threshold, sorted by p-value, to a CSV.

**Output Files:**
- `output/glmplot.pdf` — Genome-wide Manhattan plot of -log10(p)
- `output/glm_sorted.csv` — Significant SNPs (p < Bonferroni threshold), sorted by p-value

### `subsetFRQ_changes_2L.R`

**Description:** Subsets the pooled allele frequency table to the GLM-significant SNPs (from `plot_glm.R`'s output) on chromosome 2L, then visualizes per-SNP frequency trajectories across generations and the distribution of frequency change (F25-F2) by mitochondrial background.

**Method:** Joins significant SNPs to the pooled frequency table, repolarizes SNPs that are on average decreasing in Zim across replicates so Zim-specific changes are uniformly increasing, then plots per-SNP/replicate trajectories over generations and a histogram of absolute F25-F2 change by mitotype.

**Output Files:**
- `snp_trajectories.pdf` — Per-SNP, per-replicate allele frequency trajectories across generations
- `output/delta_histograms.pdf` — Histogram of |F25-F2| frequency change, by mitochondrial background

### `glm_haplofrqs.R`

**Description:** Combines the treatment GLM results with estimated haplotype frequencies for a chosen founder and generation, producing a genome-wide Manhattan plot aligned to a haplotype frequency trace for each mitochondrial background (Bei, Zim, Yak).

**Method:** Aligns GLM p-values and haplotype frequency traces on a shared cumulative genome position. Recombination-cold regions (per the Comeron 2012 map, dm6 liftover) are shown with reduced line opacity rather than removed.

**Output Files:**
- `../output/glm_haplofrqs_RR_F25.pdf` — Combined Manhattan plot and haplotype frequency trace figure
- `../output/glm_haplofrqs_RR_F25.pdf.rds` — Saved combined plot object

### `plot_haplofrqs_ext.R`

**Description:** Plots estimated haplotype frequencies using the extended founder set (OreR, DGRP-375, w1118, Bei11, ZW142). This produces plots for additional haplotype frequencies beyond OreR and 375, to test for contamination from other founder sources.

**Method:** For a chosen generation and data directory of extended haplotype frequency estimates, produces per-population frequency traces (whole genome or zoomed to a specified region) and a stacked bar plot of mean founder frequency per population; a non-trivial contribution from founders other than OreR/375 flags possible contamination.

**Usage:** `Rscript plot_haplofrqs_ext.R <data_dir> <outfile> [chosen_gen]`

**Output Files:**
- `output/stacked_haplo_F<chosen_gen>.pdf` — Stacked bar plot of mean founder frequency by population, faceted by mitochondrial background
- `output/stacked_haplo_F<chosen_gen>.rds` — Saved bar plot object

### `haplo_frqs.sh`

**Description:** Slurm job array script that submits one job per population (f1-f60) to estimate founder haplotype frequencies via `haplo_frqs.R`.

**Method:** Each array task runs `haplo_frqs.R` on the shared variant frequency table for a single population, keyed by the Slurm array task ID.

**Output Files:**
- `data/haplo_frq/<population>.tsv` — Per-population table of haplotype frequency estimates by chromosome and window

### `haplo_frqs.R`

**Description:** Estimates founder haplotype frequencies in sliding windows along each chromosome for a single population. Adapted from a script by Anthony Long ([fly_XQTL](https://github.com/tdlong/fly_XQTL/blob/main/scripts/haplotyper.limSolve.code.R)); please cite the original authors if reusing this script.

**Method:** Fits a constrained least-squares regression (sum-to-one, non-negative) of pooled SNP frequencies onto founder genotypes within each window, using distance-based weighting of SNPs and an automatically optimized weighting bandwidth (sigma).

**Output Files:**
- `data/haplo_frq/<population>.tsv` — Table of haplotype frequency estimates (population, chr, pos, NSNPs, frequencies) for the given population

## Running the Scripts

To run these scripts from the project root, open the `RRpopCage.Rproj` file in RStudio to set the working directory correctly, then:

```R
source("scripts/plot_phenos.R")
source("scripts/plot_popSize.R")
```

Ensure that the required data files are present in the `data/` directory and that the `output/` directory exists before running the scripts.

## Dependencies

The R scripts require the following packages:
- `ggplot2` — plotting
- `dplyr` — data manipulation
- `tidyr` — data tidying
- `lme4` — mixed-effects models
- `lmerTest` — ANOVA for mixed models
- `emmeans` — estimated marginal means
- `purrr` — functional programming (popSize, glm_haplofrqs, oreVcovariates)
- `broom` — extracting model summaries (popSize only)
- `patchwork` — combining plots (phenos, glm_haplofrqs, oreVcovariates)
- `knitr` and `kableExtra` — table formatting (phenos, oreVcovariates)
- `stringr` — string handling (glm_haplofrqs, oreVcovariates)
- `ggrastr` — rasterized plotting (plot_glm, glm_haplofrqs)
- `ppcor` — partial correlations (oreVcovariates only)
- `ggrepel` — outlier labels (pca only)

- `limSolve` — constrained least-squares regression (haplo_frqs.R only)

`join_vcfSync.sh`, `glm.sh`, and `haplo_frqs.sh` are Slurm batch scripts; `join_vcfSync.sh` requires `bcftools`, `glm.sh` requires `poolFreqDiff` (py27 conda environment) and R, `haplo_frqs.sh` requires R.
