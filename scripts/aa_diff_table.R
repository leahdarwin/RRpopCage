# ============================================================================
# Script: Mitochondrial amino acid difference table
#
# Description:
#   Reads the table of amino acid differences vs. the D. melanogaster
#   reference by mitochondrial gene (data/aa_diff_vs_ref_by_gene.csv, from
#   mito_translate.R) and writes a LaTeX (booktabs) table. The B, Y, and Z
#   population samples (Bei, Yak, and Zim mitochondrial backgrounds) are
#   largely duplicated rows, so only one representative sample is kept for
#   each prefix: the first sample carrying that prefix's most common
#   difference profile. Reference and founder mtDNA rows are kept as is.
#
# Output:
#   output/aa_diff_table.tex
# ============================================================================

library(dplyr)
library(knitr)
library(kableExtra)

aa = read.csv("data/aa_diff_vs_ref_by_gene.csv", row.names = 1, check.names = FALSE)
colnames(aa) = sub("^mt:", "", colnames(aa))

# Profile of each row, and its B/Y/Z prefix (NA for the reference/founder rows)
profile = apply(aa, 1, paste, collapse = ",")
prefix  = ifelse(grepl("^[BYZ][0-9]", rownames(aa)), substr(rownames(aa), 1, 1), NA)

# Keep the first sample having the most common profile within each B/Y/Z prefix
keep = rep(FALSE, nrow(aa))
for (p in c("B", "Y", "Z")) {
  idx = which(prefix == p)
  common = names(which.max(table(profile[idx])))
  keep[idx[profile[idx] == common][1]] = TRUE
}
keep[is.na(prefix)] = TRUE

tab = aa[keep, ] %>%
  tibble::rownames_to_column("Sample")

kable_out = kable(tab, format = "latex", booktabs = TRUE, row.names = FALSE,
                  caption = "Amino acid differences from the \\textit{D. melanogaster} reference (BDGP6.32) by mitochondrial gene. One representative sample is shown for each of the B, Y, and Z population groups.") %>%
  kable_styling(latex_options = c("scale_down", "hold_position"))

kable_out

dir.create("output", showWarnings = FALSE)
writeLines(as.character(kable_out), "output/aa_diff_table.tex")
