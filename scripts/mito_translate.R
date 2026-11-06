# ============================================================================
# Script: Translate mitochondrial CDS and count amino acid differences
#
# Description:
#   Takes the multiple sequence alignment of mitochondrial genomes made by
#   mito_consensus.sh and, for each protein-coding gene in the reference GTF,
#   translates every sequence with the Drosophila mitochondrial genetic code
#   (NCBI table 5; alternative start codons read as Met). CDS coordinates are
#   mapped through the alignment from the reference, so in-frame insertions
#   and deleted codons relative to the reference are counted as differences;
#   codons with N or frameshift gaps are ignored. Writes per-gene protein FASTAs,
#   per-gene pairwise AA difference matrices, a samples-by-genes summary of
#   differences vs. the reference, and a long table of every substitution.
#
# Usage:
#   Rscript scripts/mito_translate.R <msa.fa> <gtf> <chrom> <ref_name> <out_dir>
#
# Author: Leah Darwin
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
msa_file <- args[1]; gtf_file <- args[2]; chrom <- args[3]; ref_name <- args[4]; out_dir <- args[5]

read_fasta <- function(f) {
  l <- readLines(f)
  hdr <- which(startsWith(l, ">"))
  ends <- c(hdr[-1] - 1, length(l))
  seqs <- mapply(function(s, e) toupper(paste(l[(s + 1):e], collapse = "")), hdr, ends)
  names(seqs) <- sub("^>", "", sub("\\s.*", "", l[hdr]))
  strsplit(seqs, "")
}

aln <- read_fasta(msa_file)
stopifnot(ref_name %in% names(aln))
aln_mat <- do.call(rbind, aln)

# Map ungapped reference position -> alignment column
ref_cols <- which(aln_mat[ref_name, ] != "-")

# Drosophila mitochondrial genetic code (NCBI table 5)
bases <- c("T", "C", "A", "G")
codons <- apply(expand.grid(b3 = bases, b2 = bases, b1 = bases)[, 3:1], 1, paste, collapse = "")
aas <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]]
code <- setNames(aas, codons)
code["ATA"] <- "M"; code["TGA"] <- "W"; code["AGA"] <- "S"; code["AGG"] <- "S"

# Complement only: columns are already put in transcription order for "-" strand genes
complement <- function(x) {
  comp <- c(A = "T", C = "G", G = "C", T = "A", N = "N", "-" = "-")
  unname(ifelse(x %in% names(comp), comp[x], "N"))
}

# Alternative mito start codons are all read as Met; treat them as M so they are not counted as differences
start_codons <- c("ATT", "ATC", "ATA", "ATG", "GTG")

translate <- function(nt) {
  n <- length(nt) %/% 3
  vapply(seq_len(n), function(i) {
    cod <- paste(nt[(3 * i - 2):(3 * i)], collapse = "")
    if (grepl("[^ACGT]", cod)) "X" else code[[cod]]
  }, "")
}

# Number of AA differences between two protein vectors (codon-block cells; a cell may be several
# amino acids for an in-frame insertion, or "-" for a deleted codon). Cells containing X are ignored.
diff_count <- function(a, b) {
  ok <- !(grepl("X", a, fixed = TRUE) | grepl("X", b, fixed = TRUE))
  multi <- ok & (nchar(a) > 1 | nchar(b) > 1)
  single <- ok & !multi
  d <- sum(a[single] != b[single])
  for (k in which(multi)) {
    d <- d + if (nchar(a[k]) != nchar(b[k])) abs(nchar(a[k]) - nchar(b[k]))
             else sum(strsplit(a[k], "")[[1]] != strsplit(b[k], "")[[1]])
  }
  d
}

# CDS from GTF
gtf <- read.delim(gtf_file, header = FALSE, comment.char = "#", stringsAsFactors = FALSE, quote = "")
cds <- gtf[gtf$V1 == chrom & gtf$V3 == "CDS", ]
cds$gene <- sub('.*gene_name "([^"]+)".*', "\\1", cds$V9)
cds$tid <- sub('.*transcript_id "([^"]+)".*', "\\1", cds$V9)
genes <- unique(cds$gene)

samples <- rownames(aln_mat)
prot_all <- list()

for (g in genes) {
  # one transcript per gene (the longest), exons ordered by transcription direction
  sub_g <- cds[cds$gene == g, ]
  tid <- names(which.max(tapply(sub_g$V5 - sub_g$V4 + 1, sub_g$tid, sum)))
  e <- sub_g[sub_g$tid == tid, ]
  strand <- e$V7[1]
  e <- e[order(e$V4), ]
  # Alignment columns spanning each exon. Spans (not just reference bases) so that in-frame insertions
  # relative to the reference (reference-gap columns) are kept rather than silently dropped.
  cols <- unlist(lapply(seq_len(nrow(e)), function(i) ref_cols[e$V4[i]]:ref_cols[e$V5[i]]))
  if (strand == "-") cols <- rev(cols)
  is_ref <- aln_mat[ref_name, cols] != "-"
  nblk <- sum(is_ref) %/% 3                      # complete reference codons
  blk <- pmax(1L, ceiling(cumsum(is_ref) / 3))   # codon block of each column (insertions join a block)
  keep <- blk <= nblk
  cols <- cols[keep]; is_ref <- is_ref[keep]; blk <- blk[keep]
  bidx <- split(seq_along(blk), blk)

  prot <- t(sapply(samples, function(s) {
    nt <- aln_mat[s, cols]
    if (strand == "-") nt <- complement(nt)
    vapply(seq_len(nblk), function(b) {
      idx <- bidx[[b]]
      refnt <- nt[idx[is_ref[idx]]]              # bases at the 3 reference codon positions
      allnt <- nt[idx]; allnt <- allnt[allnt != "-"]
      if (all(refnt == "-")) return("-")         # whole codon deleted relative to reference
      aa <- if (!any(refnt == "-") && length(allnt) > 3 && length(allnt) %% 3 == 0) {
        translate(allnt)                          # in-frame insertion: extra amino acid(s)
      } else {
        translate(refnt)                          # otherwise reference codon positions only (gaps -> X)
      }
      if (b == 1 && paste(refnt, collapse = "") %in% start_codons) aa[1] <- "M"
      paste(aa, collapse = "")
    }, "")
  }))
  prot_all[[g]] <- prot
  safe <- gsub("[^A-Za-z0-9]", "_", g)

  # protein FASTA
  fa <- unlist(lapply(samples, function(s) c(paste0(">", s), gsub("-", "", paste(prot[s, ], collapse = ""), fixed = TRUE))))
  writeLines(fa, file.path(out_dir, paste0("protein_", safe, ".fa")))

  # pairwise AA difference matrix (ignores sites with X in either sequence)
  n <- length(samples)
  m <- matrix(0L, n, n, dimnames = list(samples, samples))
  for (i in seq_len(n)) for (j in seq_len(n)) m[i, j] <- diff_count(prot[i, ], prot[j, ])
  write.csv(m, file.path(out_dir, paste0("aa_diff_matrix_", safe, ".csv")))
}

# Summary: AA differences vs reference, samples x genes
summ <- sapply(prot_all, function(p) {
  r <- p[ref_name, ]
  apply(p, 1, function(x) diff_count(x, r))
})
write.csv(summ, file.path(out_dir, "aa_diff_vs_ref_by_gene.csv"))

# Long table of every AA substitution relative to the reference
long <- do.call(rbind, lapply(names(prot_all), function(g) {
  p <- prot_all[[g]]; r <- p[ref_name, ]
  do.call(rbind, lapply(setdiff(rownames(p), ref_name), function(s) {
    d <- which(p[s, ] != r & !grepl("X", p[s, ], fixed = TRUE) & !grepl("X", r, fixed = TRUE))
    if (length(d) == 0) return(NULL)
    data.frame(gene = g, sample = s, aa_pos = d, ref_aa = r[d], alt_aa = p[s, d])
  }))
}))
write.csv(long, file.path(out_dir, "aa_substitutions_vs_ref.csv"), row.names = FALSE)
cat("Done. Wrote outputs to", out_dir, "\n")
