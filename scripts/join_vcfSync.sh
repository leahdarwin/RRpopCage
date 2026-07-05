#!/bin/bash

# ============================================================================
# Script: VCF-to-table conversion and sync-frequency join
#
# Description:
#   Slurm batch script that subsets a VCF to a sample list and converts
#   genotypes to a CHROM/POS/REF/ALT table with a recoded numeric GT per
#   sample (1 = homozygous REF, 0 = homozygous ALT, NA = missing/
#   heterozygous), then joins that variant table with a pooled-sequencing
#   allele frequency (.frq) file on CHROM/POS, keeping only sites present
#   in both files.
#
# Author: Leah Darwin
# ============================================================================

#SBATCH -p batch
#SBATCH --mem=10G
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o %j.out
#SBATCH -e %j.err

##Computer cluster specific module loads
module load bcftools

# ----------------------------------
# Step 1: VCF -> variant table (formerly vcf2tab.sh)
# ----------------------------------
##this file is too large to upload to github, to reproduce you must call variants using short read data of the inbred lines 
VCF="data/filtered_snps.pass.vcf"
##subset of lines particular to this study that were jointly called with a larger panel of lines, use vcfsubsetNames.txt for additional lines in study
NAMES="vcfsubsetNames.txt"

##temp file — holds sample-subset VCF
TMP=$(mktemp --suffix=.vcf)
bcftools view -S "$NAMES" "$VCF" -o "$TMP"

VAR_FILE="data/filtered_var.table"

##Convert vcf to table
printf "CHROM\tPOS\tREF\tALT" > "$VAR_FILE"
bcftools query -l "$TMP" | while read S; do
    printf "\t%s" "$S" >> "$VAR_FILE"
done
printf "\n" >> "$VAR_FILE"


# Extract CHROM POS REF ALT and raw numeric GT for each sample
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$TMP" \
| awk -F'\t' 'BEGIN{OFS="\t"}
NR==1 {
    # bcftools does not output a header automatically, so we must synthesize one.
    # This is optional — you can skip if you don’t want a header.
    next
}
{
    n = split($0, a, FS)

    # Process each GT beginning in column 5
    for(i=5; i<=n; i++){
        gt = a[i]

        # Handle missing or weird
        if(gt=="./." || gt==".|." || gt=="." || gt==""){
            a[i] = "NA"
            continue
        }

        # Normalize phased -> unphased
        gsub(/\|/, "/", gt)
        split(gt, al, "/")

        # Any missing allele → NA
        if(al[1]=="." || al[2]=="."){
            a[i] = "NA"
            continue
        }

        # Both alleles REF (0/0) → 1
        if(al[1]=="0" && al[2]=="0"){
            a[i] = 1
        } else if(al[1] != al[2]){
            # Heterozygous → missing
            a[i] = "NA"
        } else {
            # Homozygous ALT (e.g. 1/1) → 0
            a[i] = 0
        }
    }

    # Print processed row
    for(i=1; i<=n; i++){
        printf "%s", a[i]
        if(i<n) printf OFS; else printf ORS
    }
}' >> "$VAR_FILE"

rm -f "$TMP"


# ----------------------------------
# Step 2: join variant table with FRQ file 
# ----------------------------------
 ##created with calfreq script (https://github.com/Yiguan/popoolation2helper)
FRQ_FILE="joined.sync.MAF01.frq"
OUT_FILE="data/var_frq_ext.tsv"

awk '
BEGIN { FS=OFS="\t" }

# ----------------------------------
# First file = FRQ
# ----------------------------------
FNR==1 && NR==1 {
    frq_header = $0
    next
}

NR==FNR {
    key = $1 FS $2
    frq[key] = $0
    next
}

# ----------------------------------
# Second file = VAR
# ----------------------------------
FNR==1 {
    var_header = $0

    # Build joined header:
    split(frq_header, fH, FS)
    split(var_header, vH, FS)

    header = "CHROM" OFS "POS"

    for (i=3; i<=length(fH); i++) header = header OFS fH[i]
    for (i=3; i<=length(vH); i++) header = header OFS vH[i]

    print header > outfile
    next
}

{
    key = $1 FS $2

    # Only print if FRQ AND VAR both have this key
    if (!(key in frq)) next

    # Split lines
    split(frq[key], fA, FS)
    split($0, vA, FS)

    # Start with CHROM POS
    out = vA[1] OFS vA[2]

    # FRQ columns 3+
    for (i=3; i<=length(fA); i++) out = out OFS fA[i]

    # VAR columns 3+
    for (i=3; i<=length(vA); i++) out = out OFS vA[i]

    print out > outfile
}
' outfile="$OUT_FILE" "$FRQ_FILE" "$VAR_FILE"
