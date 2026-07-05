# ============================================================================
# Script: Estimate haplotype frequencies from SNPs
#
# Description:
#   Estimates founder haplotype frequencies in sliding windows along each
#   chromosome for a single population, using constrained least-squares
#   (sum-to-one, non-negative) regression of pooled SNP frequencies onto
#   founder genotypes, with distance-based weighting of SNPs within each
#   window. Takes a variant frequency table, founder genotype file, output
#   path, and population name as command-line arguments.
#
# Author: Anthony Long, Leah Darwin
# Note on authorship: This script was heavily adapted from a script written by Anthony Long and coauthor Leah Darwin only changed parameters and updated the readabillity of the script for use in this experiment please cite the original authors if using this script (https://github.com/tdlong/fly_XQTL/blob/main/scripts/haplotyper.limSolve.code.R).
# ============================================================================

library(limSolve)

##choose sigma based on SNP positions
opt_sigma = function(par,tw){

	ww = exp(-unlist(tw)/(2*par^2))
	ww = ww/sum(ww)
	NumberSites = 200
	portionOfWeight = 0.50
	# the idea is to choose sigma, such that the "NumberSites" closest sites to pos
	# account for "portionOfWeight" of the weight.  In yeast I used 50:0.50
	# the trade of is bigger and bigger windows cut down resolution so this should sort of scale
	# with haplotype block size
	# maybe in flies I can use 100:0.50...
	(sum(sort(ww,decreasing=TRUE)[NumberSites:length(ww)])-portionOfWeight)^2

}

##command line args
args = commandArgs(trailingOnly = TRUE)

frqvar_file = args[1]
founder_file = args[2]
out_file = args[3]

##population names are just f1-f120
#popNames = paste(rep("f",120),1:120,sep="") 
popNames = args[4]


##founder gt names
founders = read.table(founder_file, header=FALSE, col.names=c("names"))

##read table with variants and obeserved freqs
frqvar = read.table(frqvar_file, header=TRUE, sep="\t")

# remove NA values
n_before = nrow(frqvar)
frqvar = na.omit(frqvar)
n_after = nrow(frqvar)
cat("Dropped", n_before - n_after, "rows containing NA\n")

print(head(frqvar))

##get column positions that match founder gt names
founderCols = match(founders$names, names(frqvar))

##get column positions that match pop names
popCols = match(popNames, names(frqvar))

stepSize = 10000
winSize = 100000
chrs = c("2L","2R","3L","3R","X")

##Pre-calculate the number of windows
totalWindows = 0 
for(chr_i in chrs){
	chrPOS = frqvar$POS[frqvar$CHROM == chr_i]

	print(paste("CHR:",chr_i,"minPOS:", min(chrPOS),", maxPOS:", max(chrPOS)))

	totalWindows = totalWindows + length(seq(min(chrPOS)+winSize,max(chrPOS)-winSize,stepSize))
}

print(paste("Total windows: ",totalWindows))

##pre-allocate output table
results = data.frame(
	population = rep("", totalWindows*length(popNames)),
	chr = rep("", totalWindows*length(popNames)),
	pos = rep(0, totalWindows*length(popNames)),
	NSNPs = rep(0, totalWindows*length(popNames)),
	frequencies = rep("", totalWindows*length(popNames)),
	stringsAsFactors = FALSE
)

print(dim(results))

idx = 1
skippedWindows = 0

##calculate haplotype freqs for each chromosome 
for(chr_i in chrs){
	
	##get columns for chr_i
	chrData = frqvar[frqvar$CHROM == chr_i, ]
	maxPOS = max(chrData$POS)
	minPOS = min(chrData$POS)

	##centers of each window
	winCenters = seq(minPOS + winSize, maxPOS - winSize, stepSize)

	##estimate haplotype freqs for each window 
	for(pos_i in winCenters){
		
		##subset window 
		in_window = chrData$POS > (pos_i - winSize) & chrData$POS < (pos_i + winSize)

		##if there are fewer than 100 variants in a given window
		if(sum(in_window) < 100){ 
			skippedWindows = skippedWindows+1 
			next } 
		
		predictors = chrData[in_window, founderCols]
		tw = (chrData[in_window, "POS"] - pos_i)^2

		##drop rows with missing founder frqs
		row_ok = apply(predictors, 1, function(x) sum(is.na(x)) == 0)
		predictors = predictors[row_ok,]
		tw = tw[row_ok]

		if(nrow(predictors) < 100){ 
			skippedWindows = skippedWindows+1
                        next }

		##fit sigma automatically
		opt = optimize(opt_sigma, c(20000,100000), tw)
		sigma = opt$minimum
		weights = exp(-tw /(2*sigma^2))

		##weighted subset
		keep = weights > 0.01
		A = predictors[keep,]
		if(nrow(A) < 100){skippedWindows = skippedWindows+1
                        next }

		# Sum-to-one + non-negative constraints
		d = ncol(A)
		E = t(matrix(rep(1, d)))
		F = 1
		G = diag(rep(1, d))
		H = matrix(rep(0, d))

	

		for(pop_i in popNames){
			
#			popCol = popCols[pop_i]

			chrW = chrData[in_window, ]
			chrW = chrW[row_ok, ]
			chrW = chrW[keep, ]
			popY = chrW[, pop_i]


			fit = lsei(A=A, B=popY, E=E, F=F, G=G, H=H, verbose=FALSE)
			if(fit$IsError){ skippedWindows = skippedWindows+1
                        next }

			freqs = paste(round(as.numeric(fit$X),4), collapse=";")

			results[idx, ] = c(pop_i, chr_i, pos_i, nrow(A), freqs)
			idx = idx + 1
		}
		
	}	

	write.table(results, out_file, sep="\t", quote=FALSE, row.names=FALSE)
	
}

print(paste("Total windows:",totalWindows, "Skipped windows:",skippedWindows))
