# Calculate Fst
# Input data: Matrix, SNP per row, Sample per column (with a population prefix)
# Allele calls: Major as 0, minor as 1, mixed allowed
# Populations must have a prefix which matches given population IDs
# 14 Mar 2016

# Resolve NA and -Inf issues!

# Libraries
library(gdata)

# Controls
outPrefix <- 'DEFAULTOUT'
dat <- read.table('INFILE', header=T) # Whole dataset
populations <- c('POPA','POPB','POPC','POPD') # List of population prefix IDs
allGroups <- list(c('POPA','POPB'), c('POPC','POPD'))
nameGroups <- c('GRPA','GRPB')

# Functions
# Per-population calculations
# Provide the full dataset, return a data frame of the B and H stats for each population
# Get B allele frequency
# Get H (mean of A and B allele frequency)
calcPopStats <- function(popDat, pops) {
	print('Calculating population statistics')
	# Reset A, B & H holders
	out.A <- c()
	out.B <- c()
	out.H <- c()

	for (i in pops) {
		print(paste('Calculating for ',i,sep=''))
		# Cut dataset to population
		popX <- popDat[,startsWith(names(popDat), i)]

		# Calculate B & H stats for that population
		popX.B <- rowSums(popX)/dim(popX)[2]
		popX.A <- 1 - popX.B
		popX.H <- 2*popX.A*popX.B

		# Add A, B & H stats to output dataframe
		# Format as population per column & snp/gene per row
		# NB: These are being assigned globally
		out.A <- cbind(out.A, tmp=popX.A)
		dimnames(out.A)[[2]][which(pops==i)] <- i
		out.B <- cbind(out.B, tmp=popX.B)
		dimnames(out.B)[[2]][which(pops==i)] <- i
		out.H <- cbind(out.H, tmp=popX.H)
		dimnames(out.H)[[2]][which(pops==i)] <- i
	}
	assign('out.A', out.A, envir=.GlobalEnv)
	assign('out.B', out.B, envir=.GlobalEnv)
	assign('out.H', out.H, envir=.GlobalEnv)
}

# Population group calculations
calcGroupFst <- function(pops.H, pops.A, pops.B, groups, groupNames) { # groups is a list of arrays: list(c(popA,popB), c(popC), c(popD,popE,popF))
	print('Calculating group Fst values')
	# From per-population A, B & H stats, calculate Fst.
	out.Fst <- data.frame()
	if (TRUE) { # Whole population Fst
		print('Calculating whole group Fst')
		Hs <- rowSums(pops.H)/dim(pops.H)[2]
		At <- rowSums(pops.A)/dim(pops.A)[2]
		Bt <- rowSums(pops.B)/dim(pops.B)[2]
		Ht <- 2 * At * Bt
		Fst <- (Ht - Hs) / Ht
		write.table(Fst, file=paste(outPrefix,'_whole.Fst',sep=''))
	}
	if (TRUE) { # Pre-group Fst
		for (i in groups) {
			print('Calculating Fst for subgroup')
			i.A <- pops.A[,(names(pops.A) %in% i)]
			i.B <- pops.B[,(names(pops.B) %in% i)]
			i.H <- pops.H[,(names(pops.H) %in% i)]

			Hs <- rowSums(i.H)/dim(i.H)[2]
			At <- rowSums(i.A)/dim(i.A)[2]
			Bt <- rowSums(i.B)/dim(i.B)[2]
			Ht <- 2 * At * Bt
			Fst <- (Ht - Hs) / Ht

			out.Fst <- cbind(out.Fst, groupFst=Fst)
		}
	}
	names(out.Fst) <- groupNames
	write.table(out.Fst, file=paste(outPrefix,'_groups.Fst',sep=''))
}

# Pairwise
calcPairwiseFst <- function() {
	print('Calculating pairwise Fst')
	# REF v otherPops
}

# Default run
calcPopStats(dat, populations)
calcGroupFst(out.H, out.A, out.B, allGroups, nameGroups)
