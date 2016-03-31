# Calculate Tajima's D for all populations, write output to file.
# TODO: Add plots post-output.
# Input matrix must have samples per column and SNPs per row.
# A column (dat_annot.Gene) must also exist for specifying genes, eventually I'll make it possible to just parse a given column name.
# Note that non-population columns can be included, these will be ignored.
# Individual population should be identifiable from a population-specific prefix within that column name (e.g. POPA_ID1)

# Matt Ravenhall
# 24th Feb 2016

library(gdata)

# Options
pops <- c('POPA','POPB','POPC')
minSNPs <- 4 # Ignore genes with fewer than this many SNPs
inFile <- 'FILEPATH'
geneCol <- 'dat_annot.Gene'

# Read in matrix
print('Parsing matrix')
dat.full <- read.table(inFile, header=T)

# Determine complete list of genes
print('Identifying gene names')
uniqGenes <- unique(dat.full[,geneCol])
uniqGenes <- uniqGenes[uniqGenes != '.']

# Set initial D.df
D.df <- data.frame(Gene=uniqGenes)

# Iterate through populations
print('Iterating through populations')
print('=============================================')
for (pop in pops) {

	# Set initial dT & dW arrays
	dT <- c()
	dW <- c()

	# Subset dat.full to population under consideration
	print(paste('Analysing ',pop,' (',which(pops == pop) ,' of ',length(pops),')',sep=''))
	#dat.pop <- cbind(dat.full$dat_annot.Gene,dat.full[,startsWith(names(dat.full),pop)])

	# Reset the skipLog
	skipLog <- c()

	# Iterate through genes
	print('Iterating through genes')
	print('---------------------------------------------')
	for (gene in uniqGenes) {
			# Subset to that gene
			print(paste('Subsetting for ',gene,' (',which(uniqGenes == gene),' of ',length(uniqGenes),')',sep=''))
			dat.gene <- dat.full[dat.full$dat_annot.Gene == gene,] ## TODO: Take column name from input
			# Subset to population
			dat.gene <- dat.gene[,startsWith(names(dat.gene),pop)]

			# DECISION NEEDED: Dealing with missing calls. Currently omitting those SNPs, should these be assumed as variant or reference?
			print('Omitting missing calls')
			dat.gene <- na.omit(dat.gene)

		# Avoid genes where no SNPs exist
		if (dim(dat.gene)[1] > minSNPs) {
		## Tajima's estimator
			print("Calculating Tajima's estimator, dT")
			n <- dim(dat.gene)[2] # Number of sequences/individuals

			PWdiff <- 0 # Pairwise differences between sequences, e.g. (n diff sites a&b + n diff sites a&c + n diff sites b&c)

			for (index1 in (1:n)) {
				for (index2 in (1:n)) {
					if (index2 > index1) {
						PWdiff <- PWdiff + sum(dat.gene[index1] != dat.gene[index2])
					}
				}
			}

		dT <- c(dT, PWdiff/(n*(n-1)/2)) # Number of pairwise differences (Tajima's estimator), add to population array

		## Watterson's estimator
			print("Calculating Watterson's estimator, dW")
			# Number of segregating sites in population
			# Sample sum of 0/1 calls divided by total number of samples. Segregating sites will be non-0/1. All mixed-call SNPs are ignored.
			segSites = rowSums(dat.gene) / dim(dat.gene)[2]
			segSites = segSites[!(segSites %in% c(0,1))]

#			harmonic <- 0
#			for (denominator in (1:n-1)) { # sum(1/1 + 1/2 + 1/3 + ... + 1/n-1)
#				harmonic <- harmonic + 1/denominator
#			}

			harmonic <- sum(1/(1:(n-1)))

			dW <- c(dW, length(segSites)/harmonic)  # Number of segregating sites (Watterson's estimator), add to population array

		} else {
			print(paste('Skipping ',gene,', too few SNPs',sep=''))
			#log index of skipped genes, then make that value '-' later
			skipLog <- c(skipLog, which(uniqGenes == gene))
		}
		print('---------------------------------------------')
	}
	# Note: dT & dW are arrays!
	# Add population Tajima's D array to data.frame
	print(paste("Calculating Tajima's D for ",pop,sep=''))
	D.pop = (dT - dW) / sd(dT - dW) # Tajima's D

	# Add in missing genes as NA
	if (length(skipLog) > 0) {
		for (missing in 1:length(skipLog)) {
			D.pop <- append(D.pop, NA, after=skipLog[missing]-1)
			dT <- append(dT, NA, after=skipLog[missing]-1)
			dW <- append(dW, NA, after=skipLog[missing]-1)
		}
	}
	# If
	#D.pop[is.nan(D.pop)] <- 'H'

	D.df <- data.frame(D.df, tmp=D.pop, dT=dT, dW=dW)
	names(D.df)[length(names(D.df))-2] <- pop
	names(D.df)[length(names(D.df))-1] <- paste(pop,'.dT',sep='')
	names(D.df)[length(names(D.df))] <- paste(pop,'.dW',sep='')
	print('=============================================')
}

# Once all population Tajima's D values have been calculated, print out the results
print('Printing out results')
write.table(D.df, file=paste(paste(pops,collapse='_'),'.TD',sep=''))
