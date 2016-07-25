# Calculate Fst
# Input data: Matrix, SNP per row, Sample per column (with a population prefix)
# Allele calls: Major as 0, minor as 1, mixed allowed
# Populations must have a prefix which matches given population IDs
# 14 Mar 2016

# Resolve NA and -Inf issues!

# Libraries
suppressMessages(library(gdata))
suppressMessages(library(data.table))

# Controls
if (FALSE) {
  outPrefix <- 'DEFAULTOUT'
  #dat <- read.table('INFILE', header=T) # Whole dataset
  dat <- fread('INFILE', header=T) # Whole dataset
  populations <- c('ABC','EFG','IJK') # List of population prefix IDs
  allGroups <- list(c('ABC','EFG'), c('IJK'))
  nameGroups <- c('groupA','groupB')
}

# Get A and n from input df

calcStats <- function(inData, pops) {
  # Get A and n arrays from input data
  out.A <- c()
  out.n <- c()

  for (i in pops) {
    print(paste('Calculating A and n for ',i,sep=''))
    popX <- inData[,startsWith(names(inData), i), with=F]

    out.n <- cbind(out.n, dim(popX)[2])
    out.A <- cbind(out.A, rowSums(popX)/dim(popX)[2])

  }
  list(out.A, out.n)
}

calcFstWright <- function(A, n) {
  # Calculates Fst for one SNP
  # A = array of allele A frequencies for each population
  # n = array of population sizes in same order as A

  if (length(A) != length(n)) {
    stop('Frequency array (A) and population size array (n) are different lengths.')
  }

  At <- sum(A)/length(A)
  nT <- length(n)
  nM <- mean(n)
  sum((n*((A-At)*(A-At)))/((nT-1)*nM))/(At*(1-At))
}

calcFstNei <- function(A) {
  # Calculates Fst estimation for one SNP
  # A = array of allele A frequencies for each population

  Hs = sum(2*A*(1-A)) / length(A)
  Ht = 2 * mean(A) * mean(1-A)

  1 - (Hs / Ht)
}

### Wrapper
calcFst <- function(inData, pops, groups=NULL, pairwise=T, Nei=T) {
  # inData = input SNP table, as data.frame
  # pops = array of sub-population ID prefixes
  # groups = list of groups that those sub-pops belong to

  if (pairwise) {
    # iterate pairs of populations (ie. Malawi v Tanzania)

    # Calculate pairs of populations
    pairs <- combn(pops,2)

    # For each pair
    for (i in 1:dim(pairs)[2]) {

      stats <- calcStats(pairs[inData,i])

      A <- data.frame(stats[1])
      n <- unlist(stats[2])

      if (Nei) { Fst <- apply(A, 1, calcFstNei) }
      else { Fst <- apply(A, 1, calcFstWright, n=n) } # e.g. if two, pairwise

      output <- data.frame(PopA=pairs[1,],PopB=pairs[2,],Fst,stringsAsFactors=F)

    }

    values <- combn(inData[1,],2,sum) # replace with array fst call
    output <- data.frame(PopA=pairs[1,],PopB=pairs[2,],values,stringsAsFactors=F)
    # call Fst from pairs columns
  }
}


calcFst(dat, c('wAFR','eAFR','cAFR','sASI','seASI'))
