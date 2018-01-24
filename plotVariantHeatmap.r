#!/usr/bin/env Rscript

options(scipen=99999)
SVPOP <- '/path/to/SVPop/'
variantFile='/path/to/variantFile'
annotFile='/path/to/annotation'
annot <- read.table(annotFile,header=T)

# Run as './plotVariantHeatmaps.r <chromosome> <bp0> <bp1> <title (optional)>
args <- commandArgs(trailingOnly=TRUE)
chr <- args[1]
x0 <- as.numeric(args[2])
x1 <- as.numeric(args[3])
if (length(args) == 3) {
	mainTitle <- paste(chr,':',x0,'-',x1,sep='')
} else {
	mainTitle <- args[4]
}
subsetFile <- paste('Subset_',chr,':',x0,'-',x1,'_v1.4.9.csv',sep='')
annot <- annot[annot$Chromosome == chr & annot$Start < x1 & annot$End > x0,]

# Check three args given, check first is a string, check second is less than third.
if (length(args) < 3) { stop('Invalid number of arguments provided. Please supply <chromosome> <regionStart> <regionEnd> <title (optional)>') }
if (!file.exists(subsetFile)) { cat('Subset file not found, creating.\n'); system(paste(SVPOP,' --SUBSET --variantFile=',variantFile,' --region=',chr,':',x0,'-',x1,sep='')) } # stop("Given subset file doesn't exist") }

cat('Reading in Subset file.\n')
dat <- read.csv(subsetFile)
dat <- dat[dat$Precision != 0,]
models <- c('DEL','DUP','INS','INV')

cat('Plotting png.\n')
width_in = 8
tiff(paste(mainTitle,'.tiff',sep=''), height=width_in*0.3, width=width_in, units='in', res=300, compression='lzw')
par(mfrow=c(length(models)+1,1),oma=c(0,0,2,0))
for (i in 1:(length(models)+1)) {
	if (length(models) == 1) {
		par(mar=c(0,2,0.5,2))
		showx = 's'
	} else if (i == 1) {
		par(mar=c(0,2,0.5,2))
		showx = 'n'
	} else if (i == length(models)+1) {
		par(mar=c(2,2,0,2))
		showx = 's'
		plot(1,col='white',xlim=c(x0,x1),ylim=c(0.5,-0.5),xlab='Position (bp)',ylab='',bty='n',xaxs='i',axes=F)
		for (feature in annot$Feature) {
			f = annot[annot$Feature == feature,]
			p0 = f$Start[1]
			p1 = f$End[1]
			mid = p0 + ((p1-p0)/2)
			rect(p0,0,p1,1,col='grey')
			text(mid,y=0.4,labels=feature,pos=3,offset=0,cex=0.6)
		}
		rounder <- -(nchar(x1-x0)-3)
		x0r <- round(x0,rounder)
		x1r <- round(x1,rounder)
		xtiks <- round((x1r-x0r) / 5)
		axis(1, at=c(x0,seq(x0r,x1r,xtiks),x1))
	} else {
		par(mar=c(0,2,0,2))
		showx = 'n'
	}

	if (i != (length(models)+1)) {
		d <- dat[dat$Model == models[i],]
		plot(1,col='white',xlim=c(x0,x1),ylim=c(0,dim(d)[1]),yaxt='n',xlab='',xaxt=showx,ylab='',xaxs='i')
		title(ylab=models[i], line=0)
		if (dim(d)[1] > 0) {
			for (x in 1:dim(d)[1]) {
				density <- max(d$Frequency[x], 0.01)
				colors <- c(rgb(1,0,0,density), rgb(0,0,1,density), rgb(0,1,0.2,density), rgb(1,0.5,0.1,density))
				if (models[i] == 'INS') {
					abline(v=d$Start[x], col=colors[i], lwd=1)
				} else {
					rect(d$Start[x], -10, d$End[x], dim(d)[1]*1.5, col=colors[i], border=NA)
				}
			}
		}
	}
}
title(mainTitle,outer=T)
dev.off()
cat('Done.\n')
