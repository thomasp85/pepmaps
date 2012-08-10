# TODO: Add comment
# 
# Author: Thomas
###############################################################################
### getAAtable
### loads the AAtable.csv file in inst/extdata
getAAtable <- function(){
	AAtable <- read.csv(R.home(component='library/pepmaps/extdata/AAtable.csv'), header=TRUE, as.is=TRUE)
	AAtable
}
### Qval
### Calculate Q-value
Qval <- function(pepseq, bitter=FALSE){
	AAtable <- getAAtable()
	ans <- as.data.frame(matrix(0, ncol=2, nrow=length(pepseq)))
	for (i in 1:length(pepseq)){
		ans[i, 1] <- sum(AAtable$Qval[match(strsplit(toupper(pepseq[i]), '')[[1]], AAtable$Code1)])/nchar(pepseq[i])
	}
	ans[, 2] <- ifelse(ans[,1] > 1400 & nchar(pepseq) < 7, TRUE, FALSE)
	names(ans) <- c('Q.value', 'Bitter')
	ans$Bitter <- ans$Bitter==1
	if(!bitter){
		ans <- ans[, 'Q.value', drop=FALSE]
	} else {}
	ans
}
### fragExp
### Create fragment notation from sequence
fragExp <- function(pepseq){
	AAtable <- getAAtable()
	pepseq <- toupper(pepseq)
	if(!(sum(strsplit(pepseq, '')[[1]] %in% AAtable$Code1)==nchar(pepseq))){
		stop('Invalid peptide sequence.')
	} else {}
	yion <- paste('y[', (nchar(pepseq)-1):1, ']', sep='')
	bion <- paste('b[', 1:(nchar(pepseq)-1), ']', sep='')
	pepseq <- strsplit(pepseq, '')[[1]]
	frag <- pepseq[1]
	for(i in 1:(length(pepseq)-1)){
		frag <- paste(frag, '~~', 'integral(', pepseq[i+1], ',', bion[i], ',', yion[i], ')', sep='')
	}
	frag
}
### pepMass
### create synthetic fragmentation pattern
pepMass <- function(pepseq, mono=FALSE){
	AAtable <- getAAtable()
	pepseq <- toupper(pepseq)
	if(!(sum(strsplit(pepseq, '')[[1]] %in% AAtable$Code1)==nchar(pepseq))){
		stop('Invalid peptide sequence.')
	} else {}
	if(mono){
		mass <- sapply(pepseq, function(x) sum(AAtable$Mono[match(strsplit(x, '')[[1]], AAtable$Code1)], na.rm=FALSE) + 18.01056)
	} else {
		mass <- sapply(pepseq, function(x) sum(AAtable$Avg[match(strsplit(x, '')[[1]], AAtable$Code1)], na.rm=FALSE) + 18.02)
	}
	mass
}
### peppI
### Calculates the isoelectric point of a peptide
peppI <- function(pepseq){
	AAtable <- pepmaps:::getAAtable()
	pepseq <- toupper(pepseq)
	if(!(sum(strsplit(pepseq, '')[[1]] %in% AAtable$Code1)==nchar(pepseq))){
		stop('Invalid peptide sequence.')
	} else {}
	# Sum up charges
	pepseq <- strsplit(pepseq, '')[[1]]
	charge <- function(pH, pepseq){
		ans <- list()
		ans$CTerm <- -1/(1 + 10^(AAtable$pK.cTerm[AAtable$Code1 == pepseq[1]] - pH))
		ans$Asp <- -table(pepseq)['D']/(1 + 10^(AAtable$pK.side[AAtable$Code1 == 'D'] - pH))
		ans$Glu <- -table(pepseq)['E']/(1 + 10^(AAtable$pK.side[AAtable$Code1 == 'E'] - pH))
		ans$Cys <- -table(pepseq)['C']/(1 + 10^(AAtable$pK.side[AAtable$Code1 == 'C'] - pH))
		ans$Tyr <- -table(pepseq)['Y']/(1 + 10^(AAtable$pK.side[AAtable$Code1 == 'Y'] - pH))
		ans$His <- table(pepseq)['H']/(1 + 10^(pH - AAtable$pK.side[AAtable$Code1 == 'H']))
		ans$Lys <- table(pepseq)['K']/(1 + 10^(pH - AAtable$pK.side[AAtable$Code1 == 'K']))
		ans$Arg <- table(pepseq)['R']/(1 + 10^(pH - AAtable$pK.side[AAtable$Code1 == 'R']))
		ans$NTerm <- 1/(1 + 10^(pH - AAtable$pK.nTerm[AAtable$Code1 == pepseq[length(pepseq)]]))
		ans <- do.call('c', ans)
		ans <- sum(ans, na.rm=TRUE)
		ans
	}
	# Search for pI
	pH <- 6.5
	pHprev <- 0
	pHnext <- 14
	while((pH-pHprev) > 0.01 & (pHnext-pH) > 0.01){
		Q <- charge(pH, pepseq)
		if(Q < 0){
			pHnext <- pH
			pH <- pH - (pH-pHprev)/2
		} else {
			pHprev <- pH
			pH <- pH + (pHnext-pH)/2
		}
	}
	pH
}
### fragPattern
### Calculate the theoretical fragmentation pattern by CID for a peptide
fragPattern <- function(pepseq){
	AAtable <- pepmaps:::getAAtable()
	pepseq <- toupper(pepseq)
	if(!(sum(strsplit(pepseq, '')[[1]] %in% AAtable$Code1)==nchar(pepseq))){
		stop('Invalid peptide sequence.')
	} else {}
	pepseq <- strsplit(pepseq,'')[[1]]
	bion <- data.frame(ion=paste('b', 1:(length(pepseq)-1), sep=''), mz=NA, stringsAsFactors=FALSE)
	yion <- data.frame(ion=paste('y', 1:(length(pepseq)-1), sep=''), mz=NA, stringsAsFactors=FALSE)
	for(i in 1:(length(pepseq))-1){
		bionseq <- paste(pepseq[1:i], collapse='')
		bionmz <- pepMass(bionseq)+1.0078250-18.01057
		bion$mz[i] <- bionmz
		yionseq <- paste(pepseq[(length(pepseq)+1-i):length(pepseq)], collapse='')
		yionmz <- pepMass(yionseq)+1.0078250
		yion$mz[i] <- yionmz
	}
	ans <- rbind(bion, yion)
	ans
}