### Class that merges PepID and Peaklist (and possibly additional *ID classes).
setClass(
    Class='Complist',
    representation=representation(
        raw='matrix',
		chromatograms='list',
        Sample.info='data.frame',
        Peak.info='data.frame',
        annotation='list',
        parameters='Parameters',
        IDindex='list',
        pepID='PepID',
        xcmsSet='xcmsSet',
        models='list',
        filter='list'
    ),
    prototype=prototype(
        raw=matrix(ncol=0, nrow=0),
		chromatograms=list(TIC=NULL, BPC=NULL),
        Sample.info=data.frame(),
        Peak.info=data.frame(),
        annotation=list(),
        parameters=parameters(),
        IDindex=list(),
        pepID=pepID(),
        xcmsSet=new(Class='xcmsSet'),
        models=list(),
        filter=list(Peaks=c(), Peptides=c(), Samples=c())
    ),
    validity=function(object){
        if(length(object@raw) != 0){
            if(min(object@raw) < 0){
                stop('MS data cannot be negative')
            } else if(nrow(object@Sample.info) != ncol(object@raw)){
                stop('Submitted sample information does not correspond to the number of samples in xcms data')
            } else {TRUE}
        } else {TRUE}
    }
)
### show
### Short summary of Complist object
setMethod(
    'show', 'Complist',
    function(object){
		
		# Test if object is empty
        if(length(object) == 0){
            cat('An empty Complist object\n')
		
		# Test if identification has been performed
		} else {
            cat('A Complist object with', nrow(object@Sample.info), 'samples and', nrow(object@Peak.info),'peak groups\n\n')
            cat('\t', sum(object@filter$Peptides$Matched), ' peptides matched to ', sum(object@filter$Peaks$Matched), ' peaks\n\n', sep='')
			cat('\t', sum(apply(object@filter$Peptides[,c('Matched', 'FDR')], 1, all)), ' matched peptides with FDR cutoff of ', object@pepID@FDR, '\n', sep='')
        }
    }
)
### plot
### Creates a correlation matrix plot for samples in the Complist object
setMethod(
    'plot', 'Complist',
    function(x, method='pearson', unmatched=FALSE){
        if(!unmatched){
            ggCorr(x@raw[x@filter$Peaks$Matched,], method=method)
        } else {
            ggCorr(x@raw, method=method)
        }
    }
)
### setFolder
### Sets new location of data files
setMethod(
    'setFolder', 'Complist',
    function(object, folder=choose.dir()){
		
		# Data input
		if(missing(folder)){
			if(Sys.info()["sysname"] == 'Windows'){
				folder <- choose.dir()
			} else {
				folder <- readline('Path to folder:')
			}
		}
		
		# Create and assign new filepaths
        newFP <- paste(folder, basename(object@Sample.info$Filepath), sep='/')
        object@Sample.info$Filepath <- newFP
        object@xcmsSet@filepaths <- newFP
        object
    }
)
### length
### Outputs 0 if complist object is empty and 1 if not
setMethod(
    'length', 'Complist',
    function(x){
        if(sum(length(x@raw)+length(x@Sample.info)+length(x@Peak.info)+length(x@annotation)+length(x@parameters)+length(x@IDindex)+length(x@pepID)) == 0){
            0
        } else {
            1
        }
    }
)
### getRaw
### Extract consensus matrix for the samples in either height or area
setMethod(
    'getRaw', 'Complist',
    function(object, outlier.rm=TRUE, mix.rm=TRUE, filter=NULL, unmatched=FALSE, FDR=TRUE, height=FALSE, normalize=FALSE){
		
		# Extract data in the right format
        if(height){
			ans <- groupval(object@xcmsSet, value='maxo')
		} else {
			ans <- object@raw
		}
		
		# Subsetting of data to match input
        row.names(ans) <- 1:nrow(ans)
        remR <- getRemove(object, 'Peak', outlier=outlier.rm, filter=filter, match=!unmatched, FDR=FDR)
        remC <- getRemove(object, 'Sample', outlier=outlier.rm, mix=mix.rm)
        if(length(remR) != 0){
            ans <- ans[-remR,]
        } else {}
        if(length(remC) != 0){
            ans <- ans[,-remC]
        } else {}
		if(normalize){
			ans <- t(scale(t(ans), center=F, scale=apply(ans, 1, max)))*100
		} else {}
        ans
    }
)
### getMatch
### For a given peak id (index in consensus matrix), output the lines in the peplist with the matched peptides
setMethod(
    'getMatch', 'Complist',
    function(object, ID, bestmatch=FALSE){
        ans <- list()
		
		# Extract the data from the peplist
        for(i in 1:length(ID)){
            index <- which(sapply(object@IDindex, function(x, id) id %in% x, id=as.numeric(ID[i])))
            ans[[i]] <- getPeplist(object, outlier.rm=FALSE, unmatched=TRUE)[index,]
        }
		
		# Create feedback data (nomatch feedback)
        names(ans) <- ID
        nomatch <- names(ans)[which(sapply(ans, nrow) == 0)]
        ans <- ans[which(sapply(ans, nrow) != 0)]
        if(length(nomatch) != 0){
            cat('No match found for: ', paste(nomatch, collapse=', '), '\n\n', sep='')
        } else {}
		
		# Extract best match
        if(bestmatch){
            ans <- ldply(ans, selectMatch)
            names(ans)[1] <- 'Peak.ID'
        }
        ans
    }
)
### getPeplist
### Extracts the list of peptides stored in the PepID slot
setMethod(
    'getPeplist', 'Complist',
    function(object, outlier.rm=TRUE, unmatched=FALSE, FDR=TRUE, filter=NULL){
        ans <- object@pepID@peplist
        rem <- getRemove(object, 'Peptide', outlier=outlier.rm, match=!unmatched, FDR=FDR, filter=filter)
        if(length(rem) != 0){
            ans <- ans[-rem,]
        } else {}
        ans
    }
)
### getPeakinfo
### Extract information on the peaks in the consensus matrix matched to peptides
setMethod(
    'getPeakinfo', 'Complist',
    function(object, outlier.rm=TRUE, filter=NULL, unmatched=FALSE, FDR=TRUE){
        ans <- object@Peak.info
        rem <- getRemove(object, 'Peak', outlier=outlier.rm, filter=filter, match=!unmatched, FDR=FDR)
        if(length(rem) != 0){
            ans <- ans[-rem,]
        } else {}
        ans
    }
)
### getRawID
### Extracts the raw data from the ID analysis, stored in the PepID slot
setMethod(
    'getRawID', 'Complist',
    function(object){
        ans <- object@pepID@raw
        ans
    }
)
### getPepPeakIndex
### Creates a list with indexes for matches between peaks and compounds. The format is decided with the from parameter
setMethod(
	'getPepPeakIndex', 'Complist',
	function(object, from, outlier.rm=TRUE, FDR=TRUE){
		if(from == 'Peak'){
			
			# Get vector with peak ID's
			rem <- getRemove(object, 'Peptide', outlier=outlier.rm, match=TRUE, FDR=FDR)
			if(length(rem) != 0){
				ind <- object@IDindex[-rem]
			} else {
				ind <- object@IDindex
			}
			ind <- do.call('c', ind)
			ind <- ind[!ind %in% which(object@filter$Peaks$Outlier)]
			ind <- unique(ind)
			
			# Create a list with an element for each peak ID holding the index of the matched peptide
			k <- rep(seq_along(object@IDindex), sapply(object@IDindex, length))
			ans <- sapply(ind, function(x) k[which(unlist(object@IDindex) %in% x)])
			names(ans) <- ind
			ans <- ans[order(ind)]
			ans
		} else if(from == 'Peptide'){
			
			# Trimming of the already available list in the object
			rem <- getRemove(object, 'Peptide', outlier=outlier.rm, FDR=FDR, match=TRUE)
			if(length(rem) != 0){
				ans <- object@IDindex[-rem]
			} else {
				ans <- object@IDindex
			}
			ans
		} else {
			stop('Wrong input: from must be either \'Peak\' or \'Peptide\'')
		}
	}
)
### getBestmatch
### Calculates the best match for all matched peaks in the consensus matrix
setMethod(
    'getBestmatch', 'Complist',
    function(object, outlier.rm=TRUE, filter=NULL, FDR=TRUE){
		ID <- getPepPeakIndex(object, from='Peak', outlier.rm=outlier.rm, FDR=FDR)
		ID <- lapply(ID, function(x) selectMatch(object@pepID@peplist[x, ], IDtype=object@pepID@type))
		rem <- getRemove(object, 'Peak', filter=filter)
		if(length(rem) != 0){
			ID <- ID[!names(ID) %in% as.character(rem)]
		}
        ID <- do.call('rbind', ID)
        ID
    }
)
### getFeature
### Extracts and optionally plots EIC's for peaks at given set rt and mz values 
setMethod(
    'getFeature', 'Complist',
    function(object, rt, mz, plot=TRUE){
		
		# Input check
        if(length(rt) != length(mz)){
            stop('Length of rt and mz does not match')
        } else {}
        index <- list()
        name <- list()
		
		# Iterate over the rt/mz pairs
        for(i in 1:length(rt)){
			
			# Get peaks matching the rt/mz pair
            ans <- which(object@Peak.info$rtmax > rt[i] & object@Peak.info$rtmin < rt[i] & object@Peak.info$mzmax > mz[i] & object@Peak.info$mzmin < mz[i])
            if(length(ans) != 0){
                info <- list()
                for(j in 1:length(ans)){
                    tmp <- list()
                    tmp$Index <- ans[j]
                    tmp$Info <- object@Peak.info[ans[j], 1:6]
					
					# Get index of features from the same compound
                    ann <- object@annotation$Compounds[getGroups(ans[j], object@annotation)]
                    if(length(ann) != 0){
                        tmp$Annotation <- ann
                    } else {}
					
					# Get matched peptides
                    pep <- which(sapply(object@IDindex, function(x) if(!is.null(x)) x %in% ans[j] else FALSE))
                    if(length(pep) != 0){
                        tmp$Peptide <- object@pepID@peplist[pep,]
                    } else {}
					
					# Get outlier status and intensity
                    tmp$Outlier <- object@filter$Peaks$Outlier[ans[j]]
                    tmp$Signal <- object@raw[ans[j],]
					
                    info[[j]] <- tmp
                }
				if(length(info) == 1){
					info <- info[[1]]
				} else {}
                name[[i]] <- paste(mz[i], '/', rt[i], sep='')
                index[[i]] <- info
            } else {}
        }
        name <- do.call('c', name)
        if(length(index) != 0){
			rem <- which(sapply(index, is.null))
			if(length(rem) != 0){
				index <- index[-rem]
			} else {}
			names(index) <- name
			
			# Plotting
            if(plot){
                ans <- list()
                for(i in 1:length(index)){
                    if(is.null(names(index[[i]]))){
                        df <- list()
                        for(j in 1:length(index[i])){
                            df[[j]] <- data.frame(Intensity=index[[i]][[j]]$Signal, Index=factor(index[[i]][[j]]$Index))
                        }
                        df <- do.call('rbind', df)
                        df$Feature <- names(index)[i]
                    } else {
                        df <- data.frame(Intensity=index[[i]]$Signal, Index=factor(index[[i]]$Index), Feature=names(index)[i])
                    }
                    ans[[i]] <- df
                }
                ans <- do.call('rbind', ans)
                p <- ggplot(data=ans, aes(x=Index, y=Intensity)) + theme_bw()
                p <- p + geom_boxplot() + xlab('') + facet_wrap(~Feature, scales='free_x')
                print(p)
            } else {}
            index
        } else {
            cat('No match(es) found\n')
        }
    }
)
### reAnnotate
### Run CAMERA analysis on the Complist, replacing any earlier deconvolution
setMethod(
    'reAnnotate', 'Complist',
    function(object, ...){
        an <- CAMERA::annotate(object@xcmsSet, ...)
        annotation <- conFeature(an)
        object@annotation <- annotation
		
		# Update Complist parameters
        object@parameters <- editPar(object@parameters, type='annotate', ...)
        object
    }
)
### reFindPeaks
### Repeat peak finding with new parameters, and recalculate everythin based on the old parameters
setMethod(
    'reFindPeaks', 'Complist',
    function(object, annotate=FALSE, ...){
        samples <- sampleInfo(object, outlier.rm=FALSE, mix.rm=FALSE)$Filepath
        xset <- xcmsSet(samples, method='centWave', ...)
        para <- object@parameters@parameters$retcor
        para$object <- xset
        xset <- do.call('retcor', para)
        para <- object@parameters@parameters$group
        para$object <- xset
        xset <- do.call('group', para)
        xset <- fillPeaks(xset)
        if(annotate){
            para <- object@parameters@parameters$annotate
            para$object <- xset
            xset <- do.call('annotate', para)
        } else {}
        newSinfo <- sampleInfo(object, outlier.rm=FALSE, mix.rm=FALSE)
        newSinfo <- newSinfo[,-which(names(newSinfo) %in% c('Class', 'Filepath')), drop=FALSE]
        rtwin <- object@parameters@parameters$Complist$rtwin
        mzwin <- object@parameters@parameters$Complist$mzwin
		chrom <- object@chromatograms
        object <- complist(xset, newSinfo, object@parameters, object@pepID, rtwin=rtwin, mzwin=mzwin)
        cat('Models, filters and outlier information have been cleared...\n')
		object@chromatograms <- chrom
        object@parameters <- editPar(object@parameters, type='findPeak', ...)
        object
    }
)
### reGroup
### Recalculates grouping of features based on new parameters
setMethod(
    'reGroup', 'Complist',
    function(object, annotate=FALSE, ...){
        xset <- group(object@xcmsSet, ...)
        xset <- fillPeaks(xset)
        if(annotate){
            para <- object@parameters@parameters$annotate
            para$object <- xset
            xset <- do.call('annotate', para)
        } else {}
        newSinfo <- sampleInfo(object, outlier.rm=FALSE, mix.rm=FALSE)
        newSinfo <- newSinfo[,-which(names(newSinfo) %in% c('Class', 'Filepath')), drop=FALSE]
        rtwin <- object@parameters@parameters$Complist$rtwin
        mzwin <- object@parameters@parameters$Complist$mzwin
		chrom <- object@chromatograms
        object <- complist(xset, newSinfo, object@parameters, object@pepID, rtwin=rtwin, mzwin=mzwin)
        cat('Models, filters and outlier information have been cleared...\n')
		object@chromatograms <- chrom
        object@parameters <- editPar(object@parameters, type='group', ...)
        object
    }
)
### reRTcorrect
### Recalculates retention time correction based on new parameters. Can optionally remove correction altogether
setMethod(
    'reRTcorrect', 'Complist',
    function(object, annotate=FALSE, reset=FALSE, ...){
        if(reset){
            object@xcmsSet@rt$corrected <- object@xcmsSet@rt$raw
			xset <- object@xcmsSet
        } else {
			object@xcmsSet@rt$corrected <- object@xcmsSet@rt$raw
            xset <- retcor(object@xcmsSet, method=object@parameters@parameters$retcor$method, ...)
        }
        para <- object@parameters@parameters$group
        para$object <- xset
        xset <- do.call('group', para)
        xset <- fillPeaks(xset)
        if(annotate){
            para <- object@parameters@parameters$annotate
            para$object <- xset
            xset <- do.call('annotate', para)
        } else {}
        newSinfo <- sampleInfo(object, outlier.rm=FALSE, mix.rm=FALSE)
        newSinfo <- newSinfo[,-which(names(newSinfo) %in% c('Class', 'Filepath')), drop=FALSE]
		pepID <- rescalePeplist(basename(xset@filepaths), xset@rt, object@pepID)
        rtwin <- object@parameters@parameters$Complist$rtwin
        mzwin <- object@parameters@parameters$Complist$mzwin
		chrom <- object@chromatograms
		object <- complist(xset, newSinfo, object@parameters, pepID, rtwin=rtwin, mzwin=mzwin)
        cat('Models, filters and outlier information have been cleared...\n')
		object@chromatograms <- chrom
		object@parameters <- editPar(object@parameters, type='retcor', ...)
        object
    }
)
### reMatch
### Performs a new match between peaks and peptides based on a new window
setMethod(
    'reMatch', 'Complist',
    function(object, rtwin, mzwin){
        peakID <- rep(FALSE, nrow(object@filter$Peaks))
        pep.index <- vector('list', nrow(object@filter$Peptides))
        for(i in 1:length(pep.index)){
            ind <- which(object@Peak.info$mzmin < (object@pepID@peplist$mz[i]+mzwin) & object@Peak.info$mzmax > (object@pepID@peplist$mz[i]-mzwin) & object@Peak.info$rtmin < (object@pepID@peplist$Retention.time[i]+rtwin) & object@Peak.info$rtmax > (object@pepID@peplist$Retention.time[i]-rtwin))
            if(length(ind) > 0){
                peakID[ind] <- TRUE
                pep.index[[i]] <- ind
            } else {}
        }
        object@filter$Peaks$Matched <- peakID
		object@filter$Peptides$Matched <- !sapply(pep.index, is.null)
        object@IDindex <- pep.index
        object@filter$Peptides$Outlier[] <- FALSE
        outlierID <- which(object@filter$Peaks$Outlier)
        index <- which(sapply(object@IDindex, function(x) sum(x %in% outlierID)) > 0)
        if(length(index) > 0){
            index <- unlist(index)
            object@filter$Peptides$Outlier[index] <- TRUE
        } else {}
        cat('Peptide outliers set according to peak outliers...\n')
        object@parameters <- editPar(object@parameters, type='Complist', rtwin=rtwin, mzwin=mzwin)
        object
    }
)
### reSetFDR
### Sets a new FDR cutoff value and filters the identifications
setMethod(
	'reSetFDR', 'Complist',
	function(object, FDR){
		object@pepID@FDR <- FDR
		object@filter$Peptides$FDR <- object@pepID@peplist$FDR < FDR
		object
	}
)
### findOutlier
### Performs an outlier detection based on the type (Sample, Feature, Peptide)
setMethod(
    'findOutlier', 'Complist',
    function(object, what, method, replic, qreg='linear', plot=TRUE){
        if(what == 'Sample'){
			
			# Default (and only) method for Sample outlier detection
            if(missing(method)){
                method <- 'RMD'
            } else {}
            if(method == 'RMD'){
				
				# Prepare data
                dat <- log(object@raw)
                dat[which(is.infinite(dat))] <- NA
                cr <- cor(dat, use='na.or.complete')
                for(i in 1:ncol(cr)){
                    cr[i,i] <- NA
                }
				
				# Calculate the 5 statistics
                if(missing(replic)){
                    M1 <- apply(cr, 1, function(x) mean(x, na.rm=TRUE))
                } else {
                    M1 <- rep(0, ncol(dat))
                    for(i in 1:length(M1)){
                        cl <- which(object@Sample.info[,replic] == object@Sample.info[i,replic])
                        M1[i] <- mean(cr[i, cl], na.rm=TRUE)
                    }
                }
                M2 <- rep(0, ncol(dat))
                for(i in 1:length(M2)){
                    M2[i] <- sum(is.na(dat[,i]))/nrow(dat)
                }
                M3 <- apply(dat, 2, function(x) mad(x, na.rm=TRUE))
                skew <- function(x){
                    sum((x - mean(x, na.rm=TRUE))^3, na.rm=TRUE)/(length(x) * sd(x, na.rm=TRUE)^3)
                }
                M4 <- apply(dat, 2, skew)
                kurtosis <- function(x){
                    sum((x - mean(x, na.rm=TRUE))^4, na.rm=TRUE)/(length(x) * sd(x, na.rm=TRUE)^4)-3
                }
                M5 <- apply(dat, 2, kurtosis)
                M <- data.frame(M1, M2, M3, M4, M5)
                row.names(M) <- row.names(object@Sample.info)
				
				# Perform a robust PCA by projection persuit
                cv <- covPCAproj(M, control=list(k=ncol(M), scale=function(x) mad(x)*1.483))
				
				# Calculate distances
                mdist <- mahalanobis(M, cv$center, cv$cov)
                mdist <- data.frame(Sample=names(mdist), Dist=mdist)
                mdist$Dist <- log2(mdist$Dist)
                pval <- data.frame(Confidence=factor(c('0.95', '0.99', '0.9999')), value=log2(qchisq(c(0.95,0.99,0.9999), 5)))
                lab <- mdist[which(mdist$Dist > pval$value[1]),]
				
				# Plotting
                p <- ggplot(data=mdist) + theme_bw() + xlab('') + ylab('ln(Mahalanobis distance)')
                p <- p + geom_point(aes(x=Sample, y=Dist))
                p <- p + geom_hline(data=pval, aes(yintercept=value, linetype=Confidence), colour=I('red'), size=1)
                if(nrow(lab) != 0){
                    p <- p + geom_text(data=lab, aes(x=Sample, y=Dist, label=Sample), hjust=-0.1, vjust=0, size=2.5)
                } else {}
                p <- p + theme(axis.text.x=element_text(angle=-45, hjust=0, vjust=1, size=8))
                print(p)
				invisible(lab$Sample)
            } else {}
        } else if(what == 'Peak'){
			
			# Default (and only) method for Peaks
            if(missing(method)){
                method <- 'OutlierD'
            } else {}
            if(method == 'OutlierD'){
				
				# Check input
                if(length(replic) != 2){
                    stop('Supply 2 replicates\n')
                } else {}
				
				# Create data
                if(is.numeric(replic)){
                    dat <- object@raw[,replic]
                } else {
                    if(sum(sampleNames(object, all=TRUE) %in% replic) != 2){
                        stop('Incorrect sample names provided\n')
                    } else {}
                    dat <- object@raw[,which(sampleNames(object, all=TRUE) %in% replic)]
                }
				
				# Run OutlierD
                fit <- OutlierD(dat[,1], dat[,2], k=1.5, method=qreg)
                ans <- which(fit$x$Outlier)
				
				# Create plot
                if(plot){
                    p <- qplot(fit$x$A, fit$x$M, geom='point', colour=fit$x$Outlier) + scale_colour_hue('Outlier')
                    p <- p + geom_line(aes(x=fit$x$A, y=fit$x$UB, colour=NULL))
                    p <- p + geom_line(aes(x=fit$x$A, y=fit$x$LB, colour=NULL))
                    p <- p + geom_line(aes(x=fit$x$A, y=fit$x$Q1, colour=NULL), linetype=I(3))
                    p <- p + geom_line(aes(x=fit$x$A, y=fit$x$Q3, colour=NULL), linetype=I(3))
                    p <- p + theme_bw() + xlab('A') + ylab('M')
                    print(p)
                }
                invisible(ans)
            } else {}
        } else {}
    }
)
### setOutlier
### Flag a given sample, peak or peptide as outlier
setMethod(
    'setOutlier', 'Complist',
    function(object, what, ID, value=TRUE, recursive=TRUE){
        
		# Input check
		if(!is.logical(value)){
            stop('value must be either TRUE or FALSE')
        } else {}
		
		# Samples
        if(what == 'Sample'){
			
			# Input check
            if(sum(!(ID %in% row.names(object@Sample.info))) != 0){
                ID <- ID[which(!(ID %in% row.names(object@Sample.info)))]
                stop('No samples called:', paste(ID, collapes=', '))
            } else {}
            index <- which(row.names(object@Sample.info) %in% ID)
            object@filter$Samples$Outlier[index] <- value
			
		# Peaks
        } else if(what == 'Peak'){
			
			# Input check
            if(max(ID) > nrow(object@Peak.info)){
                stop('ID must be between 1 and', nrow(object@Peak.info))
            } else {}
            object@filter$Peaks$Outlier[ID] <- value
			
			# Flag matched peptides as outliers too
            if(recursive){
                index <- which(sapply(object@IDindex, function(x) sum(x %in% ID)) > 0)
                if(length(index) > 0){
                    index <- unlist(index)
                    object@filter$Peptides$Outlier[index] <- value
                } else {}
            } else {}
			
		# Peptides
        } else if(what == 'Peptide'){
            if(is.character(ID)){
				
				# Test if input is a sequence or a name
                seqorID <- grepl('-', ID)
                index <- list()
                for(i in 1:length(seqorID)){
                    if(seqorID[i]){
                        index[[i]] <- which(object@pepID@peplist$Peptide.ID %in% ID[i])
                    } else {
                        ID1 <- object@pepID@peplist$Peptide.ID[which(object@pepID@peplist$Sequence %in% ID[i])]
                        index[[i]] <- which(object@pepID@peplist$Peptide.ID %in% ID1)
                    }
                }
                index <- do.call('c', index)
            } else {
                index <- ID
            }
            object@filter$Peptides$Outlier[index] <- value
			
			# Flag matched peaks as outliers too
            if(recursive){
                index <- unlist(object@IDindex[index])
				object@filter$Peaks$Outlier[index] <- value
            } else {}
        } else {}
        object
    }
)
### outlier
### Display the flagged outliers
setMethod(
    'outlier', 'Complist',
    function(object){
        if(sum(c(object@filter$Samples$Outlier, object@filter$Peaks$Outlier, object@filter$Peptides$Outlier)) == 0){
            cat('No outliers flagged in Complist object...\n')
        } else {
            if(sum(object@filter$Samples$Outlier) > 0){
                cat('Samples:\n\n')
                cat(row.names(object@Sample.info[object@filter$Samples$Outlier,]))
                cat('\n\n')
            } else {}
            if(sum(object@filter$Peaks$Outlier) > 0){
                cat('Peak groups:\n\n')
                print(object@Peak.info[object@filter$Peaks$Outlier,])
                cat('\n')
            } else {}
            if(sum(object@filter$Peptides$Outlier) > 0){
            cat('Peptides:\n\n')
            print(object@pepID@peplist[object@filter$Peptides$Outlier,])
            }
        }
    }
)
### editSampleinfo
### Add, delete or change the design matrix of the Complist
setMethod(
    'editSampleinfo', 'Complist',
    function(object, what='add', data, name, sort){
        if(what != 'delete'){
			
			# Input check
            if(is.data.frame(data)){
                data <- data
            } else if(is.vector(data) | is.factor(data) | is.character(data)){
                if(length(name) == 1){
                    data <- data.frame(data, stringsAsFactors=FALSE)
                    names(data) <- name
                } else {
                    stop('Dimension mismatch between data and name\n')
                }
            } else if(is.matrix(data)){
                if(length(colnames(data)) == 0){
                    if(ncol(data) == length(name)){
                        data <- data.frame(data, stringsAsFactors=FALSE)
                        names(data) <- name
                    } else {
                        stop('Dimension mismatch between data and name\n')
                    }
                } else {
                    data <- data.frame(data, stringsAsFactors=FALSE)
                }
            } else {
                stop('Wrong data input\n')
            }
            if(nrow(data) != nrow(object@Sample.info) && missing(sort)){
                stop('Wrong sample number supplied\n')
            } else {}
        } else {}
        if(what == 'add'){
			
			# Input check
            if(sum(names(object@Sample.info) %in% names(data)) > 0){
                err <- names(data)[which(names(data) %in% names(object@Sample.info))]
                stop('Sample info with name <', paste(err, collapse=', '), '> already exists\n', sep='')
            } else {}
			
			# Sort the data to match the sample info
			if(!missing(sort)){
				if(!sort %in% names(data)){
					stop('Cannot sort by', sort, '! Not part of supplied data')
				} else {}
				sortby <- data[, sort]
				data[, sort] <- NULL
				sinfo <- data.frame(object@Sample.info, data[match(rownames(object@Sample.info), sortby), ])
			} else {
				sinfo <- data.frame(object@Sample.info, data)
			}
        } else if(what == 'change'){
			
			# Input check
            if(sum(names(object@Sample.info) %in% names(data)) != ncol(data)){
                err <- names(data)[which(!(names(data) %in% names(object@Sample.info)))]
                stop('Sample info with name <', paste(err, collapse=', '), '> doesn\'t exist\n', sep='')
            } else {}
            sinfo <- object@Sample.info
            sinfo[,names(data)] <- data
        } else if(what == 'delete'){
			
			# Input check
            if(sum(names(object@Sample.info) %in% name) != length(name)){
                err <- name[which(!(name %in% names(object@Sample.info)))]
                stop('Sample info with name <', paste(err, collapse=', '), '> doesn\'t exist\n', sep='')
            } else {}
            sinfo <- object@Sample.info[,-which(names(object@Sample.info) %in% name)]
        } else{
            stop('Unknown operator <', what, '>\n', sep='')
        }
        object@Sample.info <- sinfo
        object
    }
)
### sampleInfo
### Returns the design matrix for the complist object
setMethod(
    'sampleInfo', 'Complist',
    function(object, outlier.rm=TRUE, mix.rm=TRUE, all=FALSE){
		if(all){
			outlier.rm=FALSE
			mix.rm=FALSE
		}
		rem <- getRemove(object, what='Sample', outlier=outlier.rm, mix=mix.rm)
		if(length(rem) != 0){
			sinfo <- object@Sample.info[-rem, ]
		} else {
			sinfo <- object@Sample.info
		}
        sinfo
    }
)
### sampleNames
### Returns the sample names of the Complist object
setMethod(
	'sampleNames', 'Complist',
	function(object, outlier.rm=TRUE, mix.rm=TRUE, all=FALSE){
		rownames(sampleInfo(object, outlier.rm=outlier.rm, mix.rm=mix.rm, all=all))
	}
)
### plotEIC
### Creates a plot of the EIC's associated with given peptides and optionally saves the plot
setMethod(
    'plotEIC', 'Complist',
    function(object, ID, type='Peptide', retcor=TRUE, Sample, group, outlier.rm=TRUE, mix.rm=TRUE, title=''){
        if(retcor){
            rt <- 'corrected'
        } else {
            rt <- 'raw'
        }
		
		# Input check
        if(!type %in% c('Peptide', 'Peak')){
            stop('Unknown compound type <', type, '> ...', sep='')
        } else {}
        if(!missing(group)){
            if(!group %in% names(object@Sample.info)){
                stop('Unknown Sample information <', group, '>', sep='')
            } else {}
        } else {}
		
        if(type == 'Peptide'){
			if(any(ID %in% rownames(object@pepID@peplist)[object@filter$Peptides$Outlier])){
				IDout <- ID[ID %in% rownames(object@pepID@peplist)[object@filter$Peptides$Outlier]]
				cat(paste('Warning: ', paste(IDout, collapse=', '), ' marked as outliers\n', sep=''))
			} else {}
			if(any(!ID %in% object@pepID@peplist$Peptide.ID)){
				IDout <- ID[!ID %in% object@pepID@peplist$Peptide.ID]
				cat('Warning: Invalid peptide id for: ', paste(IDout, collapse=', '), '\n', sep='')
			} else {}
			
			# Get rt mz window for all EIC's
            windows <- data.frame(matrix(NA, ncol=5, nrow=length(ID)))
            for(i in 1:length(ID)){
                index <- which(object@pepID@peplist$Peptide.ID == ID[i])
                index <- object@IDindex[[index]]
                if(is.null(index)){
                    cat('Warning: Peptide with ID: <', ID[i], '> has not been linked to a peak group.\n\n', sep='')
                    flush.console()
                    peps <- getRawID(object)[which(getRawID(object)$Peptide.ID == ID[i]),]
                    windows[i,1] <- ID[i]
                    windows[i,2:3] <- c(min(peps$rt)-5, max(peps$rt)+5)
                    windows[i,4:5] <- c(min(peps$Precursor)-0.5, max(peps$Precursor)+0.5)
                } else {
                    if(object@filter$Peaks$Outlier[index]){
                        cat('Warning: Peak group matching Peptide.ID: <', ID[i], '> has been marked as outlier.\n\n', sep='')
                    } else {}
                    windows[i,1] <- ID[i]
                    windows[i,2:3] <- c(object@Peak.info$rtmin[index]-5, object@Peak.info$rtmax[index]+5)
                    windows[i, 4:5] <- c(object@Peak.info$mzmin[index], object@Peak.info$mzmax[index])
                }
            }
            colnames(windows) <- c('ID', 'rtmin', 'rtmax', 'mzmin', 'mzmax')
			
        } else if(type == 'Peak'){
			if(any(ID < 1) | any(ID > nrow(object@raw))){
				IDout <- ID[ID < 1 | ID > nrow(object@raw)]
				cat(paste('Warning: Invalid peak id for', paste(IDout, collapse=', '), '\n', sep=''))
			} else {}
			if(any(ID %in% which(object@filter$Peaks$Outlier))){
				IDout <- ID[ID %in% which(object@filter$Peaks$Outlier)]
				cat('Warning: ', paste(IDout, collapse=', '), ' marked as outliers\n', sep='')
			} else {}
			
			windows <- data.frame(matrix(NA, ncol=5, nrow=length(ID)))
			for(i in 1:length(ID)){
				windows[i,1] <- ID[i]
				windows[i,2:3] <- c(object@Peak.info$rtmin[ID[i]]-5, object@Peak.info$rtmax[ID[i]]+5)
				windows[i, 4:5] <- c(object@Peak.info$mzmin[ID[i]], object@Peak.info$mzmax[ID[i]])
			}
			colnames(windows) <- c('ID', 'rtmin', 'rtmax', 'mzmin', 'mzmax')
		} else {}
		
		# Extract EIC's
		cat('Extracting EIC for all samples...\n\n')
		flush.console()
		if(missing(Sample)){
			Sample <- sampleNames(object, outlier.rm=outlier.rm, mix.rm=mix.rm)
		}
		eic <- getEIC(object@xcmsSet, mzrange=as.matrix(windows[,4:5], ncol=2), rtrange=as.matrix(windows[,2:3], ncol=2), sampleidx=Sample, rt=rt)
		ans <- list()
		for(i in 1:length(eic@eic)){
			eics <- eic@eic[[i]]
			eics <- lapply(seq(along=eics), function(x,dat, ID) if(nrow(dat[[x]])!=0) {data.frame(Peak.ID=ID[[x]], dat[[x]])} else {}, dat=eics, ID=windows[,1])
			eics <- do.call('rbind', eics)
			eics$Sample <- names(eic@eic)[i]
			if(!missing(group)){
				eics$Group <- sampleInfo(object, all=TRUE)[eics$Sample, group]
			} else {}
			ans[[i]] <- eics
		}
		ans <- do.call('rbind', ans)
		ans <- as.data.frame(ans)
		
		# Plotting
        p <- ggplot(data=ans, aes(x=rt, y=intensity, group=Sample)) + facet_wrap(~Peak.ID, scales='free')
        p <- p + theme_bw() + xlab('Retention time (sec)') + ylab('Intensity')
        if(missing(group)){
            p <- p + geom_line()
        } else {
            p <- p + geom_line(aes(colour=Group)) + scale_colour_hue(group)
        }
        p
    }
)
### getChrom
### Extracts the TIC and BPC for all samples. If the @chromatograms slot is empty it will get the data from the raw filles and add it to the complist object
setMethod(
    'getChrom', 'Complist',
    function(object, type, Sample, outlier.rm, mix.rm, retcor, clean, objectname){
		
		# Cleans old data
		if(!missing(clean)){
			if(clean==TRUE){
				object@chromatograms$TIC <- NULL
				object@chromatograms$BPC <- NULL
			} else {}
		} else {}
		
		# Set default
		if(missing(objectname)){
			objectname <- deparse(substitute(object))
		} else {}
		if(missing(type)){
			type <- 'TIC'
		} else {}
        if(missing(Sample)){
            Sample <- row.names(object@Sample.info)
        } else {
            fail <- !(Sample %in% row.names(object@Sample.info))
            if(sum(fail) > 0){
                stop('Unknown sample: ', paste(Sample[fail], collapse=', '))
            }
        }
		if(missing(outlier.rm)){
			outlier.rm <- TRUE
		} else {}
		if(missing(mix.rm)){
			mix.rm <- TRUE
		} else {}
		if(missing(retcor)){
			retcor <- TRUE
		}
        if(outlier.rm){
            if(all(!Sample %in% sampleNames(object, outlier.rm=outlier.rm, mix.rm=mix.rm))){
                stop('Provided samples marked as outlier')
            } else{}
        } else {}
		Sample <- Sample[Sample %in% sampleNames(object, outlier.rm=outlier.rm, mix.rm=mix.rm)]
		
		# Extracts chromatograms from raw files if missing
		if(is.null(object@chromatograms$TIC)){
			paths <- object@Sample.info$Filepath
			Samplename <- row.names(object@Sample.info)
			cat('Extracting chromatographic data: ')
			flush.console()
			TIC <- list()
			BPC <- list()
			for(i in 1:length(paths)){
				cat(Samplename[i], ' ')
				flush.console()
				xraw <- xcmsRaw(as.character(paths[i]))
				TIC[[i]] <- apply(xraw@env$profile,2,sum)
				BPC[[i]] <- apply(xraw@env$profile,2,max)
				gc()
			}
			object@chromatograms <- list(TIC=TIC, BPC=BPC)
			assign(objectname, object, envir=.GlobalEnv)
			gc()
			cat('\n')
			flush.console()
		}
		
		# Create output
        Sampleind <- which(row.names(object@Sample.info) %in% Sample)
        ans <- list()
        for(i in 1:length(Sampleind)){
            if(retcor){
                TIC <- data.frame(Time=object@xcmsSet@rt$corrected[[Sampleind[i]]], Intensity=object@chromatograms$TIC[[Sampleind[i]]], Sample=Sample[i], Type='TIC', stringsAsFactors=FALSE)
                BPC <- data.frame(Time=object@xcmsSet@rt$corrected[[Sampleind[i]]], Intensity=object@chromatograms$BPC[[Sampleind[i]]], Sample=Sample[i], Type='BPC', stringsAsFactors=FALSE)
            } else {
                TIC <- data.frame(Time=object@xcmsSet@rt$raw[[Sampleind[i]]], Intensity=object@chromatograms$TIC[[Sampleind[i]]], Sample=Sample[i], Type='TIC', stringsAsFactors=FALSE)
                BPC <- data.frame(Time=object@xcmsSet@rt$raw[[Sampleind[i]]], Intensity=object@chromatograms$BPC[[Sampleind[i]]], Sample=Sample[i], Type='BPC', stringsAsFactors=FALSE)
            }
            ans[[i]] <- rbind(TIC, BPC)
        }
        ans <- do.call('rbind', ans)
        ans <- subset(ans, ans$Type %in% type)
        ans
    }
)
### filterAnova
### Given a term in the design matrix perform ANOVA on the peaks in the design matrix and report p-value (beware of multiple comparison problem)
setMethod(
    'filterAnova', 'Complist',
    function(object, cov, p=0.05, outlier.rm=TRUE){
		
		# Input check
        if(sum(cov %in% names(object@Sample.info)) != length(cov)){
            err <- cov[which(!(cov %in% names(object@Sample.info)))]
            stop('No sample info with name <', paste(err, collapse=', '), '> defined\n')
        } else {}
		
		# ANOVA
        for(i in 1:length(cov)){
			rem <- getRemove(object, what='Sample', outlier=outlier.rm, mix=TRUE)
			if(length(rem) != 0){
				dat <- object@raw[,-rem]
            } else {
                dat <- object@raw
            }
            pval <- apply(dat, 1, function(x, cov) anova(lm(x~cov))$'Pr(>F)'[1], cov=sampleInfo(object, outlier.rm=outlier.rm, mix.rm=TRUE)[, cov[i]])
            pval <- pval < p
            object@filter$Peaks[,cov[i]] <- pval
        }
        object
    }
)
### filterNZV
### Filters peaks with near zero variance
setMethod(
	'filterNZV', 'Complist',
	function(object, freqCut=95/5, uniqueCut=10, outlier.rm=TRUE){
		NZV <- nearZeroVar(t(object@raw[, which(!object@filter$Samples$Outlier)]), freqCut=freqCut, uniqueCut=uniqueCut)
		filter <- rep(TRUE, nrow(object@raw))
		filter[NZV] <- FALSE
		object@filter$Peaks$NZV <- filter
		object
	}
)
### filterCor
### Filters peaks showing high correlation to each others
setMethod(
	'filterCor', 'Complist',
	function(object, cutoff=0.75, ...){
		Cor <- findCorrelation(cor(t(object@raw), ...), cutoff=cutoff)
		filter <- rep(TRUE, nrow(object@raw))
		filter[Cor] <- FALSE
		object@filter$Cor <- filter
		object
	}
)
### modelPCA
### Creates a PCA model based on the specification given in the parameters
setMethod(
    'modelPCA', 'Complist',
    function(object, name, outlier.rm=TRUE, mix.rm=TRUE, unmatched=FALSE, FDR=TRUE, filter=NULL, method='nipals', nPcs=2, scale='uv', center=TRUE, cv='q2'){
        
		# Extract data
		dat <- getRaw(object, outlier.rm=outlier.rm, mix.rm=mix.rm, filter=filter, FDR=FDR, unmatched=unmatched)
        pInf <- getPeakinfo(object, outlier.rm=outlier.rm, filter=filter, FDR=FDR, unmatched=unmatched)
        ID <- getBestmatch(object, outlier.rm=outlier.rm, FDR=FDR)
		ID <- ID[which(rownames(ID) %in% rownames(pInf)), ]
		newID <- data.frame(matrix(NA, nrow=nrow(pInf), ncol=ncol(ID)))
		names(newID) <- names(ID)
		ind <- match(rownames(pInf), rownames(ID))
		newID[ind[!is.na(ind)],] <- ID
		ID <- cbind(pInf, newID)
		
		# Create model
		model <- pca(t(dat), method=method, nPcs=nPcs, scale=scale, center=center, cv=cv)
        ans <- list(Type='PCA', Model=model, Var.class=ID)
        object@models[[name]] <- ans
        object
    }
)
### modelCaret
### Allows taping into the caret package's vast number of modelling tools
setMethod(
	'modelCaret', 'Complist',
	function(object, method, name, predict, outlier.rm=TRUE, mix.rm=TRUE, unmatched=FALSE, FDR=TRUE, filter=NULL, nTrain, preProcess=NULL, weights=NULL, trControl=trainControl(), tuneGrid=NULL, tuneLength=3, ...){
		
		# Extract Data
		dat <- getRaw(object, outlier.rm=outlier.rm, mix.rm=mix.rm, filter=filter, unmatched=unmatched, FDR=FDR)
		Y <- sampleInfo(object, outlier.rm=outlier.rm, mix.rm=mix.rm)[, predict]
		if(!missing(nTrain)){
			if(nTrain < 0 | nTrain > 1){
				stop('nTrain must be a proportion: 0 < nTrain < 1')
			} else {}
			trainIndex <- sample(1:ncol(dat), round(nTrain*ncol(dat)))
			dat <- dat[, trainIndex]
			testset <- list(X=t(dat[, -trainIndex]), Y=Y[-trainIndex])
			Y <- Y[trainIndex]
		} else {
			testset <- NULL
		}
		pInf <- getPeakinfo(object, outlier.rm=outlier.rm, filter=filter, unmatched=unmatched, FDR=FDR)
		ID <- getBestmatch(object, outlier.rm=outlier.rm, filter=filter, FDR=FDR)
		ID <- ID[which(rownames(ID) %in% rownames(pInf)), ]
		ID <- cbind(pInf, ID)
		
		# Create model
		model <- train(t(dat), Y, method=method, preProcess=preProcess, weights=weights, trControl=trControl, tuneGrid=tuneGrid, tuneLength=tuneLength, ...)
		ans <- list(Type=method, Model=model, Testset=testset, Var.class=ID)
		object@models[[name]] <- ans
		object
	}
)
### predictCaret
### Used to make predictions based on models created with modelCaret
setMethod(
	'predictCaret', 'Complist',
	function(object, model, newdata='Testset', type='predictions'){
		
		# Test input
		if(!model %in% names(object@models)){
			stop(paste('No model named: ', name, sep=''))
		} else {}
		caret <- object@models[[model]]$Model
		if(class(caret) != 'train'){
			stop('Model must be created using modelCaret()')
		} else {}
		
		# Perform prediction
		caret <- list(caret)
		names(caret) <- model
		if(newdata == 'none'){
			if(type=='predictions'){
				ans <- extractPrediction(caret)
			} else if(type=='probability'){
				ans <- extractProb(caret)
			}
		} else if(newdata=='Testset'){
			if(type=='predictions'){
				ans <- extractPrediction(caret, testX, testY)
			} else if(type=='probability'){
				ans <- extractProb(caret)
			}
		}
	}
)
### plotScore
### Plot scores for a given model (At the moment only PCA's
setMethod(
    'plotScore', 'Complist',
    function(object, model, title, pc1=1, pc2=2, biplot=FALSE, ncomp=max(c(pc1,pc2)), scText, scPoint, scColour, scArea, scConnect, scHide, loText, loPoint, loColour, loHide, alpha=0.95, na.rm=TRUE){
        
		# Input check
		if(!(model %in% names(object@models))){
            stop('No model called <', model, '> is defined...\n', sep='')
        } else {}
        if(!missing(scColour)){
            if(sum(scColour %in% names(object@Sample.info)) == 0){
                stop('No sample class called <', scColour, '> defined...\n')
            } else {}
        } else {}
        loRMc <- FALSE
        loRMp <- FALSE
        if(!missing(scColour) & !missing(loColour)){
            loRMc <- TRUE
            cat('Colour parameter for loadings removed. Only support for either scores or loadings.\n')
        } else {}
        if(!missing(scPoint) & !missing(loPoint)){
            if(scPoint != '.' & loPoint != '.'){
                loRMp <- TRUE
                cat('Point parameter for loadings removed. Only support for either scores or loadings.\n')
            } else {}
        } else {}
		
		# For models from modelPCA
        if(object@models[[model]]$Type == 'PCA'){
            mod <- object@models[[model]]$Model
            sco <- as.data.frame(mod@scores)[,1:ncomp]
            lo <- as.data.frame(mod@loadings)[,1:ncomp]
            int <- simconf(sco, alpha)
            angle <- seq(-pi, pi, length = 50)
            df <- data.frame(cox = sin(angle)*int[pc1,2], coy = cos(angle)*int[pc2,2])
            loadscale <- int[min(c(pc1, pc2)), 2]/max(lo[,min(c(pc1, pc2))])*1.05
            lo <- lo*loadscale
            exp <- round(100*mod@R2, 2)
            sco <- sco[,c(pc1, pc2)]
            names(sco) <- c('x','y')
            lo <- lo[, c(pc1, pc2)]
            names(lo) <- c('x','y')
            
			## Prepare scores
            if(!missing(scHide)){
                hide <- list()
                for(i in 1:length(scHide)){
                    hide[[i]] <- sampleNames(object, all=TRUE)[sampleInfo(object, all=TRUE)[,names(scHide)[i]] %in% scHide[[i]]]
                }
                hide <- unique(do.call('c', hide))
                sco <- sco[which(!(rownames(sco) %in% hide)), ]
            }
            if(!missing(scArea)){
                area <- data.frame(sco, Area=sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% row.names(sco), scArea])
                area <- dhull(area, 'x', 'y', 'Area')
                area$Area <- factor(area$Area)
            } else {}
            if(!missing(scConnect)){
                if(length(scConnect) == 2){
                    connect <- data.frame(sco, sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% row.names(sco), scConnect])
                    connect <- ddply(connect, c(scConnect), function(z) if(nrow(z)!=0) c(mean(z$x), mean(z$y)))
                    names(connect) <- c('Path', 'Group', 'x', 'y')
                    connect <- droplevels(connect[connect$Group %in% names(table(connect$Group) != 1)[table(connect$Group) != 1],])
                } else {
                    connect <- data.frame(sco, Path=sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% row.names(sco), scConnect])
                    connect <- ddply(connect, 'Path', function(z) if(nrow(z)!=0) c(mean(z$x), mean(z$y)))
                    names(connect) <- c('Path', 'x', 'y')
                    connect <- data.frame(connect, Group='A')
                }
            } else {}
            if(!missing(scText)){
                if(scText == 'Sample'){
                    sco$Text <- as.character(row.names(sco))
                } else {
                    sco$Text <- sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% row.names(sco), scText]
                }
            } else if(!missing(scPoint)){
                if(!(scPoint == '.')){
                    sco$Point <- factor(sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% row.names(sco), scPoint])
                } else {}
            } else {}
            if(!missing(scColour)){
                sco$Colour <- sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% row.names(sco), scColour]
            } else {}
            sco <- droplevels(sco)
			
            ## Create plot
            p <- ggplot(data=sco, aes(x=x, y=y))
            if(!missing(scArea)){
                p <- p + geom_polygon(data=area, aes(fill=Area), alpha=I(0.2)) + scale_fill_hue(scArea)
                p <- p + geom_path(data=area, aes(group=Area))
            } else {}
            if(!missing(scConnect)){
                if(length(scConnect) == 1){
                    p <- p + geom_path(data=connect, aes(linetype=Group), arrow=arrow(length=unit(0.1, 'inches'))) + scale_linetype(scConnect[1], breaks='A', labels='')
                } else if(length(unique(connect$Group)) > 12){
                    connect$Group1 <- 'A'
                    p <- p + geom_path(data=connect, aes(group=Group, linetype=Group1), arrow=arrow(length=unit(0.1, 'inches'))) + scale_linetype(scConnect[1], breaks='A', labels='')
                } else {
                    p <- p + geom_path(data=connect, aes(group=Group, linetype=Group), arrow=arrow(length=unit(0.1, 'inches'))) + scale_linetype(scConnect[2])
                }
            } else {}
            if(missing(scColour)){
                if(!missing(scText)){
                    p <- p + geom_text(aes(label=Text), size=I(3))
                } else if(!missing(scPoint)){
                    if(scPoint == '.'){
                        p <- p + geom_point()
                    } else {
                        p <- p + geom_point(aes(shape=Point), na.rm=TRUE) + scale_shape(scPoint)
                    }
                } else {
                    p <- p + geom_point(shape=I(''))
                }
            } else {
                if(!missing(scText)){
                    p <- p + geom_text(aes(colour=Colour, label=Text), size=I(3))
                    if(is.numeric(sco$Colour)){
                        p <- p + scale_colour_gradient(scColour, low='yellow', high='red')
                    } else {
                        p <- p + scale_colour_hue(scColour)
                    }
                } else if(!missing(scPoint)){
                    if(scPoint == '.'){
                        p <- p + geom_point(aes(colour=Colour))
                    } else {
                        p <- p + geom_point(aes(colour=Colour, shape=Point)) + scale_shape(scPoint)
                    }
                    if(is.numeric(sco$Colour)){
                        p <- p + scale_colour_gradient(scColour, low='yellow', high='red')
                    } else {
                        p <- p + scale_colour_hue(scColour)
                    }
                } else {
                    p <- p + geom_point(shape=I(''))
                }
            }
            p <- p + xlab(paste('\nPC', pc1, ' (', exp[pc1], ' %)', sep='')) + ylab(paste('PC', pc2, ' (', exp[pc2], ' %)', sep='')) + theme_bw()
            p <- p + geom_path(aes(x=cox, y=coy, label=NULL), data=df, linetype=I(2) , colour=I('red'))
            
			## Add loadings
            if (biplot != FALSE){
				
                ### Prepare loadings
                if(!missing(loHide)){
                    hide <- list()
                    for(i in 1:length(loHide)){
                        hide[[i]] <- rownames(object@models[[model]]$Var.class)[object@models[[model]]$Var.class[,names(loHide)[i]] %in% loHide[[i]]]
                    }
                    hide <- unique(do.call('c', hide))
                } else {}
                if(!missing(loText)){
                    lo$Text <- object@models[[model]]$Var.class[,loText]
                } else if(!missing(loPoint) & !loRMp){
                    if(!(loPoint == '.')){
                        lo$Point <- factor(object@models[[model]]$Var.class[,loPoint])
                    } else {}
                } else {
                    lo$Text <- rownames(lo)
                }
                if(!missing(loColour)){
                    lo$Colour <- object@models[[model]]$Var.class[,loColour]
                } else {}
                if(na.rm){
                    lo <- lo[complete.cases(lo),]
                } else {}
                if(!missing(loHide)){
                    lo <- droplevels(lo[which(!(row.names(lo) %in% hide)),])
                } else {}
				
                ### Plot loadings
                if(missing(loColour) | loRMc){
                    if(!missing(loText)){
                        p <- p + geom_text(data=lo, aes(label=Text), size=I(2))
                    } else if(!missing(loPoint) & !loRMp){
                        if(loPoint == '.'){
                            p <- p + geom_point(data=lo)
                        } else {
                            p <- p + geom_point(data=lo, aes(shape=Point)) + scale_shape(loPoint)
                        }
                    } else {
                        p <- p + geom_text(data=lo, aes(label=Text), size=I(2))
                    }
                } else {
                    if(!missing(loText)){
                        p <- p + geom_text(data=lo, aes(colour=Colour, label=Text), size=I(2))
                    } else if(!missing(loPoint) & !loRMp){
                        if(loPoint == '.'){
                            p <- p + geom_point(data=lo, aes(colour=Colour))
                        } else{
                            p <- p + geom_point(data=lo, aes(colour=Colour, shape=Point)) + scale_shape(loPoint)
                        }
                    } else {
                        p <- p + geom_text(data=lo, aes(colour=Colour, label=Text), size=I(2))
                    }
                    if(is.numeric(lo$Colour)){
                        p <- p + scale_colour_gradient(loColour, low='yellow', high='red')
                    } else {
                        p <- p + scale_colour_hue(loColour)
                    }
                }
            } else {}
        } else {}
	if(!missing(title)){
            p <- p + ggtitle(title)
        } else {}
        p
    }
)
### plotStat
### Plot statistics for a model (only PCA supported)
setMethod(
    'plotStat', 'Complist',
    function(object, model, title, nPcs){
        if(object@models[[model]]$Type == 'PCA'){
            mod <- object@models[[model]]$Model
            if(missing(nPcs)){
                nPcs <- mod@nPcs
            } else {}
            df <- as.data.frame(matrix(0, ncol=3, nrow=3*nPcs))
            names(df) <- c('PC', 'value', 'type')
            df$PC <- rep(1:nPcs, 3)
            df$value <- c(mod@R2[1:nPcs], mod@R2cum[1:nPcs], mod@cvstat[1:nPcs])
            df$type <- rep(c('R2', 'R2cum', 'Q2'), each=nPcs)
            p <- ggplot(data=df, aes(x=factor(PC), y=value, fill=type)) + geom_histogram(position='dodge', colour=I('black')) + theme_bw()
            p <- p + scale_x_discrete('\nPrincipal component', breaks= 1:nPcs, labels=colnames(mod@scores)[1:nPcs]) + scale_y_continuous('', limits=c(0,1)) + scale_fill_hue('Model statistic')
        }
		if(!missing(title)){
            p <- p + ggtitle(title)
        } else {}
        p
    }
)
### plotEval
### Creates a scatterplot of Hottelling T2 and DmodX for the given PCA model
setMethod(
    'plotEval', 'Complist',
    function(object, model, title, scText, scPoint, scColour, scHide){
        mod <- object@models[[model]]$Model
        s0<-sqrt(sum(residuals(mod)^2)/((mod@nObs-mod@nPcs-1)*(mod@nVar-mod@nPcs)))
		df <- data.frame(matrix(0, ncol=3, nrow=nrow(mod@scores)))
		names(df) <- c('DmodX', 'T2', 'Sample')
		df$Sample <- rownames(mod@scores)
		
		# Prepares plot info
		if(!missing(scHide)){
            hide <- list()
            for(i in 1:length(scHide)){
                hide[[i]] <- sampleNames(object, all=TRUE)[sampleInfo(object, all=TRUE)[,names(scHide)[i]] %in% scHide[[i]]]
            }
            hide <- unique(do.call('c', hide))
            df <- df[which(!(df$Sample %in% hide)), ]
        } else {}
        if(!missing(scText)){
            if(scText == 'Sample'){
                df$Text <- as.character(df$Sample)
            } else {
                df$Text <- sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% df$Sample, scText]
            }
        } else if(!missing(scPoint)){
            if(!(scPoint == '.')){
                df$Point <- factor(sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% df$Sample, scPoint])
            } else {}
        } else {}
        if(!missing(scColour)){
            df$Colour <- sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% df$Sample, scColour]
        } else {}
        df <- droplevels(df)
		
		# Calculates statistics
		df$DmodX <- sqrt(apply(residuals(mod)^2, 1, sum)/(mod@nVar-mod@nPcs))/s0
		df$T2 <- apply(mod@scores^2/mod@sDev^2, 1, sum)
		T2conf <- mod@nPcs*(mod@nObs-1)/(mod@nObs-mod@nPcs)*qf(0.95,mod@nPcs,mod@nObs-mod@nPcs)
		DmodXconf <- sd(df$DmodX)*qf(0.95,mod@nPcs,mod@nObs-mod@nPcs)
		
		# Creates plot
        p <- ggplot(data=df, aes(x=T2, y=DmodX)) + theme_bw()
        p <- p + geom_vline(xintercept=T2conf, linetype=2, colour='red') + geom_hline(yintercept=DmodXconf, linetype=2, colour='red')
        if(!missing(scText)){
            if(missing(scColour)){
                p <- p + geom_text(aes(label=Text), size=I(3))
            } else {
                p <- p + geom_text(aes(label=Text, colour=Colour), size=I(3))
                if(is.numeric(df$Colour)){
                    p <- p + scale_colour_gradient(scColour, low='yellow', high='red')
                } else {
                    p <- p + scale_colour_hue(scColour)
                }
            }
        } else if(!missing(scPoint)){
            if(scPoint == '.'){
                if(missing(scColour)){
                    p <- p + geom_point()
                } else {
                    p <- p + geom_point(aes(colour=Colour)) + scale_shape(scPoint)
                    if(is.numeric(df$Colour)){
                        p <- p + scale_colour_gradient(scColour, low='yellow', high='red')
                    } else {
                        p <- p + scale_colour_hue(scColour)
                    }
                }
            } else {
                if(missing(scColour)){
                    p <- p + geom_point(aes(shape=Point)) + scale_shape(scPoint)
                } else {
                    p <- p + geom_point(aes(shape=Point, colour=Colour)) + scale_shape(scPoint)
                    if(is.numeric(df$Colour)){
                        p <- p + scale_colour_gradient(scColour, low='yellow', high='red')
                    } else {
                        p <- p + scale_colour_hue(scColour)
                    }
                }
            }
        }
        if(!missing(title)){
            p <- p + ggtitle(title)
        } else {}
		p
    }
)
### plotMDS
### Calculates a two dimensional MDS and creates a scatterplot
setMethod(
	'plotMDS', 'Complist',
	function(object, outlier.rm=TRUE, mix.rm=TRUE, filter=NULL, unmatched=FALSE, FDR=TRUE, scText, scPoint, scColour, scArea, scConnect, scHide, title){
		
		# Calculates data
		raw <- getRaw(object, outlier.rm=outlier.rm, mix.rm=mix.rm, filter=filter, unmatched=unmatched, FDR=FDR)
		scaled <- apply(raw, 2, scale)
		mds <- cmdscale(dist(t(scaled)), eig=TRUE)
		explained <- round(c(mds$eig[1]/sum(mds$eig)*100, mds$eig[2]/sum(mds$eig)*100))
		sco <- as.data.frame(mds$points)
		names(sco) <- c('x', 'y')
		
		# Preparing score info
		if(!missing(scHide)){
			hide <- list()
			for(i in 1:length(scHide)){
				hide[[i]] <- sampleNames(object, all=TRUE)[sampleInfo(object, all=TRUE)[,names(scHide)[i]] %in% scHide[[i]]]
			}
			hide <- unique(do.call('c', hide))
			sco <- sco[which(!(rownames(sco) %in% hide)), ]
		}
		if(!missing(scArea)){
			area <- data.frame(sco, Area=sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% rownames(sco), scArea])
			area <- dhull(area, 'x', 'y', 'Area')
			area$Area <- factor(area$Area)
		} else {}
		if(!missing(scConnect)){
			if(length(scConnect) == 2){
				connect <- data.frame(sco, sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% rownames(sco), scConnect])
				connect <- ddply(connect, c(scConnect), function(z) if(nrow(z)!=0) c(mean(z$x), mean(z$y)))
				names(connect) <- c('Path', 'Group', 'x', 'y')
				connect <- droplevels(connect[connect$Group %in% names(table(connect$Group) != 1)[table(connect$Group) != 1],])
			} else {
				connect <- data.frame(sco, Path=sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% rownames(sco), scConnect])
				connect <- ddply(connect, 'Path', function(z) if(nrow(z)!=0) c(mean(z$x), mean(z$y)))
				names(connect) <- c('Path', 'x', 'y')
				connect <- data.frame(connect, Group='A')
			}
		} else {}
		if(!missing(scText)){
			if(scText == 'Sample'){
				sco$Text <- as.character(rownames(sco))
			} else {
				sco$Text <- sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% rownames(sco), scText]
			}
		} else if(!missing(scPoint)){
			if(!(scPoint == '.')){
				sco$Point <- factor(sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% rownames(sco), scPoint])
			} else {}
		} else {}
		if(!missing(scColour)){
			sco$Colour <- sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% rownames(sco), scColour]
		} else {}
		
		# Creates plot
		p <- ggplot(data=sco, aes(x=x, y=y))
		if(!missing(scArea)){
			p <- p + geom_polygon(data=area, aes(fill=Area), alpha=I(0.2)) + scale_fill_hue(scArea)
			p <- p + geom_path(data=area, aes(group=Area))
		} else {}
		if(!missing(scConnect)){
			if(length(scConnect) == 1){
				p <- p + geom_path(data=connect, aes(linetype=Group), arrow=arrow(length=unit(0.1, 'inches'))) + scale_linetype(scConnect[1], breaks='A', labels='')
			} else if(length(unique(connect$Group)) > 12){
				connect$Group1 <- 'A'
				p <- p + geom_path(data=connect, aes(group=Group, linetype=Group1), arrow=arrow(length=unit(0.1, 'inches'))) + scale_linetype(scConnect[1], breaks='A', labels='')
			} else {
				p <- p + geom_path(data=connect, aes(group=Group, linetype=Group), arrow=arrow(length=unit(0.1, 'inches'))) + scale_linetype(scConnect[2])
			}
		} else {}
		if(missing(scColour)){
			if(!missing(scText)){
				p <- p + geom_text(aes(label=Text), size=I(3))
			} else if(!missing(scPoint)){
				if(scPoint == '.'){
					p <- p + geom_point()
				} else {
					p <- p + geom_point(aes(shape=Point), na.rm=TRUE) + scale_shape(scPoint)
				}
			} else {
				p <- p + geom_point(shape=I(''))
			}
		} else {
			if(!missing(scText)){
				p <- p + geom_text(aes(colour=Colour, label=Text), size=I(3))
				if(is.numeric(sco$Colour)){
					p <- p + scale_colour_gradient(scColour, low='yellow', high='red')
				} else {
					p <- p + scale_colour_hue(scColour)
				}
			} else if(!missing(scPoint)){
				if(scPoint == '.'){
					p <- p + geom_point(aes(colour=Colour))
				} else {
					p <- p + geom_point(aes(colour=Colour, shape=Point)) + scale_shape(scPoint)
				}
				if(is.numeric(sco$Colour)){
					p <- p + scale_colour_gradient(scColour, low='yellow', high='red')
				} else {
					p <- p + scale_colour_hue(scColour)
				}
			} else {
				p <- p + geom_point(shape=I(''))
			}
		}
		p <- p + xlab(paste('Eig%: ', explained[1], sep='')) + ylab(paste('Eig%: ', explained[2], sep='')) + theme_bw()
		if(!missing(title)){
			p <- p + ggtitle(title)
		} else {}
		p		
	}
)
### plotCont
### Creates a contribution plot for a selected subset of samples in a PCA model
setMethod(
    'plotCont', 'Complist',
    function(object, model, title, groupA, groupB, continuousVar=FALSE, loColour, loHide, sort=FALSE){
        if(object@models[[model]]$Type == 'PCA'){
			
			# Extract info on Group A
            A <- list()
            for(i in 1:length(groupA)){
                if(names(groupA)[i] == 'Sample'){
                    A[[i]] <- sampleNames(object, all=TRUE)[sampleNames(object, all=TRUE) %in% groupA[[i]]]
                } else {
                    A[[i]] <- sampleNames(object, all=TRUE)[sampleInfo(object, all=TRUE)[,names(groupA)[i]] %in% groupA[[i]]]
                }
            }
            A <- unique(do.call('c', A))
            vec <- object@models[[model]]$Model@scores[which(row.names(object@models[[model]]$Model@scores) %in% A),]
            if(!is.null(dim(vec))){
				vec <- apply(vec, 2, mean)
            } else {}
			
			# Extracts information on Group B
            if(!missing(groupB)){
                B <- list()
                for(i in 1:length(groupB)){
                    if(names(groupB)[i] == 'Sample'){
                        B[[i]] <- sampleNames(object, all=TRUE)[sampleNames(object, all=TRUE) %in% groupB[[i]]]
                    } else {
                        B[[i]] <- sampleNames(object, all=TRUE)[sampleInfo(object, all=TRUE)[,names(groupB)[i]] %in% groupB[[i]]]
                    }    
                }
                B <- unique(do.call('c', B))
                diffvec <- object@models[[model]]$Model@scores[which(row.names(object@models[[model]]$Model@scores) %in% B),]
                if(!is.null(dim(diffvec))){
                    diffvec <- apply(diffvec, 2, mean)
                } else {}
                vec <- vec-diffvec
            } else {}
			
			# Calculates contribution
            mat <- sqrt(solve(diag(object@models[[model]]$Model@sDev))) %*% t(as.matrix(object@models[[model]]$Model@loadings))
            res <- apply(vec*mat, 2, sum)
			
			# Prepares plot info
            if(continuousVar){
                res <- data.frame(variable=as.numeric(names(res)), value=res)
            } else {
                res <- data.frame(variable=factor(names(res), levels=names(res)), value=res)
                if(!missing(loColour)){
                    res$Colour <- object@models[[model]]$Var.class[,loColour]
                } else {}
            }
            if(!missing(loHide)){
                hide <- list()
                for(i in 1:length(loHide)){
                    hide[[i]] <- rownames(object@models[[model]]$Var.class)[object@models[[model]]$Var.class[,names(loHide)[i]] %in% loHide[[i]]]
                }
                hide <- unique(do.call('c', hide))
                res <- droplevels(res[which(!(res$variable %in% hide)),])
            } else {}
			if(sort){
				res <- transform(res, variable=reorder(variable, -value))
			} else {}
			
			# Creates plot
            p <- ggplot(data=res, aes(x=variable, y=value)) + theme_bw() + theme(axis.text.x = element_text(size=8, angle=90, hjust=1, vjust=0.5))
            p <- p + xlab('') + ylab('')
            if(continuousVar){
                p <- p + geom_area(fill=I('grey'))
            } else {
                if(!missing(loColour)){
                    p <- p + geom_bar(aes(fill=Colour))
                    if(is.numeric(res$Colour)){
                        p <- p + scale_fill_gradient(loColour)
                    } else {
                        p <- p + scale_fill_hue(loColour)
                    }
                } else {
                    p <- p + geom_bar(fill=I('grey'))
                }
				p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
            }
        } else {}
        if(!missing(title)){
            p <- p + ggtitle(title)
        } else {}
        p
    }
)
### plotData
### Creates a plot of the underlying data in a PCA model
setMethod(
    'plotData', 'Complist',
    function(object, model, title, samples, groupA, mean=FALSE, scale=FALSE, continuousVar=FALSE, loColour, loHide){
        
		# Input check
		if(missing(samples)){
            samples <- groupA
        } else {}
        if(object@models[[model]]$Type == 'PCA'){
			
			# Input check
            samp <- list()
            for(i in 1:length(samples)){
                if(names(samples)[i] == 'Sample'){
                    samp[[i]] <- sampleNames(object, all=TRUE)[sampleNames(object, all=TRUE) %in% samples[[i]]]
                } else {
                    samp[[i]] <- sampleNames(object, all=TRUE)[sampleInfo(object, all=TRUE)[,names(samples)[i]] %in% samples[[i]]]
                }
            }
			
			# Create data
            samp <- unique(do.call('c', samp))
            data <- data.frame(object@models[[model]]$Model@completeObs)
            if(scale){
                data <- data.frame(scale(data, center=FALSE))
            } else {}
            if(!missing(loHide)){
                hide <- list()
                for(i in 1:length(loHide)){
                    hide[[i]] <- rownames(object@models[[model]]$Var.class)[object@models[[model]]$Var.class[,names(loHide)[i]] %in% loHide[[i]]]
                }
                hide <- unique(do.call('c', hide))
                data <- droplevels(data[which(!(row.names(data) %in% hide)),])
            } else {}
            avg <- apply(data, 2, mean)
            avg <- data.frame(variable=names(avg), value=avg, zero=rep(0, length(avg)), fill='Global Average')
            data <- data[which(rownames(data) %in% samp), ]
            if(mean){
                data <- apply(data, 2, mean)
                data <- data.frame(variable=names(data), value=data, group='Mean')
            } else {
                data <- data.frame(Sample=row.names(data), data)
                data <- melt(data, id=1)
            }
            avg$variable <- sub('X', '', avg$variable)
            data$variable <- sub('X', '', data$variable)
            if(!missing(loHide)){
                data <- droplevels(data[which(!(data$variable %in% hide)), ])
                avg <- droplevels(avg[which(!(avg$variable %in% hide)), ])
            }
            if(!mean){
                avg <- ddply(data.frame(Sample=unique(data$Sample)), .(Sample), function(x, dat) data.frame(dat, Sample=x), dat=avg)
            }
            if(!missing(loColour)){
                colour <- data.frame(variable=row.names(object@models[[model]]$Var.class), Colour=object@models[[model]]$Var.class[,loColour])
                data <- merge(data, colour, all.x=TRUE)
            } else {}
            if(continuousVar){
                data$variable <- as.numeric(as.character(data$variable))
                avg$variable <- as.numeric(as.character(data$variable))
            } else {}
			
			# Create plot
            p <- ggplot(data=data, aes(x=variable, y=value)) + theme_bw() + theme(axis.text.x = element_text(size=8, angle=90, hjust=1, vjust=0.5))
            if(continuousVar){
                if(mean){
                    p <- p + geom_line(aes(linetype=group)) + geom_area(data=avg, aes(fill=fill), alpha=I(0.5))
                    p <- p + scale_linetype('') + scale_fill_hue('')
                } else {
                    p <- p + geom_line() + geom_area(data=avg, aes(fill=fill), alpha=I(0.5)) + facet_wrap(~Sample) + scale_fill_hue('')
                }
            } else {
                if(!missing(loColour)){
                    p <- p + geom_bar(aes(fill=Colour), position='dodge')
                    if(is.numeric(data$Colour)){
                        p <- p + scale_fill_gradient(loColour)
                    } else {
                        p <- p + scale_fill_hue(loColour)
                    }
                } else {
                    p <- p + geom_bar(fill='grey')
                }
                p <- p + geom_point(data=avg, aes(shape=fill))
                p <- p + geom_linerange(data=avg, aes(ymin=zero, ymax=value))
                p <- p + scale_shape('')
                if(!mean){
                    p <- p + facet_wrap(~Sample)
                } else {}
				p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
            }
            p <- p + xlab('') + ylab('')
            if(!missing(title)){
                p <- p + ggtitle(title)
            } else {}
        } else {}
        p
    }
)
### plotDendro
### Creates a dendrogram plot with optional factor info in a facet below
setMethod(
    'plotDendro', 'Complist',
    function(object, model, compound, title, method='ward', distance='euclidean', outlier.rm=TRUE, mix.rm=TRUE, unmatched=FALSE, FDR=TRUE, what, factor, cutHeight, cutGroup, filter=NULL){
        if(missing(model)){
			
			# Extract data
            data <- getRaw(object, outlier.rm=outlier.rm, mix.rm=mix.rm, filter=filter, unmatched=unmatched, FDR=FDR)
			pInf <- getPeakinfo(object, outlier.rm=outlier.rm, filter=filter, unmatched=unmatched, FDR=FDR)
			ID <- getBestmatch(object, outlier.rm=outlier.rm, FDR=FDR)
			ID <- ID[which(rownames(ID) %in% rownames(pInf)), ]
			ID <- cbind(pInf, ID)
            if(what == 'Sample'){
				
				# Samples
				## Create clustering
                hc <- hclust(dist(t(data), method=distance), method=method)
                den <- dendro_data(hc)
				
				## Cut the tree
                if(!missing(cutHeight)){
                    Colour <- cutree(hc, h=cutHeight)
                    den$label$Colour <- factor(Colour[match(den$label$text, names(Colour))])
                    den$segments$Colour <- NA
                    den$segments$Colour[which(den$segments$yend == 0)] <- den$label$Colour
                    den1 <- den
                    den1$segments <- den1$segments[!is.na(den1$segments$Colour), ]
                } else if(!missing(cutGroup)){
                    Colour <- cutree(hc, k=cutGroup)
                    den$label$Colour <- factor(Colour[match(den$label$text, names(Colour))])
                    den$segments$Colour <- NA
                    den$segments$Colour[which(den$segments$yend == 0)] <- den$label$Colour
                    den1 <- den
                    den1$segments <- den1$segments[!is.na(den1$segments$Colour), ]
                } else {}
                den$label <- merge(den$label, data.frame(text=sampleNames(object, all=TRUE), sampleInfo(object, all=TRUE)), sort=FALSE)
				
				## Create factor data
                if(!missing(factor)){
                    fac <- recplot(sampleInfo(object)[,factor], size=max(den$segments$y))
                    fac[[1]] <- data.frame(fac[[1]], den$label)
                    fac[[2]]$x <- max(fac[[1]]$right)+1
                    lev <- c('Dendrogram', factor)
                    fac[[1]]$plot <- factor(factor, levels=lev)
                    fac[[2]]$plot <- factor(factor, levels=lev)
                    den$segments$plot <- factor('Dendrogram', levels=lev)
                    den$label$plot <- factor('Dendrogram', levels=lev)
                } else {}
				
            } else if(what == 'Peak'){
				
				# Peaks
				## Create clustering
                hc <- hclust(dist(data, method=distance), method=method)
                den <- dendro_data(hc)
				
				## Cut the tree
                if(!missing(cutHeight)){
                    Colour <- cutree(hc, h=cutHeight)
                    den$label$Colour <- factor(Colour[match(den$label$text, names(Colour))])
                    den$segments$Colour <- NA
                    den$segments$Colour[which(den$segments$yend == 0)] <- den$label$Colour
                    den1 <- den
                    den1$segments <- den1$segments[!is.na(den1$segments$Colour), ]
                } else if(!missing(cutGroup)){
                    Colour <- cutree(hc, k=cutGroup)
                    den$label$Colour <- factor(Colour[match(den$label$text, names(Colour))])
                    den$segments$Colour <- NA
                    den$segments$Colour[which(den$segments$yend == 0)] <- den$label$Colour
                    den1 <- den
                    den1$segments <- den1$segments[!is.na(den1$segments$Colour), ]
                } else {}
                den$label <- merge(den$label, data.frame(text=row.names(ID), ID), sort=FALSE)
				
				## Create factor data
                if(!missing(factor)){
                    fac <- recplot(ID[,factor], size=max(den$segments$y))
                    fac[[1]] <- data.frame(fac[[1]], den$label)
                    fac[[2]]$x <- max(fac[[1]]$right)+1
                    lev <- c('Dendrogram', factor)
                    fac[[1]]$plot <- factor(factor, levels=lev)
                    fac[[2]]$plot <- factor(factor, levels=lev)
                    den$segments$plot <- factor('Dendrogram', levels=lev)
                    den$label$plot <- factor('Dendrogram', levels=lev)
                } else {}
            } else {
                stop('No plotDendro method for <', what, '>\n', sep='')
            }
        } else if(object@models[[model]]$Type == 'PCA'){
			
			# Dendrogram based on PCA model
            if(what == 'Score' | what == 'Sample'){
				
				## For Scores
				### Create clustering
                data <- object@models[[model]]$Model@scores
                hc <- hclust(dist(data, method=distance), method=method)
                den <- dendro_data(hc)
				
				### Cut tree
                if(!missing(cutHeight)){
                    Colour <- cutree(hc, h=cutHeight)
                    den$label$Colour <- factor(Colour[match(den$label$text, names(Colour))])
                    den$segments$Colour <- NA
                    den$segments$Colour[which(den$segments$yend == 0)] <- den$label$Colour
                    den1 <- den
                    den1$segments <- den1$segments[!is.na(den1$segments$Colour), ]
                } else if(!missing(cutGroup)){
                    Colour <- cutree(hc, k=cutGroup)
                    den$label$Colour <- factor(Colour[match(den$label$text, names(Colour))])
                    den$segments$Colour <- NA
                    den$segments$Colour[which(den$segments$yend == 0)] <- den$label$Colour
                    den1 <- den
                    den1$segments <- den1$segments[!is.na(den1$segments$Colour), ]
                } else {}
				
				### Create factor data
                if(!missing(factor)){
                    fac <- recplot(sampleInfo(object, all=TRUE)[sampleNames(object, all=TRUE) %in% row.names(data), factor], size=max(den$segments$y))
                    fac[[1]] <- data.frame(fac[[1]], den$label)
                    fac[[2]]$x <- max(fac[[1]]$right)+1
                    lev <- c('Dendrogram', factor)
                    fac[[1]]$plot <- factor(factor, levels=lev)
                    fac[[2]]$plot <- factor(factor, levels=lev)
                    den$segments$plot <- factor('Dendrogram', levels=lev)
                    den$label$plot <- factor('Dendrogram', levels=lev)
                } else {}
            } else if(what == 'Loading'){
				
				## For loadings
				### Create clustering
                data <- object@models[[model]]$Model@loadings
                hc <- hclust(dist(data, method=distance), method=method)
                den <- dendro_data(hc)
				
				### Cut the tree
                if(!missing(cutHeight)){
                    Colour <- cutree(hc, h=cutHeight)
                    den$label$Colour <- factor(Colour[match(den$label$text, names(Colour))])
                    den$segments$Colour <- NA
                    den$segments$Colour[which(den$segments$yend == 0)] <- den$label$Colour
                    den1 <- den
                    den1$segments <- den1$segments[!is.na(den1$segments$Colour), ]
                } else if(!missing(cutGroup)){
                    Colour <- cutree(hc, k=cutGroup)
                    den$label$Colour <- factor(Colour[match(den$label$text, names(Colour))])
                    den$segments$Colour <- NA
                    den$segments$Colour[which(den$segments$yend == 0)] <- den$label$Colour
                    den1 <- den
                    den1$segments <- den1$segments[!is.na(den1$segments$Colour), ]
                } else {}
				
				### Create factor info
                if(!missing(factor)){
                    fac <- recplot(object@models[[model]]$Var.class[,factor], size=max(den$segments$y))
                    fac[[1]] <- data.frame(fac[[1]], den$label)
                    fac[[2]]$x <- max(fac[[1]]$right)+1
                    lev <- c('Dendrogram', factor)
                    fac[[1]]$plot <- factor(factor, levels=lev)
                    fac[[2]]$plot <- factor(factor, levels=lev)
                    den$segments$plot <- factor('Dendrogram', levels=lev)
                    den$label$plot <- factor('Dendrogram', levels=lev)
                } else {}
            }
        }
		
		# Create plot
        p <- ggplot(data=segment(den)) + theme_bw()
        p <- p + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
        if(!missing(cutHeight) | !missing(cutGroup)){
            p <- p + geom_segment(data=segment(den1), aes(x=x, y=y, xend=xend, yend=yend, colour=factor(Colour)), size=2)
        } else {}
        if(!missing(factor)){
            if(!missing(cutHeight) | !missing(cutGroup)){
                p <- p + geom_rect(data=fac[[1]], aes(xmin=left, xmax=right, ymin=down, ymax=up, fill=Colour), colour=I('black'))
            } else {
                p <- p + geom_rect(data=fac[[1]], aes(xmin=left, xmax=right, ymin=down, ymax=up), colour=I('black'), fill=I('grey'))
            }
            p <- p + facet_grid(plot~., scales='free', height=c(4,1))
            p <- p + scale_y_continuous(breaks=c(fac[[2]]$y, seq(0, max(den$segments$y), by=inter(max(den$segments$y)))), labels=c(as.character(fac[[2]]$Levels), as.character(seq(0, max(den$segments$y), by=inter(max(den$segments$y))))))
        }
        p <- p + scale_x_continuous(breaks=den$label$x, labels=den$label$label) + theme(axis.text.x=element_text(angle=-90, hjust=0))
        p <- p + xlab('') + ylab('') + theme(legend.position='none')
        if(!missing(title)){
            p <- p + ggtitle(title)
        } else {}
        p
    }
)
### plotHeat
### Creates a heatplot based on ggHeat().
setMethod(
    'plotHeat', 'Complist',
    function(object, pep.info, sample.info, sample.name, outlier.rm=TRUE, mix.rm=TRUE, filter=NULL, FDR=TRUE, ...){
		
		# Input check
        if(!missing(pep.info)){
            if(sum(!pep.info %in% names(object@pepID@peplist)) != 0){
                stop('Unknown rowfactor(s). Specify names found in the peplist...')
            } else {}
        } else {}
        if(!missing(sample.info)){
            if(sum(!sample.info %in% names(object@Sample.info)) != 0){
                stop('Unknown colfactor(s). Specify names found in the sampleInfo...')
            } else {}
        } else {}
		
		# Creates data
		dat <- getRaw(object, outlier.rm=outlier.rm, mix.rm=mix.rm, filter=filter, FDR=FDR)
		if(!missing(sample.name)){
			colnames(dat) <- as.vector(sampleInfo(object, outlier.rm=outlier.rm, mix.rm=mix.rm)[, sample.name])
		} else {}
        int <- t(scale(t(dat), center=F, scale=apply(dat, 1, max)))*100
		int[is.na(int)] <- 0
        if(!missing(pep.info)){
            rinf <- getBestmatch(object, outlier.rm=outlier.rm, filter=filter, FDR=FDR)
            rinf <- rinf[,pep.info, drop=FALSE]
        } else {}
        if(!missing(sample.info)){
            sinf <- sampleInfo(object, outlier.rm=outlier.rm, mix.rm=mix.rm)
            sinf <- sinf[,sample.info, drop=FALSE]
        } else {}
		
		# Call to ggHeat
        if(missing(pep.info)){
            if(missing(sample.info)){
                ggHeat(int, dataname='Normalized\nIntensity', hideticks='rows', ...)
            } else {
                ggHeat(int, dataname='Normalized\nIntensity', hideticks='rows', colfactor=sinf, ...)
            }
        } else {
            if(missing(sample.info)){
                ggHeat(int, dataname='Normalized\nIntensity', hideticks='rows', rowfactor=rinf, ...)
            } else {
                ggHeat(int, dataname='Normalized\nIntensity', hideticks='rows', colfactor=sinf, rowfactor=rinf, ...)
            }
        }
    }
)
### plotCoverage
### Plots the coverage of the proteins in the database
setMethod(
	'plotCoverage', 'Complist',
	function(object, Sample, outlier.rm=TRUE, mix.rm=TRUE, filter=NULL, FDR=TRUE, useFDR=FALSE, protein, poswin, anova){
		if(missing(protein)){
			protein <- names(object@pepID@database)
		} else {}
		if(missing(Sample)){
			if(missing(anova)){
				raw <- data.frame(getBestmatch(object, outlier.rm=outlier.rm, filter=filter, FDR=FDR))
			} else {
				factor <- sampleInfo(object, outlier.rm=outlier.rm, mix.rm=mix.rm)[, anova]
				response <- as.list(data.frame(t(getRaw(object, outlier.rm=outlier.rm, mix.rm=mix.rm, unmatched=FALSE, FDR=FDR))))
				aov <- mapply(function(factor, response) anova(lm(response~factor))[['Pr(>F)']][1], response=response, MoreArgs=list(factor=factor))
				raw <- data.frame(getBestmatch(object, outlier.rm=outlier.rm, filter=filter, FDR=FDR), aov=aov)
			}
		} else {
			raw <- data.frame(getBestmatch(object, outlier.rm=outlier.rm, filter=filter, FDR=FDR), intensity=getRaw(object, outlier.rm=outlier.rm, filter=filter, FDR=FDR, unmatched=FALSE, normalize=T)[, Sample])
		}
		raw <- raw[which(raw$Protein %in% protein), ]
		ambig <- which(grepl('[0-9]/[0-9]', raw$Peptide.name))
		if(length(ambig) != 0){
			ambigRaw <- raw[ambig, ]
			raw <- raw[-ambig, ]
			ambigStart <- strsplit(sub(".*\\(([0-9/]+)\\-[0-9/]+\\)", '\\1', ambigRaw$Peptide.name, perl=T), '/')
			ambigEnd <- strsplit(sub(".*\\([0-9/]+\\-([0-9/]+)\\)", '\\1', ambigRaw$Peptide.name, perl=T), '/')
			ambigRaw <- ambigRaw[rep.int(1:nrow(ambigRaw), times=sapply(ambigStart, length)), ]
			ambigRaw$start <- as.numeric(unlist(ambigStart)) - 0.5
			ambigRaw$end <- as.numeric(unlist(ambigEnd)) + 0.5
		} else {}
		raw$start <- as.numeric(sub(".*\\(([0-9]+)\\-[0-9]+\\)", '\\1', raw$Peptide.name, perl=T)) - 0.5
		raw$end <- as.numeric(sub(".*\\([0-9]+\\-([0-9]+)\\)", '\\1', raw$Peptide.name, perl=T)) + 0.5
		if(length(ambig) != 0){
			raw <- rbind(raw, ambigRaw)
		} else {}
		raw <- ddply(raw, .(Protein), function(x) data.frame(x[order(x$start, x$end), ], yPos=1:nrow(x)))
		protein <- unique(raw$Protein)
		seq <- list()
		for(i in 1:length(protein)){
			res <- strsplit(as.character(object@pepID@database[which(names(object@pepID@database) %in% protein[i])]), '')[[1]]
			present <- as.vector(unlist((mapply('seq', raw$start[raw$Protein==protein[i]]+0.5, raw$end[raw$Protein==protein[i]]-0.5))))
			posFDR <- present
			for(j in 1:length(posFDR)){
				posFDR[j] <- min(raw$FDR[raw$Protein==protein[i] & raw$start < posFDR[j] & raw$end > posFDR[j]])
			}
			seq[[i]] <- data.frame(Protein=protein[i], x=1:length(res), y=-nrow(raw[raw$Protein==protein[i], ])*0.01, res=res, present=1:length(res) %in% present, FDR=1)
			seq[[i]]$FDR[present] <- posFDR
		}
		seq <- do.call('rbind', seq)
		if(!missing(poswin)){
			raw <- raw[raw$end >= poswin[1] & raw$start <= poswin[2], ]
			seq <- seq[seq$x >= poswin[1] & seq$x <= poswin[2], ]
			seq <- seq[seq$Protein %in% raw$Protein, ]
		} else {}
		seq1 <- seq[seq$present, ]
		seqScale <- diff(c(min(seq$x), max(seq$x)))/dev.size('cm')[1]
		if(seqScale < 10){
			if(seqScale > 4){
				seqScale <- 10/seqScale
			} else {
				seqScale <- 10/4
			}
		} else {
			seqScale <- 1
		}
		p <- ggplot(data=raw)
		p <- p + geom_text(data=seq, aes(x=x, y=y, label=res), colour='grey60', fontface=2, size=seqScale, hjust=0.5, vjust=0.5)
		if(useFDR & !FDR){
			p <- p + geom_text(data=seq1, aes(x=x, y=y, label=res, alpha=FDR), colour='red', fontface=2, size=seqScale, hjust=0.5, vjust=0.5, show_guide=FALSE)
			if(missing(Sample)){
				if(missing(anova)){
					p <- p + geom_segment(aes(x=start, xend=end, y=yPos, yend=yPos, alpha=FDR, size=FDR))
				} else {
					p <- p + geom_segment(aes(x=start, xend=end, y=yPos, yend=yPos, colour=aov, alpha=FDR, size=FDR))
					p <- p + scale_colour_gradient(paste('Anova on\n', anova, sep=''), limits=c(0,1), trans='sqrt', guide=guide_colourbar())
				}
			} else {
				p <- p + geom_segment(aes(x=start, xend=end, y=yPos, yend=yPos, colour=intensity, alpha=FDR, size=FDR))
				p <- p + scale_colour_gradient('Normalized\nintensity', low='yellow', high='red', guide=guide_colourbar())
			}
			p <- p + scale_alpha(limits=c(0, 1), range=c(1, 0))
			p <- p + scale_size(limits=c(0, 1), range=c(0.5, 0), trans='sqrt')
		} else {
			p <- p + geom_text(data=seq1, aes(x=x, y=y, label=res), colour='red', fontface=2, size=seqScale, hjust=0.5, vjust=0.5)
			if(missing(Sample)){
				if(missing(anova)){
					p <- p + geom_segment(aes(x=start, xend=end, y=yPos, yend=yPos))
				} else {
					p <- p + geom_segment(aes(x=start, xend=end, y=yPos, yend=yPos, colour=aov))
					p <- p + scale_colour_gradient(paste('Anova on\n', anova, sep=''), limits=c(0,1), trans='sqrt', guide=guide_colourbar())
				}
			} else {
				p <- p + geom_segment(aes(x=start, xend=end, y=yPos, yend=yPos, colour=intensity))
				p <- p + scale_colour_gradient('Normalized\nintensity', low='yellow', high='red', guide=guide_colourbar())
			}
		}
		p <- p + facet_grid(Protein~., space='free', scales='free')
		p <- p + theme_bw() + xlab('Position of residue')
		p <- p + scale_y_continuous('', breaks=NULL)
		p
	}
)
### plotSample
### Plots an overview of a sample either as a heatplot type or in 3d (requires rgl)
setMethod(
	'plotSample', 'Complist',
	function(object, Sample, ddd=FALSE, chrom='TIC', spec='sum', rtwin, mzwin, precursors=FALSE, nPeptides, filter=1000){
		if(!missing(mzwin)){
			if(diff(mzwin) < 100){
				pStep <- 0.1
			} else {
				pStep <- 1
			}
		} else {
			pStep <- 1
		}
		if(ddd){
			require(rgl)
			file <- as.character(sampleInfo(object, all=TRUE)$Filepath[sampleNames(object, all=TRUE) == Sample])
			para <- list(object=xcmsRaw(file, profstep=pStep))
			if(!missing(rtwin)){
				para$rtrange <- rtwin
			} else {}
			if(!missing(mzwin)){
				para$mzrange <- mzwin
			} else {}
			do.call('plotSurf', para)
		} else {
			file <- as.character(sampleInfo(object, all=TRUE)$Filepath[sampleNames(object, all=TRUE) == Sample])
			if(precursors){
				xRaw <- xcmsRaw(file, profstep=pStep, includeMSn=TRUE)
				prec <- data.frame(mz=xRaw@msnPrecursorMz, Time=xRaw@scantime[xRaw@msnPrecursorScan], shape='1', gridcol=factor(1, levels=c(1,2)), gridrow=factor(2, levels=c(1,2)))
			} else {
				xRaw <- xcmsRaw(file, profstep=pStep)
			}
			if(missing(rtwin)){
				rtwin <- range(xRaw@scantime)
			} else {}
			if(missing(mzwin)){
				mzwin <- range(xRaw@mzrange)
			} else {}
			prof <- xRaw@env$profile
			rts <- xRaw@scantime
			mzs <- seq(xRaw@mzrange[1], xRaw@mzrange[2], by=profStep(xRaw))
			prof <- prof[mzs >= mzwin[1] & mzs <= mzwin[2], rts >= rtwin[1] & rts <= rtwin[2]]
			rts <- rts[rts >= rtwin[1] & rts <= rtwin[2]]
			mzs <- mzs[mzs >= mzwin[1] & mzs <= mzwin[2]]
			rtScale <- 50 * (dev.size('cm')[1]-2) * ifelse(chrom=='none', 1, 0.75)
			mzScale <- 50 * dev.size('cm')[2] * ifelse(spec=='none', 1, 0.75)
			txScale <- (dev.size('cm')[1]-2) * 0.25 * 0.375
			if(length(rts) < rtScale & diff(mzwin) < mzScale){
				sScale <- max(c(length(rts), diff(mzwin)))/max(c(rtScale, mzScale))
			} else {
				sScale <- 1
			}
			colnames(prof) <- rts
			rownames(prof) <- mzs
			if(precursors){
				prec <- prec[prec$mz >= mzwin[1] & prec$mz <= mzwin[2] & prec$Time >=rtwin[1] & prec$Time <= rtwin[2], ]
			} else {}
			if(chrom != 'none'){
				if(chrom == 'TIC'){
					chromdat <- apply(prof,2,sum)
				} else if(chrom=='BPC'){
					chromdat <- apply(prof,2,max)
				} else {}
				chromdat <- data.frame(Intensity=chromdat, Time=rts, gridcol=factor(1, levels=c(1,2)), gridrow=factor(1, levels=c(1,2)))
			} else {}
			if(spec != 'none'){
				if(spec == 'sum'){
					specdat <- apply(prof,1,sum)
				} else if(spec == 'max'){
					specdat <- apply(prof,1,max)
				} else {}
				specdat <- data.frame(zero=0, Intensity=specdat, mz=mzs, gridcol=factor(2, levels=c(1,2)), gridrow=factor(2, levels=c(1,2)))
			} else {}
			prof <- melt(prof)
			names(prof) <- c('mz', 'Time', 'Intensity')
			prof <- prof[prof$Intensity > filter, ]
			prof$gridcol <- factor(1, levels=c(1,2))
			prof$gridrow <- factor(2, levels=c(1,2))
			if(spec != 'none' & chrom != 'none'){
				stat <- paste(round(min(prof$Time), digits=1), ' - ', round(max(prof$Time), digits=1), '\n', min(prof$mz), ' - ', max(prof$mz), '\n', sprintf('%.3g', max(prof$Intensity)), '\n', filter, '\n\n', profStep(xRaw), '\n~ ', round(diff(range(rts))/length(rts), digits=2), sep='')
				stat <- data.frame(label=stat, x=mean(range(specdat$Intensity))*1.55, y=max(chromdat$Intensity), gridcol=factor(2, levels=c(1,2)), gridrow=factor(1, levels=c(1,2)))
				units <- data.frame(label='sec\nDa\ncount\ncount\n\nDa\nscan/sec', x=mean(range(specdat$Intensity))*1.6, y=max(chromdat$Intensity), gridcol=factor(2, levels=c(1,2)), gridrow=factor(1, levels=c(1,2)))
				text <- data.frame(label='Retention time:\nMass/charge:\nMax Intensity:\nIntensity cutoff:\n\nBinwidth:\nScanspeed', x=0, y=max(chromdat$Intensity), gridcol=factor(2, levels=c(1,2)), gridrow=factor(1, levels=c(1,2)))
				if(precursors){
					text$label <- paste(text$label, '\n\n# Precursors', sep='')
					stat$label <- paste(stat$label, '\n\n', nrow(prec), sep='')
					units$label <- paste(units$label, '\n\n', sep='')
				}
				back <- data.frame(x=c(min(specdat$Intensity), min(specdat$Intensity), max(specdat$Intensity), max(specdat$Intensity)), y=c(min(chromdat$Intensity), max(chromdat$Intensity), max(chromdat$Intensity), min(chromdat$Intensity)), gridcol=factor(2, levels=c(1,2)), gridrow=factor(1, levels=c(1,2)))
			}
			p <- ggplot()
			p <- p + geom_point(data=prof, aes(x=Time, y=mz, colour=Intensity, size=log(Intensity)))
			p <- p + scale_size(range=c(0.05, 0.3)/sScale, guide='none')
			p <- p + scale_colour_gradient(low='yellow', high='red', guide=guide_colorbar(), trans='log10')
			if(chrom != 'none'){
				p <- p + geom_line(data=chromdat, aes(x=Time, y=Intensity))
			} else {}
			if(spec != 'none'){
				p <- p + geom_segment(data=specdat, aes(x=zero, xend=Intensity, y=mz, yend=mz))
			} else {}
			if(chrom != 'none'){
				if(spec != 'none'){
					p <- p + facet_grid(gridrow~gridcol, scales='free', height=c(1, 3), width=c(3, 1))
					p <- p + geom_polygon(data=back, aes(x=x, y=y), fill='white', colour='white', size=10)
					p <- p + geom_text(data=text, aes(x=x, y=y, label=label), hjust=0, vjust=1, size=txScale)
					p <- p + geom_text(data=stat, aes(x=x, y=y, label=label), hjust=1, vjust=1, size=txScale)
					p <- p + geom_text(data=units, aes(x=x, y=y, label=label), hjust=0, vjust=1, size=txScale)	
				} else {
					p <- p + facet_grid(gridrow~., scales='free', height=c(1, 3))
				}
			} else {
				if(spec != 'none'){
					p <- p + facet_grid(.~gridcol, scales='free', width=c(3, 1))
				} else {}
			}
			if(precursors){
				p <- p + geom_point(data=prec, aes(x=Time, y=mz, shape=shape))
				p <- p + scale_shape('', solid=FALSE, labels='Precursors')
			} else {}
			p <- p + scale_x_continuous(expand=c(0.05, 0)) + scale_y_continuous(expand=c(0.05, 0))
			p <- p + theme_bw() + theme(strip.background=element_blank(), strip.text=element_blank())
			p <- p + theme(axis.text.x=element_text(size=5, angle=315, hjust=0, vjust=1), axis.text.y=element_text(size=5))
			p
		}
	}
)
### plotChromatogram
### Plots the selected chromatograms in overlay and optionally labels the top intensity peptides
setMethod(
    'plotChromatogram', 'Complist',
    function(object, Sample, retcor, type, colour, outlier.rm, mix.rm, filter=NULL, FDR=TRUE, rtwin, title, nPeptides){
		
		# Remember name of object in case extracted chromatograms must be assigned
		objectname <- deparse(substitute(object))
		
		# Set default
		if(missing(retcor)){
			retcor <- TRUE
		} else {}
		if(missing(type)){
			type <- c('TIC', 'BPC')
		} else {}
		if(missing(outlier.rm)){
			outlier.rm <- TRUE
		} else {}
		if(missing(mix.rm)){
			mix.rm <- TRUE
		} else {}
		if(missing(colour)){
			colour <- NULL
		}
		
		# Extract chromatograms
		if(missing(Sample)){
			Sample <- sampleNames(object, outlier.rm=outlier.rm, mix.rm=mix.rm)
		} else {}
		chroms <- getChrom(object, Sample=Sample, retcor=retcor, type=type, outlier.rm=outlier.rm, mix.rm=mix.rm, objectname=objectname)
		if(!is.null(colour)){
			if(colour == 'Sample'){
				chroms$Colour <- chroms$Sample
			} else {
				if(!colour %in% names(sampleInfo(object))){
					stop('Wrong colour argument: Not part of Sample info...')
				} else {}
				col <- sampleInfo(object, all=TRUE)[, colour, drop=FALSE]
				names(col) <- 'Colour'
				col$Sample <- rownames(col)
				chroms <- merge(chroms, col, all.x=TRUE, sort=FALSE)
			}
		}
		
		# Labelling of hign intensity peaks
		if(!missing(nPeptides)){
			
			# Find top peaks
			raw <- getRaw(object, outlier.rm = outlier.rm, mix.rm = mix.rm, filter=filter, unmatched=FALSE, FDR=FDR)
			peakInfo <- getPeakinfo(object, outlier.rm = outlier.rm, filter=filter, unmatched=FALSE, FDR=FDR)
			pepInfo <- getBestmatch(object, outlier.rm = outlier.rm, filter=filter, FDR=FDR)
			raw <- raw[, which(colnames(raw) %in% Sample)]
			if (!missing(rtwin)) {
				ind <- which(peakInfo$rtmin > rtwin[1] & peakInfo$rtmax < rtwin[2])
				peakInfo <- peakInfo[ind, ]
				raw <- raw[ind, ]
				pepInfo <- pepInfo[ind, ]
			} else {}
			if (length(Sample) > 1) {
				raw <- apply(raw, 1, max)
			} else {}
			peakInfo <- peakInfo[order(raw, decreasing = TRUE), ]
			pepInfo <- pepInfo[order(raw, decreasing = TRUE), ]
			raw <- raw[order(raw, decreasing = TRUE)]
			
			if(nPeptides > nrow(pepInfo)){
				nPeptides <- nrow(pepInfo)
				cat('Only ', nPeptides, ' peptides in plotregion\n', sep='')
			} else {}
			peakInfo <- peakInfo[1:nPeptides, ]
			pepInfo <- pepInfo[1:nPeptides, ]
			raw <- raw[1:nPeptides]
			
			# Create annotation data
			annotateTIC <- cbind(pepInfo, peakInfo, raw)
			annotateTIC$Intensity <- NA
			annotateBPC <- annotateTIC
			annotateTIC$Type <- 'TIC'
			annotateBPC$Type <- 'BPC'
			window <- (max(chroms$Time) - min(chroms$Time))/500
			for (i in 1:nrow(annotateTIC)) {
				subset <- chroms[which(chroms$Type == 'TIC' & chroms$Time > annotateTIC[i, 'rtmed'] - window & chroms$Time < annotateTIC[i, 'rtmed'] + window), ]
				annotateTIC[i, 'Intensity'] <- max(subset$Intensity)
				subset <- chroms[which(chroms$Type == 'BPC' & chroms$Time > annotateBPC[i, 'rtmed'] - window & chroms$Time < annotateBPC[i, 'rtmed'] + window), ]
				annotateBPC[i, 'Intensity'] <- max(subset$Intensity)
			}
			pretty <- annotateTIC[, c('Peptide.ID', 'Peptide.name', 'Sequence', 'FDR')]
			pretty <- pretty[order(pretty$Peptide.ID), ]
			annotate <- rbind(annotateTIC, annotateBPC)
			annotate <- annotate[which(annotate$Type %in% type), ]
		} else {}
		
		# Create plot
		if(is.null(colour)){
			p <- ggplot(data=chroms, aes(x=Time, y=Intensity, group=Sample))
		} else {
			p <- ggplot(data=chroms, aes(x=Time, y=Intensity, group=Sample, colour=Colour))
			p <- p + scale_colour_hue(colour)
		}
        p <- p + geom_line() + facet_grid(Type~., scales='free')
        p <- p + theme_bw()
        if(!missing(rtwin)){
            p <- p + coord_cartesian(xlim=rtwin)
        } else {}
        if(!missing(title)){
            p <- p + ggtitle(title)
        } else {}
		if(!missing(nPeptides)){
			p <- p + geom_text(data = annotate, aes(x = rtmed, y = Intensity, group = NULL, colour = NULL, label = Peptide.ID), size=3, vjust=-0.5, show_guide=FALSE)
			print(pretty, row.names=F)
		}
        p
    }
)
### plotDetection
### Plots EICs on common scale for high quality peptide identification and splits them whether on whether the peak has been detected
setMethod(
	'plotDetection', 'Complist',
	function(object){
		peptides <- getPeplist(object, FDR=TRUE, unmatched=TRUE)
		eics <- getEIC(object@xcmsSet, mzrange=matrix(c(peptides$mz-0.1, peptides$mz+0.1), ncol=2), rtrange=matrix(c(peptides$Retention.time-10, peptides$Retention.time+10), ncol=2), sampleidx=sampnames(object@xcmsSet), rt='corrected')
		data <- list()
		for(i in 1:length(eics@eic[[1]])){
			TMPdata <- adply(1:length(eics@eic), 1, function(x, eic, name) data.frame(Sample=name[x], eic[[x]][[i]]), eic=eics@eic, name=names(eics@eic))
			TMPdata[, 'Peak'] <- as.character(i)
			data[[i]] <- TMPdata
		}
		data <- do.call('rbind', data)
		data$Sample.Peak <- paste(data$Sample, ': ', data$Peak, sep='')
		data$Matched <- rep(FALSE, nrow(data))
		data$Matched[which(data$Peak %in% as.character(which(peptides$Peptide.ID %in% getPeplist(object, FDR=TRUE, unmatched=FALSE)$Peptide.ID)))] <- TRUE
		qplot(rt, intensity, data=data, geom='line', colour=Peak, group=Sample.Peak) + theme_bw() + guides(colour='none') + facet_grid(Matched~.)
	}
)
### plotRetcor
### Plots deviation profile, before and after views of the chromatograms to evaluate the effect of the retention time correction
setMethod(
    'plotRetcor', 'Complist',
    function(object, outlier.rm, mix.rm, rtwin, type, title){
		
		# Get object name in case of chromatogram assignement during extraction
		objectname <- deparse(substitute(object))
		
		# Set defaults
		if(missing(outlier.rm)){
			outlier.rm <- TRUE
		} else {}
		if(missing(mix.rm)){
			mix.rm <- TRUE
		} else {}
		if(missing(type)){
			type <- 'BPC'
		} else {}
		
		# Extract deviation profile
        devi <- list()
        for(i in 1:length(object@xcmsSet@rt$raw)){
            ans <- data.frame(
                Time=object@xcmsSet@rt$corrected[[i]],
                Intensity=object@xcmsSet@rt$raw[[i]] - object@xcmsSet@rt$corrected[[i]],
                Sample=rownames(object@xcmsSet@phenoData)[i],
                Type='RT Dev.'
            )
            devi[[i]] <- ans
        }
        devi <- do.call('rbind', devi)
        devi <- devi[devi$Sample %in% sampleNames(object, outlier.rm=outlier.rm, mix.rm=mix.rm), ]
		
		# Extract chromatograms
        rawRT <- getChrom(object, type=type, outlier.rm=outlier.rm, mix.rm=mix.rm, retcor=FALSE, objectname=objectname)
        rawRT$Type <- 'Raw'
        corRT <- split(rawRT, rawRT$Sample)
        sInd <- match(names(corRT), rownames(object@xcmsSet@phenoData))
        for(i in 1:length(corRT)){
            corRT[[i]] <- corRT[[i]][order(corRT[[i]]$Time), ]
            corRT[[i]]$Time <- object@xcmsSet@rt$corrected[[sInd[i]]]
        }
        corRT <- do.call('rbind', corRT)
        corRT$Type <- 'Corrected'
        ans <- rbind(devi, rawRT, corRT)
		ans <- droplevels(ans)
		
		# Create plot
        p <- ggplot(data=ans, aes(x=Time, y=Intensity, group=Sample, colour=Sample))
        p <- p + geom_line() + facet_grid(Type~., scales='free', height=c(1,2,2))
        p <- p + theme_bw()
        if(!missing(rtwin)){
            p <- p + coord_cartesian(xlim=rtwin)
        } else {}
        if(!missing(title)){
            p <- p + ggtitle(title)
        } else {}
        p
    }
)
### getCont
### Extracts the contributions plottet by plotCont
setMethod(
    'getCont', 'Complist',
    function(object, model, groupA, groupB, loHide){
        if(object@models[[model]]$Type == 'PCA'){
			
			# Extracts data for group A
            A <- list()
            for(i in 1:length(groupA)){
                if(names(groupA)[i] == 'Sample'){
                    A[[i]] <- sampleNames(object, all=TRUE)[sampleNames(object, all=TRUE) %in% groupA[[i]]]
                } else {
                    A[[i]] <- sampleNames(object, all=TRUE)[sampleInfo(object, all=TRUE)[,names(groupA)[i]] %in% groupA[[i]]]
                }
            }
            A <- unique(do.call('c', A))
            vec <- object@models[[model]]$Model@scores[which(row.names(object@models[[model]]$Model@scores) %in% A),]
            if(!is.null(dim(vec))){
				vec <- apply(vec, 2, mean)
            } else {}
			
			# Extract data for group B
            if(!missing(groupB)){
                B <- list()
                for(i in 1:length(groupB)){
                    if(names(groupB)[i] == 'Sample'){
                        B[[i]] <- sampleNames(object, all=TRUE)[sampleNames(object, all=TRUE) %in% groupB[[i]]]
                    } else {
                        B[[i]] <- sampleNames(object, all=TRUE)[sampleInfo(object, all=TRUE)[,names(groupB)[i]] %in% groupB[[i]]]
                    }    
                }
                B <- unique(do.call('c', B))
                diffvec <- object@models[[model]]$Model@scores[which(row.names(object@models[[model]]$Model@scores) %in% B),]
                if(!is.null(dim(diffvec))){
                    diffvec <- apply(diffvec, 2, mean)
                } else {}
                vec <- vec-diffvec
            } else {}
			
			# Calculate contributions
            mat <- sqrt(solve(diag(object@models[[model]]$Model@sDev))) %*% t(as.matrix(object@models[[model]]$Model@loadings))
            res <- apply(vec*mat, 2, sum)
            res <- data.frame(variable=factor(names(res), levels=names(res)), Contribution=res)
			
			# Hides selcted variables
            if(!missing(loHide)){
                hide <- list()
                for(i in 1:length(loHide)){
                    hide[[i]] <- row.names(object@models[[model]]$Var.class)[object@models[[model]]$Var.class[,names(loHide)[i]] %in% loHide[[i]]]
                }
                hide <- unique(do.call('c', hide))
                res <- droplevels(res[which(!(res$variable %in% hide)),])
            } else {}
            info <- object@models[[model]]$Var.class[row.names(object@models[[model]]$Var.class) %in% res$variable, ]
            res <- data.frame(res, info)
            res$variable <- NULL
        }
        res
    }
)
### getRemove
### Extract index of samples, peptides or peaks to remove given a set of filterering options
setMethod(
    'getRemove', 'Complist',
    function(object, what, outlier=FALSE, mix=FALSE, filter, match=FALSE, annotation=FALSE, index, FDR=FALSE){
		
		# For Samples
        if(what == 'Sample'){
            ans <- list()
            if(outlier){
                ans$o <- object@filter$Samples$Outlier
            } else {}
			if(mix){
				ans$m <- object@filter$Samples$Mix
			}
            if(!missing(index)){
                ans$i2 <- rep(FALSE, nrow(object@Sample.info))
				if(is.numeric(index)){
					ans$i2[index] <- TRUE
				} else {
					ans$i2[index %in% sampleNames(object)] <- TRUE
				}
            } else {}
            ans <- do.call('cbind', ans)
            if(!is.null(ans)){
                ans <- apply(ans, 1, any)
                ans <- which(ans)
            } else {
                ans <- integer()
            }
		
		# For peaks
        } else if(what == 'Peak'){
            ans <- list()
            if(outlier){
                ans$o <- object@filter$Peaks$Outlier
            } else {}
            if(!missing(filter)){
                ans$f <- object@filter$Peaks[,which(names(object@filter$Peaks) %in% filter)]
            } else {}
            if(match){
                ans$m <- !object@filter$Peaks$Matched
            } else {}
			if(FDR){
				ans$f2 <- rep(TRUE, nrow(object@raw))
				ans$f2[do.call('c', object@IDindex[object@filter$Peptides$FDR])] <- FALSE
			} else {}
            if(annotation){
                ans$a <- laply(object@annotation$Features, length) == 0
            } else {}
            if(!missing(index)){
                ans$i2 <- rep(FALSE, nrow(object@raw))
                ans$i2[index] <- TRUE
            } else {}
            ans <- do.call('cbind', ans)
            if(!is.null(ans)){
                ans <- apply(ans, 1, any)
                ans <- which(ans)
            } else {
                ans <- integer()
            }
			
		# For peptides
        } else if(what == 'Peptide'){
            ans <- list()
            if(outlier){
                ans$o <- object@filter$Peptides$Outlier
            } else {}
            if(match){
                ans$m <- !object@filter$Peptides$Matched
            } else {}
            if(!missing(index)){
                ans$i2 <- rep(FALSE, nrow(object@filter$Peptides))
                ans$i2[index] <- TRUE
            } else {}
			if(!missing(filter)){
				ans$f <- object@filter$Peptides[,which(names(object@filter$Peptides) %in% filter)]
			} else {}
			if(FDR){
				ans$f2 <- !object@filter$Peptides$FDR
			}
            ans <- do.call('cbind', ans)
            if(!is.null(ans)){
                ans <- apply(ans, 1, any)
                ans <- which(ans)
            } else {
                ans <- integer()
            }
        } else {
            stop('Unknown input <', what, '>\n', sep='')
        }
        ans
    }
)
### intensityReport
### Extract intensity (normalized by default) and id for matches and writes the result to csv files
setMethod(
	'intensityReport', 'Complist',
	function(object, folder=NULL, outlier.rm=TRUE, mix.rm=TRUE, unmatched=FALSE, FDR=TRUE, normalize=TRUE, height=FALSE){
		
		# Extracts data
		dat <- getRaw(object, outlier.rm=outlier.rm, mix.rm=mix.rm, unmatched=unmatched, FDR=FDR, normalize=normalize, height=height)
		ID <- getPepPeakIndex(object, from='Peak', outlier.rm=outlier.rm, FDR=FDR)
		ID <- lapply(ID, function(x) object@pepID@peplist[x, ])
			
		# Formatting for xlsx files
		if(is.null(folder)){
			filename <- 'Intensity.xlsx'
		} else {
			dir.create(file.path(getwd(), folder), showWarnings=FALSE)
			filename <- file.path(folder, 'Intensity.xlsx')
		}
		ID1 <- lapply(ID, function(x) selectMatch(x, IDtype=object@pepID@type))
		ID1 <- do.call('rbind', ID1)
		dat <- data.frame(ID1, dat, check.names=FALSE)
		ID <- data.frame(Peak.ID=rep(seq(along=ID), sapply(ID, nrow)), do.call('rbind', ID))
		rownames(dat) <- 1:nrow(dat)
		wb <- createWorkbook()
		dat_sheet <- createSheet(wb, 'Data')
		ID_sheet <- createSheet(wb, 'ID')
		addDataFrame(dat, dat_sheet)
		addDataFrame(ID, ID_sheet, row.names=FALSE)
		saveWorkbook(wb, filename)
	}
)
### peakReport
### Creates a pdf file with EICs for all peaks identified as peptides
setMethod(
	'peakReport', 'Complist',
	function(object, folder=NULL, outlier.rm=TRUE, mix.rm=TRUE, FDR=TRUE, group=NULL, retcor=TRUE){
		if(retcor){
			rt <- 'corrected'
		} else {
			rt <- 'raw'
		}
		if(is.null(folder)){
			filename <- 'Peakplot.pdf'
		} else {
			dir.create(file.path(getwd(), folder), showWarnings=FALSE)
			filename <- file.path(folder, 'Peakplot.pdf')
		}
		
		peptides <- getPeplist(object, outlier.rm=outlier.rm, unmatched=FALSE, FDR=FDR)
		peaks <- object@Peak.info[do.call('c', object@IDindex[as.numeric(rownames(peptides))]),]
		peptides <- peptides[rep(1:nrow(peptides), sapply(object@IDindex[as.numeric(rownames(peptides))], length)), ]
			
		windows <- data.frame(peptides$Peptide.ID, peaks[, c('rtmin', 'rtmax', 'mzmin', 'mzmax')])
		windows[, 2] <- windows[, 2] - 5
		windows[, 3] <- windows[, 3] + 5
		names(windows) <- c('ID', 'rtmin', 'rtmax', 'mzmin', 'mzmax')
		
		# Extract EIC's
		cat('Extracting EIC for all samples...\n\n')
		flush.console()
		samples <- sampleNames(object, outlier.rm=outlier.rm, mix.rm=mix.rm)
		eic <- getEIC(object@xcmsSet, mzrange=as.matrix(windows[,4:5], ncol=2), rtrange=as.matrix(windows[,2:3], ncol=2), sampleidx=samples, rt=rt)
		ans <- list()
		for(i in 1:length(eic@eic)){
			eics <- eic@eic[[i]]
			eics <- lapply(seq(along=eics), function(x,dat, ID) if(nrow(dat[[x]])!=0) {data.frame(Peak.ID=ID[[x]], dat[[x]])} else {}, dat=eics, ID=windows[,1])
			eics <- do.call('rbind', eics)
			eics$Sample <- names(eic@eic)[i]
			if(!is.null(group)){
				eics$Group <- object@Sample.info[i,group]
			} else {}
			ans[[i]] <- eics
			rm(eics)
		}
		rm(eic)
		gc()
		ans <- do.call('rbind', ans)
		
		# Creating pdf
		cat('\nPlotting...\n')
		flush.console()
		subs <- split(unique(ans$Peak.ID), rep(1:ceiling(length(unique(ans$Peak.ID))/12), each=12)[1:length(unique(ans$Peak.ID))])
		pdf(filename, width=9, height=12)
		for(i in 1:length(subs)){
			subdat <- subset(ans, ans[,1] %in% subs[[i]])
			p <- ggplot(data=subdat, aes(x=rt, y=intensity, group=Sample)) + facet_wrap(~Peak.ID, scales='free', ncol=3, nrow=4)
			p <- p + theme_bw() + xlab('Retention time (sec)') + ylab('Intensity') + ggtitle('Extracted Ion Chromatograms for all\nmatched peptides.\n')
			if(is.null(group)){
				p <- p + geom_line()
			} else {
				p <- p + geom_line(aes(colour=Group)) + scale_colour_hue(group)
			}
			gc()
			print(p)
		}
		invisible(dev.off())
	}
)
### pcaReport
### Creates a 2-page pdf with a simple PCA model
setMethod(
	'pcaReport', 'Complist',
	function(object, folder=NULL, outlier.rm=TRUE, mix.rm=TRUE, unmatched=FALSE, FDR=TRUE, filter=NULL, group=NULL, alpha=0.95){
		
		# Create model
		object <- modelPCA(object, 'PCA', outlier.rm=outlier.rm, mix.rm=mix.rm, unmatched=unmatched, FDR=FDR, filter=filter, nPcs=10)
		
		# Setting filename
		if(is.null(folder)){
			filename <- 'PCA.pdf'
		} else {
			dir.create(file.path(getwd(), folder), showWarnings=FALSE)
			filename <- file.path(folder, 'PCA.pdf')
		}
		
		# Creating plots
		if(is.null(group)){
			score1 <- plotScore(object, 'PCA', 'Scoreplot for PC1 and PC2', 1, 2, scText='Sample', alpha=alpha)
			score2 <- plotScore(object, 'PCA', 'Scoreplot for PC3 and PC4', 3, 4, scText='Sample', alpha=alpha)
		} else {
			score1 <- plotScore(object, 'PCA', 'Scoreplot for PC1 and PC2', 1, 2, scText='Sample', scColour=group, alpha=alpha)
			score2 <- plotScore(object, 'PCA', 'Scoreplot for PC3 and PC4', 3, 4, scText='Sample', scColour=group, alpha=alpha)
		}
		scoregrid <- plotmatrix(as.data.frame(object@models$PCA$Model@scores)) + theme_bw() + ggtitle('Overview of the 10 first PCs')
		stat <- plotStat(object, 'PCA', 'Statistics for the model')
		
		# Creating pdf
		pdf(filename, width=9, height=12)
		vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(2,1)))
		print(score1, vp=vplayout(1,1))
		print(score2, vp=vplayout(2,1))
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(2,1)))
		print(scoregrid, vp=vplayout(1,1))
		print(stat, vp=vplayout(2,1))
		invisible(dev.off())
	}
)
### sampleReport
### Creates a pdf file with sample plots for each sample
setMethod(
		'sampleReport', 'Complist',
		function(object, folder=NULL, outlier.rm=TRUE, mix.rm=TRUE, precursors=FALSE, chrom='TIC', spec='sum'){
			samples <- sampleNames(object, outlier.rm=outlier.rm, mix.rm=mix.rm)
			if(is.null(folder)){
				filename <- 'Samples.pdf'
			} else {
				dir.create(file.path(getwd(), folder), showWarnings=FALSE)
				filename <- file.path(folder, 'Samples.pdf')
			}
			
			# Create plot
			pdf(filename, height=9, width=12)
			for(i in 1:length(samples)){
				p <- plotSample(object, samples[i], chrom=chrom, spec=spec, precursors=precursors)
				p <- p + ggtitle(paste('Raw data overview of ', samples[i], sep=''))
				print(p)
			}
			invisible(dev.off())
		}
)
### coverageReport
### Creates a pdf file with coverage plots for each sample
setMethod(
	'coverageReport', 'Complist',
	function(object, folder=NULL, outlier.rm=TRUE, mix.rm=TRUE, FDR=TRUE, useFDR=FALSE, anova){
		if(!missing(anova)){
			if(is.null(folder)){
				filename <- 'Coverage_anova.pdf'
			} else {
				dir.create(file.path(getwd(), folder), showWarnings=FALSE)
				filename <- file.path(folder, 'Coverage_anova.pdf')
			}
			pdf(filename, height=12, width=9)
			for(i in 1:length(anova)){
				p <- plotCoverage(object, outlier.rm=outlier.rm, mix.rm=mix.rm, FDR=FDR, useFDR=useFDR, anova=anova[i])
				print(p)
			}
		} else {
			samples <- sampleNames(object, outlier.rm=outlier.rm, mix.rm=mix.rm)
			if(is.null(folder)){
				filename <- 'Coverage.pdf'
			} else {
				dir.create(file.path(getwd(), folder), showWarnings=FALSE)
				filename <- file.path(folder, 'Coverage.pdf')
			}
			
			# Create plot
			pdf(filename, height=12, width=9)
			for(i in 1:length(samples)){
				p <- plotCoverage(object, samples[i], outlier.rm=outlier.rm, FDR=FDR, useFDR=useFDR)
				p <- p + ggtitle(paste('Protein coverage for ', samples[i], sep=''))
				print(p)
			}
		}
		invisible(dev.off())
	}
)
### chromReport
### Creates a pdf file with chromatograms for each sample
setMethod(
	'chromReport', 'Complist',
	function(object, folder=NULL, outlier.rm, mix.rm){
		
		# Get object name in case of chromatogram assignement during extraction
		objectname <- deparse(substitute(object))
		
		# Set defaults
		if(missing(outlier.rm)){
			outlier.rm=TRUE
		} else {}
		if(missing(mix.rm)){
			mix.rm=TRUE
		} else {}
		if(is.null(folder)){
			filename <- 'Chromatograms.pdf'
		} else {
			dir.create(file.path(getwd(), folder), showWarnings=FALSE)
			filename <- file.path(folder, 'Chromatograms.pdf')
		}
		
		# Extracts data
		chrom <- getChrom(object, type=c('TIC', 'BPC'), outlier.rm=outlier.rm, mix.rm=mix.rm, objectname=objectname)
		
		# Create pdf
		pdf(filename, width=9, height=12)
		for(i in 1:length(unique(chrom$Sample))){
			sname <- unique(chrom$Sample)[i]
			schrom <- subset(chrom, chrom$Sample == sname)
			p <- ggplot(data=schrom, aes(x=Time, y=Intensity, group=Sample)) + facet_grid(Type~., scales='free')
			p <- p + geom_line() + theme_bw() + ggtitle(paste('MS chromatograms for sample:', sname, sep=' ')) + xlab('Retention time (sec)')
			print(p)
		}
		invisible(dev.off())
	}
)
### groupReport
### Creates an xls document containing sheets for each group extracted with plotHeat() and the arguments rGroup and/or cGroup. The sheets contain information on the peptides in each group
setMethod(
	'groupReport', 'Complist',
	function(object, group, folder){
		for(i in 1:length(group)){
			if(names(group)[i] == 'Rows'){
				data <- data.frame(Group = group[[i]], Peak = as.numeric(names(group[[i]])))
				data <- dlply(data, .(Group), function(x) getMatch(object, ID=x$Peak))
				wb <- createWorkbook()
				filename <- paste(names(group)[i], '.xlsx', sep='')
				if(!missing(folder)){
					dir.create(file.path(getwd(), folder), showWarnings=FALSE)
					filename <- file.path(folder, filename)
				}
				for(j in names(data)){
					group <- data[[j]]
					ans <- list()
					for(k in names(group)){
						ans[[k]] <- data.frame(Peak=k, group[[k]])
					}
					ans <- do.call('rbind', ans)
					sheet <- createSheet(wb, paste('Group ', j, sep=''))
					addDataFrame(ans, sheet, row.names=FALSE)
				}
				saveWorkbook(wb, filename)
			} else if(names(group)[i] == 'Columns'){
				data <- data.frame(Group = group[[i]], Sample = names(group[[i]]))
				data <- dlply(data, .(Group), function(x) sampleInfo(object, all=T)[x$Sample,])
				wb <- createWorkbook()
				filename <- paste(names(group)[i], '.xlsx', sep='')
				if(!missing(folder)){
					dir.create(file.path(getwd(), folder), showWarnings=FALSE)
					filename <- file.path(folder, filename)
				}
				for(j in names(data)){
					ans <- data[[j]]
					ans <- do.call('rbind', ans)
					sheet <- createSheet(wb, paste('Group ', j, sep=''))
					addDataFrame(ans, sheet, row.names=FALSE)
				}
				saveWorkbook(wb, filename)
			} else {}
		}
	}
)
### mixReport
### Creates a 2-page pdf with a simple PCA model
setMethod(
		'mixReport', 'Complist',
		function(object, outlier.rm=TRUE, unmatched=TRUE, FDR=FALSE, alpha=0.95){
			
			# Create model
			mod <- modelPCA(object, 'PCA', outlier.rm=outlier.rm, mix.rm=FALSE, unmatched=unmatched, FDR=FDR, nPcs=4)@models[['PCA']]$Model
			
			sco <- as.data.frame(mod@scores)[,1:4]
			
			info <- sampleInfo(object, outlier.rm=outlier.rm, mix.rm=FALSE)
			sco1 <- sco[!info$Master.mix, ]
			sco2 <- sco[info$Master.mix, ]
			sco1 <- data.frame(Sample=rep(rownames(info)[!info$Master.mix], 2), plot=rep(c('PC1 vs PC2', 'PC3 vs PC4'), each=nrow(sco1)), x=c(sco1[,1], sco1[,3]), y=c(sco1[,2], sco1[,4]))
			sco2 <- data.frame(Sample=rep(rownames(info)[info$Master.mix], 2), plot=rep(c('PC1 vs PC2', 'PC3 vs PC4'), each=nrow(sco2)), x=c(sco2[,1], sco2[,3]), y=c(sco2[,2], sco2[,4]))
			
			int <- pepmaps:::simconf(sco, alpha)
			angle <- seq(-pi, pi, length = 50)
			df <- data.frame(co1 = sin(angle)*int[1,2], co2 = cos(angle)*int[2,2], co3 = sin(angle)*int[3,2], co4 = cos(angle)*int[4,2])
			df <- data.frame(plot=rep(c('PC1 vs PC2', 'PC3 vs PC4'), each=nrow(df)), x=c(df[,1], df[,3]), y=c(df[,2], df[,4]))
			
			p <- ggplot(data=sco1, aes(x=x, y=y)) + theme_bw()
			p <- p + geom_point(colour=I('grey')) + geom_text(aes(label=Sample), data=sco2, size=3)
			p <- p + geom_path(data=df, colour=I('red'), linetype=I(2))
			p <- p + facet_wrap(~plot, ncol=2, scales='free') + xlab('') + ylab('')
			p
		}
)
### updateComplist
### Updates a Complist object created in an earlier version
updateComplist <- function(object){
	if(class(object) == 'Complist'){
		new(Class='Complist', object)
	} else {
		stop('Argument must be an old Complist object')
	}
}
### complist
### Constructor for complist objects
complist <- function(data, Sample.info, para, PepID, rtwin=0, mzwin=0){
    if(missing(data)){
        new(
            Class='Complist'
        )
    } else {
		
		# Input check
        if(missing(para)){
            para <- parameters()
        } else if(class(para) != 'Parameters'){
            stop('Apply parameters as a Parameters object')
        } else {}
        if(class(data)[[1]] == 'xsAnnotate'){
            pdata <- data@xcmsSet
            annotation <- conFeature(data)
        } else if(class(data)[[1]] == 'xcmsSet'){
            pdata <- data
            annotation <- list()
        } else {
            stop('Peaklist objects can only be created from xcmsSet or xsAnnotate objects')
        }
		
		# Extracts information from the xcmsSet object
        rawdata <- groupval(pdata, value='into')
        peak.info <- data.frame(pdata@groups)
        if(class(data)[[1]] == 'xsAnnotate'){
            pinfo <- getPeaklist(data)[,c('isotopes', 'adduct', 'pcgroup')]
            peak.info <- data.frame(peak.info, pinfo)
        } else {}
        sinfo <- data.frame(pdata@phenoData, I(pdata@filepaths))
        names(sinfo) <- c('Class', 'Filepath')
		
		# Combines supplied sample info with the one extracted from the xcmsSet
        if(!missing(Sample.info)){
            if(nrow(sinfo) != nrow(Sample.info)){
                stop('Submitted sample information does not correspond to the number of samples in xcms data')
            } else {
                sinfo <- data.frame(sinfo, Sample.info)
            }
        } else {}
		
		# In case no PepID object is provided
        if(missing(PepID)){
            pep.index <- list()
            PepID <- pepID()
			filter=list()
		
		# Creates link between ID and peakgroups
        } else {
            peakID <- rep(FALSE, nrow(peak.info))
            pep.index <- vector('list', nrow(PepID@peplist))
            for(i in 1:nrow(PepID@peplist)){
                ind <- which(peak.info$mzmin < (PepID@peplist$mz[i]+mzwin) & peak.info$mzmax > (PepID@peplist$mz[i]-mzwin) & peak.info$rtmin < (PepID@peplist$Retention.time[i]+rtwin) & peak.info$rtmax > (PepID@peplist$Retention.time[i]-rtwin))
                if(length(ind) > 0){
                    peakID[ind] <- TRUE
                    pep.index[[i]] <- ind
                } else {}
            }
			filter <- list(Peaks=data.frame(Outlier=FALSE, Matched=peakID), Peptides=data.frame(Outlier=FALSE, FDR=PepID@peplist$FDR < PepID@FDR, Matched=!sapply(pep.index, is.null)), Samples=data.frame(Outlier=FALSE, Mix=Sample.info$Master.mix))
        }
		
		# Creates object
        new(
            Class='Complist',
            raw=rawdata,
            Sample.info=sinfo,
            Peak.info=peak.info,
            annotation=annotation,
            IDindex=pep.index,
            pepID=PepID,
            parameters=para,
            xcmsSet=pdata,
            models=list(),
            filter=filter
        )
    }
}