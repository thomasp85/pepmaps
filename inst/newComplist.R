# TODO: Add comment
# 
# Author: Thomas
###############################################################################
setGeneric(
		'combine',
		def=function(...){standardGeneric('combine')}
)
setGeneric(
		'consensus',
		def=function(object, ...){standardGeneric('consensus')}
)
### Class that merges PepID and Peaklist (and possibly additional *ID classes).
setClass(
		Class='Complist',
		representation=representation(
				raw='FeatureSet',
				chromatograms='list',
				Sample.info='data.frame',
				annotation='list',
				parameters='Parameters',
				peakID='data.frame',
				IDindex='list',
				pepID='PepID',
				xcmsSet='character',
				outlier='list',
				models='list',
				filter='list'
		),
		prototype=prototype(
				raw=featureSet(),
				chromatograms=list(TIC=NULL, BPC=NULL),
				Sample.info=data.frame(),
				annotation=list(),
				parameters=parameters(),
				peakID=data.frame(),
				IDindex=list(),
				pepID=pepID(),
				xcmsSet=character(),
				outlier=list(),
				models=list(),
				filter=list()
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
setClass(
		Class='FeatureSet',
		representation=representation(
				data='list',
				samples='vector',
				peak.info='matrix'
		),
		prototype=prototype(
				data=list(into=matrix(), maxo=matrix()),
				samples=vector(),
				peak.info=matrix()
		),
		validity=function(object){
			if(nrow(object@data$into) != nrow(object@data$maxo)){
				stop('Dimensions does not match between the two raw matrices')
			} else if(ncol(object@data$into) != ncol(object@data$maxo)){
				stop('Dimensions does not match between the two raw matrices')
			} else if(length(object@data) != 2){
				stop('Wrong raw data format. The list must have length == 2')
			} else if(all(names(object@data) %in% c('into', 'maxo'))){
				stop('Invalid names for the raw data. Must be \'into\' and \'maxo\'')
			} else if(length(object@samples) != ncol(object@data$into)){
				stop('Number of samples does not match dimensions of raw data')
			} else if(nrow(object@peak.info) != nrow(object@data$into)){
				stop('Dimensions of peak.info does not match raw data')
			} else {
				TRUE
			}
		}
)
setMethod(
	'length', 'FeatureSet',
	function(x){
		if(length(x@samples) == 0){
			0
		} else {
			length(x@samples)
		}
	}
)
setMethod(
	'show', 'FeatureSet',
	function(object){
		if(length(object) == 0){
			cat('An empty FeatureSet object\n')
		} else {
			cat('A FeatureSet object with ', length(object), ' samples\n', sep='')
			cat('\n')
			cat('The set contains:\t', nrow(object@peak.info), ' peak groups\n', sep='')
			cat('\tin the range:\t', min(object@peak.info[, 'rtmin']), ' - ', max(object@peak.info[, 'rtmax']), ' sec\n', sep='')
			cat('\t\t\t and:\t', min(object@peak.info[, 'mzmin']), ' - ', max(object@peak.info[, 'mzmax']), ' m/z\n', sep='')
			cat('\n')
		}
	}
)
setMethod(
	'dim', 'FeatureSet',
	function(x){
		dim(x@data$into)
	}
)
setMethod(
	'as.matrix', 'FeatureSet',
	function(x, type='into'){
		if(!all(type %in% c('into', 'maxo'))){
			stop('Wrong type argument')
		} else {}
		x@data[[type]]
	}
)
setMethod(
		'as.data.frame', 'FeatureSet',
		function(x, type='into'){
			if(!all(type %in% c('into', 'maxo'))){
				stop('Wrong type argument')
			} else {}
			data.frame(x@data[[type]], check.names=FALSE, row.names=1:nrow(x@data[[type]]))
		}
)
featureSet <- function(object){
	if(missing(object)){
		new(Class='FeatureSet')
	} else {
		if(class(object) == 'xcmsSet'){
			data <- list(into=groupval(object, value='into'), maxo=groupval(object, value='maxo'))
			samples <- rownames(object@phenoData)
			peak.info <- object@groups[, colnames(object@groups) %in% c('mzmed', 'mzmin', 'mzmax', 'rtmed', 'rtmin', 'rtmax')]
		} else if(class(object) == 'xsAnnotate'){
			data <- list(into=groupval(object@xcmsSet, value='into'), maxo=groupval(object@xcmsSet, value='maxo'))
			samples <- rownames(object@xcmsSet@phenoData)
			peak.info <- object@xcmsSet@groups[, colnames(object@xcmsSet@groups) %in% c('mzmed', 'mzmin', 'mzmax', 'rtmed', 'rtmin', 'rtmax')]
		} else {
			stop('Can only create FeatureSet from xcmsSet or xsAnnotate objects')
		}
		new(
			Class='FeatureSet',
			data=data,
			samples=samples,
			peak.info=peak.info
		)
	}
}
setClass(
		Class='Spectrum',
		representation=representation(
				mz='vector',
				int='vector'
		),
		prototype=prototype(
				mz=numeric(),
				int=numeric()
		),
		validity=function(object){
			if(!any(class(object@mz) %in% c('integer', 'numeric'))){
				stop('mz values must be numeric')
			} else if(!any(class(object@int) %in% c('integer', 'numeric'))){
				stop('Intensity values must be numeric')
			} else if(length(object@mz) != length(object@int)){
				stop('Intensity and mz values must have same length')
			} else {
				TRUE
			}
		}
)
setMethod(
		'initialize', 'Spectrum',
		function(.Object, mz, int){
			if(length(mz) != 0){
				names(mz) <- NULL
				names(int) <- NULL
				.Object@mz <- mz
				.Object@int <- int
			} else {}
			validObject(.Object)
			return(.Object)
		}
)
setMethod(
	'as.matrix', 'Spectrum',
	function(x){
		matrix(c(x@mz, x@int), ncol=2, dimnames=list(NULL, c('mz', 'intensity')))
	}
)
setMethod(
		'as.data.frame', 'Spectrum',
		function(x){
			data.frame(mz=x@mz, intensity=x@int)
		}
)
setMethod(
		'as.list', 'Spectrum',
		function(x){
			list(mz=x@mz, intensity=x@int)
		}
)
setMethod(
	'length', 'Spectrum',
	function(x){
		length(x@mz)
	}
)
setMethod(
	'show', 'Spectrum',
	function(object){
		if(length(object) != 0){
			cat('A Spectrum object\n')
			cat('\n')
			cat('mz range:\n')
			cat('\t', min(object@mz), ' - ', max(object@mz), '\n', sep='')
			cat('\n')
			cat('Intensity:\n')
			print(summary(object@int))
		} else {
			cat('An empty Spectrum object\n')
		}
	}
)
setMethod(
	'plot', 'Spectrum',
	function(x){
		x <- as.data.frame(x)
		qplot(x=mz, xend=mz, y=0, yend=intensity, data=x, geom='segment') + xlab('m/z') + ylab('intensity') + theme_bw()
	}
)
setMethod(
	'consensus', 'Spectrum',
	function(object, ...){
		spectra <- list(object, ...)
		if(length(spectra) == 1){
			spectra[[1]]
		} else {
			nSpectra <- length(spectra)
			specid <- as.list(rep(1:nSpectra, sapply(spectra, length)))
			spectra <- do.call('rbind', lapply(spectra, as.matrix))
			
			tollevels <- c(0.01, 0.02, 0.03, 0.04)
			
			for(i in tollevels){
				specid <- specid[order(spectra[, 'mz'])]
				spectra <- spectra[order(spectra[, 'mz']), ]
				groups <- diff(diff(spectra[,'mz']) < i)
				if(all(groups == 0)){
					next
				} else {}
				if(rev(groups[groups!=0])[1] == 1){
					groups <- c(groups, -1)
				} else {}
				groupInt <- which(groups != 0) + 1
				if(groups[groups!=0][1] == -1){
					groupInt <- c(1, groupInt)
				} else {}
				groups <- lapply(split(groupInt, rep(1:(length(groupInt)/2), each=2)), function(x) x[1]:x[2])
				spectraTMP1 <- spectra[-unlist(groups), ]
				specidTMP1 <- specid[-unlist(groups)]
				spectraTMP2 <- lapply(groups, function(x, spectra) c(weighted.mean(spectra[x, 'mz'], spectra[x, 'intensity']), sum(spectra[x, 'intensity'])), spectra=spectra)
				specidTMP2 <- lapply(split(specid[unlist(groups)], rep(1:length(groups), sapply(groups, length))), unlist)
				spectra <- rbind(spectraTMP1, do.call('rbind', spectraTMP2))
				specid <- c(specidTMP1, specidTMP2)
			}
			
			specid <- specid[order(spectra[, 'mz'])]
			spectra <- spectra[order(spectra[, 'mz']), ]
			
			scale <- 0.95 + 0.05*(1 + sapply(specid, function(x, nSpectra) length(unique(x))/nSpectra, nSpectra=nSpectra))^5
			spectra[, 'intensity'] <- spectra[, 'intensity']*scale
			
			keep <- rep(TRUE, nrow(spectra))
			for(i in 1:nrow(spectra)){
				mz <- spectra[i, 'mz']
				ind <- which(spectra[, 'mz'] > mz-50 & spectra[, 'mz'] < mz+50)
				cutoff <- sort(spectra[ind, 'intensity'], decreasing=T)[5]
				if(is.na(cutoff)){
					next
				} else {
					remove <- ind[which(spectra[ind, 'intensity'] < cutoff)]
					keep[remove] <- FALSE
				}
			}
			spectra <- spectra[keep, ]
			spectrum(mz=spectra[, 'mz'], int=spectra[, 'intensity'])
		}
	}
)
setMethod(
	'combine', 'Spectrum',
	function(...){
		ans <- do.call('rbind', lapply(list(...), as.matrix))
		ans <- ans[order(ans[, 'mz']), ]
		spectrum(int=ans[, 'intensity'], mz=ans[, 'mz'])
	}
)
spectrum <- function(int=numeric(), mz=numeric()){
	new(
		Class='Spectrum',
		mz=mz,
		int=int
	)
}
setClass(
		Class='MSnSet',
		representation=representation(
				consensus='Spectrum',
				raw='list',
				similarity='numeric',
				samples='vector',
				acquisitionNum='vector',
				charge='vector',
				precursorMz='vector',
				rt='vector'
		),
		prototype=prototype(
				consensus=spectrum(),
				raw=list(),
				similarity=numeric(),
				samples=character(),
				acquisitionNum=numeric(),
				charge=numeric(),
				precursorMz=numeric(),
				rt=numeric()
		),
		validity=function(object){
			if(length(object) != 0){
				if(length(unique(c(length(object@raw), length(object@samples), length(object@acquisitionNum), length(object@charge), length(object@precursorMz), length(object@rt)))) != 1){
					stop('Length of elements must match')
				} else if(!all(sapply(object@raw, class) == 'Spectrum')){
					stop('Raw data must be a list of Spectrum objects')
#				} else if(!any(sapply(object@similarity, class) %in% c('integer', 'numeric'))){
#					stop('Similarity must be numeric')
				} else if(!any(class(object@acquisitionNum) %in% c('integer', 'numeric'))){
					stop('Acquisition number must be numeric')
				} else if(!any(class(object@charge) %in% c('integer', 'numeric'))){
					stop('Charge must be numeric')
				} else if(!any(class(object@precursorMz) %in% c('integer', 'numeric'))){
					stop('Precursor mz must be numeric')
				} else if(!any(class(object@rt) %in% c('integer', 'numeric'))){
					stop('Retention time must be numeric')
				} else {
					TRUE
				}
			} else {
				TRUE
			}
		}
)
setMethod(
		'initialize', 'MSnSet',
		function(.Object, raw, samples, acquisitionNum, charge, precursorMz, rt){
			if(length(raw) != 0){
#				.Object@consensus <- do.call('consensus', raw)
#				.Object@similarity <- sapply(raw, function(x, consensus) similarity(x, consensus), consensus=.Object@consensus)
				.Object@raw <- raw
				.Object@samples <- samples
				.Object@acquisitionNum <- acquisitionNum
				.Object@charge <- charge
				.Object@precursorMz <- precursorMz
				.Object@rt <- rt
			} else {}
			validObject(.Object)
			return(.Object)
		}
)
setMethod(
	'show', 'MSnSet',
	function(object){
		if(length(object) == 0){
			cat('An empty MSnSet object\n')
		} else {
			cat('An MSnSet object containing ', length(object), ' samples\n', sep='')
			cat('\n')
			cat('Mean Precursor mz:\t\t', mean(object@precursorMz), ' m/z\n', sep='')
			cat('Mean Retention time:\t', mean(object@rt), ' sec\n', sep='')
			cat('Median charge:\t\t\t', median(object@charge), '\n', sep='')
		}
	}
)
setMethod(
	'length', 'MSnSet',
	function(x){
		if(length(x@raw) == 0){
			0
		} else {
			length(x@raw)
		}
	}
)
setMethod(
	'plot', 'MSnSet',
	function(x){
		plot(x@consensus)
	}
)
setMethod(
		'combine', 'MSnSet',
		function(...){
			spectra <- list(...)[sapply(list(...), length) != 0]
			raw <- do.call('c', lapply(spectra, function(x) x@raw))
			samples <- do.call('c', lapply(spectra, function(x) x@samples))
			acquisitionNum <- do.call('c', lapply(spectra, function(x) x@acquisitionNum))
			charge <- do.call('c', lapply(spectra, function(x) x@charge))
			precursorMz <- do.call('c', lapply(spectra, function(x) x@precursorMz))
			rt <- do.call('c', lapply(spectra, function(x) x@rt))
			new(
				Class='MSnSet',
				raw=raw,
				samples=samples,
				acquisitionNum=acquisitionNum,
				charge=charge,
				precursorMz=precursorMz,
				rt=rt
			)
		}
)
setMethod(
	'consensus', 'MSnSet',
	function(object){
		if(length(object@raw) != 0){
			object@consensus <- do.call('consensus', object@raw)
		}
		object
	}
)
similarity <- function(s1, s2){
	if(class(s1) != 'Spectrum' | class(s2) != 'Spectrum'){
		stop('Supplied data must be Spectrum objects')
	} else {}
	k <- round(15 * max(s1@mz, s2@mz)/1000)
	s1 <- as.matrix(s1)
	s1 <- s1[order(s1[,'intensity'], decreasing=T)[1:k],]
	s2 <- as.matrix(s2)
	s2 <- s2[order(s2[,'intensity'], decreasing=T)[1:k],]
	cMZ <- c(s1[, 'mz'], s2[, 'mz'])
	0
}
createMSnSet <- function(xcmsRaw, msnIndex, specList){
	if(missing(specList) & !missing(xcmsRaw)){
		if(class(xcmsRaw) != 'xcmsRaw'){
			stop('xcmsRaw must be an xcmsRaw object')
		} else if(!any(class(msnIndex) %in% c('integer', 'numeric'))){
			stop('Index must be numeric')
		} else {}
		if(length(msnIndex) == 1){
			raw <- list(getMSnScan(object=xcmsRaw, scan=msnIndex))
		} else {
			raw <- lapply(msnIndex, function(x, object) getMSnScan(object=xcmsRaw, scan=x), object=xcmsRaw)
		}
		raw <- lapply(raw, function(x) spectrum(mz=x[,'mz'], int=x[,'intensity']))
		samples <- rep(xcmsRaw@filepath, length(msnIndex))
		acquisitionNum <- xcmsRaw@msnAcquisitionNum[msnIndex]
		charge <- xcmsRaw@msnPrecursorCharge[msnIndex]
		precursorMz <- xcmsRaw@msnPrecursorMz[msnIndex]
		rt <- xcmsRaw@msnRt[msnIndex]
	} else if(missing(xcmsRaw)){
		raw=list()
		samples=character()
		acquisitionNum=numeric()
		charge=numeric()
		precursorMz=numeric()
		rt=numeric()
	} else {
		if(!all(sapply(specList, class) == 'Spectrum')){
			stop('spectra must be provided as a list of Spectrum objects')
		} else {}
		raw <- specList
		samples <- rep(NA, length(msnIndex))
		acquisitionNum <- rep(NA, length(msnIndex))
		charge <- rep(NA, length(msnIndex))
		precursorMz <- rep(NA, length(msnIndex))
		rt <- rep(NA, length(msnIndex))
	}
	new(
		Class='MSnSet',
		raw=raw,
		samples=samples,
		acquisitionNum=acquisitionNum,
		charge=charge,
		precursorMz=precursorMz,
		rt=rt
	)
}
getMSnSet <- function(object){
	msnSets <- lapply(seq_along(object@groupidx), function(x) createMSnSet())
	for(i in unique(object@peaks[, 'sample'])){
		filename <- object@filepaths[i]
		cat(basename(filename), ': Reading file, ', sep='')
		flush.console()
		raw <- xcmsRaw(filename, includeMSn=T)
		cat('extracting spectra, ', sep='')
		flush.console()
		ind <- mapply(function(rtmin, rtmax, mzmin, mzmax, rt, mz) which(rt>rtmin & rt<rtmax & mz>round(mzmin, digits=2)-0.01 & mz<round(mzmax, digits=2)+0.01), rtmin=object@peaks[object@peaks[, 'sample'] == i,'rtmin'], rtmax=object@peaks[object@peaks[, 'sample'] == i,'rtmax'], mzmin=object@peaks[object@peaks[, 'sample'] == i,'mzmin'], mzmax=object@peaks[object@peaks[, 'sample'] == i,'mzmax'], MoreArgs=list(rt=raw@msnRt, mz=raw@msnPrecursorMz))
		group <- rep(seq_along(object@groupidx), sapply(object@groupidx, length))[match(which(object@peaks[, 'sample'] == i), unlist(object@groupidx))]
		real <- which(sapply(ind, length) > 0 & !is.na(group))
		ind <- ind[real]
		group <- group[real]
		ind <- split(unlist(ind), rep(match(group, unique(group)), sapply(ind, length)))
		group <- unique(group)
		tmpSet <- lapply(ind, function(x, raw) createMSnSet(raw, x), raw=raw)
		cat('combining spectra', sep='')
		flush.console()
		msnSets[group] <- mapply(combine, msnSets[group], tmpSet)
		cat('\n')
	}
	cat('\n')
	msnSets
}