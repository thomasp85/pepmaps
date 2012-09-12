# List of wrapper functions

### createComplist
### Wrapper function to create xcmsSet and PepID and combine the information into a Complist object
createComplist <- function(mzXML, Sample.info, ID='MSGF+', IDdir, database='milk_prot', cache=FALSE, par='standard', sep='\t', dec='.', annotate=FALSE, retcor=TRUE, useMGF=FALSE){
    first <- proc.time()
	
	# Input check
    if(missing(mzXML)){
		if(Sys.info()["sysname"] == 'Windows'){
			readline('\nChoose directory of mzXML data: <Press Return>')
			path <- choose.dir()
		} else {
			path <- readline('Path to directory with mzXML data:')
		}
        mzXML <- list.files(path, pattern='*.mzXML', full.names=TRUE, ignore.case=TRUE)
    } else {}
    cat('\n', dirname(mzXML[1]), '\n\n')
    for(i in 1:length(mzXML)){
        cat(paste(i, ': ', basename(mzXML[i]), '\n', sep=''))
    }
    cat('\n')
    flush.console()
    mix <- readline('Type the numbers of the samples that are master mixtures (space separated): ')
    if(mix != ''){
        mix <- as.numeric(unlist(strsplit(mix, ' ')))
        mmix <- rep(FALSE, length(mzXML))
        mmix[mix] <- TRUE
        if(missing(Sample.info)){
            Sample.info <- data.frame(Master.mix=mmix)
        } else {
            Sample.info$Master.mix <- mmix
        }
        cat('The following samples has been set as master mixtures:\n\n', paste(basename(mzXML[mix]), collapse='\n'), '\n\n', sep='')
    } else {
        mix <- c()
        cat('No master mixtures in experiment...\n')
    }
    flush.console()
	
	# Extracts parameters
    par <- parameters(par)
	
	# Creates PepID object
	cat('\nAnalysing peptide ID data...\n\n')
	flush.console()
	
	## Creates analysis from scratch using MSGF+
	if(ID == 'MSGF+'){
		if(cache){
			if(file.exists(R.home(component='library/pepmaps/msgfdb_res.cache.RData'))){
				ID <- local(get(load(R.home(component='library/pepmaps/msgfdb_res.cache.RData'))))
			} else {
				stop('No cache on drive...')
			}
			if(!all(unique(ID@raw$SpecFile) %in% mzXML)){
				cat('Warning: Samples in cache does not match samples in analysis\n')
				flush.console()
			} else {}
		} else {
			type <- ID
			para <- c(type=type, directory=dirname(mzXML[1]), database=database, useMGF=useMGF, par@parameters$MSGFplus)
			ID <- do.call('pepID', para)
			save(ID, file=R.home(component='library/pepmaps/msgfdb_res.cache.RData'))
		}
		
	## Load result file from MassAI/Crosslink
	} else {
	    if(missing(IDdir)){
			if(Sys.info()["sysname"] == 'Windows'){
				readline('\nChoose location of ID file(s): <Press Return>')
				path <- choose.dir()
			} else {
				path <- readline('Path to directory with ID file(s):')
			}
	        massai <- list.files(path, pattern='*.txt', full.names=TRUE, recursive=TRUE, ignore.case=TRUE)
	    } else {
	        massai <- list.files(IDdir, pattern='*.txt', full.names=TRUE, recursive=TRUE, ignore.case=TRUE)
	    }
	    type <- ID
	    ID <- pepID(type, path=massai[1], sep=sep, dec=dec, ...)
		if(length(massai) > 1){
			for(i in 2:length(massai)){
				ID2 <- pepID(type, path=massai[i], sep=sep, dec=dec, ...)
				ID <- mergePepID(ID, ID2)
			}
		} else {}
		par@parameters$MSGFplus <- NULL
	}
    if(par@Type == 'LC-MS/MS'){
		
		# Run xcms pipeline to create xcmsSet
        cat('Running xcms analysis...\n\n')
		cat('Detecting peaks...\n')
        flush.console()
        xset <- do.call('xcmsSet', c(list(files=mzXML), par@parameters$findPeak))
        if(retcor){
			cat('\nCorrecting shifts in retention time...\n\n')
			flush.console()
            if(length(mix) != 0){
                par@parameters$retcor$center <- mix[round(median(1:length(mix)))]
            } else {}
            xset <- do.call('retcor', c(list(object=xset), par@parameters$retcor))
        } else {}
		cat('\nGrouping peaks between samples...\n\n')
		flush.console()
        xset <- do.call('group', c(list(object=xset), par@parameters$group))
		cat('\nDetecting signal in samples missing from peak groups...\n\n')
		flush.console()
        xset <- fillPeaks(xset)
        if(annotate){
            cat('\nRunning annotation with CAMERA...\n\n')
            flush.console()
            gc(verbose=FALSE)
            xset <- do.call(function(...) CAMERA::annotate(...), c(list(object=xset), par@parameters$annotate))
        } else {}
        cat('\nCreating Complist...\n\n')
        flush.console()
		
		# Creating complist object
        if(missing(Sample.info)){
            comp <- complist(xset, para=par, PepID=ID, rtwin=par@parameters$Complist$rtwin, mzwin=par@parameters$Complist$mzwin)
        } else {
            comp <- complist(xset, Sample.info=Sample.info, para=par, PepID=ID, rtwin=par@parameters$Complist$rtwin, mzwin=par@parameters$Complist$mzwin)
        }
    }
    last <- proc.time()
    cat('Done!\n\tProcess lasted ', (last-first)[3], ' sec. (', (last-first)[3]/60, ' min.)\n', sep='')
    comp
}
### pepReport
### Report function for Complist
pepReport <- function(object, group=NULL, folder=NULL, outlier.rm=TRUE, mix.rm=TRUE){
	
	# Input check
    if(class(object) != 'Complist'){
        stop('Wrong object type. Provide a Complist object...\n')
    } else {}
    if(sum(names(object@peakID) %in% 'Peptides') != 1){
        stop('Incorect peptide data...\n')
    } else {}
    if(!missing(folder)){
        if(!is.character(folder) | length(folder) != 1){
            stop('Wrong dir input. Provide directory name...')
        } else {
            dir.create(file.path(getwd(), folder), showWarnings=FALSE)
        }
    } else {}
	
	# Create reports
    intensityReport(object, folder=folder, outlier.rm=outlier.rm, mix.rm=mix.rm)
	coverageReport(object, folder=folder, outlier.rm=outlier.rm, mix.rm=mix.rm)
    chromReport(object, folder=folder, outlier.rm=outlier.rm, mix.rm=mix.rm)
    pcaReport(object, folder=folder, outlier.rm=outlier.rm, mix.rm=mix.rm, group=group)
	gc()
	peakReport(object, folder=folder, outlier.rm=outlier.rm, mix.rm=mix.rm, group=group)
	
	# Feedback
    cat('\nDone...\n')
    if(outlier.rm){
        out <- row.names(object@Sample.info)[object@outlier$Sample]
        if(length(out) > 0){
            cat('\nOutlier samples removed:', paste(out, collapse=', ', '\n', sep=' '))
        } else {}
    } else {}
}