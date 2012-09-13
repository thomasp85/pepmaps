### Class to handle peptide ID results from MassAI or Crosswork
### Contains both raw and truncated data
### Type must be either 'MassAI' or 'Crosswork'
setClass(
    Class = 'PepID',
    representation = representation(
        type='character',
        raw='data.frame',
        peplist='data.frame',
		database='AAStringSet',
        npep='numeric',
        nsample='numeric'
    ),
    validity = function(object){
        if(length(object@type) != 0){
            if(object@type %in% c('MassAI', 'Crosswork', 'MSGF+') == FALSE){
                stop('Invalid type argument. Only MSGF+, MassAI and Crosswork supported')
            } else {}
        } else {}
        return(TRUE)
    }
)
### Initialize for PepID
### Handles creation of truncated list as well as the nsample and npep arguments
### Calls validObject
setMethod(
    'initialize', 'PepID',
    function(.Object, type, raw, database){
        if(length(type) == 0){
            .Object@type <- character()
            .Object@raw <- data.frame()
            .Object@peplist <- data.frame()
			.Object@database <- Biostrings:::AAStringSet()
            .Object@npep <- numeric()
            .Object@nsample <- numeric()
        } else {
            .Object@type <- type
            .Object@raw <- raw
			.Object@database <- database
            if(type == 'MassAI'){
                .Object@nsample <- length(unique(raw$File))
                .Object@npep <- length(unique(paste(raw$ProteinAID, raw$PeptideAID, sep='-')))
                peplist <- ddply(raw, c(.(ProteinAID), .(PeptideAID)), function(x) data.frame(paste(x$Protein.A[1], ' (', x$First.residue[1], '-', x$Last.residue[1], ')', sep=''), paste(x$Peptide.A[1]), mean(x$Precursor), mean(x$Charge), mean(x$Parent.mass), mean(x$Mass.error), x$Protein.A[1], if(any(!is.na(x$Score))) max(x$Score, na.rm=TRUE) else NA, mean(x$Retention), if('Lev2' %in% x$Notes){'Lev2'}else{'Lev1'}))
				peplist <- peplist[, -1]
                names(peplist) <- c('Peptide.ID', 'Peptide.name', 'Sequence', 'Retention.time', 'mz', 'Charge', 'Mass', 'Mass.error', 'Protein', 'Score', 'Notes')
                peplist$Peptide.name <- sub('>', '', peplist$Peptide.name)
                peplist$Protein <- sub('>', '', peplist$Protein)
                peplist$Notes <- as.character(peplist$Notes)
                peplist$Peptide.ID <- paste(peplist$Protein.ID, peplist$Peptide.ID, sep='-')
                .Object@peplist <- peplist
            } else if(type == 'Crosswork'){
                .Object@nsample <- length(unique(raw$Run))
                pepname <- paste(raw$Protein, ' (', raw$Residue, '-', raw$Residue+nchar(as.character(raw$Sequence))-1, ')', sep='')
                .Object@npep <- length(unique(pepname))
                peplist <- data.frame(pepname, raw)
                peplist <- ddply(peplist, .(pepname), function(x) data.frame(x$Sequence[1], mean(x$retention), mean(x$Precursor), x$Charge[1], x$Precursor[1]*x$Charge[1]-x$Charge[1], mean(x$mass.error), x$Protein[1], max(x$Matches)))
                peplist <- data.frame(1:nrow(peplist), peplist)
                names(peplist) <- c('Peptide.ID', 'Peptide.name', 'Sequence', 'Retention.time', 'mz', 'Charge', 'Mass', 'Mass.error', 'Protein', 'Notes')
                peplist$Notes <- as.character(peplist$Notes)
                .Object@peplist <- peplist
            } else if(type == 'MSGF+'){
				.Object@nsample <- length(unique(raw$SpecFile))
				.Object@npep <- length(unique(raw$Peptide.ID))
				cPep <- function(x){
					seq <- strsplit(as.character(x$Peptide[1]), '.', fixed=T)[[1]][2]
					ind <- gregexpr(seq, as.character(database[names(database) == x$Protein[1]]))[[1]]
					name <- paste(x$Protein[1], ' (', paste(ind, collapse='/'), '-', paste(ind + attr(ind, 'match.length') - 1, collapse='/'), ')', sep='')
					rt <- mean(x$rt)
					mz <- mean(x$Precursor)
					charge <- mean(x$Charge)
					mass <- pepMass(seq, mono=TRUE)
					m.err <- mean(x$'PrecursorError(Da)')
					prot <- x$Protein[1]
					fdr <- min(x$QValue)
					data.frame(name, seq, rt, mz, charge, mass, m.err, prot, fdr)
				}
				peplist <- ddply(raw, .(Peptide.ID), cPep)
				names(peplist) <- c('Peptide.ID', 'Peptide.name', 'Sequence', 'Retention.time', 'mz', 'Charge', 'Mass', 'Mass.error', 'Protein', 'FDR')
				.Object@peplist <- peplist
			} else {}
            .Object@peplist <- data.frame(.Object@peplist, Qval(as.character(.Object@peplist$Sequence)), Length=nchar(as.character(.Object@peplist$Sequence)))
        }
        validObject(.Object)
        return(.Object)
    }
)
### Show for PepID
### Very short summary of the PepID object
setMethod(
    'show', 'PepID',
    function(object){
        if(length(object) == 0){
            cat('An empty PepID object.\n')
        } else {
            cat('A PepID object created from a', object@type, 'analysis.\n\n')
            cat(object@npep, 'different peptides detected in', object@nsample, 'samples.\n')
        }
    }
)
### length for PepID
### Very short summary of the PepID object
setMethod(
    'length', 'PepID',
    function(x){
        if(length(x@type)+length(x@raw)+length(x@peplist)+length(x@npep)+length(x@nsample) == 0){
            ans <- 0
        } else {
            ans <- 1
        }
        ans
    }
)
### getPeplist
setMethod(
    'getPeplist', 'PepID',
    function(object){
        ans <- object@peplist
        ans
    }
)
### getRawID
setMethod(
    'getRawID', 'PepID',
    function(object){
        ans <- object@raw
        ans
    }
)
### mergePepID method for PepID
setMethod(
    'mergePepID', signature=c(object1='PepID', object2='PepID'),
    function(object1, object2){
        if(object1@type != object2@type){
            stop('Cannot merge PepIDs from different analyses\n')
        } else {}
        prob <- names(which(table(rbind(object1@peplist, object2@peplist)$Peptide.name) > 1))
        ans1 <- object1@peplist[object1@peplist$Peptide.name %in% prob,c('Peptide.name','Protein.ID', 'Peptide.ID')]
        ans2 <- object2@peplist[object2@peplist$Peptide.name %in% prob,c('Peptide.name','Protein.ID', 'Peptide.ID')]
        names(ans1) <- c('Name', 'ProA', 'PepA')
        names(ans2) <- c('Name', 'ProB', 'PepB')
        ans <- merge(ans1, ans2)
        mx <- ddply(rbind(object1@raw, object2@raw), .(ProteinID), function(x) max(x$PeptideID))
        addNewID <- function(data, mx){
            mx <- mx[which(mx$ProteinID == data$ProA[1]), 2]
            newID <- (mx+1):(mx+nrow(data))
            ans <- data.frame(data, newPep=newID)
            ans
        }
        ans <- ddply(ans, .(ProA), function(x, mx) addNewID(x, mx), mx=mx)
        raw1 <- object1@raw
        raw2 <- object2@raw
        for(i in 1:nrow(ans)){
            pID <- ans$newPep[i]
            raw1$PeptideID[which(raw1$PeptideID == ans$PepA[i] & raw1$ProteinID == ans$ProA[i])] <- pID
            raw2$PeptideID[which(raw2$PeptideID == ans$PepB[i] & raw2$ProteinID == ans$ProB[i])] <- pID
        }
        mx <- ddply(rbind(raw1, raw2), .(ProteinID), function(x) max(x$PeptideID))
        addNewID <- function(data, mx){
            mx <- mx[which(mx$ProteinID == data$Protein.ID[1]), 2]
            newID <- (mx+1):(mx+nrow(data))
            ans <- data.frame(data, newPep=newID)
            ans
        }
        findDup <- function(x){
            if(nrow(x) > 1){
                if(length(unique(x$Peptide.name)) > 1 ){
                    x[1, c('Protein.ID','Peptide.ID')]
                } else {}
            } else {}
        }
        ans <- ddply(rbind(object1@peplist, object2@peplist), c(.(Protein.ID), .(Peptide.ID)), findDup)
        ans <- ddply(ans, .(Protein.ID), function(x, mx) addNewID(x, mx), mx=mx)
        for(i in 1:nrow(ans)){
            pID <- ans$newPep[i]
            raw2$PeptideID[which(raw2$PeptideID == ans$Peptide.ID[i] & raw2$ProteinID == ans$Protein.ID[i])] <- pID
        }
        raw <- rbind(raw1, raw2)
        new(Class='PepID', type=object1@type, raw=raw)
    }
)
### evalID for PepID
### Creates bootstrap estimation of number of identified peptides as a function of number of samples
setMethod(
    'evalID', 'PepID',
    function(object){
        if(object@type == 'MassAI'){
            data <- object@raw
            set1 <- dlply(data, .(File), function(x) x$Peptide)
            set2 <- dlply(data[which(data$Notes == 'Lev2'),], .(File), function(x) x$Peptide)
            boot1 <- matrix(NA, ncol=length(set1), nrow=1000)
            cat('Sampling... Lev1+2')
            flush.console()
            for(i in 1:ncol(boot1)){
                for(j in 1:nrow(boot1)){
                    index <- sample(1:ncol(boot1), i)
                    boot1[j,i] <- length(unique(unlist(set1[index])))
                }
            }
            cat(' Done.')
            flush.console()
            boot1 <- data.frame(boot1)
            names(boot1) <- 1:ncol(boot1)
            boot1 <- data.frame(boot.no=1:nrow(boot1), boot1)
            boot1 <- melt(boot1, id=1, variable_name='sample.no')
            boot1$sample.no <- sub('X', '', boot1$sample.no)
            boot1$sample.no <- as.numeric(boot1$sample.no)
            boot1 <- data.frame(boot1, quality='Mass and fragment')
            boot2 <- matrix(NA, ncol=length(set2), nrow=1000)
            cat(' Lev2')
            flush.console()
            for(i in 1:ncol(boot2)){
                for(j in 1:nrow(boot2)){
                    index <- sample(1:ncol(boot2), i)
                    boot2[j,i] <- length(unique(unlist(set2[index])))
                }
            }
            cat(' Done.\n')
            flush.console()
            boot2 <- data.frame(boot2)
            names(boot2) <- 1:ncol(boot2)
            boot2 <- data.frame(boot.no=1:nrow(boot2), boot2)
            boot2 <- melt(boot2, id=1, variable_name='sample.no')
            boot2$sample.no <- sub('X', '', boot2$sample.no)
            boot2$sample.no <- as.numeric(boot2$sample.no)
            boot2 <- data.frame(boot2, quality='Fragment')
            boot <- rbind(boot1, boot2)
        } else if(object@type == 'Crosswork'){
            data <- object@raw
            set1 <- dlply(data, .(Run), function(x) x$Sequence)
            boot1 <- matrix(NA, ncol=length(set1), nrow=1000)
            for(i in 1:ncol(boot1)){
                for(j in 1:nrow(boot1)){
                    index <- sample(1:ncol(boot1), i)
                    boot1[j,i] <- length(unique(unlist(set1[index])))
                }
            }
            boot1 <- data.frame(boot1)
            names(boot1) <- 1:ncol(boot1)
            boot1 <- data.frame(boot.no=1:nrow(boot1), boot1)
            boot1 <- melt(boot1, id=1, variable_name='sample.no')
            boot1$sample.no <- sub('X', '', boot1$sample.no)
            boot1$sample.no <- as.numeric(boot1$sample.no)
            boot <- data.frame(boot1, quality='Mass and fragment')
        } else {stop('Only MassAI and Crosswork analysis supported')}
        p <- ggplot(boot, aes(x=sample.no, y=value)) + facet_grid(quality~., scales='free')
        p <- p + geom_point(alpha=I(0.2))
        if(object@nsample>10){
            p <- p + geom_smooth()
        } else {}
        p <- p + xlab('Number Of Samples') + ylab('Number Of Identified Peptides') + theme_bw()
        p
    }
)

### Run MSGFDB through R
MSGFplus <- function(file, database, tolerance, tda=TRUE, instrument, protease, lengthRange, chargeRange, verbose=FALSE) {
	if(missing(file)){
		cat('Choose spectrum file...\n')
		flush.console()
		file <- file.choose()
	} else {}
	if(Sys.info()["sysname"] == 'Windows'){
		file <- paste('\"', file, '\"', sep='')
	} else {
		file <- gsub(' ', '\\ ', file, fixed=T)
	}
	call <- paste('-s ', file, sep='')
	
	if(missing(database)){
		cat('Choose database file...\n')
		flush.console()
		database <- file.choose()
	} else {}
	if(Sys.info()["sysname"] == 'Windows'){
		database <- paste('\"', database, '\"', sep='')
	} else {
		database <- gsub(' ', '\\ ', database, fixed=T)
	}
	call <- paste(call, ' -d ', database, sep='')
	
	if(missing(tolerance)){
		stop('Supply a mass tolerance...')
	} else {}
	call <- paste(call, ' -t ', tolerance, sep='')
	
	tmp <- R.home(component='library/pepmaps/msgfplus.cache.mzid')
	call <- paste(call, ' -o ', tmp, sep='')
	
	if(tda){
		call <- paste(call, ' -tda 1', sep='')
	} else {
		call <- paste(call, ' -tda 0', sep='')
	}
	
	if(!missing(instrument)){
		if(instrument == 'TOF'){
			call <- paste(call, ' -inst 2', sep='')
		} else if(instrument == 'LowLTQ'){
			call <- paste(call, ' -inst 0', sep='')
		} else if(instrument == 'HighLTQ'){
			call <- paste(call, ' -inst 1', sep='')
		} else {}
	} else {}
	
	if(!missing(protease)){
		if(protease == 'None'){
			call <- paste(call, ' -e 0', sep='')
		} else if(protease == 'Trypsin'){
			call <- paste(call, ' -e 1', sep='')
		} else if(protease == 'Chymotrypsin'){
			call <- paste(call, ' -e 2', sep='')
		} else if(protease == 'Lys-C'){
			call <- paste(call, ' -e 3', sep='')
		} else if(protease == 'Lys-N'){
			call <- paste(call, ' -e 4', sep='')
		} else if(protease == 'Glu-C'){
			call <- paste(call, ' -e 5', sep='')
		} else if(protease == 'Arg-C'){
			call <- paste(call, ' -e 6', sep='')
		} else if(protease == 'Asp-N'){
			call <- paste(call, ' -e 7', sep='')
		} else if(protease == 'alphaLP'){
			call <- paste(call, ' -e 8', sep='')
		} else if(protease == 'Unknown'){
			call <- paste(call, ' -e 9', sep='')
		} else {}
	} else {}
	
	if(!missing(lengthRange)){
		if(length(lengthRange) != 2){
			stop('Please provide a vector of length 2 for the lengthRange...')
		} else {}
		call <- paste(call, ' -minLength ', lengthRange[1], ' -maxLength ', lengthRange[2], sep='')
	} else {}
	
	if(!missing(chargeRange)){
		if(length(chargeRange) != 2){
			stop('Please provide a vector of length 2 for the chargeRange...')
		} else {}
		call <- paste(call, ' -minCharge ', chargeRange[1], ' -maxCharge ', chargeRange[2], sep='')
	} else {}
	
	call <- paste('java -Xmx1500M -jar ', R.home(component='library/pepmaps/java/MSGFplus.jar'), ' ', call, sep='')
	
	unlink(paste(tmp, '*', sep=''))
	
	system(call, ignore.stderr=TRUE, ignore.stdout=!verbose)
	cat('Importing results...')
	flush.console()
	callConv <- paste('java -Xmx1500M -cp ', R.home(component='library/pepmaps/java/MSGFplus.jar'), ' edu.ucsd.msjava.ui.MzIDToTsv -i ', tmp, ' -o ', paste(tmp, '.tsv', sep=''), ' -unroll 1')
	system(callConv)
	
	if(length(scan(paste(tmp, '.tsv', sep=''), skip=1, nlines=1, what='character', quiet=T)) == 0){
		warning(paste('No peptides detected in ', basename(file), sep=''))
		unlink(paste(tmp, '*', sep=''))
	} else {
		ans <- read.table(paste(tmp, '.tsv', sep=''), sep='\t')
		names(ans) <- scan(paste(tmp, '.tsv', sep=''), nlines=1, what=character(), quiet=TRUE)
		cat('DONE\n')
		flush.console()
		unlink(paste(tmp, '*', sep=''))
		ans
	}
}
collateMSGFplus <- function(directory, database, useML=FALSE, ...){
	if(missing(directory)){
		if(Sys.info()["sysname"] == 'Windows'){
			cat('Choose directory of .mzXML files to analyse...\n')
			flush.console()
			directory <- choose.dir()
		} else {
			directory <- readline('Path to directory with .mzXML files to analyse:')
		}
	} else {}
	if(useML){
		files <- list.files(directory, pattern='*.mzML', ignore.case=TRUE, full.names=TRUE)
	} else {
		files <- list.files(directory, pattern='*.mzXML', ignore.case=TRUE, full.names=TRUE)
	}
	
	if(length(files) == 0){
		stop('Directory doesn\'t contain any data files')
	} else {
		cat(paste('Read ', length(files), ' data files...\n', sep=''))
		flush.console()
	}
	if(missing(database)){
		cat('Choose database file...\n')
		flush.console()
		database <- file.choose()
	} else {}
	cat(paste('Database contains ', length(read.AAStringSet(database)), ' sequences...\n', sep=''))
	flush.console()
	res <- list()
	for(i in 1:length(files)){
		cat(paste(basename(files[i]), ' ', sep=''))
		flush.console()
		res[[i]] <- MSGFplus(files[i], database, ...)
	}
	cat('\n')
	flush.console()
	res <- do.call('rbind', res)
	identifier <- unique(res$Peptide)
	identifier <- data.frame(Peptide=identifier, Peptide.ID=1:length(identifier))
	res <- merge(res, identifier, all.x=TRUE, sort=FALSE)
	names(res) <- sub('#', '', names(res))
	if(useML){
		res$SpecFile <- sub('.mzML$', '.mzXML', res$SpecFile)
	}
	res
}
### Constructor for PepID objects
### Ask for location of MassAI/Crosswork results and creates the PepID object accordingly
pepID <- function(type, path=file.choose(), sep='\t', dec='.', directory, database, useML=FALSE, ...){
    if(missing(type)){
        new(Class='PepID', type=character(), raw=data.frame())
    } else {
        if(type == 'MassAI'){
            data <- read.table(path, skip=2, sep=sep, dec=dec, col.names=scan(path, sep=sep, nlines=1, what='character', quiet=TRUE), stringsAsFactors=FALSE)
            data <- data[!is.na(data$PeptideAID),]
			db <- Biostrings:::AAStringSet()
        } else if (type == 'Crosswork'){
            data <- read.table(path, sep=sep, dec=dec, header=TRUE, row.names=NULL, stringsAsFactors=FALSE)
            names(data) <- scan(path, sep=sep, nlines=1, what='character', quiet=TRUE)
            data <- data[,-ncol(data)]
			db <- Biostrings:::AAStringSet()
        } else if(type == 'MSGF+'){
			if(missing(directory)){
				if(Sys.info()["sysname"] == 'Windows'){
					cat('Choose the folder containing the mzXML files to analyse\n')
					flush.console()
					directory <- choose.dir()					
				} else {
					path <- readline('Path to directory with .mzXML files to analyse:')
				}
			} else {}
			if(missing(database)){
				cat('Choose database file (fasta format required)\n')
				flush.console()
				database <- file.choose()
			} else if(basename(database) == database){
				database <- paste(R.home(component='library/pepmaps/extdata/'), database, '.fasta', sep='')
			} else {}
			data <- collateMSGFplus(directory=directory, database=database, useML=useML, ...)
			data$rt <- NA
			datafiles <- list.files(directory, pattern='*.mzXML', ignore.case=TRUE, full.names=TRUE)
			for(i in 1:length(datafiles)){
				raw <- mzR::header(mzR::openMSfile(datafiles[i]))
				ind <- which(data$SpecFile == basename(datafiles)[i])
				scan <- data$Scan[ind]
				rt <- raw$retentionTime[match(scan, raw$acquisitionNum)]
				data$rt[ind] <- rt
			}
			db <- read.AAStringSet(database)
		} else {stop('Only MSGF+, MassAI and Crosswork supported')}
        new(Class='PepID', type=type, raw=data, database=db)
    }
}