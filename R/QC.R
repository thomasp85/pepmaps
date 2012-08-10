### Class that contains process control history and new entries for the history
setClass(
    Class='QC',
    representation=representation(
        path='character',
        history='list',
        newEntry='list',
        oldSamples='data.frame',
        newSamples='data.frame',
        event='data.frame',
        statistics='list',
        oldStatistics='list'
    ),
    prototype=prototype(
        path=character(),
        history=list(),
        newEntry=list(),
        oldSamples=data.frame(),
        newSamples=data.frame(),
        event=data.frame(),
        statistics=list(),
        oldStatistics=list()
    )
)
## Show method for QC
setMethod(
    'show', 'QC',
    function(object){
        if(length(object) == 0){
            cat('An empty QC object...\n')
        } else {
            nsamp <- nrow(object@oldSamples)
            if(length(object@newSamples) != 0){
                nsamp <- nsamp + nrow(object@newSamples)
            } else {}
            cat('A QC object with ', nsamp, ' samples.\n', sep='')
        }
    }
)
## length method for QC
setMethod(
    'length', 'QC',
    function(x){
        if(length(x@history)+length(x@newEntry) == 0){
            return(0)
        } else {
            nsamp <- nrow(x@oldSamples)
            if(length(x@newSamples) != 0){
                nsamp <- nsamp + nrow(x@newSamples)
            } else {}
            return(nsamp)
        }
    }
)
## getStatQC
setMethod(
    'getStatQC', 'QC',
    function(object, statistics, where='new'){
        if(where == 'new'){
            ans <- findName(object@newEntry, statistics)
            if(is.list(ans) & !is.data.frame(ans)){
                ans <- llply(ans, function(x, data) data.frame(Stat=x,data), data=object@newSamples)
            } else {
                ans <- data.frame(Stat=ans, object@newSamples)
            }
        } else if(where == 'old'){
            ans <- findName(object@history, statistics)
            if(is.list(ans) & !is.data.frame(ans)){
                ans <- llply(ans, function(x, data) data.frame(Stat=x,data), data=object@oldSamples)
            } else {
                ans <- data.frame(Stat=ans, object@oldSamples)
            }
        }
        ans
    }
)
## modelQC method
setMethod(
    'modelQC', 'QC',
    function(object, saveOld=TRUE){
        if(length(object@statistics) != 0 & saveOld){
            if(length(object@oldStatistics) == 0){
                object@oldStatistics <- list(object@statistics)
            } else {
                object@oldStatistics <- c(list(object@statistics), object@oldStatistics)
            }
        } else {}
        statQC <- function(x){
            x <- x[!x$Outlier,]
            if(ncol(x) > 5){
                ans <- list()
                for(i in 1:(ncol(x)-4)){
                    avg <- mean(x[,i])
                    UCL <- avg + sd(x[,i])*3
                    LCL <- avg - sd(x[,i])*3
                    ans[[i]] <- data.frame(avg, UCL, LCL, stat=names(x)[i])
                }
                ans <- do.call('rbind', ans)
            } else {
                avg <- mean(x[,1])
                UCL <- avg + sd(x[,1])*3
                LCL <- avg - sd(x[,1])*3
                ans <- data.frame(avg, UCL, LCL)
            }
            ans
        }
        ans <- list()
        for(i in 1:length(object@history)){
            stat <- getStatQC(object, names(object@history)[i], 'old')
            ans[[i]] <- llply(stat, statQC)
        }
        names(ans) <- names(object@history)
        object@statistics <- ans
        object
    }
)
## plotQC method for QC object
setMethod(
    'plotQC', 'QC',
    function(object, statistics, highlight=c('new', 'outlier', 'rules'), event=FALSE){
        if(length(object) == 0){
            stop('Cannot plot an empty WC object...\n')
        } else {}
        if(length(object@newEntry) == 0){
            highlight[highlight == 'new'] <- 'none'
        } else {}
        if(missing(statistics)){
            for(i in 1:length(object@history)){
                p <- plotQC(object, statistics=names(object@history)[i], highlight=highlight, event=event)
                print(p)
                readline('Press return for next plot...')
            }
        } else {
            if(statistics %in% names(object@history)){
                ansold <- getStatQC(object, statistics, where='old')
                ansold <- lapply(seq(along=ansold), function(x, data) data.frame(data[[x]], Statistics=names(data)[x], collection='old'), data=ansold)
                if(length(object@newEntry) != 0){
                    ansnew <- getStatQC(object, statistics)
                    ansnew <- lapply(seq(along=ansnew), function(x, data) data.frame(data[[x]], Statistics=names(data)[x], collection='new'), data=ansnew)
                    ans <- bind(ansold, ansnew)
                } else {
                    ans <- ansold
                }
                if(any('new' %in% highlight)){
                    p <- ggplot(data=ans[[1]], aes(x=Date, y=Stat, colour=factor(Run.order), shape=collection, size=collection))
                } else {
                    p <- ggplot(data=ans[[1]], aes(x=Date, y=Stat, colour=factor(Run.order)))
                }
                if(statistics != 'Calibrant'){
                    for(i in seq(along=ans)){
                        if(ans[[i]]$Statistics[1] != 'RS'){
                            stat <- data.frame(value=as.numeric(object@statistics[[statistics]][[i]]), stat=c('avg','cl','cl'), Statistics=ans[[i]]$Statistics[1])
                            p <- p + geom_hline(data=stat, aes(yintercept=value, linetype=stat), colour=I('red'))
                            p <- p + geom_point(data=ans[[i]])
                        } else {
                            stat <- object@statistics[[statistics]][[i]]
                            stat <- do.call('rbind', apply(stat, 1, function(x) data.frame(value=as.numeric(x[1:3]), stat=c('avg','cl','cl'))))
                            stat$Statistics <- ans[[i]]$Statistics[1]
                            p <- p + geom_hline(data=stat, aes(yintercept=value, linetype=stat), colour=I('red'))
                            p <- p + geom_linerange(data=ans[[i]], aes(x=Date, ymin=Stat.25., ymax=Stat.75.), legend=FALSE, position=position_dodge(width=0.5))
                        }
                        if(event){
                            if(length(object@event) != 0){
                                eve <- data.frame(object@event, Statistics=ans[[i]]$Statistics[1], y=mean(stat$value)+(max(stat$value)-mean(stat$value))*1.2)
                                p <- p + geom_vline(data=eve, aes(xintercept=Date), linetype=I(5), legend=FALSE)
                                if(i == 1){
                                    altstat <- stat[which(stat$stat == 'avg'),'value']+(max(stat$value)-stat[which(stat$stat == 'avg'),'value'])*1.2
                                    eve$y <- max(max(eve$y), altstat)
                                    eve <- ddply(eve, c(.(Date), .(Statistics), .(y)), function(x) paste('   ', paste(x$Event, collapse='\n   '), sep=''))
                                    names(eve)[which(names(eve) == 'V1')] <- 'Event'
                                    p <- p + geom_text(data=eve, aes(x=Date, y=y, label=Event, colour=NULL), size=I(2), hjust=0, legend=FALSE)
                                } else {}
                            } else {}
                        } else {}
                    }
                    p <- p + facet_grid(Statistics~., scales='free')
                } else {
                    ans <- melt(ans, measure=1:3)
                    stat <- melt(object@statistics$Calibrant, measure=1:3)
                    names(stat) <- c('variable', 'stat', 'value', 'Statistics')
                    ans$variable <- factor(sub('Stat.', '', ans$variable), levels=c('Presence', 'Median', 'Sd'))
                    stat$variable <- factor(sub('Stat.', '', stat$variable), levels=c('Presence', 'Median', 'Sd'))
                    stat$Statistics <- factor(stat$Statistics)
                    levels(stat$stat) <- c(levels(stat$stat), 'cl')
                    stat$stat[stat$stat == 'UCL'] <- 'cl'
                    stat$stat[stat$stat == 'LCL'] <- 'cl'
                    levels(stat$stat)[2:3] <- NA
                    if(event){
                        comb <- ddply(ans, c(.(Statistics), .(variable)), function(x) mean(x$value)+(max(x$value)-mean(x$value))*1.2)
                        names(comb)[3] <- 'y'
                        eve <- list()
                        for(i in 1:nrow(object@event)){
                            eve[[i]] <- data.frame(object@event[i,], comb, row.names=row.names(comb))
                        }
                        eve <- do.call('rbind', eve)
                        lab <- eve[which(eve$variable == 'Presence'),]
                        lab <- ddply(lab, c(.(Date), .(Statistics), .(y), .(variable)), function(x) paste('   ', paste(x$Event, collapse='\n   '), sep=''))
                        names(lab)[which(names(lab) == 'V1')] <- 'Event'
                        altstat <- stat[which(stat$variable == 'Presence'),]
                        altstat <- ddply(altstat, .(Statistics), function(x) x[which(x$stat == 'avg'),'value']+(max(x$value)-x[which(x$stat == 'avg'),'value'])*1.2)
                        lab$y <- max(max(lab$y), max(altstat$V1))
                        p <- p + geom_vline(data=eve, aes(xintercept=Date), linetype=I(5), legend=FALSE)
                        p <- p + geom_text(data=lab, aes(x=Date, y=y, label=Event, colour=NULL), size=I(2), hjust=0, legend=FALSE)
                    } else {}
                    p <- p + geom_hline(data=stat, aes(yintercept=value, linetype=stat), colour=I('red'))
                    p <- p + geom_point(data=ans, aes(y=value))
                    p <- p + facet_grid(variable~Statistics, scale='free')
                }
                p <- p + opts(title=paste('Control chart for variables related to: ', statistics, '\n', sep=''))
            } else {
                ansold <- getStatQC(object, statistics, where='old')
                if(is.null(ansold)){
                    stop('Unknown statistic..\n')
                } else {}
                ansold <- data.frame(ansold, collection='old')
                if(length(object@newEntry) != 0){
                    ansnew <- getStatQC(object, statistics)
                    ansnew <- data.frame(ansnew, collection='new')
                    ans <- rbind(ansold, ansnew)
                } else {
                    ans <- ansold
                }
                if(any('new' %in% highlight)){
                    p <- ggplot(data=ans, aes(x=Date, y=Stat, colour=factor(Run.order), shape=collection, size=collection))
                } else {
                    p <- ggplot(data=ans, aes(x=Date, y=Stat, colour=factor(Run.order)))
                }
                if(statistics != 'RS'){
                    stat <- findName(object@statistics, statistics)
                    stat <- data.frame(value=as.numeric(stat), stat=c('avg','cl','cl'))
                    p <- p + geom_hline(data=stat, aes(yintercept=value, linetype=stat), colour=I('red'))
                    p <- p + geom_point()
                } else {
                    stat <- findName(object@statistics, statistics)
                    stat <- do.call('rbind', apply(stat, 1, function(x) data.frame(value=as.numeric(x[1:3]), stat=c('avg','cl','cl'))))
                    p <- p + geom_hline(data=stat, aes(yintercept=value, linetype=stat), colour=I('red'))
                    p <- p + geom_linerange(aes(x=Date, ymin=Stat.25., ymax=Stat.75.), position=position_dodge(width=0.5))
                }
                if(event){
                    if(length(object@event) != 0){
                        eve <- data.frame(object@event, y=mean(stat$value)+(max(stat$value)-mean(stat$value))*1.2)
                        altstat <- stat[which(stat$stat == 'avg'),'value']+(max(stat$value)-stat[which(stat$stat == 'avg'),'value'])*1.2
                        eve$y <- max(max(eve$y), altstat)
                        eve <- ddply(eve, c(.(Date), .(y)), function(x) paste('   ', paste(x$Event, collapse='\n   '), sep=''))
                        names(eve)[which(names(eve) == 'V1')] <- 'Event'
                        p <- p + geom_vline(data=eve, aes(xintercept=Date), linetype=I(5), legend=FALSE)
                        p <- p + geom_text(data=eve, aes(x=Date, y=y, label=Event, colour=NULL), size=I(2), hjust=0, legend=FALSE)
                    } else {}
                } else {}
                p <- p + opts(title=paste('Control chart for: ', statistics, '\n', sep=''))
            }
            p <- p + theme_bw() + ylab('') + xlab('\nDate')
            p <- p + scale_colour_hue('Run order')
            p <- p + scale_linetype(legend=FALSE)
            p <- p + scale_y_continuous(expand=c(0.25,0))
            p <- p + opts(axis.text.x=theme_text(angle=45, hjust=1, vjust=1))
            if(any('new' %in% highlight)){
                p <- p + scale_shape('') + scale_size_manual('', values=c(2,4))
            } else {}
            p
        }
    }
)
## saveQC method for QC
setMethod(
    'saveQC', 'QC',
    function(object, newLocation=FALSE){
        if(length(object@path) == 0 | newLocation){
            name <- readline('\nType the name of the object: ')
			if(Sys.info()["sysname"] == 'Windows'){
	            readline('\nChoose save location: <Press Return>')
	            path <- choose.dir()				
			} else {
				path <- readline('Path to save location:')
			}
            filepath <- paste(path, '\\', name, '.RData', sep='')
            object@path <- filepath
            filelist <- read.csv(R.home(component='library/pepmaps/extdata/QCfiles.csv'), header=TRUE, stringsAsFactors=FALSE)
            if(name %in% filelist$Name){
                again <- TRUE
                while(again){
                    say <- paste('The name is already in use..\n\n', paste(filelist$Name, collapse='\n'), '\n\nChoose a name not on the above list: ', sep='')
                    name <- readline(say)
                    again <- name %in% filelist$Name
                }
            }
            filelist <- rbind(filelist, c(name, filepath))
            names(filelist) <- c('Name', 'Path')
            write.csv(filelist, file=R.home(component='library/pepmaps/extdata/QCfiles.csv'), row.names=FALSE)
        } else {
            filepath <- object@path
        }
        save(object, file=filepath)
    }
)
## loadQC
loadQC <- function(name){
    if(missing(name)){
        name <- readline('Type the name of the QC object to load: ')
    } else {}
    filelist <- read.csv(R.home(component='library/pepmaps/extdata/QCfiles.csv'), header=TRUE, stringsAsFactors=FALSE)
    if(!name %in% filelist$Name){
        again <- TRUE
        while(again){
            say <- paste('No QC object called <', name, '> is registered\n\n', paste(filelist$Name, collapse='\n'), '\n\nChoose a name on the above list (or type <Q> to abort): ', sep='')
            name <- readline(say)
            if(name == 'Q'){
                stop()
            } else {}
            again <- !name %in% filelist$Name
        }
    }
    path <- filelist$Path[which(filelist$Name == name)]
    object <- local(get(load(path)))
}
## addQCevent for QC
addQCevent <- function(Date, type){
    if(missing(Date)){
        Date <- readline('Type the date of the event (yyyy-mm-dd): ')
    } else {}
    Date <- as.Date(Date, format='%Y-%m-%d')
    event <- c('Service', 'Column change', 'Precolumn change', 'Filter change', 'New eluent', 'New Calibrant A', 'New Calibrant B', 'Cleaning of ion source')
    if(missing(type)){
        cat('Select the type of event:\n\n')
        cat('1: Service\n')
        cat('2: Column change\n')
        cat('3: Precolumn change\n')
        cat('4: Filter change\n')
        cat('5: New eluent\n')
        cat('6: New Calibrant A\n')
        cat('7: New Calibrant B\n')
        cat('8: Cleaning of ion source\n')
        cat('\n')
        type <- readline('Type the number: ')
    } else {
        if(!is.numeric(type)){
            stop('Wrong type argument. Must be a number')
        } else {
            if(!type %in% 1:8){
                stop('Wrong type argument. Must be between 1 and 8')
            } else {}
        }
    }
    if(type == 1){
        com <- readline('Any comments regarding service?: ')
    } else {
        com <- ''
    }
    ans <- data.frame(Date=Date, Event=event[as.numeric(type)], Comment=com)
    if(type %in% c('1', '5', '6', '7', '8')){
        files <- QCfile()
        if(is.null(files)){
            stop('No QC objects to add event to...')
        } else {}
        for(i in 1:nrow(files)){
            temp <- loadQC(files$Name[i])
            temp@event <- rbind(ans, temp@event)
            saveQC(temp)
        }
        cat('Event added to: ', paste(files$Name, collapse=', '), sep='')
    } else {
        files <- QCfile()
        cat('Which QC object does the event apply to:\n\n')
        for(i in 1:nrow(files)){
            cat(i, ': ', files$Name[i], '\n', sep='')
        }
        cat('\n')
        fileind <- readline('Type the number(s) (space separated): ')
        fileind <- as.numeric(unlist(strsplit(fileind, ' ')))
        files <- files[fileind,]
        for(i in 1:nrow(files)){
            temp <- loadQC(files$Name[i])
            temp@event <- rbind(ans, temp@event)
            saveQC(temp)
        }
        cat('Event added to: ', paste(files$Name, collapse=', '), sep='')
    }
    cat('\n')
}
## addQCfile for QC
addQCfile <- function(){
    readline('Choose the file to add: <Press Return>')
    path <- file.choose()
    type <- basename(path)
    if(!grepl('.rdata', type, ignore.case=TRUE)){
        stop('Unknown fileformat...\n')
    } else {}
    object = local(get(load(path)))
    if(class(object) != 'QC'){
        stop('Selected file is not a QC object...\n')
    }
    object@path <- path
    name <- readline('\nType the name of the object: ')
    filelist <- read.csv(R.home(component='library/pepmaps/extdata/QCfiles.csv'), header=TRUE, stringsAsFactors=FALSE)
    if(name %in% filelist$Name){
        again <- TRUE
        while(again){
            say <- paste('The name is already in use..\n\n', paste(filelist$Name, collapse='\n'), '\n\nChoose a name not on the above list: ', sep='')
            name <- readline(say)
            again <- name %in% filelist$Name
        }
    }
    filelist <- rbind(filelist, c(name, path))
    names(filelist) <- c('Name', 'Path')
    write.csv(filelist, file=R.home(component='library/pepmaps/extdata/QCfiles.csv'), row.names=FALSE)
    save(object, file=path)
}
## removeQCfile for QC
removeQCfile <- function(name){
    if(missing(name)){
        name <- readline('Type the name of the QC object to delete: ')
    } else {}
    filelist <- read.csv(R.home(component='library/pepmaps/extdata/QCfiles.csv'), header=TRUE, stringsAsFactors=FALSE)
    if(!name %in% filelist$Name){
        again <- TRUE
        while(again){
            say <- paste('No QC object called <', name, '> is registered\n\n', paste(filelist$Name, collapse='\n'), '\n\nChoose a name on the above list (or type <Q> to abort): ', sep='')
            name <- readline(say)
            if(name == 'Q'){
                stop()
            } else {}
            again <- !name %in% filelist$Name
        }
    }
    index <- which(filelist$Name == name)
    filelist <- filelist[-index,]
    names(filelist) <- c('Name', 'Path')
    write.csv(filelist, file=R.home(component='library/pepmaps/extdata/QCfiles.csv'), row.names=FALSE)
}
## show QCfile
QCfile <- function(){
    filelist <- read.csv(R.home(component='library/pepmaps/extdata/QCfiles.csv'), header=TRUE, stringsAsFactors=FALSE)
    if(nrow(filelist) == 0){
        cat('No QC objects defined...\n')
    } else {
        filelist
    }
}
## setOutlier method for QC
setMethod(
    'setOutlier', 'QC',
    function(object, ID, value=TRUE, where='new'){
        if(where == 'new'){
            if(max(ID) > nrow(object@newSamples)){
                stop('ID extends range of new samples...\n')
            } else {}
            object@newSamples$Outlier[ID] <- value
        } else if(where == 'old'){
            if(max(ID) > nrow(object@oldSamples)){
                stop('ID extends range of old samples...\n')
            } else {}
            object@oldSamples$Outlier[ID] <- value
        } else {
            stop('Unknown <where> argument...\n')
        }
        object
    }
)
## creator function for QC
createQC <- function(mzXML, ID, Date, Runorder, calibrant=c(121.0509, 922.0098)){
    if(missing(mzXML)){
		if(Sys.info()["sysname"] == 'Windows'){
			readline('\nChoose directory containing only control samples: <Press Return>')
			path <- choose.dir()
		} else {
			path <- readline('Path to directory containing only control samples:')
		}
        mzXML <- list.files(path, pattern='*.mzXML', full.names=TRUE, recursive=TRUE, ignore.case=TRUE)
        if(length(mzXML) == 0){
            stop('Directory does not contain any compatible files...\n')
        } else {}
    } else {}
    if(length(mzXML) < 25){
        warning('A new QC set should not be based on less than 25 repetitions...\n')
    } else {}
    if(missing(ID)){
        readline('\nChoose file with MassAI result from controls: <Press Return>')
        ID <- file.choose()
    } else {}
    pID <- pepID('MassAI', path=ID)
    pID <- dlply(pID@raw, .(File))
    names(pID) <- sub('.mgf', '', basename(names(pID)), ignore.case=TRUE)
    name <- sub('.mzXML', '', basename(mzXML), ignore.case=TRUE)
    name <- sub('.xml', '', name, ignore.case=TRUE)
    if(sum(name %in% names(pID)) != length(name)){
        stop('\nData and ID files does not match...\n')
    } else if(length(mzXML) != length(pID)){
        warning('\nNumber of data files does not correspond to the number of samples in ID file...\n')
        pID <- pID[which(names(pID) %in% name)]
    } else {}
    pID <- pID[match(name, names(pID))]
    if(missing(Date)){
        anadate <- readline('\nType the run-date for each sample (yyyy-mm-dd)(space seperated): ')
        anadate <- unlist(strsplit(anadate, ' '))
        anadate <- as.Date(anadate, format='%Y-%m-%d')
        cat('\n')
    } else {
        anadate <- Date
    }    
    if(length(anadate != 1)){
        if(length(mzXML) != length(anadate)){
            stop('Provided dates does not match the number of samples...\n')
        } else {}
    } else {
        warning('A new QC set should be based on runs from more than one date...\n')
    }
    if(missing(Runorder)){
        dates <- unique(anadate)
        if(length(dates) == 1){
            cat(paste(name,collapse='\n'))
            cat('\n')
            flush.console()
            anaorder <- readline('\nType run-order of the experiments (space seperated): ')
            anaorder <- as.numeric(unlist(strsplit(anaorder, ' ')))
            if(length(anaorder) != length(name)){
                retry <- TRUE
                while(retry){
                    cat('\nNumber of run-order entries does not match the number of experiments...\n\n')
                    cat(paste(name,collapse='\n'))
                    cat('\n')
                    flush.console()
                    anaorder <- readline('\nType run-order of the experiments (space seperated): ')
                    anaorder <- as.numeric(unlist(strsplit(anaorder, ' ')))
                    retry <- length(anaorder) != length(name)
                }
            } else {}
        } else {
            anaorder <- rep(NA, length(name))
            for(i in 1:length(dates)){
                index <- which(anadate %in% dates[i])
                cat(paste(name[index],collapse='\n'))
                cat('\n')
                flush.console()
                say <- paste('\nType run-order of the experiments from ', dates[i], ' (space seperated)(leave blank if experiments are listed in order): ', sep='')
                anaorder1 <- readline(say)
                cat('\n')
                if(anaorder1 == ''){
                    anaorder1 <- 1:length(index)
                } else {
                    anaorder1 <- as.numeric(unlist(strsplit(anaorder1, ' ')))
                }
                if(length(anaorder1) != length(index)){
                    retry <- TRUE
                    while(retry){
                        cat('\nNumber of run-order entries does not match the number of experiments...\n\n')
                        cat(paste(name[index],collapse='\n'))
                        cat('\n')
                        flush.console()
                        anaorder1 <- readline(say)
                        cat('\n')
                        if(anaorder1 == ''){
                            anaorder1 <- 1:length(index)
                        } else {
                            anaorder1 <- as.numeric(unlist(strsplit(anaorder1, ' ')))
                        }
                        retry <- length(anaorder1) != length(index)
                    }
                } else {}
                anaorder[index] <- anaorder1
            }
        }
    } else {
        if(length(Runorder) != length(name)){
            retry <- TRUE
            while(retry){
                cat('\nNumber of run-order entries does not match the number of experiments...\n\n')
                cat(paste(name,collapse='\n'))
                cat('\n')
                flush.console()
                anaorder <- readline('\nType run-order of the experiments (space seperated): ')
                anaorder <- as.numeric(unlist(strsplit(anaorder, ' ')))
                retry <- length(anaorder) != length(name)
            }
        } else {
            anaorder <- Runorder
        }
    }
    dateorder <- sort(unique(anadate))
    ans <- list()
    for(i in 1:length(dateorder)){
        cat('Analysing experiment: ', paste(dateorder[i]), '\n', sep='')
        flush.console()
        index <- which(anadate %in% dateorder[i])
        ans[[i]] <- calcQC(mzXML[index], ID=pID[index], Date=dateorder[i], Runorder=anaorder[index], calibrant=calibrant)
    }
    ans <- do.call('bind', ans)
    experiment <- ans$ExperimentInfo
    ans$ExperimentInfo <- NULL
    experiment <- data.frame(experiment, Outlier=FALSE)
    experiment$Date <- as.Date(experiment$Date)
    object <- new(
        Class='QC',
        history=ans,
        oldSamples=experiment
    )
    object <- modelQC(object)
    object
}
## update function for QC
updateQC <- function(name, mzXML, ID, Date, Runorder, calibrant=c(121.0509, 922.0098)){
    if(missing(name)){
        name <- readline('Type the name of the QC object to load: ')
    } else {}
    filelist <- read.csv(R.home(component='library/pepmaps/extdata/QCfiles.csv'), header=TRUE, stringsAsFactors=FALSE)
    if(!name %in% filelist$Name){
        again <- TRUE
        while(again){
            say <- paste('No QC object called <', name, '> is registered\n\n', paste(filelist$Name, collapse='\n'), '\n\nChoose a name on the above list (or type <Q> to abort): ', sep='')
            name <- readline(say)
            if(name == 'Q'){
                stop()
            } else {}
            again <- !name %in% filelist$Name
        }
    }
    if(missing(mzXML)){
		if(Sys.info()["sysname"] == 'Windows'){
			readline('\nChoose directory containing only control samples: <Press Return>')
			path <- choose.dir()
		} else {
			path <- readline('Path to directory containing only control samples:')
		}
        mzXML <- list.files(path, pattern='*.mzXML', full.names=TRUE, recursive=TRUE, ignore.case=TRUE)
        if(length(mzXML) == 0){
            stop('Directory does not contain any compatible files...\n')
        } else {}
    } else {}
    if(missing(ID)){
        readline('\nChoose file with MassAI result from controls: <Press Return>')
        ID <- file.choose()
    } else {}
    pID <- pepID('MassAI', path=ID)
    pID <- dlply(pID@raw, .(File))
    names(pID) <- sub('.mgf', '', basename(names(pID)), ignore.case=TRUE)
    name1 <- sub('.mzXML', '', basename(mzXML), ignore.case=TRUE)
    name1 <- sub('.xml', '', name1, ignore.case=TRUE)
    if(sum(name1 %in% names(pID)) != length(name1)){
        stop('\nData and ID files does not match...\n')
    } else if(length(mzXML) != length(pID)){
        warning('\nNumber of data files does not correspond to the number of samples in ID file...\n')
        pID <- pID[which(names(pID) %in% name1)]
    } else {}
        if(missing(Date)){
        anadate <- readline('\nType the date of the experiment (yyyy-mm-dd): ')
        anadate <- as.Date(anadate, format='%Y-%m-%d')
        cat('\n')
    } else {
        anadate <- Date
    }
    if(missing(Runorder)){
        cat(paste(name1,collapse='\n'))
        cat('\n')
        flush.console()
        anaorder <- readline('\nType run-order of the experiments (space seperated)(leave blank if experiments are listed in order): ')
        if(anaorder == ''){
            anaorder <- 1:length(name1)
        } else {
            anaorder <- as.numeric(unlist(strsplit(anaorder, ' ')))
        }
    } else {
        anaorder <- Runorder
    }
    if(length(anaorder) != length(name1)){
        retry <- TRUE
        while(retry){
            cat('\nNumber of run-order entries does not match the number of experiments...\n\n')
            cat(paste(name1,collapse='\n'))
            cat('\n')
            flush.console()
            anaorder <- readline('\nType run-order of the experiments (space seperated): ')
            anaorder <- as.numeric(unlist(strsplit(anaorder, ' ')))
            retry <- length(anaorder) != length(name1)
        }
    } else {}
    ans <- calcQC(mzXML=mzXML, ID=pID, Date=anadate, Runorder=anaorder, calibrant=calibrant)
    path <- filelist$Path[which(filelist$Name == name)]
    object <- local(get(load(path)))
    if(length(object@newEntry) != 0){
        object@history <- do.call('bind', list(object@history, object@newEntry))
        object@oldSamples <- rbind(object@oldSamples, object@newSamples)
    } else {}
    experiment <- ans$ExperimentInfo
    ans$ExperimentInfo <- NULL
    experiment <- data.frame(experiment, Outlier=FALSE)
    experiment$Date <- as.Date(experiment$Date)
    object@newEntry <- ans
    object@newSamples <- experiment
    object
}
## calculate QC parameters for a set of controls
calcQC <- function(mzXML, ID, Date, Runorder, calibrant=c(121.0509, 922.0098)){
    if(missing(mzXML)){
		if(Sys.info()["sysname"] == 'Windows'){
			readline('\nChoose directory containing only control samples: <Press Return>')
			path <- choose.dir()
		} else {
			path <- readline('Path to directory containing only control samples:')
		}
        mzXML <- list.files(path, pattern='*.mzXML', full.names=TRUE, recursive=TRUE, ignore.case=TRUE)
        if(length(mzXML) == 0){
            stop('Directory does not contain any compatible files...\n')
        } else {}
    } else {}
    if(missing(ID)){
        readline('\nChoose file with MassAI result from controls: <Press Return>')
        ID <- file.choose()
        pID <- pepID('MassAI', path=ID)
        pID <- dlply(pID@raw, .(File))
    } else if(is.character(ID)){
        pID <- pepID('MassAI', path=ID)
        pID <- dlply(pID@raw, .(File))
    } else {
        pID <- ID
    }
    names(pID) <- sub('.mgf', '', basename(names(pID)), ignore.case=TRUE)
    name <- sub('.mzXML', '', basename(mzXML), ignore.case=TRUE)
    name <- sub('.xml', '', name, ignore.case=TRUE)
    if(sum(name %in% names(pID)) != length(name)){
        stop('\nData and ID files does not match...\n')
    } else if(length(mzXML) != length(pID)){
        warning('\nNumber of data files does not correspond to the number of samples in ID file...\n')
    } else {}
    if(missing(Date)){
        anadate <- readline('\nType the date of the experiment (yyyy-mm-dd): ')
        anadate <- as.Date(anadate, format='%Y-%m-%d')
        cat('\n')
    } else {
        anadate <- Date
    }
    if(missing(Runorder)){
        cat(paste(name,collapse='\n'))
        cat('\n')
        flush.console()
        anaorder <- readline('\nType run-order of the experiments (space seperated)(leave blank if experiments are listed in order): ')
        if(anaorder == ''){
            anaorder <- 1:length(name)
        } else {
            anaorder <- as.numeric(unlist(strsplit(anaorder, ' ')))
        }
    } else {
        anaorder <- Runorder
    }
    if(length(anaorder) != length(name)){
        retry <- TRUE
        while(retry){
            cat('\nNumber of run-order entries does not match the number of experiments...\n\n')
            cat(paste(name,collapse='\n'))
            cat('\n')
            flush.console()
            anaorder <- readline('\nType run-order of the experiments (space seperated): ')
            anaorder <- as.numeric(unlist(strsplit(anaorder, ' ')))
            retry <- length(anaorder) != length(name)
        }
    } else {}
    ans <- list()
    for(i in 1:length(mzXML)){
        pID1 <- pID[[which(names(pID) == name[i])]]
        cat(paste('\nAnalysing ', name[i], '...\n'))
        flush.console()
        raw <- xcmsRaw(mzXML[i], includeMSn=T)
        peakID <- findPeaks.centWavePrec(raw, ppm=25, peakwidth=c(5, 25), snthresh=8, tolerance=1, fitgauss=TRUE, noise=0, verbose.columns=TRUE, combine.peaks=FALSE)
        uorder <- order(peakID[, "into"], decreasing = TRUE)
        pm <- as.matrix(peakID[, c("mzmin", "mzmax", "rtmin", "rtmax"), drop = FALSE])
        uindex <- xcms:::rectUnique(pm, uorder, -0.001, ydiff = -1e-05)
        peakIDc <- peakID[uindex, , drop = FALSE]
        peak <- findPeaks.centWave(raw, ppm=25, peakwidth=c(5, 25), snthresh=0, prefilter=c(1,1000), fitgauss=TRUE, noise=0, verbose.columns=TRUE)
        RS <- quantile(unique(pID1[,c('Retention', 'Precursor')])$Retention)
        RS <- RS[c(2,4)]
        ans$Chromatography$RS <- rbind(ans$Chromatography$RS, RS)
        ans$Chromatography$Peak.width <- c(ans$Chromatography$Peak.width, median(peak[,'sigma']*2*sqrt(2*ln(2)), na.rm=T))
        ans$Chromatography$Peak.width.var <- c(ans$Chromatography$Peak.width.var, summary(peak[,'sigma']*2*sqrt(2*ln(2)))[5]-summary(peak[,'sigma']*2*sqrt(2*ln(2)))[2])
        ans$Chromatography$Pep.pr.min <- c(ans$Chromatography$Pep.pr.min, nrow(unique(pID1[which(pID1$Retention > RS[1] & pID1$Retention < RS[2]),c('Retention', 'Precursor')]))/(RS[2]/60-RS[1]/60))
        ans$Ion.Source$Percent1H <- c(ans$Ion.Source$Percent1H, sum(pID1$Charge == 1)/length(pID1$Charge)*100)
        ans$Ion.Source$Percent3H <- c(ans$Ion.Source$Percent3H, sum(pID1$Charge == 3)/length(pID1$Charge)*100)
        ans$Ion.Source$Percent4H <- c(ans$Ion.Source$Percent4H, sum(pID1$Charge == 4)/length(pID1$Charge)*100)
        tic <- apply(raw@env$profile,2,sum)
        tic <- tic[which(raw@scantime > RS[1] & raw@scantime < RS[2])]
        tic <- tic[1:(length(tic)-1)]/tic[2:length(tic)]
        ans$Ion.Source$TIC.jump <- c(ans$Ion.Source$TIC.jump, sum(tic < 0.1 & tic > 10, na.rm=T))
        ans$Ion.Source$Prec.mz <- c(ans$Ion.Source$Prec.mz, median(pID1$Precursor))
        RSpeak <- peakIDc[which(peakIDc[,'rtmax'] > RS[1] & peakIDc[,'rtmin'] < RS[2]),]
        ans$MS1$Pep.intens <- c(ans$MS1$Pep.intens, log(median(RSpeak[,'maxo'])))
        dyn <- quantile(RSpeak[,'maxo'], probs=c(0.05,0.95), type=3)
        ans$MS1$Pep.dyn <- c(ans$MS1$Pep.dyn, log(dyn[2]/dyn[1]))
        scans <- which(raw@scantime > RS[1] & raw@scantime < RS[2])
        spec <- alply(scans, 1, function(x, r) getScan(r, x), r=raw)
        ans$MS1$MS1.SN <- c(ans$MS1$MS1.SN, median(laply(spec, function(x) max(x[,'intensity'])/median(x[,'intensity']))))
        spec1 <- do.call('rbind', spec)
        ans$MS1$MS1.SN.Total <- c(ans$MS1$MS1.SN.Total, max(spec1[,'intensity'])/median(spec1[,'intensity']))
        ans$MS1$Med.TIC <- c(ans$MS1$Med.TIC, log(median(apply(raw@env$profile,2,sum)[which(raw@scantime > RS[1] & raw@scantime < RS[2])])/1000))
        ans$Dynamic.Sampling$Trigger <- c(ans$Dynamic.Sampling$Trigger, median(peakID[,'precursorInt']/peakID[,'maxo'], na.rm=T))
        ans$Dynamic.Sampling$TriggerLow <- c(ans$Dynamic.Sampling$TriggerLow, median(peakID[which(peakID[,'into'] < mean(peakID[,'into'])),'precursorInt']/peakID[which(peakID[,'into'] < mean(peakID[,'into'])),'maxo'], na.rm=T))
        ans$Dynamic.Sampling$Oversamp <- c(ans$Dynamic.Sampling$Oversamp, sum(table(peakID[,'maxo']) > 1)/length(table(peakID[,'maxo'])))
        ans$Dynamic.Sampling$N.MS1 <- c(ans$Dynamic.Sampling$N.MS1, sum(raw@scantime > RS[1] & raw@scantime < RS[2]))
        ans$Dynamic.Sampling$N.MS2 <- c(ans$Dynamic.Sampling$N.MS2, sum(raw@msnRt > RS[1] & raw@msnRt < RS[2]))
        sind <- raw@msnScanindex
        ans$MS2$N.peak.MS2 <- c(ans$MS2$N.peak.MS2, median(sind[2:length(sind)]-sind[1:length(sind)-1]))
        top25 <- quantile(peak[,'maxo'], probs=c(0.25, 0.75), type=3)
        ans$MS2$N.top.PepID <- c(ans$MS2$N.top.PepID, sum(peakIDc[,'maxo'] > top25[2]))
        ans$MS2$N.low.PepID <- c(ans$MS2$N.low.PepID, sum(peakIDc[,'maxo'] < top25[1]))
        scans <- which(raw@msnAcquisitionNum %in% pID1$Scan.number)
        spec <- alply(scans, 1, function(x, r) getMSnScan(r, x), r=raw)
        ans$MS2$MS2.SN <- c(ans$MS2$MS2.SN, median(laply(spec, function(x) max(x[,'intensity'])/median(x[,'intensity']))))
        spec1 <- do.call('rbind', spec)
        ans$MS2$MS2.SN.Total <- c(ans$MS2$MS2.SN.Total, max(spec1[,'intensity'])/median(spec1[,'intensity']))
        ans$Identification$Lev2.frac <- c(ans$Identification$Lev2.frac, sum(pID1$Notes == 'Lev2')/length(pID1$Notes))
        ans$Identification$N.Pep <- c(ans$Identification$N.Pep, length(unique(pID1$Peptide)))
        ans$Identification$ID.score <- c(ans$Identification$ID.score, median(pID1$Score, na.rm=T))
        ans$Identification$Mass.error <- c(ans$Identification$Mass.error, median(pID1$Mass.error, na.rm=T))
        for(j in 1:length(calibrant)){
            mzwin <- calibrant[j] + c(-0.0005, 0.0005)
            cal <- xcms:::rawMZ(raw, mzrange=mzwin)
            exist <- which(round(cal, digits=4) == calibrant[j])
            ans$Calibrant[[paste(calibrant[j])]]$Presence <- c(ans$Calibrant[[paste(calibrant[j])]]$Presence, length(exist)/length(cal))
            intcal <- xcms:::rawEIC(raw, mzrange=mzwin)
            ans$Calibrant[[paste(calibrant[j])]]$Median <- c(ans$Calibrant[[paste(calibrant[j])]]$Median, median(intcal$intensity[exist]))
            ans$Calibrant[[paste(calibrant[j])]]$Sd <- c(ans$Calibrant[[paste(calibrant[j])]]$Sd, sd(intcal$intensity[exist]))
        }
    }
    if(length(mzXML) > 1){
        rtc <- retcorRaw(mzXML, anaorder, response=10)
        minlength <- min(sapply(rtc$original, length))
        rtcdev <- list()
        for(i in 1:length(rtc$original)){
            rtcdev[[i]] <- (rtc$original[[i]] - rtc$corrected[[i]])[1:minlength]
        }
        rtcdev <- do.call('cbind', rtcdev)
        ans$Chromatography$rtdevmedianmax <- rep(max(apply(rtcdev, 2, function(x) median(abs(x)))), length(mzXML))
        rtcdev <- apply(rtcdev, 1, function(x) diff(range(x)))
        ans$Chromatography$rtdevmax <- rep(max(rtcdev), length(mzXML))
        ans$Chromatography$rtdevmaxtime <- rep(rtc$original[[1]][which.max(rtcdev)], length(mzXML))    
    } else {
        ans$Chromatography$rtdevmedianmax <- rep(NA, length(mzXML))
        ans$Chromatography$rtdevmax <- rep(NA, length(mzXML))
        ans$Chromatography$rtdevmaxtime <- rep(NA, length(mzXML))    
    }
    ans$ExperimentInfo <- data.frame(Sample.names=name, Run.order=anaorder, Date=anadate)
    names(ans$Chromatography$RS) <- NULL
    row.names(ans$Chromatography$RS) <- NULL
    names(ans$Chromatography$Peak.width.var) <- NULL
    names(ans$Chromatography$Pep.pr.min) <- NULL
    names(ans$MS1$Pep.dyn) <- NULL
    invisible(ans)
}