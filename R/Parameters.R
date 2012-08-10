### Class for run parameter information
setClass(
    Class='Parameters',
    representation=representation(
        Name='character',
        Type='character',
        parameters='list'
    )
)
### Show method for Parameters
setMethod(
    'show', 'Parameters',
    function(object){
        if(length(object@Name) == 0){
            cat('An empty Parameters object\n')
        } else {
            cat(object@Name, '\n\n')
            cat('A parameter set of the type:',object@Type, '\n\n')
            cat('The set contains parameters for the following:\n\n')
            for (i in 1:length(names(object@parameters))){
                cat('\t', names(object@parameters)[i], '\n')
            }
        }
    }
)
### length method for Parameters
setMethod(
    'length', 'Parameters',
    function(x){
        if(sum(length(x@Name)+length(x@Type)+length(x@parameters)) == 0){
            0
        } else {
            1
        }
    }
)
### editPar for parameters
setMethod(
    'editPar', 'Parameters',
    function(object, type, ...){
        changes <- list(...)
        if(length(changes) != 0){
            para <- object@parameters
            if(!type %in% names(para)){
                stop('Unknown type...')
            } else {}
            for(i in 1:length(changes)){
                name <- names(changes)[i]
                para[[type]][[name]] <- changes[[i]]
            }
            if(!grepl('edit', object@Name)){
                object@Name <- paste(object@Name, '- edit', sep=' ')
            } else {}
            object@parameters <- para
        } else {}
        object
    }
)
### Constructor for Parameters
parameters <- function(Name){
    source(R.home(component='library/pepmaps/extdata/run_parameters.R'))
    if(missing(Name)){
        new(Class='Parameters', Name=character(), Type=character(), parameters=list())
    } else {
        if((Name %in% names(param)) == FALSE){
            stop('\n\nNo parameter set with name: <', Name, '> defined.\nUse setParameter() to define new set.')
        } else {
            param <- param[[which(names(param) == Name)]]
        }
        new(Class='Parameters', Name=Name, Type=param$Type, parameters=param[-1])
    }
}
### Function to create new parameter set
setParameters <- function(name, assist=TRUE){
    load()
    if(name %in% names(parameters)){
        stop('Parameter set with name: <', name, '> already set!\n\nUse editParameters() to edit or delete it.', sep='')
    } else {}
    if(assist){
        cat('Creating parameter set: ', name, '\n\n')
        type <- readline('Type:\n
                         Which type of parameter set should be created?\n\n
                         (only <LC-MS/MS supported)\t\t\t')
        while((type %in% c('LC-MS/MS')) == FALSE){
            type <- readline('Unknown type!\n\n
                             Only <LC-MS/MS> available at the moment')
        }
        if(type == 'LC-MS/MS'){
            parameters <- list(
                Type='LC-MS/MS',
                XCMS=list(
                    findPeak='centWave',
                    ppm=5,
                    peakwidth=c(7, 25),
                    snthresh=10,
                    prefilter=c(0, 0),
                    fwhm=NULL,
                    retcor='obiwarp',
                    plottype='deviation'
                ),
                CAMERA=list(
                    polarity='positive'
                ),
                Complist=list(
                    mzwin=0,
                    rtwin=0
                )
            )
            cat('Parameters used by XCMS:\n\n')
            cat('Which method should be used for peak picking?\n\n
                centWave:\tDetects peaks using continous wavelet transform\n
                \ton regions of interest. (default)
                matchedFilter:\tRuns filter matching using a Gaussian\n
                \tmodel peak on the MS data.')
            parameters$XCMS$findPeak <- readline('')
            while((parameters$XCMS$findPeak %in% c('centWave','matchedFilter')) == FALSE){
                parameters$XCMS$findPeak <- readline('Wrong input (either centWave or matchedFilter):')
            }
            if(parameters$XCMS$findPeak == 'centWave'){
                parameters$XCMS$ppm <- readline('ppm: The mz accuracy of the equipment. Should be set generously high (default=25)')
            }
        }
    }
}