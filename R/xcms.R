### xcms expansion
setMethod(
    'plotFragment', 'xcmsRaw',
    function(object, ID, pepseq, tolerance=0.2, precursor.rm=FALSE){
        if(length(object@msnLevel) == 0){
            stop('No Fragment information in xcmsRaw object')
        } else if(!(ID %in% object@msnAcquisitionNum)){
            stop('No fragment at scan number: ', ID)
        } else if(!missing(pepseq)){
            if(!(sum(strsplit(tolower(pepseq), '')[[1]] %in% c('g','a','s','p','v','t','c','l','i','n','d','q','k','e','m','h','f','r','y','w'))==nchar(pepseq))){
                stop('Invalid peptide sequence.')
            } else {}
        } else {}
        ID <- which(object@msnAcquisitionNum == ID)
        index <- c(object@msnScanindex[ID], object@msnScanindex[ID+1]-1)
        data <- data.frame('mz'=object@env$msnMz[index[1]:index[2]], 'Intensity'=object@env$msnIntensity[index[1]:index[2]], 'ymin'=0)
        precurmz <- data$mz[which(abs(data$mz-(object@msnPrecursorMz[ID])) == min(abs(data$mz-(object@msnPrecursorMz[ID]))))]
        precurint <- max(data$Intensity[which(data$mz > precurmz-0.1 & data$mz < precurmz+0.1)])
        precur <- data.frame('mz'=precurmz, 'Intensity'=precurint)
		if(precursor.rm){
			data <- data[-which.min(abs(data$mz-(object@msnPrecursorMz[ID]))),]
		} else {}
		info <- paste('Precursor info:\nIntensity: ', object@msnPrecursorIntensity[ID], '\nRetention time: ', object@msnRt[ID], '\nm/z: ', object@msnPrecursorMz[ID], '\nCharge: ', object@msnPrecursorCharge[ID], '\nCollision energy: ', object@msnCollisionEnergy[ID], sep='')
		info <- data.frame('mz'=min(data$mz), 'y'=max(data$Intensity), 'label'=info)
        p <- ggplot(data=data) + geom_linerange(aes(x=mz, ymax=Intensity, ymin=ymin))
        p <- p + geom_text(aes(x=mz, y=y, label=label), data=info, hjust=0, vjust=1, size=3, alpha=0.5)
        if(!precursor.rm){
			p <- p + geom_point(aes(x=mz, y=Intensity, colour=I('red')), size=4, data=precur) + scale_colour_hue('', breaks='red', labels='Precursor')
		}
        p <- p + theme_bw(base_size=10) + xlab('\nm/z') + ylab('Intensity')
        if(!missing(pepseq)){
            p <- p + ggtitle(parse(text=fragExp(pepseq)))
            ionlab <- data.frame(data, ion=NA)
            ionlist <- fragPattern(pepseq)
            for(i in 1:nrow(ionlist)){
                ind <- which(ionlab$mz < ionlist$mz[i]+tolerance & ionlab$mz > ionlist$mz[i]-tolerance)
                if(length(ind) == 1){
                    ionlab$ion[ind] <- ionlist$ion[i]
                } else if(length(ind) > 1){
                    ind <- ind[which(abs(ionlab$mz[ind]-ionlist$mz[i]) == min(abs(ionlab$mz[ind]-ionlist$mz[i])))]
                    ionlab$ion[ind] <- ionlist$ion[i]
                } else {}
            }
            ionlab <- ionlab[!is.na(ionlab$ion),]
            if(nrow(ionlab) != 0){
                cat(nrow(ionlab), 'fragment match\n')
                p <- p + geom_text(data=ionlab, aes(x=mz, y=Intensity, label=ion), colour=I('blue'))
            } else {
                cat('No fragment match\n')
            }
        } else {}
        p
    }
)
### getMSnScan for xcmsRaw
setMethod(
    'getMSnScan', 'xcmsRaw',
    function(object, scan){
        idx <- seq(object@msnScanindex[scan] + 1, min(object@msnScanindex[scan + 1], length(object@env$msnMz), na.rm = TRUE))
        points <- cbind(mz = object@env$msnMz[idx], intensity = object@env$msnIntensity[idx])
        invisible(points)
    }
)
### New findmzROI method
findmzROIprec <- function(object, scanrange = scanrange, dev = dev, tolerance = tolerance, noise = noise){
    if(length(object@msnPrecursorScan) == 0){
        stop('MSn data required...\n')
    } else {}
    precursorindex <- which(object@msnPrecursorScan %in% scanrange[1]:scanrange[2])
    perc <- round(seq(min(precursorindex), max(precursorindex), length=10))
    cat('% finished: 0')
    flush.console()
    prog <- 1
    ans <- list()
    for(i in precursorindex){
        if(which(i <= perc)[1] > prog){
            cat(' ')
            cat(prog*10)
            flush.console()
            prog <- prog+1
        } else {}
        if(!is.na(object@msnPrecursorScan[i])){
            mzdev <- object@msnPrecursorMz[i]*dev
            mzwin <- c(object@msnPrecursorMz[i]-mzdev, object@msnPrecursorMz[i]+mzdev)
            searchback <- TRUE
            miss <- 0
            scmin <- object@msnPrecursorScan[i]
            while(searchback){
                sc <- getScan(object, scmin, mzrange=mzwin)
                if(nrow(sc) == 0){
                    miss <- miss + 1
                } else if(max(sc[,'intensity']) < noise){
                    miss <- miss + 1
                } else {}
                scmin <- scmin-1
                if(miss > tolerance){
                    searchback <- FALSE
                } else {}
                if(scmin < 1){
                    scmin <- 1
                    searchback <- FALSE
                }
            }
            searchforward <- TRUE
            miss <- 0
            scmax <- object@msnPrecursorScan[i]
            while(searchforward){
                sc <- getScan(object, scmax, mzrange=mzwin)
                if(nrow(sc) == 0){
                    miss <- miss + 1
                } else if(max(sc[,'intensity']) < noise){
                    miss <- miss + 1
                } else {}
                scmax <- scmax+1
                if(miss > tolerance){
                    searchforward <- FALSE
                } else {}
                if(scmax > length(object@scantime)){
                    scmax <- scmax-1
                    searchforward <- FALSE
                } else {}
            }
            sclength <- length(scmin:scmax)
            intensity <- sum(do.call('rbind', alply(scmin:scmax, 1, function(x, object) getScan(object, x, mzrange=mzwin), object=object))[,2])
            precursorInt <- max(getScan(object, object@msnPrecursorScan[i], mzrange=mzwin)[,'intensity'])
            if(length(precursorInt) == 0){
                precursorInt <- NA
            } else {}
            ans[[i]] <- list(mz=object@msnPrecursorMz[i], mzmin=mzwin[1], mzmax=mzwin[2], scmin=scmin, scmax=scmax, length=sclength, intensity=intensity, MS2scan=i, precursorInt=precursorInt)
        }
    }
    cat(' 100')
    flush.console()
    ans <- ans[sapply(ans, function(x) if(is.null(x)) FALSE else x$intensity > 0)]
    cat('\n', length(ans), ' Precursor ROIs detected\n')
    flush.console()
    ans
}
### New findPeaks method
findPeaks.centWavePrec <- function (object, ...){
    .local <- function (object, ppm = 25, peakwidth = c(20, 50), 
        snthresh = 10, tolerance = 1, mzCenterFun = "wMean", 
        integrate = 1, mzdiff = -0.001, fitgauss = FALSE, scanrange = numeric(), 
        noise = 0, sleep = 0, verbose.columns = FALSE, combine.peaks = TRUE) 
    {
        if(!xcms:::isCentroided(object)){
            warning("It looks like this file is in profile mode. centWave can process only centroid mode data !\n")
        } else {}
        mzCenterFun <- paste("mzCenter", mzCenterFun, sep = ".")
        if(!exists(mzCenterFun, mode = "function", envir=loadNamespace('xcms'))){
            stop("Error: >", mzCenterFun, "< not defined ! \n")
        } else {}
        scanrange.old <- scanrange
        if(length(scanrange) < 2){ 
            scanrange <- c(1, length(object@scantime))
        } else {
            scanrange <- range(scanrange)
        }
        scanrange[1] <- max(1, scanrange[1])
        scanrange[2] <- min(length(object@scantime), scanrange[2])
        if(!(identical(scanrange.old, scanrange)) && (length(scanrange.old) > 0)){
            cat("Warning: scanrange was adjusted to ", scanrange, "\n")
        } else {}
        basenames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "intb", "maxo", "sn")
        verbosenames <- c("egauss", "mu", "sigma", "h", "f", "dppm", "scale", "scpos", "scmin", "scmax", "lmin", "lmax", "MS2scan", "precursorInt")
        scalerange <- round((peakwidth/mean(diff(object@scantime)))/2)
        if(length(scalerange) > 1){ 
            scales <- seq(from = scalerange[1], to = scalerange[2], by = 2)
        } else {
            scales <- scalerange
        }
        dev <- ppm * 1e-06
        minPeakWidth <- scales[1]
        noiserange <- c(minPeakWidth * 3, max(scales) * 3)
        maxGaussOverlap <- 0.5
        minPtsAboveBaseLine <- max(4, minPeakWidth - 2)
        minCentroids <- minPtsAboveBaseLine
        scRangeTol <- maxDescOutlier <- floor(minPeakWidth/2)
        peaklist <- list()
        cat("\n Detecting mass traces at", ppm, "ppm ... \n")
        flush.console()
        featlist <- findmzROIprec(object, scanrange = scanrange, dev = dev, tolerance = tolerance, noise = noise)
        scantime <- object@scantime
        Nscantime <- length(scantime)
        lf <- length(featlist)
        if(lf == 0){
            cat("No ROIs found ! \n")
            if(verbose.columns){
                nopeaks <- new("xcmsPeaks", matrix(nrow = 0, ncol = length(basenames) + length(verbosenames)))
                colnames(nopeaks) <- c(basenames, verbosenames)
            } else {
                nopeaks <- new("xcmsPeaks", matrix(nrow = 0, ncol = length(basenames)))
                colnames(nopeaks) <- c(basenames)
            }
            return(invisible(nopeaks))
        } else {}
        cat("\n Detecting chromatographic peaks ... \n % finished: ")
        lp <- -1
        for(f in 1:lf){
            perc <- round((f/lf) * 100)
            if((perc%%10 == 0) && (perc != lp)){
                cat(perc, " ", sep = "")
                lp <- perc
            } else {}
            flush.console()
            feat <- featlist[[f]]
            N <- feat$scmax - feat$scmin + 1
            peaks <- peakinfo <- NULL
            mzrange <- c(feat$mzmin, feat$mzmax)
            sccenter <- feat$scmin[1] + floor(N/2) - 1
            scrange <- c(feat$scmin, feat$scmax)
            sr <- c(max(scanrange[1], scrange[1] - max(noiserange)), min(scanrange[2], scrange[2] + max(noiserange)))
            eic <- rawEIC(object, mzrange = mzrange, scanrange = sr)
            d <- eic$intensity
            td <- sr[1]:sr[2]
            scan.range <- c(sr[1], sr[2])
            mzROI.EIC <- rawEIC(object, mzrange = mzrange, scanrange = scrange)
	    omz <- xcms:::rawMZ(object, mzrange = mzrange, scanrange = scrange)
            od <- mzROI.EIC$intensity
            otd <- mzROI.EIC$scan
            if(all(od == 0)){
                stop("centWave: debug me: (all(od == 0))?\n")
            } else {}
            ftd <- max(td[1], scrange[1] - scRangeTol):min(td[length(td)], scrange[2] + scRangeTol)
            fd <- d[match(ftd, td)]
            if(N >= 10 * minPeakWidth){ 
                noised <- rawEIC(object, mzrange = mzrange, scanrange = scanrange)$intensity
            } else {
                noised <- d
            }
            noise <- xcms:::estimateChromNoise(noised, c(0.05, 0.95), minPts = 3 * minPeakWidth)
            if(!xcms:::continuousPtsAboveThreshold(fd, threshold = noise, num = minPtsAboveBaseLine)){
                next
            } else {}
            lnoise <- xcms:::getLocalNoiseEstimate(d, td, ftd, noiserange, Nscantime)
            baseline <- max(1, min(lnoise[1], noise))
            sdnoise <- max(1, lnoise[2])
            sdthr <- sdnoise * snthresh
            if(!(any(fd - baseline >= sdthr))){
                next
            } else {}
	    singlezero <- which(aaply(1:length(d), 1, function(x,d){
    		if(d[x] == 0){
        	    if(d[max(x-1, 1)] != 0 & d[min(x+1, length(d))] != 0){
            		TRUE
        	    } else {
                        FALSE
                    }
    		} else {
                    FALSE
                }
	    }, d=d))
	    d[singlezero] <- mean(c(d[singlezero-1], d[singlezero+1]))
            wCoefs <- xcms:::MSW.cwt(d, scales = scales, wavelet = "mexh")
            if(!(!is.null(dim(wCoefs)) && any(wCoefs - baseline >= sdthr))){
                next
            } else {}
            if(td[length(td)] == Nscantime){
                wCoefs[nrow(wCoefs), ] <- wCoefs[nrow(wCoefs) - 1, ] * 0.99
            } else {}
            localMax <- xcms:::MSW.getLocalMaximumCWT(wCoefs)
            rL <- xcms:::MSW.getRidge(localMax)
            wpeaks <- sapply(rL, function(x){
                w <- min(1:length(x), ncol(wCoefs))
                any(wCoefs[x, w] - baseline >= sdthr)
            })
            if(any(wpeaks)){
                wpeaksidx <- which(wpeaks)
                for(p in 1:length(wpeaksidx)){
                    opp <- rL[[wpeaksidx[p]]]
                    pp <- unique(opp)
                    if(length(pp) >= 1){
                        dv <- td[pp] %in% ftd
                        if(any(dv)){
                            if(any(d[pp[dv]] - baseline >= sdthr)){
                                inti <- numeric(length(opp))
                                irange = rep(ceiling(scales[1]/2), length(opp))
                                for(k in 1:length(opp)){
                                    kpos <- opp[k]
                                    r1 <- ifelse(kpos - irange[k] > 1, kpos - irange[k], 1)
                                    r2 <- ifelse(kpos + irange[k] < length(d), kpos + irange[k], length(d))
                                    inti[k] <- sum(d[r1:r2])
                                }
                                maxpi <- which.max(inti)
                                if(length(maxpi) > 1){
                                    m <- wCoefs[opp[maxpi], maxpi]
                                    bestcol <- which(m == max(m), arr = T)[2]
                                    best.scale.nr <- maxpi[bestcol]
                                } else {
                                    best.scale.nr <- maxpi
                                }
                                best.scale <- scales[best.scale.nr]
                                best.scale.pos <- opp[best.scale.nr]
                                pprange <- min(pp):max(pp)
                                lwpos <- max(1, best.scale.pos - best.scale)
                                rwpos <- min(best.scale.pos + best.scale, length(td))
                                p1 <- match(td[lwpos], otd)[1]
                                p2 <- match(td[rwpos], otd)
                                p2 <- p2[length(p2)]
                                if (is.na(p1)){
                                    p1 <- 1
                                } else {}
                                if (is.na(p2)){
                                    p2 <- N
                                } else {}
                                mz.value <- omz[p1:p2]
                                if(all(mz.value == 0)){
                                    next
                                } else {}
                                mz.int <- od[p1:p2]
                                maxint <- max(mz.int)
                                mzrange <- range(mz.value[mz.value!=0])
                                mzmean <- do.call(mzCenterFun, list(mz = mz.value[mz.value!=0], intensity = mz.int[mz.value!=0]), envir=loadNamespace('xcms'))
                                dppm <- NA
                                if(verbose.columns){ 
                                    if (length(mz.value) >= (minCentroids + 1)){
                                        dppm <- round(min(xcms:::running(abs(diff(mz.value))/(mzrange[2] * 1e-06), fun = max, width = minCentroids)))
                                    } else {}
                                } else {
                                    dppm <- round((mzrange[2] - mzrange[1])/(mzrange[2] * 1e-06))
                                }
                                peaks <- rbind(peaks, c(mzmean, mzrange, NA, NA, NA, NA, NA, maxint, round((maxint - baseline)/sdnoise), NA, NA, NA, NA, f, dppm, best.scale, td[best.scale.pos], td[lwpos], td[rwpos], NA, NA, feat$MS2scan, feat$precursorInt))
                                peakinfo <- rbind(peakinfo, c(best.scale, best.scale.nr, best.scale.pos, lwpos, rwpos))
                            } else {}
                        } else {}
                    } else {}
                }
            } else {}
            if(!is.null(peaks)){
                colnames(peaks) <- c(basenames, verbosenames)
                colnames(peakinfo) <- c("scale", "scaleNr", "scpos", "scmin", "scmax")
                for(p in 1:dim(peaks)[1]){
                    if(integrate == 1){
                        lm <- xcms:::descendMin(wCoefs[, peakinfo[p, "scaleNr"]], istart = peakinfo[p, "scpos"])
                        gap <- all(d[lm[1]:lm[2]] == 0)
                        if((lm[1] == lm[2]) || gap){
                            lm <- xcms:::descendMinTol(d, startpos = c(peakinfo[p, "scmin"], peakinfo[p, "scmax"]), maxDescOutlier)
                        } else {}
                    } else {
                        lm <- descendMinTol(d, startpos = c(peakinfo[p, "scmin"], peakinfo[p, "scmax"]), maxDescOutlier)
                    }
                    pd <- d[lm[1]:lm[2]]
                    np <- length(pd)
                    lm.l <- xcms:::findEqualGreaterUnsorted(pd, 1)
                    lm.l <- max(1, lm.l - 1)
                    lm.r <- xcms:::findEqualGreaterUnsorted(rev(pd), 1)
                    lm.r <- max(1, lm.r - 1)
                    lm <- lm + c(lm.l - 1, -(lm.r - 1))
                    peakrange <- td[lm]
                    peaks[p, "rtmin"] <- scantime[peakrange[1]]
                    peaks[p, "rtmax"] <- scantime[peakrange[2]]
                    peaks[p, "maxo"] <- max(d[lm[1]:lm[2]])
                    pwid <- (scantime[peakrange[2]] - scantime[peakrange[1]])/(peakrange[2] - peakrange[1])
                    if(is.na(pwid)){
                        pwid <- 1
                    } else {}
                    peaks[p, "into"] <- pwid * sum(d[lm[1]:lm[2]])
                    db <- d[lm[1]:lm[2]] - baseline
                    peaks[p, "intb"] <- pwid * sum(db[db > 0])
                    peaks[p, "lmin"] <- lm[1]
                    peaks[p, "lmax"] <- lm[2]
                    if(fitgauss){
                        md <- max(d[lm[1]:lm[2]])
                        d1 <- d[lm[1]:lm[2]]/md
                        pgauss <- xcms:::fitGauss(td[lm[1]:lm[2]], d[lm[1]:lm[2]], pgauss = list(mu = peaks[p, "scpos"], sigma = peaks[p, "scmax"] - peaks[p, "scmin"], h = peaks[p, "maxo"]))
                        rtime <- peaks[p, "scpos"]
                        if(!any(is.na(pgauss)) && all(pgauss > 0)){
                            gtime <- td[match(round(pgauss$mu), td)]
                            if(!is.na(gtime)){
                                rtime <- gtime
                                peaks[p, "mu"] <- pgauss$mu
                                peaks[p, "sigma"] <- pgauss$sigma
                                peaks[p, "h"] <- pgauss$h
                                peaks[p, "egauss"] <- sqrt((1/length(td[lm[1]:lm[2]])) * sum(((d1 - xcms:::gauss(td[lm[1]:lm[2]], pgauss$h/md, pgauss$mu, pgauss$sigma))^2)))
                            } else {}
                        } else {}
                        peaks[p, "rt"] <- scantime[rtime]
                        if(peaks[p, "rt"] < peaks[p, "rtmin"]){
                            peaks[p, "rt"] <- scantime[peaks[p, "scpos"]]
                        } else {}
                    } else {
                        peaks[p, "rt"] <- scantime[peaks[p, "scpos"]]
                    }
                }
                peaks <- joinOverlappingPeaks(td, d, otd, omz, od, scantime, scan.range, peaks, maxGaussOverlap, mzCenterFun = mzCenterFun)
            } else {}
            if((sleep > 0) && (!is.null(peaks))){
                tdp <- scantime[td]
                trange <- range(tdp)
                egauss <- paste(round(peaks[, "egauss"], 3), collapse = ", ")
                cdppm <- paste(peaks[, "dppm"], collapse = ", ")
                csn <- paste(peaks[, "sn"], collapse = ", ")
                par(bg = "white")
                l <- layout(matrix(c(1, 2, 3), nr = 3, nc = 1, byrow = T), heights = c(0.5, 0.75, 2))
                par(mar = c(2, 4, 4, 2) + 0.1)
                plotRaw(object, mzrange = mzrange, rtrange = trange, log = TRUE, title = "")
                title(main = paste(f, ": ", round(mzrange[1], 4), " - ", round(mzrange[2], 4), " m/z , dppm=", cdppm, ", EGauss=", egauss, ",  S/N =", csn, sep = ""))
                par(mar = c(1, 4, 1, 2) + 0.1)
                image(y = scales[1:(dim(wCoefs)[2])], z = wCoefs, col = terrain.colors(256), xaxt = "n", ylab = "CWT coeff.")
                par(mar = c(4, 4, 1, 2) + 0.1)
                plot(tdp, d, ylab = "Intensity", xlab = "Scan Time")
                lines(tdp, d, lty = 2)
                lines(scantime[otd], od, lty = 2, col = "blue")
                abline(h = baseline, col = "green")
                bwh <- length(sr[1]:sr[2]) - length(baseline)
                if(xcms:::odd(bwh)){
                    bwh1 <- floor(bwh/2)
                    bwh2 <- bwh1 + 1
                } else {
                    bwh1 <- bwh2 <- bwh/2
                }
                if(any(!is.na(peaks[, "scpos"]))){
                    abline(v = scantime[na.omit(peaks[(peaks[, "scpos"] > 0), "scpos"])], col = "red")
                } else {}
                abline(v = na.omit(c(peaks[, "rtmin"], peaks[, "rtmax"])), col = "green", lwd = 1)
                if(fitgauss){
                    fitted.peaks <- which(!is.na(peaks[, "mu"]))
                    for(p in fitted.peaks){
                        yg <- xcms:::gauss(tdp, peaks[p, "h"], peaks[p, "rt"], peaks[p, "sigma"])
                        lines(tdp, yg, col = "blue")
                    }
                } else {}
                Sys.sleep(sleep)
            } else {}
            if(!is.null(peaks)){
                peaklist[[length(peaklist) + 1]] <- peaks
            } else {}
        }
        if(length(peaklist) == 0){
            cat("\nNo peaks found !\n")
            if(verbose.columns){
                nopeaks <- new("xcmsPeaks", matrix(nrow = 0, ncol = length(basenames) + length(verbosenames)))
                colnames(nopeaks) <- c(basenames, verbosenames)
            } else {
                nopeaks <- new("xcmsPeaks", matrix(nrow = 0, ncol = length(basenames)))
                colnames(nopeaks) <- c(basenames)
            }
            return(invisible(nopeaks))
        } else {}
        p <- do.call(rbind, peaklist)
        if(!verbose.columns){
            p <- p[, basenames, drop = FALSE]
        } else {}
        if(combine.peaks){
            uorder <- order(p[, "into"], decreasing = TRUE)
            pm <- as.matrix(p[, c("mzmin", "mzmax", "rtmin", "rtmax"), drop = FALSE])
            uindex <- xcms:::rectUnique(pm, uorder, mzdiff, ydiff = -1e-05)
            pr <- p[uindex, , drop = FALSE]
        } else {
            pr <- p
        }
        cat("\n", dim(pr)[1], " Peaks.\n")
        invisible(new("xcmsPeaks", pr))
    }
    .local(object, ...)
}
### retcor without peak detection
retcorRaw <- function (object, ...){
    .local <- function (filepath, runorder,
        profStep = 1, center = NULL, response = 1, 
        distFunc = "cor_opt", gapInit = NULL, gapExtend = NULL, 
        factorDiag = 2, factorGap = 1, localAlignment = 0, initPenalty = 0) 
    {
        if(is.null(gapInit)){
            if(distFunc == "cor"){
                gapInit = 0.3
            } else {}
            if(distFunc == "cor_opt"){
                gapInit = 0.3
            } else {}
            if(distFunc == "cov"){
                gapInit = 0
            } else {}
            if(distFunc == "euc"){
                gapInit = 0.9
            } else {}
            if(distFunc == "prd"){
                gapInit = 0
            } else {}
        } else {}
        if(is.null(gapExtend)){
            if(distFunc == "cor"){
                gapExtend = 2.4
            } else {}
            if(distFunc == "cor_opt"){
                gapExtend = 2.4
            } else {}
            if(distFunc == "cov"){
                gapExtend = 11.7
            } else {}
            if(distFunc == "euc"){
                gapExtend = 1.8
            } else {}
            if(distFunc == "prd"){
                gapExtend = 7.8
            } else {}
        } else {}
        N <- length(filepath)
        filenames <- basename(filepath)
        rtimecor <- vector("list", N)
        center <- which(runorder == quantile(runorder, probs=0.5, type=3))
        cat("center sample: ", filenames[center], "\nProcessing: ")
        idx <- which(seq(1, N) != center)
        obj1 <- xcmsRaw(filepath[center], profmethod = "bin", profstep = 0)
        corr <- list()
        orig <- list()
        for(si in 1:length(idx)){
            s <- idx[si]
            cat(filenames[s], " ")
            flush.console()
            xcms:::profStepPad(obj1) <- profStep
            obj2 <- xcmsRaw(filepath[s], profmethod = "bin", profstep = 0)
            xcms:::profStepPad(obj2) <- profStep
            mzmin <- min(obj1@mzrange[1], obj2@mzrange[1])
            mzmax <- max(obj1@mzrange[2], obj2@mzrange[2])
            mz <- seq(mzmin, mzmax, by = profStep)
            mz <- as.double(mz)
            mzval <- length(mz)
            scantime1 <- obj1@scantime
            scantime2 <- obj2@scantime
            mstdiff <- median(c(diff(scantime1), diff(scantime2)))
            rtup1 <- c(1:length(scantime1))
            rtup2 <- c(1:length(scantime2))
            mst1 <- which(diff(scantime1) > 5 * mstdiff)[1]
            if(!is.na(mst1)){
                rtup1 <- which(rtup1 <= mst1)
                cat("Found gaps: cut scantime-vector at ", scantime1[mst1], "seconds", "\n")
            } else {}
            mst2 <- which(diff(scantime2) > 5 * mstdiff)[1]
            if(!is.na(mst2)){
                rtup2 <- which(rtup2 <= mst2)
                cat("Found gaps: cut scantime-vector at ", scantime2[mst2], "seconds", "\n")
            } else {}
            scantime1 <- scantime1[rtup1]
            scantime2 <- scantime2[rtup2]
            rtmaxdiff <- abs(diff(c(scantime1[length(scantime1)], scantime2[length(scantime2)])))
            if(rtmaxdiff > (5 * mstdiff)){
                rtmax <- min(scantime1[length(scantime1)], scantime2[length(scantime2)])
                rtup1 <- which(scantime1 <= rtmax)
                rtup2 <- which(scantime2 <= rtmax)
            } else {}
            scantime1 <- scantime1[rtup1]
            scantime2 <- scantime2[rtup2]
            valscantime1 <- length(scantime1)
            valscantime2 <- length(scantime2)
            if(length(obj1@scantime) > valscantime1){
                obj1@env$profile <- obj1@env$profile[, -c((valscantime1 + 1):length(obj1@scantime))]
            } else {}
            if(length(obj2@scantime) > valscantime2){
                obj2@env$profile <- obj2@env$profile[, -c((valscantime2 + 1):length(obj2@scantime))]
            } else {}
            if(mzmin < obj1@mzrange[1]){
                seqlen <- length(seq(mzmin, obj1@mzrange[1], profStep)) - 1
                x <- matrix(0, seqlen, dim(obj1@env$profile)[2])
                obj1@env$profile <- rbind(x, obj1@env$profile)
            } else {}
            if(mzmax > obj1@mzrange[2]){
                seqlen <- length(seq(obj1@mzrange[2], mzmax, profStep)) - 1
                x <- matrix(0, seqlen, dim(obj1@env$profile)[2])
                obj1@env$profile <- rbind(obj1@env$profile, x)
            } else {}
            if(mzmin < obj2@mzrange[1]){
                seqlen <- length(seq(mzmin, obj2@mzrange[1], profStep)) - 1
                x <- matrix(0, seqlen, dim(obj2@env$profile)[2])
                obj2@env$profile <- rbind(x, obj2@env$profile)
            } else {}
            if(mzmax > obj2@mzrange[2]){
                seqlen <- length(seq(obj2@mzrange[2], mzmax, profStep)) - 1
                x <- matrix(0, seqlen, dim(obj2@env$profile)[2])
                obj2@env$profile <- rbind(obj2@env$profile, x)
            } else {}
            intensity1 <- obj1@env$profile
            intensity2 <- obj2@env$profile
            if((mzval * valscantime1 != length(intensity1)) || (mzval * valscantime2 != length(intensity2))){
                stop("Dimensions of profile matrices do not match !\n")
            } else {}
            rtimecor[[s]] <- .Call("R_set_from_xcms", valscantime1, 
                scantime1, mzval, mz, intensity1, valscantime2, 
                scantime2, mzval, mz, intensity2, response, distFunc, 
                gapInit, gapExtend, factorDiag, factorGap, localAlignment, 
                initPenalty, PACKAGE='xcms')
            if(length(obj2@scantime) > valscantime2){
                corr[[s]] <- c(rtimecor[[s]], obj2@scantime[(max(rtup2) + 1):length(obj2@scantime)])
            } else {
                corr[[s]] <- rtimecor[[s]]
            }
            orig[[s]] <- obj2@scantime
            rm(obj2)
            gc()
        }
        corr[[center]] <- obj1@scantime
        orig[[center]] <- obj1@scantime
        names(corr) <- filenames
        names(orig) <- filenames
        ans <- list(original=orig, corrected=corr)
        cat("\n")
        invisible(ans)
    }
    .local(object, ...)
}
#####################
joinOverlappingPeaks <- function(td, d, otd, omz, od, scantime, scan.range, peaks, maxGaussOverlap = 0.5, mzCenterFun){
    gausspeaksidx <- which(!is.na(peaks[, "mu"]))
    Ngp <- length(gausspeaksidx)
    if(Ngp == 0){
        return(peaks)
    } else {}
    newpeaks <- NULL
    gpeaks <- peaks[gausspeaksidx, , drop = FALSE]
    if(dim(peaks)[1] - Ngp > 0){
        notgausspeaks <- peaks[-gausspeaksidx, , drop = FALSE]
    } else {}
    if(Ngp > 1){
        comb <- which(upper.tri(matrix(0, Ngp, Ngp)), arr = TRUE)
        overlap <- rep(FALSE, dim(comb)[1])
        for (k in 1:dim(comb)[1]) {
            p1 <- comb[k, 1]
            p2 <- comb[k, 2]
            overlap[k] <- xcms:::gaussCoverage(xlim = scan.range, h1 = gpeaks[p1, "h"], mu1 = gpeaks[p1, "mu"], s1 = gpeaks[p1, "sigma"], h2 = gpeaks[p2, "h"], mu2 = gpeaks[p2, "mu"], s2 = gpeaks[p2, "sigma"]) >= maxGaussOverlap
        }
    } else {
        overlap <- FALSE
    }
    if(any(overlap) && (Ngp > 1)){
        jlist <- list()
        if(length(which(overlap)) > 1){
            gm <- comb[overlap, ]
            cc <- list()
            cc[[1]] <- gm[1, ]
            for(j in 2:dim(gm)[1]){
                ccl <- unlist(cc)
                nl <- sapply(cc, function(x) length(x))
                ccidx <- rep(1:length(nl), nl)
                idx <- match(gm[j, ], ccl)
                if(any(!is.na(idx))){
                  pos <- ccidx[idx[which(!is.na(idx))[1]]]
                  cc[[pos]] <- c(cc[[pos]], gm[j, ])
                } else {
                    cc[[length(cc) + 1]] <- gm[j, ]
                }
            }
            ccn <- list()
            lcc <- length(cc)
            ins <- rep(FALSE, lcc)
            if(lcc > 1){
                jcomb <- which(upper.tri(matrix(0, lcc, lcc)), arr = TRUE)
                for(j in 1:dim(jcomb)[1]){
                    j1 <- jcomb[j, 1]
                    j2 <- jcomb[j, 2]
                    if(any(cc[[j1]] %in% cc[[j2]])){
                        ccn[[length(ccn) + 1]] <- unique(c(cc[[j1]], 
                        cc[[j2]]))
                    } else {
                        if(!ins[j1]){
                            ccn[[length(ccn) + 1]] <- unique(cc[[j1]])
                            ins[j1] <- TRUE
                        } else {}
                        if(!ins[j2]){
                            ccn[[length(ccn) + 1]] <- unique(cc[[j2]])
                            ins[j2] <- TRUE
                        } else {}
                    }
                }
            } else {
                ccn <- cc
            }
            size <- sapply(ccn, function(x) length(x))
            s2idx <- which(size >= 2)
            if(length(s2idx) > 0){
                for(j in 1:length(s2idx)){
                    pgroup <- unique(ccn[[s2idx[j]]])
                    jlist[[j]] <- pgroup
                }
            } else {
                stop("(length(s2idx) = 0) ?!?")
            }
        } else {
            jlist[[1]] <- comb[overlap, ]
        }
        for(j in 1:length(jlist)){
            jidx <- jlist[[j]]
            newpeak <- gpeaks[jidx[1], , drop = FALSE]
            newmin <- min(gpeaks[jidx, "lmin"])
            newmax <- max(gpeaks[jidx, "lmax"])
            newpeak[1, "scpos"] <- -1
            newpeak[1, "scmin"] <- -1
            newpeak[1, "scmax"] <- -1
            newpeak[1, "scale"] <- -1
            newpeak[1, "maxo"] <- max(gpeaks[jidx, "maxo"])
            newpeak[1, "sn"] <- max(gpeaks[jidx, "sn"])
            newpeak[1, "lmin"] <- newmin
            newpeak[1, "lmax"] <- newmax
            newpeak[1, "rtmin"] <- scantime[td[newmin]]
            newpeak[1, "rtmax"] <- scantime[td[newmax]]
            newpeak[1, "rt"] <- weighted.mean(gpeaks[jidx, "rt"], w = gpeaks[jidx, "maxo"])
            p1 <- match(td[newmin], otd)[1]
            p2 <- match(td[newmax], otd)
            p2 <- p2[length(p2)]
            if(is.na(p1)){
                p1 <- 1
            } else {}
            if(is.na(p2)){
                p2 <- length(omz)
            } else {}
            mz.value <- omz[p1:p2]
            mz.int <- od[p1:p2]
            mz.int <- mz.int[mz.value != 0]
            mz.value <- mz.value[mz.value != 0]
            mzmean <- do.call(mzCenterFun, list(mz = mz.value, intensity = mz.int), envir=loadNamespace('xcms'))
            mzrange <- range(mz.value)
            newpeak[1, "mz"] <- mzmean
            newpeak[1, c("mzmin", "mzmax")] <- mzrange
            md <- max(d[newmin:newmax])
            d1 <- d[newmin:newmax]/md
            pgauss <- xcms:::fitGauss(td[newmin:newmax], d[newmin:newmax], pgauss = list(mu = td[newmin] + (td[newmax] - td[newmin])/2, sigma = td[newmax] - td[newmin], h = max(gpeaks[jidx, "h"])))
            if(!any(is.na(pgauss)) && all(pgauss > 0)){
                newpeak[1, "mu"] <- pgauss$mu
                newpeak[1, "sigma"] <- pgauss$sigma
                newpeak[1, "h"] <- pgauss$h
                newpeak[1, "egauss"] <- sqrt((1/length(td[newmin:newmax])) * sum(((d1 - xcms:::gauss(td[newmin:newmax], pgauss$h/md, pgauss$mu, pgauss$sigma))^2)))
            } else {
                newpeak[1, "mu"] <- NA
                newpeak[1, "sigma"] <- NA
                newpeak[1, "h"] <- NA
                newpeak[1, "egauss"] <- NA
            }
            if(is.null(newpeaks)){
                newpeaks <- newpeak
            } else {
                newpeaks <- rbind(newpeaks, newpeak)
            }
        }
        jp <- unique(unlist(jlist))
        if (dim(peaks)[1] - length(jp) > 0){
            newpeaks <- rbind(newpeaks, gpeaks[-jp, ])
        } else {}
    } else {
        newpeaks <- gpeaks
    }
    grt.min <- newpeaks[, "rtmin"]
    grt.max <- newpeaks[, "rtmax"]
    if(nrow(peaks) - Ngp > 0){
        for(k in 1:nrow(notgausspeaks)){
            if(!any((notgausspeaks[k, "rtmin"] >= grt.min) & (notgausspeaks[k, "rtmax"] <= grt.max))){
                newpeaks <- rbind(newpeaks, notgausspeaks[k, ])
            } else {}
        }
    } else {}
    rownames(newpeaks) <- NULL
    newpeaks
}
### evaluate parameters
evalPar <- function(directory, par, nSample=5){
	if(missing(directory)){
		if(Sys.info()["sysname"] == 'Windows'){
			readline('\nChoose directory of mzXML data: <Press Return>')
			directory <- choose.dir()
		} else {
			directory <- readline('Path to directory with mzXML data:')
		}
	} else {}
	mzXML <- list.files(directory, pattern='*.mzXML', full.names=TRUE, ignore.case=TRUE)
	mzXML <- sample(mzXML, nSample)
	cat('\nSamples to analyse\n\n')
	cat(paste(mzXML, collapse='\n'))
	cat('\n\n')
	flush.console()
	param <- parameters(par)
	nROI <- list()
	for(i in 1:length(mzXML)){
		raw <- xcmsRaw(mzXML[i])
		scanrange <- c(1, length(raw@scantime))
		dev <- param@parameters$findPeak$ppm * 1e-06
		scalerange <- round((param@parameters$findPeak$peakwidth/mean(diff(raw@scantime)))/2)
		if(length(scalerange) > 1){
			scales <- seq(from=scalerange[1], to=scalerange[2], by=2)
		} else {
			scales <- scalerange
		}
		minPeakWidth <- scales[1]
		minCentroids <- max(4, minPeakWidth-2)
		featlist <- xcms:::findmzROI(raw, scanrange=scanrange, dev=dev, minCentroids=minCentroids, prefilter=param@parameters$findPeak$prefilter, noise=0)
		nROI[[i]] <- length(featlist)
	}
	nROI <- do.call('c', nROI)
	p <- qplot('Samples', nROI, geom=c('boxplot', 'point')) + theme_bw() + scale_x_discrete('') + scale_y_continuous('Number of ROIs')
	print(p)
	cat('\nMean: ')
	cat(mean(nROI))
	cat('\n')
}