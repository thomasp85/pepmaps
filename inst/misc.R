### Collection of different scripts

### Plot TIC overlay for all .mzXML in a folder
files <- list.files(choose.dir(), pattern='*.mzXML', ignore.case=TRUE, full.names=TRUE)
ans <- list()
for(i in 1:length(files)){
    raw <- xcmsRaw(files[i])
    TIC <- data.frame(Intensity=raw@tic, Time=raw@scantime, Sample=basename(files[i]))
    ans[[i]] <- TIC
}
ans <- do.call('rbind', ans)
qplot(Time, Intensity, data=ans, group=Sample, geom='line')

### Evaluate MassAI samples
samples <- as.character(unique(pID@raw$File))
ans <- data.frame(Sample=basename(as.character(samples)), N.peptides=NA, N.unique=NA, Lev2ID=NA)
collect <- as.character(unique(pID@raw$Peptide))
for(i in 1:length(samples)){
    sub <- subset(pID@raw, as.character(pID@raw$File) == samples[i])
    other <- subset(pID@raw, as.character(pID@raw$File) != samples[i])
    ans$N.peptides[i] <- length(unique(sub$Peptide))
    ans$N.unique[i] <- sum(!(unique(sub$Peptide) %in% other$Peptide))
    ans$Lev2ID[i] <- sum(sub$Notes == 'Lev2')
    collect <- collect[which(collect %in% as.character(sub$Peptide))]
}
### convert factors in data frame to character
factostr <- function(x){
    if(is.factor(x)){
        x <- as.character(x)
    } else {}
    x
}
ans <- apply(test, 2, factostr)
### combine eic list from getEIC to one data frame
ans <- list()
for(i in 1: length(eic@eic)){
    eics <- eic@eic[[i]]
    eics <- adply(seq(along=eics), 1, function(x,ID) if(nrow(ID[[x]])!=0) {data.frame(Peak.ID=x, ID[[x]])} else {}, ID=eics)[,-1]
    eics <- data.frame(eics, Sample=names(eic@eic)[i])
    ans[[i]] <- eics
}
ans <- do.call('rbind', ans)
### get desired sequence
x <- 11503
y <- TRUE
z <- 1
while(y){
    y <- x > z
    z <- z*10
}
l <- 11
z <- z/1000
y <- 0
if(length(seq(0 , x, by=z)) < 10){
    ans <- z
} else {
    while(l > 10){
        y <- y+2
        l <- length(seq(0, x, by=z*y))
        ans <- z*y
    }
}

ms <- list.files(choose.dir(), full.names=TRUE)
xset <- xcmsSet(ms, method='centWave', ppm=25, peakwidth=c(5,25), snthresh=5, prefilter=c(1,1000))
xset <- retcor(xset, method='obiwarp', response=10, localAlignment=1)
xset <- group(xset, bw=1.5, minfrac=0.1)
xset <- fillPeaks(xset)
xset <- annotate(xset, polarity='positive')