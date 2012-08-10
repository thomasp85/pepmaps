\name{testPar}
\alias{testPar}
\title{Test Parameters For Peak Detection}
\description{
This function presents plots for evaluating the quality of the parameters used by \code{xcmsSet} from \pkg{xcms}.
}
\usage{
testPar(xcmsRaw = file.choose(), ppm = 25, peakwidth = c(5, 25), snthresh = 5, prefilter = c(3, 100))
}
\arguments{
  \item{xcmsRaw}{An .mzXML file representing the data the parameters should be used for.}
  \item{ppm}{The sensitivity in the mz direction.}
  \item{peakwidth}{A vector with minumum and maximum peak width in seconds.}
  \item{snthresh}{The signal to noise threshold.}
  \item{prefilter}{Prefiltering parameters.}
}
\details{
This function runs \code{xcmsRaw} and \code{findPeaks.centWave} and returns a 5*5 grid plot of 25 randomly picked peaks from the run overlayed with their signal to noise.
Additionally it creates .png files with plots of chromatograms for each 100 mz slice of the dataset with detected peak location. The .png files are writen to the working directory.\cr
\cr
\bold{ppm}\cr
The \sQuote{ppm} parameter controls the width of the window used for finding regions of interest in the mz data.
This parameter is not linked to the manufacturers ppm claims and should usually be quite a bit higher (around 25 for an LC-qTOF).\cr
\cr
\bold{snthresh}\cr
The signal to noise is calculated by \eqn{(maxo - baseline)/sd}, where \verb{maxo} is the maximum intensity of the peak,
\verb{baseline} is the mean of the data without the 10 \% highest and lowest data points and \verb{sd} is the standard deviation of the data
without the 10 \% highest and lowest data points. This value should be kept low (i.e. 5) as it is not very robust and exclusion of real peaks
are likely to occur at higher values.\cr
\cr
\bold{prefilter}\cr
The prefilter is a vector \code{c(k, I)} giving the number of peaks (k) with intensity above I that a given region of interest should contain
for it to be kept for further analysis.
}
\value{
Creates .png files in the working directory of name: mzXXX-YYY.png with XXX and YYY being the mz window of the plot.\cr
\cr
Additionally a 5*5 grid plot is produced showing 25 randomly picked peaks and a 5*5 matrix giving the index of the peaks is written.
}
\author{Thomas Dybdal Pedersen \email{thomasp85@gmail.com}}
\note{
A high number of false positives is unavoidable using current peak detection methods.
}
\seealso{\code{\link{xcmsSet}, \code{\link{findPeaks.centWave}}}
