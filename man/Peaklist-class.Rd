\name{Peaklist-class}
\Rdversion{1.1}
\docType{class}
\alias{Peaklist-class}
\alias{show,Peaklist-method}

\title{Class \code{"Peaklist"}}
\description{
This class contains information from either \code{peakTable()} or \code{getPeaklist()} from \pkg{xcms} and \pkg{CAMERA} respectively.
}
\section{Objects from the Class}{
Objects are created by the constructor \code{peaklist()}. Objects can also be created by calls of the form \code{new("Peaklist", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{raw}:}{A matrix containing the raw intensities of each peak group for all samples.}
    \item{\code{Sample.info}:}{A data frame containing information about the samples.}
    \item{\code{Peak.info}:}{A data frame containing information on each peak group}
    \item{\code{annotated}:}{Logical indicating whether CAMERA analysis has been carried out.}
    \item{\code{parameters}:}{A Parameters object with the analysis parameters used.}
  }
}
\author{
Thomas Dybdal Pedersen
}
\seealso{
\code{\link{peakTable}}, \code{\link{getPeaklist}}
}