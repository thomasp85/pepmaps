\name{Complist-class}
\Rdversion{1.1}
\docType{class}
\alias{Complist-class}
\alias{getIntensity,Complist-method}
\alias{plotPeak,Complist-method}
\alias{show,Complist-method}

\title{Class \code{"Complist"}}
\description{This class combines \code{Peaklist} object with different types of ID lists (only \code{\linkS4class{PepID}} objects currently supported), matching ID to the quantified peakgroups}
\section{Objects from the Class}{Objects can be created using the \code{\link{complist}} constructor which automatically links ID's to peakgroups.
Furthermore the wrapper function \code{\link{createComplist}} will also output a \code{Complist} object.
Objects can also be created by calls of the form \code{new("Complist", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{peakID}:}{A logical data frame giving the peakgroups that have been matched to ID}
    \item{\code{IDindex}:}{A list giving the mathc between ID and peakgroups}
    \item{\code{pepID}:}{The \code{PepID} object used}
    \item{\code{raw}:}{A matrix containing the intensities of each sample for all peakgroups}
    \item{\code{Sample.info}:}{A data frame containing information on each sample}
    \item{\code{Peak.info}:}{A data frame containing chromatographic information for each peak group}
    \item{\code{annotated}:}{TRUE/FALSE giving whether CAMERA analysis has been performed on data}
    \item{\code{parameters}:}{A \code{Parameters} object containing all different analysis parameters used}
  }
}
\section{Extends}{
Class \code{"\linkS4class{Peaklist}"}, directly.
}
\section{Methods}{
  \describe{
    \item{\link{getIntensity}}{\code{signature(object = "Complist")}: Extract intensity for all ID'ed peak groups or a specific ID}
    \item{\link{plotPeak}}{\code{signature(object = "Complist")}: Plot EIC for a given ID}
	 }
}
\author{Thomas Dybdal Pedersen}
\seealso{\code{\link{complist}}}