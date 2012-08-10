\name{PepID-class}
\Rdversion{1.1}
\docType{class}
\alias{PepID-class}
\alias{initialize,PepID-method}
\alias{show,PepID-method}

\title{Class \code{"PepID"}}
\description{
This class contains output from peptide indentfication by either Crosswork or MassAI.
}
\section{Objects from the Class}{
Objects can be created using the \code{pepID} constructor. Objects can also be created by calls of the form \code{new("PepID", type, raw)}.
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{A character string giving the type of object(either 'Crosswork' or 'MassAI')}
    \item{\code{raw}:}{A data frame containing the raw results from the peptide identification}
    \item{\code{peplist}:}{A data frame containing an overview of each unique identified peptide}
    \item{\code{npep}:}{The number of unique identified peptides}
    \item{\code{nsample}:}{The number of samples used in the analysis}
  }
}
\section{Methods}{
  \describe{
    \item{evalID}{\code{signature(object = "PepID")}: Bootstrap estimation of the number of peptides identified as a function of the number of samples used.}
	 }
}
\references{
Download \href{http://proteomics.dk/}{MassAI or Crosswork}
}
\author{
Thomas Dybdal Pedersen
}
\seealso{
\code{\link{pepID}}
}