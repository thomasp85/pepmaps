\name{Parameters-class}
\Rdversion{1.1}
\docType{class}
\alias{Parameters-class}
\alias{show,Parameters-method}

\title{Class \code{"Parameters"}}
\description{This class holds information about the analysis parameters used by \code{createComplist}.}
\section{Objects from the Class}{
Objects are created by calls to the function \code{parameters()}. Objects can also be created by calls of the form \code{new("Parameters", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{Name}:}{String giving the name of the parameter set}
    \item{\code{Type}:}{String giving the type of parameter set (only LC-MS/MS supported at the moment).}
    \item{\code{parameters}:}{List with elements corresponding to the specific type of parameter set, containing the parameters.}
  }
}
\author{
Thomas Dybdal Pedersen
}
\seealso{
\code{\link{parameters}}, \code{\link{createComplist}}
}
