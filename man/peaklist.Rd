\name{peaklist}
\alias{peaklist}
\title{
Constructor For Peaklist Objects
}
\description{
This function constructs a \code{Peaklist} object by taking either a \code{xcmsSet}
or a \code{xAnnotate} object, a data frame giving additional sample information and
a \code{Parameters} object containing the parameters used in the analysis.
}
\usage{
peaklist(data, Sample.info, parameters)
}
\arguments{
  \item{data}{
Either a \code{xcmsSet} or a \code{xAnnotate} object with analysed mz data.
}
  \item{Sample.info}{
A data frame containing additional information about the samples. The sample order should correspond to the order in the data object.
}
  \item{parameters}{
A \code{Parameters} object that holds information on the analytical parameters.
}
}
\value{
Returns a \code{Peaklist} object.
}
\author{
Thomas Dybdal Pedersen
}
\seealso{
\code{\linkS4class{Peaklist}}, \code{\linkS4class{Parameters}}, \code{\linkS4class{xcmsSet}}, \code{\linkS4class{xAnnotate}}
}