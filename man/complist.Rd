\name{complist}
\alias{complist}
\title{Merge Peaklist And Compound ID}
\description{This function creates a \code{\linkS4class{Complist}} object by automatically matching a \code{\linkS4class{Peaklist}} object to a \code{\linkS4class{PepID}} object.
A wider matching window can optionally be defined.
}
\usage{
complist(Peaklist, PepID, rtwin = 0, mzwin = 0)
}
\arguments{
  \item{Peaklist}{The \code{\linkS4class{Peaklist}} object containing peak groups from a xcms analysis}
  \item{PepID}{The \code{\linkS4class{PepID}} object with identified peptides}
  \item{rtwin}{Extension of the matching window in the retention time direction}
  \item{mzwin}{Extension of the matching window in the mz direction}
}
\details{
The constructor compares the mean retention time and mz values of the ID object to the minimum and maximum retention time and mz values for the peak groups. 
One ID can be matched to several peak groups and one peak group can be matched to several ID's if the complexity of the sample is higher than the resolving power of the LC-MS.
}
\value{
A \code{\linkS4class{Complist}} object.
}
\author{
Thomas Dybdal Pedersen
}
\seealso{
\code{\linkS4class{Complist}}, \code{\link{getIntensity}}
}