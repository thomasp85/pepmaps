\name{pepID}
\alias{pepID}
\title{
Constructor For PepID Objects
}
\description{
This function handles the creation of a \code{PepID} object by automatic formating the result from either MassAI or Crosswork.
}
\usage{
pepID(type, path = file.choose(), sep = "\t", dec = ".")
}
\arguments{
  \item{type}{
Either 'MassAI' or 'Crosswork' to indicate the formating necessary
}
  \item{path}{
The location of the result file. By default the location can be navigated to.
}
  \item{sep}{
The seperator used in the result file.
}
  \item{dec}{
The decimal seperator in the result file
}
}
\value{
A \code{PepID} object.
}
\author{
Thomas Dybdal Pedersen
}
\seealso{
\code{\linkS4class{PepID}}
}