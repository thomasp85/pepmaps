\name{parameters}
\alias{parameters}
\title{
Constructor For \code{Parameters} Objects.
}
\description{
This function creates a \code{Parameters} object by looking for parameters sets matching the name of the name given in the function call.
}
\usage{
parameters(Name)
}
\arguments{
  \item{Name}{String with the name of the set to create. Must match an already specified parameter set.}
}
\details{
This function looks for sets found in \\library\\pepmaps\\extdata\\run_parameters.R that match the name given in the function call and convert this set to a \code{Parameters} object.
}
\value{
A \code{Parameters} object.
}
\author{
Thomas Dybdal Pedersen
}
\note{
Currently new parameter sets have to be added manually. Add, edit and delete functions for new sets are in the work.
}
\seealso{
\code{\linkS4class{Parameters}}
}