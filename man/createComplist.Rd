\name{createComplist}
\alias{createComplist}
\title{
Wrapper Function To Create A Complist Object from Scratch
}
\description{
This function guides the user through the creation of a \code{\linkS4class{Complist}} object.
It handles all \pkg{xcms} and \pkg{CAMERA} function calls as well as the \pkg{pepmaps} calls.
}
\usage{
createComplist(mzXML, Sample.info, ID = "MassAI", par = "standard", sep = "\t", dec = ".", annotate = TRUE)
}
\arguments{
  \item{mzXML}{
The path to the folder containing the .mzXML files used for the analysis. If missing the function lets you choose a folder.
}
  \item{Sample.info}{
A data frame containing information about the samples arranged in alphabetical order to match the order of the .mzXML files. (optional)
}
  \item{ID}{
The type of identification carried out (only 'MassAI' supported at the moment).
}
  \item{par}{
The name of the parameter set used for the analysis.
}
  \item{sep}{
The seperator used in the ID file.
}
  \item{dec}{
The decimal seperator used in the ID file
}
  \item{annotate}{
Logical. Should CAMERA annotation be carried out.
}
}
\details{
This function conveniently handles all function calls necessary to create a \code{Complist} object. In it's barebone form only \code{'par'} should be specified 
(except when using default values) making the function call: \code{ans <- createComplist(par='name')} sufficient. Under these circumstances the function directs
the user to choosing the location of the .mzXML and ID files and outputs a \code{Complist} object.
}
\value{
A \code{Complist} object.
}
\author{
Thomas Dybdal Pedersen
}
\seealso{
\code{\listS4class{Complist}}, \code{\link{complist}}, \code{\link{peaklist}}, \code{\link{pepID}}, \code{\linkS4class{Parameters}}
}