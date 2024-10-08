% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Read_NPX_data.R
\name{read_NPX}
\alias{read_NPX}
\title{Function to read NPX data into long format}
\usage{
read_NPX(filename)
}
\arguments{
\item{filename}{Path to Olink Software output file.}
}
\value{
A "tibble" in long format. Columns include:
\itemize{
\item{SampleID:} Sample ID
\item{Index:} Index
\item{OlinkID:} Olink ID
\item{UniProt:} UniProt ID
\item{Assay:} Protein symbol
\item{MissingFreq:} Proportion of sample below LOD
\item{Panel_Version:} Panel Version
\item{PlateID:} Plate ID
\item{QC_Warning:} QC Warning Status
\item{LOD:} Limit of detection
\item{NPX:} Normalized Protein Expression
}
Additional columns may be present or missing depending on the platform
}
\description{
Imports an NPX or QUANT file exported from Olink Software.
No alterations to the output format is allowed.
}
\examples{
\donttest{
file <- system.file("extdata", "Example_NPX_Data.csv", package = "OlinkAnalyze")
read_NPX(file)
}
}
\keyword{NPX}
