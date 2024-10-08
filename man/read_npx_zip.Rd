% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_npx_zip.R
\name{read_npx_zip}
\alias{read_npx_zip}
\title{Helper function to read in Olink Explore zip csv files}
\usage{
read_npx_zip(filename)
}
\arguments{
\item{filename}{Path to Olink Software output zip file.}
}
\value{
A "tibble" in long format. Some of the columns are:
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
Helper function to read in Olink Explore zip csv files
}
\examples{
\donttest{
try({ # May fail if dependencies are not installed
file <- system.file("extdata", "Example_NPX_Data.csv", package = "OlinkAnalyze")
read_NPX(file)
})

}

}
\keyword{NPX}
