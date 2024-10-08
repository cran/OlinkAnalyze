% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/olink_normalization.R
\name{norm_internal_reference_median}
\alias{norm_internal_reference_median}
\title{Internal reference median normalization function}
\usage{
norm_internal_reference_median(
  ref_df,
  ref_samples,
  ref_name,
  ref_cols,
  reference_medians
)
}
\arguments{
\item{ref_df}{The reference dataset to be used in normalization (required).}

\item{ref_samples}{Character vector of sample identifiers to be used for
adjustment factor calculation in the reference dataset (required).}

\item{ref_name}{Project name of the reference dataset (required).}

\item{ref_cols}{Named list of column names in the reference dataset
(required).}

\item{reference_medians}{Dataset with columns "OlinkID" and "Reference_NPX"
(required). Used for reference median normalization.}
}
\value{
Tibble or ArrowObject with the normalized dataset.
}
\description{
Internal reference median normalization function
}
\author{
Klev Diamanti
}
