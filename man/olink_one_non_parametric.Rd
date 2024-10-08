% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Olink_one_non_parametric.R
\name{olink_one_non_parametric}
\alias{olink_one_non_parametric}
\title{Function which performs a Kruskal-Wallis Test or Friedman Test per protein}
\usage{
olink_one_non_parametric(
  df,
  variable,
  dependence = FALSE,
  subject = NULL,
  verbose = TRUE
)
}
\arguments{
\item{df}{NPX or Quantified_value data frame in long format with at least protein name (Assay), OlinkID, UniProt, Panel and a factor with at least 3 levels.}

\item{variable}{Single character value.}

\item{dependence}{Boolean. Default: FALSE. When the groups are independent, the kruskal-Wallis will run, when the groups are dependent, the Friedman test will run.}

\item{subject}{Group information for the repeated measurement. If (dependence = TRUE), this parameter need to be specified.}

\item{verbose}{Boolean. Default: True. If information about removed samples, factor conversion and final model formula is to be printed to the console.}
}
\value{
A tibble containing the Kruskal-Wallis Test or Friedman Test results for every protein.

Columns include:
\itemize{
\item{Assay:} "character" Protein symbol
\item{OlinkID:} "character" Olink specific ID
\item{UniProt:} "character" UniProt ID
\item{Panel:} "character" Name of Olink Panel
\item{term:} "character" term in model
\item{df:} "numeric" degrees of freedom
\item{method:} "character" which method was used
\item{statistic:} "named numeric" the value of the test statistic with a name describing it
\item{p.value:} "numeric" p-value for the test
\item{Adjusted_pval:} "numeric" adjusted p-value for the test (Benjamini&Hochberg)
\item{Threshold:} "character" if adjusted p-value is significant or not (< 0.05)
}
}
\description{
Performs an Kruskal-Wallis Test for each assay (by OlinkID) in every panel using stats::kruskal.test.
Performs an Friedman Test for each assay (by OlinkID) in every panel using rstatix::friedman_test.
The function handles factor variable. \cr\cr
Samples that have no variable information or missing factor levels are automatically removed from the analysis (specified in a message if verbose = TRUE).
Character columns in the input dataframe are automatically converted to factors (specified in a message if verbose = T).
Numerical variables are not converted to factors.
If a numerical variable is to be used as a factor, this conversion needs to be done on the dataframe before the function call. \cr\cr
Inference is specified in a message if verbose = TRUE. \cr
The formula notation of the final model is specified in a message if verbose = TRUE. \cr\cr
Adjusted p-values are calculated by stats::p.adjust according to the Benjamini & Hochberg (1995) method (“fdr”).
The threshold is determined by logic evaluation of Adjusted_pval < 0.05.
}
\examples{
\donttest{

library(dplyr)

# One-way Kruskal-Wallis Test
try({ # May fail if dependencies are not installed
kruskal_results <- olink_one_non_parametric(df = npx_data1,
                                            variable = "Site")
})

#Friedman Test
friedman_results <- olink_one_non_parametric(df = npx_data1,
                                             variable = "Time",
                                             subject = "Subject",
                                             dependence = TRUE)}

}
