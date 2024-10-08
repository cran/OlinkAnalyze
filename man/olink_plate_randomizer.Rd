% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Olink_plate_randomizer.R
\name{olink_plate_randomizer}
\alias{olink_plate_randomizer}
\title{Randomly assign samples to plates}
\usage{
olink_plate_randomizer(
  Manifest,
  PlateSize = 96,
  Product,
  SubjectColumn,
  iterations = 500,
  available.spots,
  num_ctrl = 8,
  rand_ctrl = FALSE,
  seed
)
}
\arguments{
\item{Manifest}{tibble/data frame in long format containing all sample ID's. Sample ID column must be named SampleID.}

\item{PlateSize}{Integer. Either 96 or 48. 96 is default.}

\item{Product}{String. Name of Olink product used to set PlateSize if not provided. Optional.}

\item{SubjectColumn}{(Optional) Column name of the subject ID column. Cannot contain missing values. If provided, subjects are kept on the same plate. This argument is used for longitudinal studies and must be a separate column from the SampleID column.}

\item{iterations}{Number of iterations for fitting subjects on the same plate.}

\item{available.spots}{Numeric. Number of wells available on each plate. Maximum 40 for T48 and 88 for T96. Takes a vector equal to the number of plates to be used indicating the number of wells available on each plate.}

\item{num_ctrl}{Numeric. Number of controls on each plate (default = 8)}

\item{rand_ctrl}{Logical. Whether controls are added to be randomized across the plate (default = FALSE)}

\item{seed}{Seed to set. Highly recommend setting this for reproducibility.}
}
\value{
A "tibble" including SampleID, SubjectID etc. assigned to well positions.
Columns include same columns as Manifest with additional columns:
\itemize{
\item{plate:} Plate number
\item{column:} Column on the plate
\item{row:} Row on the plate
\item{well:} Well location on the plate
}
}
\description{
Generates a scheme for how to plate samples with an option to keep subjects on the same plate.
}
\details{
Variables of interest should if possible be randomized across plates to avoid confounding with potential plate effects. In the case of multiple samples per subject (e.g. in longitudinal studies), Olink recommends keeping each subject on the same plate. This can be achieved using the SubjectColumn argument.
}
\examples{
\donttest{
#Generate randomization scheme using complete randomization
randomized.manifest_a <- olink_plate_randomizer(manifest, seed=12345)

#Generate randomization scheme that keeps subjects on the same plate (for longitudinal studies)
randomized.manifest_b <- olink_plate_randomizer(manifest,SubjectColumn="SubjectID",
                                                        available.spots=c(88,88), seed=12345)

#Visualize the generated plate layouts
olink_displayPlateLayout(randomized.manifest_a, fill.color = 'Site')
olink_displayPlateLayout(randomized.manifest_a, fill.color = 'SubjectID')
olink_displayPlateLayout(randomized.manifest_b, fill.color = 'Site')
olink_displayPlateLayout(randomized.manifest_b, fill.color = 'SubjectID')

#Validate that sites are properly randomized
olink_displayPlateDistributions(randomized.manifest_a, fill.color = 'Site')
olink_displayPlateDistributions(randomized.manifest_b, fill.color = 'Site')
}


}
\seealso{
\itemize{
\item{\code{\link[OlinkAnalyze:olink_displayPlateLayout]{olink_displayPlateLayout()}} for visualizing the generated plate layouts}
\item{\code{\link[OlinkAnalyze:olink_displayPlateDistributions]{olink_displayPlateDistributions()}} for validating that sites are properly randomized}
}
}
\keyword{plates}
\keyword{randomized}
