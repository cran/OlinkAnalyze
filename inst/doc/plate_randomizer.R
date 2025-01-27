## ----echo=FALSE,message=FALSE, eval=FALSE-------------------------------------
# library(OlinkAnalyze)

## ----include=FALSE------------------------------------------------------------
library(OlinkAnalyze)
data(manifest)

## ----echo=FALSE,message=FALSE-------------------------------------------------
library(dplyr)
library(ggplot2)

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(warning = FALSE,
                      fig.height=3,
                      fig.width=6,
                      message = FALSE) 

## ----manifest, echo = F, message=FALSE----------------------------------------
 manifest %>% head(10)  %>% 
  knitr::kable(format = "html",
               booktabs = T,
               linesep = "",
               digits = 4,
               longtable = T,
               row.names = F,
               caption = "Table 1. Example of Sample manifest",
               label = "manifest")

## ----complete_randomization, message=FALSE, warning=FALSE,results='hide'------
randomized.manifest <- olink_plate_randomizer(manifest, seed=123456)

## ----subjects_randomization, message=FALSE,results='hide'---------------------
randomized.manifest <- olink_plate_randomizer(manifest,
                                                SubjectColumn="SubjectID", 
                                                available.spots=c(48,48,42), 
                                                iterations = 500,
                                                seed=123456)

## ----randomize_controls, message = FALSE, eval = FALSE------------------------
# randomized_manifest <- olink_plate_randomizer(manifest,
#                                               Product = "Explore HT",
#                                               SubjectColumn = "SampleID",
#                                               num_ctrl = 10,
#                                               rand_ctrl = TRUE)

## ----rand_ctrl_code, echo=FALSE, results='hide', message=FALSE, warning=FALSE----
randomized_manifest <- olink_plate_randomizer(manifest, 
                                              Product = "Explore HT",
                                              SubjectColumn = "SampleID",
                                              num_ctrl = 10, 
                                              rand_ctrl = TRUE)

## ----randomized_controls_table, echo=FALSE, message=FALSE, warning=FALSE------
randomized_manifest %>% 
  head(10)  %>% 
  knitr::kable(format = "html",
               booktabs = T,
               linesep = "",
               digits = 4,
               longtable = T,
               row.names = F,
               caption = "Table 2. Example of Randomized Sample manifest",
               label = "manifest")


## -----------------------------------------------------------------------------
randomized.manifest_c <- olink_plate_randomizer(manifest, study = "Site")
olink_displayPlateLayout(randomized.manifest_c, fill.color = "Site")

## ----figure_captions, message=FALSE, echo=FALSE-------------------------------
fcap1 <- "Figure 1. Randomized samples in 96 well plate format, colored by Subject ID."
fcap2 <- "Figure 2. Randomized samples in 96 well plate format with labeled wells."
fcap3 <- "Figure 3. Distribution of Subject ID across randomized plates."
fcap4 <- "Figure 4. Distribution of Site across randomized plates."


## ----fig.height = 6, fig.width = 7.5, fig.align = "center", fig.cap= fcap1----
olink_displayPlateLayout(randomized.manifest, fill.color = 'SubjectID', include.label = FALSE)

## ----fig.height = 6, fig.width = 7.5, fig.align = "center", fig.cap= fcap2----
olink_displayPlateLayout(randomized.manifest, fill.color = 'SubjectID', include.label = TRUE)

## ----fig.align = "center", fig.cap=fcap3--------------------------------------
olink_displayPlateDistributions(randomized.manifest, fill.color = 'SubjectID')

## ----fig.align = "center", fig.cap= fcap4-------------------------------------
olink_displayPlateDistributions(randomized.manifest, fill.color = 'Site')

## ----eval = FALSE-------------------------------------------------------------
# library(writexl)
# write_xlsx(randomized.manifest,"randomized.manifest.xlsx")

