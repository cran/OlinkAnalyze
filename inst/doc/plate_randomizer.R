## ---- echo=FALSE,message=FALSE------------------------------------------------
library(OlinkAnalyze)
library(dplyr)
library(ggplot2)

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, message = FALSE) 

## ----Outlier_data, include=FALSE----------------------------------------------
library(OlinkAnalyze)
data(manifest)

## ----manifest, echo = F, message=FALSE----------------------------------------
 manifest |> head(10)  %>% 
  knitr::kable(format = "html",
               booktabs = T,
               linesep = "",
               digits = 4,
               longtable = T,
               row.names = F,
               caption = " Sample manifest.",
               label = "manifest")

## ----complete_randomization, message=FALSE, warning=FALSE,results='hide'------
randomized.manifest_a <- olink_plate_randomizer(manifest, seed=123456)

## ----subjects_randomization, message=FALSE,results='hide'---------------------
randomized.manifest_b <- olink_plate_randomizer(manifest,SubjectColumn="SubjectID", 
                                                available.spots=c(48,48,42), 
                                                iterations = 500,
                                                seed=123456)

## ---- fig.height = 8, fig.width = 8, fig.align = "center"---------------------
olink_displayPlateLayout(randomized.manifest_b, fill.color = 'SubjectID', include.label = FALSE)

## ----fig.height = 8, fig.width = 8, fig.align = "center"----------------------
olink_displayPlateLayout(randomized.manifest_b, fill.color = 'SubjectID', include.label = TRUE)

## ----fig.height = 8, fig.width = 8, fig.align = "center"----------------------
olink_displayPlateDistributions(randomized.manifest_b, fill.color = 'SubjectID')

## ----fig.height = 8, fig.width = 8, fig.align = "center"----------------------
olink_displayPlateDistributions(randomized.manifest_b, fill.color = 'Site')

## ---- eval = FALSE------------------------------------------------------------
#  library(writexl)
#  write_xlsx(randomized.manifest_b,"randomized.manifest_b.xlsx")

