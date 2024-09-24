## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy = FALSE,
  tidy.opts = list(width.cutoff = 95),
  fig.width = 6,
  fig.height = 3,
  message = FALSE,
  warning = FALSE,
  time_it = TRUE,
  fig.align = "center"
)

library(OlinkAnalyze)
library(dplyr)


## ----echo=FALSE---------------------------------------------------------------
set.seed(1234)
table1<-npx_data1 |> 
  head() |> 
  dplyr::select(-c(Index, MissingFreq, Panel_Version, QC_Warning, Subject, Treatment, Site, Time, Project, Panel, PlateID)) |>
  dplyr::mutate(Count = round(NPX * (100+sample(seq(-5,15), size = 1)))) |> 
  dplyr::mutate(SampleType = "SAMPLE") |> 
  dplyr::mutate(Normalization = "Plate control") |> 
  dplyr::mutate(NPX = round(NPX,digits = 2)) |> 
  dplyr::mutate(LOD = round(LOD, digits = 2)) |> 
  dplyr::mutate(PCNormalizedNPX = NPX) |> 
  dplyr::mutate(PCNormalizedLOD = LOD) |> 
  dplyr::select(SampleID, SampleType, OlinkID, UniProt, Assay, Count, NPX,  PCNormalizedNPX,Normalization, LOD, PCNormalizedLOD) 

table1 |> 
  knitr::kable(caption = "Example results from Plate Control Normalized Project") |> 
  kableExtra::kable_styling(font_size = 10)

table1 |>  
  dplyr::mutate(Normalization = "Intensity") |> 
  dplyr::mutate(NPX = round(NPX + 4.16,digits = 2)) |> 
  dplyr::mutate(LOD = round(LOD + 4.16, digits = 2)) |> 
  dplyr::select(SampleID, SampleType, OlinkID, UniProt, Assay, Count, NPX,  PCNormalizedNPX,Normalization, LOD, PCNormalizedLOD) |> 
  knitr::kable(caption = "Example results from Intensity Normalized Project") |> 
  kableExtra::kable_styling(font_size = 10)


## ----dataset_generation, eval = FALSE, message=FALSE, warning=FALSE-----------
#  explore_npx <- read_NPX("~/Explore_NPX_file.parquet")

## ----NCLOD_example, eval = FALSE, message=FALSE, warning=FALSE----------------
#  # Integrating negative control LOD for intensity normalized data
#  explore_npx <- read_NPX("Path_to/Explore_NPX_file.parquet")
#  olink_lod(explore_npx, lod_method = "NCLOD")

## ----FixedLOD, eval = FALSE, message=FALSE, warning=FALSE---------------------
#  # Reading in Fixed LOD file path into R environment
#  fixedLOD_filepath <- "Path_to/ExploreHT_fixedLOD.csv"
#  
#  # Integrating Fixed LOD for intensity normalized data
#  explore_npx <- read_NPX("~/Explore_NPX_file.parquet")
#  olink_lod(explore_npx, lod_file_path = fixedLOD_filepath, lod_method = "FixedLOD")

## ----echo=FALSE---------------------------------------------------------------
table1 |>  
  dplyr::mutate(Normalization = "Intensity") |> 
  dplyr::mutate(PCNormalizedNPX = round(NPX,digits = 2)) |> 
  dplyr::mutate(PCNormalizedLOD = round(LOD, digits = 2)) |> 
  dplyr::mutate(NPX = round(NPX + 4.16, digits = 2))|> 
  dplyr::mutate(LOD = round(LOD + 4.16, digits = 2)) |>
  dplyr::rename(FixedLOD = LOD,
                FixedPCNormalizedLOD = PCNormalizedLOD) |> 
  dplyr::mutate(NCLOD = FixedLOD - 2.34,
                NCPCNormalizedLOD = FixedPCNormalizedLOD - 2.34) |> 
  dplyr::select(SampleID, SampleType, OlinkID, UniProt, Assay, Count, NPX, Normalization, PCNormalizedNPX, FixedLOD, FixedPCNormalizedLOD, NCLOD, NCPCNormalizedLOD) |> 
  knitr::kable(caption = "Example results using both LOD calculation methods") |> 
  kableExtra::kable_styling(font_size = 10)

## ----explore_npx_export, eval = FALSE, message=FALSE, warning=FALSE-----------
#  # Exporting Olink Explore data with LOD information as a parquet file
#  explore_npx <- read_NPX("Path_to/Explore_NPX_file.parquet")
#  
#  explore_npx_NC_LOD <- explore_npx %>%
#    olink_lod(lod_method = "NCLOD")
#  
#  # Add metadata for export
#  df <- explore_npx_NC_LOD |>
#    arrow::as_arrow_table()
#  
#  df$metadata$FileVersion <- "NA"
#  df$metadata$ExploreVersion <- "NA"
#  df$metadata$ProjectName <- "NA"
#  df$metadata$SampleMatrix <- "NA"
#  df$metadata$DataFileType <- "Olink Analyze Export File"
#  df$metadata$ProductType <- "ExploreHT" # "ExploreHT" or "Explore3072"
#  df$metadata$Product <- "ExploreHT" # "ExploreHT" or "Explore3072"
#  
#  arrow::write_parquet(x = df, sink = "path_to_output.parquet")

