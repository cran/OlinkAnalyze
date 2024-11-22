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

## ----echo=FALSE, eval = TRUE--------------------------------------------------
library(OlinkAnalyze)
library(dplyr)
library(stringr)
library(ggplot2)
library(kableExtra)

## ----brnrtab, eval=TRUE, message=FALSE, echo = FALSE--------------------------
data.frame(Platform = c("Target 96",
                        paste0("Explore 384: \n",
                               "Cardiometabolic, Inflammation, ",
                               "Neurology, and Oncology"),
                        paste0("Explore 384: \n",
                               "Cardiometabolic II, Inflammation II,",
                               "Neurology II, and Oncology II"),
                        "Explore HT",
                        "Explore 3072 to Explore HT"),
           BridgingSamples = c("8-16",
                               "8-16",
                               "16-24",
                               "16-32",
                               "40-64")) |>
  kbl(booktabs = TRUE,
      digits = 2,
      caption = "Recommended number of bridging samples for Olink platforms") |>
  kable_styling(bootstrap_options = "striped",
                full_width = FALSE,
                position = "center",
                latex_options = "HOLD_position")


## ----fig.cap= fcap, eval = TRUE, echo = FALSE, out.width="50%"----------------
knitr::include_graphics(normalizePath("../man/figures/Bridging_schematic.png"),
                        error = FALSE)
fcap <- "Schematic of Explore 3072 to Explore HT Bridging Workflow"

## ----message=FALSE, eval=FALSE, echo = TRUE-----------------------------------
#  # Note: Explore 3072 NPX files can be CSV or parquet.
#  data_explore3072 <- read_NPX("~/NPX_Explore3072_location.parquet")
#  data_exploreht <- read_NPX("~/NPX_ExploreHT_location.parquet")

## ----echo=TRUE, eval = FALSE--------------------------------------------------
#  data_explore3072_samples <- data_explore3072 |>
#    dplyr::filter(SampleType == "SAMPLE") |>
#    dplyr::distinct(SampleID) |>
#    dplyr::pull()
#  
#  data_exploreht_samples <- data_exploreht |>
#    dplyr::filter(SampleType == "SAMPLE") |>
#    dplyr::distinct(SampleID) |>
#    dplyr::pull()
#  
#  overlapping_samples <- unique(intersect(data_explore3072_samples,
#                                          data_exploreht_samples))
#  # Note that if `SampleType` is not is input data:
#  # stringr::str_detect can be used to exclude control samples based on SampleID.

## ----echo=FALSE---------------------------------------------------------------
try(
  readRDS(normalizePath("../man/figures/overlapping_samples_table.rds")) |> 
    kableExtra::kbl(booktabs = TRUE,
                    digits = 2,
                    caption = "Overlapping bridging samples") |>
    kableExtra::kable_styling(bootstrap_options = "striped",
                              full_width = FALSE,
                              position = "center",
                              latex_options = "HOLD_position")
)

## ----include=FALSE------------------------------------------------------------
f3 <- paste0("PCA plot prior to bridging for Explore 3072 and Explore HT data.",
             " Bridge samples are indicated by color.",
             " PCA plots can be helpful in assessing",
             " if any bridge samples were outliers in one of the platforms.")

## ----eval = FALSE-------------------------------------------------------------
#  #### Extract bridging samples
#  
#  data_explore3072_before_br <- data_explore3072 |>
#    dplyr::filter(SampleType == "SAMPLE") |>
#    # Note that if `SampleType` is not is input data,
#    # stringr::str_detect can be used to exclude control samples
#    #  based on naming convention.
#    dplyr::mutate(Type = if_else(SampleID %in% overlapping_samples,
#                                 paste0("Explore 3072 Bridge"),
#                                 paste0("Explore 3072 Sample")))
#  
#  data_exploreht_before_br <- data_exploreht |>
#    dplyr::filter(SampleType == "SAMPLE") |>
#    dplyr::mutate(Type = if_else(SampleID %in% overlapping_samples,
#                                 paste0("Explore HT Bridge"),
#                                 paste0("Explore HT Sample")))
#  
#  ### PCA plot
#  pca_E3072 <- OlinkAnalyze::olink_pca_plot(df = data_explore3072_before_br,
#                                           color_g = "Type",
#                                           quiet = TRUE)
#  pca_EHT <- OlinkAnalyze::olink_pca_plot(df = data_exploreht_before_br,
#                                          color_g = "Type",
#                                          quiet = TRUE)

## ----echo=FALSE, fig.cap=f3, fig.height= 8, fig.width= 6----------------------
knitr::include_graphics(normalizePath("../man/figures/PCA_btw_product_before.png"),
                        error = FALSE)


## ----eval = FALSE-------------------------------------------------------------
#  # Find shared samples
#  npx_ht <- data_exploreht |>
#    dplyr::mutate(Project = "data1")
#  npx_3072 <- data_explore3072 |>
#    dplyr::mutate(Project = "data2")
#  
#  npx_br_data <- olink_normalization(df1 = npx_ht,
#                                     df2 = npx_3072,
#                                     overlapping_samples_df1 =
#                                       overlapping_samples,
#                                     df1_project_nr = "Explore HT",
#                                     df2_project_nr = "Explore 3072",
#                                     reference_project = "Explore HT")

## ----echo=FALSE, fig.cap=fcap, out.width="50%"--------------------------------
knitr::include_graphics(normalizePath("../man/figures/assay_bridgeability.jpg"), 
                        error = FALSE)
fcap <- paste("Criteria to determine the bridging recommendation for an assay.",
"The assessment of linearity ensures bridging between signal in both platforms",
"or noise in both platforms (but not between signal and noise).",
"Similar NPX ranges and sufficient counts provide additional insight into",
"an assay's bridgeability.",
"Distribution shape is assessed to determine recommended bridging method.", 
sep = " ")

## ----echo=FALSE---------------------------------------------------------------
try( 
  readRDS(normalizePath("../man/figures/bridging_results.rds")) |> 
    kableExtra::kbl(booktabs = TRUE,
        digits = 1,
        caption = "Table 4. First 5 rows of combined datasets after bridging.") |>
    kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE, font_size = 10, 
                  position = "center", latex_options = "HOLD_position") |> 
    kableExtra::scroll_box(width = "100%")
)

## ----include=FALSE------------------------------------------------------------
f8 <- "Combined PCA of sample controls from both platforms prior to normalization."
f9 <- "Combined PCA of bridging samples from both platforms prior to normalization."
f10 <- "Combined PCA of sample controls from both platforms after normalization."
f11 <- "Combined PCA of bridging samples from both platforms after normalization."

## ----pca_pre_sc, echo=TRUE, eval = FALSE--------------------------------------
#  ## Before Bridging
#  npx_br_data |>
#    dplyr::filter(SampleType == "SAMPLE_CONTROL") |>
#    dplyr::mutate(OlinkID = paste0(OlinkID, "_", OlinkID_E3072)) |>
#    dplyr:::mutate(SampleID = paste0(Project, SampleID)) |>
#    OlinkAnalyze::olink_pca_plot(color_g = "Project")

## ----pca_pre_sc_fig, echo=FALSE, fig.cap=f8, message=FALSE--------------------
## Before Bridging
knitr::include_graphics(normalizePath("../man/figures/SCs_pre_bridging.png"), 
                        error = FALSE)

## ----pca_pre_bridge, echo=TRUE, eval=FALSE------------------------------------
#  ## Before Bridging
#  npx_br_data |>
#    dplyr::filter(SampleType == "SAMPLE") |>
#    dplyr::filter(SampleID %in% overlapping_samples) |>
#    dplyr::mutate(OlinkID = paste0(OlinkID, "_", OlinkID_E3072)) |>
#    dplyr:::mutate(SampleID = paste0(Project, SampleID)) |>
#    OlinkAnalyze::olink_pca_plot(color_g = "Project")
#  
#  

## ----echo=FALSE, fig.cap=f9---------------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/bridges_pre_bridging.png"),
                        error = FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  ## After bridging PCA
#  
#  ### Keep the data following BridgingRecommendation
#  npx_after_br_reco <- npx_br_data |>
#    dplyr::filter(BridgingRecommendation != "NotBridgeable") |>
#    dplyr::mutate(NPX = case_when(
#      BridgingRecommendation == "MedianCentering" ~ MedianCenteredNPX,
#      BridgingRecommendation == "QuantileSmoothing" ~ QSNormalizedNPX,
#      .default = NPX)) |>
#    dplyr::filter(AssayType == "assay") |>
#    dplyr::mutate(OlinkID = paste0(OlinkID, "_", OlinkID_E3072))
#  

## ----pca_post_SC, eval=FALSE, echo = TRUE-------------------------------------
#  
#  ### Generate unique SampleIDs
#  npx_after_br_final <- npx_after_br_reco |>
#    dplyr:::mutate(SampleID = paste0(Project, SampleID))
#  
#  ### PCA plot of the data from SCs
#  npx_after_br_final |>
#      dplyr::filter(SampleType == "SAMPLE_CONTROL") |>
#      OlinkAnalyze::olink_pca_plot(color_g = "Project")
#  

## ----echo=FALSE, fig.cap=f10--------------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/SCs_post_bridging.png"),
                        error = FALSE)

## ----echo=TRUE, eval = FALSE--------------------------------------------------
#  ### PCA plot of the data from bridging samples
#  npx_after_br_reco |>
#    dplyr::filter(SampleType == "SAMPLE") |>
#    dplyr::filter(SampleID %in% overlapping_samples) |>
#    dplyr:::mutate(SampleID = paste0(Project, SampleID)) |>
#    OlinkAnalyze::olink_pca_plot(color_g = "Project")

## ----echo=FALSE, fig.cap=f11--------------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/bridges_post_bridging.png"), 
                        error = FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  df <- npx_br_data |>
#      dplyr::filter(Project == "Explore_3072") |>
#      arrow::as_arrow_table()
#  
#  df$metadata$FileVersion <- "NA"
#  df$metadata$ExploreVersion <- "NA"
#  df$metadata$ProjectName <- "NA"
#  df$metadata$SampleMatrix <- "NA"
#  df$metadata$DataFileType <- "Olink Analyze Export File"
#  df$metadata$ProductType <- "Explore3072"
#  df$metadata$Product <- "Explore3072"
#  arrow::write_parquet(x = df, sink = "path_to_output.parquet")

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  # Option 1: Exclude non bridgeable assays from both products
#  npx_recommended <- npx_after_br_final |>
#    dplyr::mutate(NPX_original = NPX) |>
#    dplyr::filter(BridgingRecommendation != "NotBridgeable") |>
#    dplyr::mutate(NPX = case_when(
#      BridgingRecommendation == "MedianCentering" ~ MedianCenteredNPX,
#      BridgingRecommendation == "QuantileSmoothing" ~ QSNormalizedNPX,
#      .default = NPX)) |>
#    dplyr::mutate(OlinkID_HT = OlinkID) |>
#    dplyr::mutate(OlinkID = paste0(OlinkID, "_", OlinkID_E3072))
#  
#  # Option 2: Analyze non bridgeable assays separately
#  npx_recommended <- npx_after_br_final |>
#    dplyr::mutate(NPX_original = NPX) |>
#    dplyr::mutate(NPX = case_when(
#      BridgingRecommendation == "MedianCentering" ~ MedianCenteredNPX,
#      BridgingRecommendation == "QuantileSmoothing" ~ QSNormalizedNPX,
#      .default = NPX)) |>
#    dplyr::mutate(OlinkID_HT = OlinkID) |>
#    dplyr::mutate(OlinkID = ifelse(BridgingRecommendation != "NotBridgeable",
#                                   paste0(OlinkID, "_", OlinkID_E3072),
#                                   # Concatenated OlinkID for bridgeable Assays
#                                   ifelse(Project == "Explore HT",
#                                          # replace with HT project name as set in function
#                                          OlinkID,
#                                          OlinkID_E3072))

