## ----Outlier_data, include=FALSE----------------------------------------------
library(OlinkAnalyze)
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
knitr::opts_chunk$set(fig.pos = 'H')
knitr::opts_chunk$set(fig.align="center")
knitr::opts_knit$set(eval.after = "fig.cap")
# Create Dataset with outliers
outlier_data <- npx_data1 |> 
  dplyr::mutate(NPX = ifelse(SampleID == "A25", NPX + 4, NPX)) |> 
  dplyr::mutate(NPX = ifelse(SampleID == "A52", NPX - 4, NPX)) |> 
  dplyr::filter(!stringr::str_detect(SampleID, "CONTROL"))

group_data <- npx_data1 |> 
  dplyr::mutate(NPX = ifelse(Site == "Site_D", NPX + 3, NPX)) |> 
  dplyr::filter(!stringr::str_detect(SampleID, "CONTROL"))

 

## ----dataset_generation, eval = FALSE, message=FALSE, warning=FALSE-----------
# # Create Datasets with outliers
# outlier_data <- npx_data1 |>
#   dplyr::mutate(NPX = ifelse(SampleID == "A25", NPX + 4, NPX)) |>
#   dplyr::mutate(NPX = ifelse(SampleID == "A52", NPX - 4, NPX)) |>
#   dplyr::filter(!stringr::str_detect(SampleID, "CONTROL"))
# 
# group_data <- npx_data1 |>
#   dplyr::mutate(NPX = ifelse(Site == "Site_D", NPX + 3, NPX)) |>
#   dplyr::filter(!stringr::str_detect(SampleID, "CONTROL"))
# 

## ----Outlier_example_code, eval = FALSE---------------------------------------
# p1<- outlier_data |> olink_pca_plot(label_samples = T, quiet = T)
# p2<- group_data |> olink_pca_plot(color_g = "Site", quiet = T)
# ggpubr::ggarrange(p1[[1]], p2[[1]], nrow = 2, labels = "AUTO")

## ----Outlier_Example, echo=FALSE, fig.cap=fcap--------------------------------
knitr::include_graphics(normalizePath("../man/figures/PCA_Outlier_Fig1.png"),error = FALSE)
fcap <- "**Figure 1** **A.** PCA can be used to identify individual outlier samples as shown by samples A25 and A52. **B.** PCA can be used to identify difference in groups as seen by Site_D samples. This is not suggesting Site_D is an outlier, but rather that there may be a global difference between sites."

## ----PCA_treatment, eval=FALSE------------------------------------------------
# OlinkAnalyze::npx_data1 |>
#   dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter duplicate SampleIDs
#   olink_pca_plot(color_g = "Treatment")

## ----PCA_Panel, eval=FALSE----------------------------------------------------
# OlinkAnalyze::npx_data2 |>
#   dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter out control SampleIDs
#   olink_pca_plot(byPanel = TRUE) # Specify by panel

## ----PCA_Panel_fig, echo=FALSE------------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/PCA_Panel.png"), error = FALSE)

## ----PCA_object, eval = FALSE-------------------------------------------------
# pca_plots<-OlinkAnalyze::npx_data2|> # Save the PCA plot to a variable
#   dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter duplicate SampleIDs
#   olink_pca_plot(byPanel = TRUE, quiet = TRUE) # By panel
# # quiet argument suppresses export
# 
# pca_plots[[1]] #Cardiometabolic PCA
# pca_plots[[2]] #Inflammation PCA

## ----Outlier_PCA, echo = FALSE------------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/Outlier_PCA.png"), error = FALSE)

## ----eval = FALSE-------------------------------------------------------------
# outlier_data |>
#   dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter duplicate SampleIDs
#   olink_pca_plot(label_samples = TRUE)

## ----Label_samples, echo = FALSE----------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/label_samples_pca.png"), error = FALSE)

## ----eval = FALSE-------------------------------------------------------------
# outlier_data |>
#   dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter duplicate SampleIDs
#   olink_pca_plot(outlierDefX = 3, outlierDefY = 3,
#                  outlierLines = TRUE, label_outliers = TRUE)

## ----outlier_line_pca, echo = FALSE-------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/outlier_line_pca.png"), error = FALSE)

## ----eval = FALSE-------------------------------------------------------------
# outlier_data |>
#   dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter duplicate SampleIDs
#   olink_pca_plot(outlierDefX = 3, outlierDefY = 3,
#                  outlierLines = FALSE, label_outliers = TRUE)

## -----------------------------------------------------------------------------
outliers_pca_labeled <- outlier_data |> 
  dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter duplicate SampleIDs
  olink_pca_plot(outlierDefX = 3, outlierDefY = 3, outlierLines = FALSE, 
                 label_outliers = TRUE, quiet = TRUE) 

outliers_pca_labeled[[1]]$data |> 
  dplyr::filter(Outlier == TRUE) |> 
  dplyr::select(SampleID) |> 
  dplyr::distinct()


## ----eval = FALSE-------------------------------------------------------------
# outlier_data |>
#   dplyr::filter(SampleID %in% c("A25", "A52", "A1", "A2", "A3", "A5", "A15", "A16", "A18", "A19", "A20"))|>
#   olink_dist_plot()

## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/dist_boxplot.png"), error = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# group_data |>
#   dplyr::filter(Site %in% c("Site_A", "Site_D")) |> # Only visualizing 2 sites to see all samples
#   olink_dist_plot(color_g = "Site")

## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/site_boxplot.png"), error = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# # Calculate SampleID Median NPX
# median_NPX<-group_data |>
#   dplyr::group_by(SampleID) |>
#   dplyr::summarise(Median_NPX = median(NPX))
# 
# # Adjust by sample median
# adjusted_data <- group_data |>
#   dplyr::inner_join(median_NPX, by = "SampleID")|>
#   dplyr::mutate(NPX = NPX - Median_NPX)
# 
# adjusted_data|>
#   dplyr::filter(Site %in% c("Site_A", "Site_D")) |> # Only visualizing 2 sites to see all samples
#   olink_dist_plot(color_g = "Site")

## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/sample_med_boxplot.png"), error = FALSE)

## ----eval = FALSE-------------------------------------------------------------
# outlier_data |>
#   olink_qc_plot()

## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/qc_plot.png"), error = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# group_data |>
#   olink_qc_plot(color_g = "Site")

## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/qc_site_plot.png"), error = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# outlier_data |>
#   olink_qc_plot(median_outlierDef = 2, IQR_outlierDef = 4,
#                  outlierLines = TRUE, label_outliers = TRUE)

## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/qc_label_plot.png"), error = FALSE)

## ----eval = FALSE-------------------------------------------------------------
# outlier_data |>
#   olink_qc_plot(median_outlierDef = 2, IQR_outlierDef = 4,
#                  outlierLines = FALSE, label_outliers = TRUE)

## -----------------------------------------------------------------------------
outliers_qc_labeled <- outlier_data |> 
  olink_qc_plot(median_outlierDef = 2, IQR_outlierDef = 4, 
                 outlierLines = FALSE, label_outliers = TRUE) 

outliers_qc_labeled$data |> 
  dplyr::filter(Outlier == TRUE) |> 
  dplyr::select(SampleID) |> 
  dplyr::distinct()


