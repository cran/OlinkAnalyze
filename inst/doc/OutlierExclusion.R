## ----Outlier_data, include=FALSE----------------------------------------------
library(OlinkAnalyze)
knitr::opts_chunk$set(fig.height=3,fig.width=6)
knitr::opts_chunk$set(warning = F)
knitr::opts_chunk$set(message = F)
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
#  # Create Datasets with outliers
#  outlier_data <- npx_data1 |>
#    dplyr::mutate(NPX = ifelse(SampleID == "A25", NPX + 4, NPX)) |>
#    dplyr::mutate(NPX = ifelse(SampleID == "A52", NPX - 4, NPX)) |>
#    dplyr::filter(!stringr::str_detect(SampleID, "CONTROL"))
#  
#  group_data <- npx_data1 |>
#    dplyr::mutate(NPX = ifelse(Site == "Site_D", NPX + 3, NPX)) |>
#    dplyr::filter(!stringr::str_detect(SampleID, "CONTROL"))
#  

## ----Outlier_Example, fig.cap=fcap--------------------------------------------
p1<- outlier_data |> olink_pca_plot(label_samples = T, quiet = T)
p2<- group_data |> olink_pca_plot(color_g = "Site", quiet = T)
ggpubr::ggarrange(p1[[1]], p2[[1]], nrow = 2, labels = "AUTO")
fcap <- "**Figure 1** **A.** PCA can be used to identify individual outlier samples as shown by samples A25 and A52. **B.** PCA can be used to identify difference in groups as seen by Site_D samples. This is not suggesting Site_D is an outlier, but rather that there may be a global difference between sites."

## -----------------------------------------------------------------------------
OlinkAnalyze::npx_data1 |> 
  dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter duplicate SampleIDs
  olink_pca_plot(color_g = "Treatment")

## -----------------------------------------------------------------------------
OlinkAnalyze::npx_data2 |> 
  dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter out control SampleIDs
  olink_pca_plot(byPanel = TRUE) # Specify by panel

## -----------------------------------------------------------------------------
pca_plots<-OlinkAnalyze::npx_data2|> # Save the PCA plot to a variable
  dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter duplicate SampleIDs
  olink_pca_plot(byPanel = TRUE, quiet = TRUE) # By panel, do not print ggarranged plot 
pca_plots[[1]] #Cardiometabolic PCA
pca_plots[[2]] #Inflammation PCA

## -----------------------------------------------------------------------------
outlier_data |> 
  dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter duplicate SampleIDs
  olink_pca_plot(byPanel = TRUE) 

## -----------------------------------------------------------------------------
outlier_data |> 
  dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter duplicate SampleIDs
  olink_pca_plot(label_samples = TRUE) 

## -----------------------------------------------------------------------------
outlier_data |> 
  dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter duplicate SampleIDs
  olink_pca_plot(outlierDefX = 3, outlierDefY = 3, 
                 outlierLines = TRUE, label_outliers = TRUE) 

## -----------------------------------------------------------------------------
outlier_data |> 
  dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter duplicate SampleIDs
  olink_pca_plot(outlierDefX = 3, outlierDefY = 3, 
                 outlierLines = FALSE, label_outliers = TRUE) 

## -----------------------------------------------------------------------------
outliers_pca_labeled <- outlier_data |> 
  dplyr::filter(stringr::str_detect(SampleID, "CONTROL", negate = T)) |> # Filter duplicate SampleIDs
  olink_pca_plot(outlierDefX = 3, outlierDefY = 3, outlierLines = FALSE, 
                 label_outliers = TRUE, quiet = TRUE) 

outliers_pca_labeled[[1]]$data |> 
  dplyr::filter(Outlier == TRUE) |> 
  dplyr::select(SampleID) |> 
  dplyr::distinct()


## -----------------------------------------------------------------------------
outlier_data |> 
  dplyr::filter(SampleID %in% c("A25", "A52", "A1", "A2", "A3", "A5", "A15", "A16", "A18", "A19", "A20"))|> 
  olink_dist_plot()

## -----------------------------------------------------------------------------
group_data |> 
  dplyr::filter(Site != "Site_E") |> # Site E filtered out so that all samples can be seen
  olink_dist_plot(color_g = "Site")

## -----------------------------------------------------------------------------
# Calculate SampleID Median NPX
median_NPX<-group_data |> 
  dplyr::group_by(SampleID) |> 
  dplyr::summarise(Median_NPX = median(NPX)) 

# Adjust by sample median
adjusted_data <- group_data |> 
  dplyr::inner_join(median_NPX, by = "SampleID")|> 
  dplyr::mutate(NPX = NPX - Median_NPX) 

adjusted_data|> 
  dplyr::filter(Site != "Site_E") |>  #Filter out Site E to see all the samples
  olink_dist_plot(color_g = "Site")

## -----------------------------------------------------------------------------
outlier_data |> 
  olink_qc_plot()

## -----------------------------------------------------------------------------
group_data |> 
  olink_qc_plot(color_g = "Site")

## -----------------------------------------------------------------------------
outlier_data |> 
  olink_qc_plot(median_outlierDef = 2, IQR_outlierDef = 4, 
                 outlierLines = TRUE, label_outliers = TRUE) 

## -----------------------------------------------------------------------------
outlier_data |> 
  olink_qc_plot(median_outlierDef = 2, IQR_outlierDef = 4, 
                 outlierLines = FALSE, label_outliers = TRUE) 

## -----------------------------------------------------------------------------
outliers_qc_labeled <- outlier_data |> 
  olink_qc_plot(median_outlierDef = 2, IQR_outlierDef = 4, 
                 outlierLines = FALSE, label_outliers = TRUE) 

outliers_qc_labeled$data |> 
  dplyr::filter(Outlier == TRUE) |> 
  dplyr::select(SampleID) |> 
  dplyr::distinct()


