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

## ----echo=FALSE---------------------------------------------------------------
library(OlinkAnalyze)
library(dplyr)
library(stringr)
library(ggplot2)
library(kableExtra)

## ----brnrtab, message=FALSE, echo=FALSE---------------------------------------
data.frame(Platform = c("Target 96",
                         "Explore 384 Cardiometabolic, Inflammation, Neurology, and Oncology",
                        "Explore 384 Cardiometabolic II, Inflammation II, Neurology II, and Oncology II"),
           BridgingSamples = c("8-16",
                               "8-16",
                              "16-24")) %>%
  kbl(booktabs = TRUE,
      digits = 2,
      caption = "Table 1. Recommended number of bridging samples for Olink platforms") %>%
  kable_styling(bootstrap_options = "striped",
                full_width = FALSE,
                position = "center",
                latex_options = "HOLD_position")


## ----bridge_sample_selection_example, echo=T, eval=T--------------------------
bridge_Samples<-olink_bridgeselector(df = npx_data1,
                     sampleMissingFreq = 0.1,
                     n = 16)

## ----bridge_sample_selection, echo=FALSE--------------------------------------
olink_bridgeselector(df = npx_data1,
                     sampleMissingFreq = 0.1,
                     n = 16) %>%
  kableExtra::kbl(booktabs = TRUE,
      digits = 2,
      caption = "Table 2. Selected Bridging Samples") %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE,
                position = "center", latex_options = "HOLD_position")
  

## ----captions, echo=FALSE-----------------------------------------------------
f1<- "Figure 1. PCA plot of bridging samples and other samples in npx_data1. Control samples are excluded from the PCA plot."

f2 <- "Figure 2. Density plot of NPX distribution in both datasets before bridging."

f3 <- "Figure 3. PCA plot of both datasets before bridging."

f4 <- "Figure 4. Density plot of NPX distribution in both datasets after bridging."

f5 <- "Figure 5. Histogram of adjustment factors in normalized data from Project \"data2\"."

f6 <- "Figure 6. Violin plot of CHL1 in both datasets prior to bridging. Bridge samples are indicated by black points."

f7 <- "Figure 7. Density plot of inter-project CV before and after bridging."

f8 <- "Figure 8. PCA plot of both datasets after bridging."

## ----bridge_sample_selection_example_pca, echo=T, eval=T, fig.cap= f1---------
npx_data1 %>% 
  filter(!str_detect(SampleID, 'CONT')) %>%
  mutate(Bridge = ifelse(SampleID %in% bridge_Samples$SampleID, "Bridge", "Sample")) %>% 
  olink_pca_plot(color_g = "Bridge")


## ----message=FALSE, eval=FALSE, echo = TRUE-----------------------------------
#  data1 <- read_NPX("~/NPX_file1_location.xlsx")
#  data2 <- read_NPX("~/NPX_file2_location.xlsx")

## ----eval= FALSE--------------------------------------------------------------
#  data.frame(SampleID = intersect(npx_data1$SampleID, npx_data2$SampleID)) %>%
#    dplyr::filter(!stringr::str_detect(SampleID, "CONTROL_SAMPLE"))

## ----echo=FALSE---------------------------------------------------------------

data.frame(SampleID = intersect(npx_data1$SampleID, npx_data2$SampleID)) %>%
  dplyr::filter(!stringr::str_detect(SampleID, "CONTROL_SAMPLE")) %>% #Remove control samples
  kableExtra::kbl(booktabs = TRUE,
      digits = 2,
      caption = "Table 3. Overlapping Bridge Samples") %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE,
                position = "center", latex_options = "HOLD_position")

## ----check, message=FALSE, fig.cap= f2----------------------------------------
# Load datasets
npx_1 <- npx_data1 %>%
  mutate(Project = "data1")
npx_2 <- npx_data2 %>%
  mutate(Project = "data2")

npx_df <- bind_rows(npx_1, npx_2)


# Plot NPX density before bridging normalization
npx_df %>%
  mutate(Panel = gsub("Olink ", "", Panel)) %>%
  ggplot(aes(x = NPX, fill = Project)) +
  geom_density(alpha = 0.4) +
  facet_grid(~Panel) +
  olink_fill_discrete(coloroption = c("red", "darkblue")) +
  set_plot_theme() +
  ggtitle("Before bridging normalization: NPX distribution") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "top")

## ----pca1, message=FALSE, fig.cap=f3------------------------------------------
#### Extract bridging samples

overlapping_samples <-  data.frame(SampleID = intersect(npx_1$SampleID, npx_2$SampleID)) %>%  
  filter(!str_detect(SampleID, "CONTROL_SAMPLE")) %>% #Remove control samples
  pull(SampleID)

npx_before_br <- npx_data1 %>%
  dplyr::filter(!str_detect(SampleID, "CONTROL_SAMPLE")) %>% #Remove control samples
  dplyr::mutate(Type = if_else(SampleID %in% overlapping_samples,
                        paste0("data1 Bridge"),
                        paste0("data1 Sample"))) %>%
  rbind({
    npx_data2 %>%
      filter(!str_detect(SampleID, "CONTROL_SAMPLE")) %>% #Remove control samples %>% 
      mutate(Type = if_else(SampleID %in% overlapping_samples,
                            paste0("data2 Bridge"),
                            paste0("data2 Sample"))) %>%
      mutate(SampleID = if_else(SampleID %in% overlapping_samples,
                                paste0(SampleID, "_new"),
                                SampleID))
  })

### PCA plot
OlinkAnalyze::olink_pca_plot(df          = npx_before_br,
                             color_g     = "Type",
                             byPanel     = TRUE)

## ----bridging, message=FALSE--------------------------------------------------
# Find shared samples
npx_1 <- npx_data1 %>%
  mutate(Project = "data1") 
npx_2 <- npx_data2 %>%
  mutate(Project = "data2")

overlap_samples <-data.frame(SampleID = intersect(npx_1$SampleID, npx_2$SampleID)) %>%
  filter(!str_detect(SampleID, "CONTROL_SAMPLE")) %>% #Remove control samples
  pull(SampleID)

overlap_samples_list <- list("DF1" = overlap_samples,
                             "DF2" = overlap_samples)

# Perform Bridging normalization
npx_br_data <- olink_normalization_bridge(project_1_df = npx_1,
                                          project_2_df = npx_2,
                                          bridge_samples = overlap_samples_list,
                                          project_1_name = "data1",
                                          project_2_name = "data2",
                                          project_ref_name = "data1")



## ----norm_data_table , echo = FALSE-------------------------------------------
npx_br_data %>% 
  head(10) %>% 
  kableExtra::kbl(booktabs = TRUE,
      digits = 1,
      caption = "Table 4. First 10 rows of combined datasets after bridging.") %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE, font_size = 10, 
                position = "center", latex_options = "HOLD_position") %>% 
  kableExtra::scroll_box(width = "100%")

## ----densitybr, message=FALSE, fig.cap= f4------------------------------------
# Plot NPX density after bridging normalization

npx_br_data %>%
  mutate(Panel = gsub("Olink ", "", Panel)) %>%
  ggplot2::ggplot(ggplot2::aes(x = NPX, fill = Project)) +
  ggplot2::geom_density(alpha = 0.4) +
  ggplot2::facet_grid(~Panel) +
  olink_fill_discrete(coloroption = c("red", "darkblue")) +
  set_plot_theme() +
  ggplot2::ggtitle("After bridging normalization: NPX distribution") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "top")


## ----dist_adj_fct, message=FALSE, echo = FALSE, fig.cap = f5------------------
npx_br_data %>% 
  dplyr::filter(Project == "data2") %>%  # Only looking at Project 2 since project 1 is unadjusted
  dplyr::select(OlinkID, Adj_factor) %>% 
  dplyr::distinct() %>% 
  ggplot2::ggplot(ggplot2::aes(x = Adj_factor)) +
    ggplot2::geom_histogram() +
  set_plot_theme()


## ----voilin plot, fig.cap = f6------------------------------------------------
# Bridge sample data
bridge_samples <- npx_1 %>%
  rbind(npx_2) %>%
  filter(SampleID %in% overlapping_samples) %>%
  filter(Assay == "CHL1") %>%
  mutate(Assay_OID = paste(Assay, OlinkID, sep = "\n"))

# Generate violin plot for CHL1
npx_data1 %>%
  mutate(Project = "data1") %>%
  bind_rows({
    npx_data2 %>%
      mutate(Project = "data2")
  }) %>%
  filter(Assay == "CHL1") %>%
  filter(!str_detect(SampleID, "CONTROL*.")) %>%
  mutate(Assay_OID = paste(Assay, OlinkID, sep = "\n")) %>% 
  ggplot2::ggplot(aes(Project, NPX)) +
  ggplot2::geom_violin(aes(fill = Project)) +
  geom_point(data = bridge_samples, position = position_jitter(0.1)) +
  theme(legend.position = "none") +
  set_plot_theme() +
  facet_wrap(. ~ Assay_OID, scales='free_y')

## ----CV_calculation, fig.cap= f7----------------------------------------------

explore_cv <- function(npx, na.rm = F) {
  sqrt(exp((log(2) * sd(npx, na.rm = na.rm))^2) - 1)*100
}

t96_cv <- function(NPX, na.rm = T) {
  100*sd(2^NPX)/mean(2^NPX)
}

tech <- "Target"

cv_before <- npx_1 %>% 
  rbind(npx_2) %>% 
  filter(str_detect(SampleID,"CONTROL*.")) %>%
  filter(NPX > LOD) %>%
  group_by(OlinkID) %>%
  mutate(CV = ifelse(tech=='Explore',explore_cv(NPX), t96_cv(NPX))) %>%
  ungroup() %>%
  distinct(OlinkID,CV)


cv_after <- npx_br_data %>%
  filter(str_detect(SampleID, "CONTROL")) %>%
  filter(NPX > LOD) %>%
  group_by(OlinkID) %>%
  mutate(CV = ifelse(tech=='Explore',explore_cv(NPX), t96_cv(NPX))) %>%
  ungroup() %>%
  distinct(OlinkID,CV)

cv_before %>%
  mutate(Analysis = "Before") %>%
  rbind((cv_after %>%
           mutate(Analysis = "After"))) %>%
  ggplot2::ggplot(ggplot2::aes(x = CV, fill = Analysis)) +
  ggplot2::geom_density(alpha = 0.7) +
  set_plot_theme() +
  olink_fill_discrete()+
  ggplot2::theme(text = ggplot2::element_text(size = 20)) + ggplot2::xlim(-50,400)



## ----pca2, message=FALSE, fig.cap= f8-----------------------------------------
## After bridging

### Generate unique SampleIDs

npx_after_br <- npx_br_data %>%
  dplyr::mutate(Type = ifelse(SampleID %in% overlapping_samples, 
                              paste(Project, "Bridge"),
                              paste(Project, "Sample"))) %>%
  dplyr:::mutate(SampleID = paste0(Project, PlateID, SampleID))

### PCA plot
OlinkAnalyze::olink_pca_plot(df          = npx_after_br,
                             color_g     = "Type",
                             byPanel     = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  new_normalized_data <- npx_br_data %>%
#    dplyr::filter(Project == "data2") %>%
#    dplyr::select(-Project, -Adj_factor) %>%
#    write.table(, file = "New_Normalized_NPX_data.csv", sep = ";")
#  

