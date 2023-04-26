## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy = FALSE,
  tidy.opts = list(width.cutoff = 95),
  fig.width = 10,
  message = FALSE,
  warning = FALSE,
  time_it = TRUE,
  fig.width = 7
)

## ---- echo=FALSE--------------------------------------------------------------
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
      caption = "Recommended number of bridging samples for Olink platforms") %>%
  kable_styling(bootstrap_options = "striped",
                full_width = FALSE,
                position = "center",
                latex_options = "HOLD_position")


## ----bridge_sample_selection, echo=FALSE--------------------------------------
olink_bridgeselector(df = npx_data1,
                     sampleMissingFreq = 0.1,
                     n = 16) %>%
  kableExtra::kbl(booktabs = TRUE,
      digits = 2,
      caption = "Selected Bridging Samples") %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE,
                position = "center", latex_options = "HOLD_position")
  

## ----message=FALSE, eval=FALSE, echo = TRUE-----------------------------------
#  data1 <- read_NPX("~/NPX_file1_location.xlsx")
#  data2 <- read_NPX("~/NPX_file2_location.xlsx")

## ---- echo=FALSE--------------------------------------------------------------

data.frame(SampleID = intersect(npx_data1$SampleID, npx_data2$SampleID)) %>%
  dplyr::filter(!stringr::str_detect(SampleID, "CONTROL_SAMPLE")) %>% #Remove control samples
  kableExtra::kbl(booktabs = TRUE,
      digits = 2,
      caption = "Overlapping Bridging Samples") %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE,
                position = "center", latex_options = "HOLD_position")

## ----check, message=FALSE-----------------------------------------------------
# Load datasets
npx_1 <- npx_data1 %>%
  mutate(dataset = "data1")
npx_2 <- npx_data2 %>%
  mutate(dataset = "data2")

npx_df <- bind_rows(npx_1, npx_2)


# Plot NPX density before bridging normalization
npx_df %>%
  mutate(Panel = gsub("Olink ", "", Panel)) %>%
  ggplot(aes(x = NPX, fill = dataset)) +
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

## ----pca1, message=FALSE, fig.cap="PCA plot of combined datasets before bridging"----
## before bridging

#### Extract bridging samples

overlapping_samples <-  data.frame(SampleID = intersect(npx_1$SampleID, npx_2$SampleID)) %>%  
  filter(!str_detect(SampleID, "CONTROL_SAMPLE")) %>% #Remove control samples
  pull(SampleID)


npx_before_br <- npx_data1 %>%
  dplyr::filter(!str_detect(SampleID, "CONTROL_SAMPLE")) %>% #Remove control samples
  dplyr::mutate(Type = if_else(SampleID %in% overlapping_samples,
                        paste0("20200001 Bridge"),
                        paste0("20200001 Sample"))) %>%
  rbind({
    npx_data2 %>%
      filter(!str_detect(SampleID, "CONTROL_SAMPLE")) %>% #Remove control samples %>% 
      mutate(Type = if_else(SampleID %in% overlapping_samples,
                            paste0("20200002 Bridge"),
                            paste0("20200002 Sample"))) %>%
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
  mutate(dataset = "data1") 
npx_2 <- npx_data2 %>%
  mutate(dataset = "data2")

overlap_samples <-data.frame(SampleID = intersect(npx_1$SampleID, npx_2$SampleID)) %>%
  filter(!str_detect(SampleID, "CONTROL_SAMPLE")) %>% #Remove control samples
  pull(SampleID)

overlap_samples_list <- list("DF1" = overlap_samples,
                             "DF2" = overlap_samples)

# Perform Bridging normalization
npx_br_data <- olink_normalization_bridge(project_1_df = npx_1,
                                          project_2_df = npx_2,
                                          bridge_samples = overlap_samples_list,
                                          project_1_name = "20200001",
                                          project_2_name = "20200002",
                                          project_ref_name = "20200001")
dplyr::glimpse(npx_br_data)


## ----densitybr, message=FALSE-------------------------------------------------
# Plot NPX density after bridging normalization

npx_br_data %>%
  mutate(Panel = gsub("Olink ", "", Panel)) %>%
  ggplot2::ggplot(ggplot2::aes(x = NPX, fill = dataset)) +
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


## ----dist_adj_fct, message=FALSE, echo = FALSE--------------------------------
npx_br_data %>% 
  dplyr::filter(Project == "20200002") %>%  # Only looking at Project 2 since project 1 is unadjusted
  dplyr::select(OlinkID, Adj_factor) %>% 
  dplyr::distinct() %>% 
  ggplot2::ggplot(ggplot2::aes(x = Adj_factor)) +
    ggplot2::geom_histogram() +
  set_plot_theme()


## ----voilin plot--------------------------------------------------------------
# Bridge sample data
bridge_samples <- npx_1 %>%
  rbind(npx_2) %>%
  filter(SampleID %in% overlapping_samples) %>%
  filter(Assay == "CHL1") %>%
  mutate(Assay_OID = paste(Assay, OlinkID, sep = "\n"))

# Generate violin plot for CHL1
npx_data1 %>%
  mutate(Project = "20200001") %>%
  bind_rows({
    npx_data2 %>%
      mutate(Project = "20200002")
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

## ----CV_calculation-----------------------------------------------------------

explore_cv <- function(npx, na.rm = F) {
  sqrt(exp((log(2) * sd(npx, na.rm = na.rm))^2) - 1)*100
}

t96_cv <- function(NPX, na.rm = T) {
  100*sd(2^NPX)/mean(2^NPX)
}

tech <- "Explore"

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



## ----pca2, message=FALSE, fig.cap="PCA plot of combined datasets after bridging"----
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

## ---- eval = FALSE------------------------------------------------------------
#  new_normalized_data <- npx_br_data %>%
#    dplyr::filter(Project == "20200002") %>%
#    dplyr::select(-Project, -Adj_factor) %>%
#    write.table(, file = "New_Normalized_NPX_data.csv", sep = ";")
#  

