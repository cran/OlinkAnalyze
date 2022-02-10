## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("remotes")

## ---- eval=FALSE--------------------------------------------------------------
#  remotes::install_github(repo ='Olink-Proteomics/OlinkRPackage/OlinkAnalyze', ref = "main", build_vignettes = TRUE)

## ----message=FALSE, warning=FALSE---------------------------------------------
# Load OlinkAnalyze
library(OlinkAnalyze)

# Load other libraries used in Vignette
library(dplyr)
library(ggplot2)
library(stringr)

## ----message=FALSE, eval=FALSE------------------------------------------------
#  data <- read_NPX("~/NPX_file_location.xlsx")

## ----message=FALSE, eval=FALSE------------------------------------------------
#  olink_plate_randomizer(manifest,
#                         SubjectColumn ="SampleID",
#                         seed=111)

## ----message=FALSE, eval=FALSE------------------------------------------------
#  # Select overlapping samples
#  olink_bridgeselector(df = npx_data1,
#                       sampleMissingFreq = 0.1,
#                       n = 8)

## ----message=FALSE, eval=FALSE------------------------------------------------
#  # Find overlapping samples
#  overlap_samples <- intersect(npx_data1$SampleID, npx_data2$SampleID) %>%
#    data.frame() %>%
#    filter(!str_detect(., 'CONTROL_SAMPLE')) %>% #Remove control samples
#    pull(.)
#  # Perform Bridging normalization
#  olink_normalization(df1 = npx_data1,
#                      df2 = npx_data2,
#                      overlapping_samples_df1 = overlap_samples,
#                      df1_project_nr = '20200001',
#                      df2_project_nr = '20200002',
#                      reference_project = '20200001')

## ----message=FALSE, eval=FALSE------------------------------------------------
#  olink_ttest(df = npx_data1,
#              variable = 'Treatment')

## ----message=FALSE, eval=FALSE------------------------------------------------
#  # One-way ANOVA, no covariates
#  anova_results_oneway <- olink_anova(df = npx_data1,
#                                      variable = 'Site')
#  # Two-way ANOVA, no covariates
#  anova_results_twoway <- olink_anova(df = npx_data1,
#                                      variable = c('Site', 'Time'))
#  # One-way ANOVA, Treatment as covariates
#  anova_results_oneway <- olink_anova(df = npx_data1,
#                                      variable = 'Site',
#                                      covariates = 'Treatment')

## ----message=FALSE, eval=FALSE------------------------------------------------
#  # calculate the p-value for the ANOVA
#  anova_results_oneway <- olink_anova(df = npx_data1,
#                                      variable = 'Site')
#  # extracting the significant proteins
#  anova_results_oneway_significant <- anova_results_oneway %>%
#    filter(Threshold == 'Significant') %>%
#    pull(OlinkID)
#  anova_posthoc_oneway_results <- olink_anova_posthoc(df = npx_data1,
#                                                      olinkid_list = anova_results_oneway_significant,
#                                                      variable = 'Site',
#                                                      effect = 'Site')

## ----message=FALSE, eval=FALSE------------------------------------------------
#  # Linear mixed model with one variable.
#  lmer_results_oneway <- olink_lmer(df = npx_data1,
#                                    variable = 'Site',
#                                    random = 'Subject')
#  # Linear mixed model with two variables.
#  lmer_results_twoway <- olink_lmer(df = npx_data1,
#                                    variable = c('Site', 'Treatment'),
#                                    random = 'Subject')

## ----message=FALSE, eval=FALSE------------------------------------------------
#  # Linear mixed model with two variables.
#  lmer_results_twoway <- olink_lmer(df = npx_data1,
#                                    variable = c('Site', 'Treatment'),
#                                    random = 'Subject')
#  # extracting the significant proteins
#  lmer_results_twoway_significant <- lmer_results_twoway %>%
#    filter(Threshold == 'Significant', term == 'Treatment') %>%
#    pull(OlinkID)
#  # performing post-hoc analysis
#  lmer_posthoc_twoway_results <- olink_lmer_posthoc(df = npx_data1,
#                                                    olinkid_list = lmer_results_twoway_significant,
#                                                    variable = c('Site', 'Treatment'),
#                                                    random = 'Subject',
#                                                    effect = 'Treatment')

## ----message=FALSE, fig.width=8, fig.height=4---------------------------------
plot <- npx_data1 %>%
  na.omit() %>% # removing missing values which exists for Site
  olink_boxplot(variable = "Site", 
                olinkid_list = c("OID01216", "OID01217"),
                number_of_proteins_per_plot  = 2)
plot[[1]]

## ----message=FALSE, fig.width=8, fig.height=4---------------------------------
npx_data1 %>% 
  filter(Panel == 'Olink Cardiometabolic') %>% # For this example only plotting one panel.
  olink_dist_plot() +
  theme(axis.text.x = element_blank()) # Due to the number of samples one can remove the text or rotate it

## ----message=FALSE, fig.width=8, fig.height=4---------------------------------
plot <- olink_lmer_plot(df = npx_data1, 
                        olinkid_list = c("OID01216", "OID01217"), 
                        variable = c('Site', 'Treatment'), 
                        x_axis_variable =  'Site',
                        col_variable = 'Treatment',
                        random = 'Subject')
plot[[1]]

## ----message=FALSE, fig.width=8, fig.height=4---------------------------------
npx_data1 %>% 
  filter(!str_detect(SampleID, 'CONTROL_SAMPLE')) %>% 
  olink_pca_plot(df = .,
                 color_g = "QC_Warning", byPanel = TRUE)  

## ----message = FALSE----------------------------------------------------------
npx_data <- npx_data1 %>%
    mutate(SampleID = paste(SampleID, "_", Index, sep = ""))
g <- olink_pca_plot(df=npx_data, color_g = "QC_Warning",
                    outlierDefX = 2.5, outlierDefY = 4, byPanel = TRUE, quiet = TRUE)
lapply(g, function(x){x$data}) %>%
  bind_rows() %>%
  filter(Outlier == 1) %>% 
  select(SampleID, Outlier, Panel)

## ----message=FALSE, fig.width=8, fig.height=4---------------------------------
npx_data1 %>% 
  filter(!str_detect(SampleID, 'CONTROL_SAMPLE'),
         Panel == 'Olink Inflammation') %>% 
  olink_qc_plot(color_g = "QC_Warning")   

## ----message = FALSE----------------------------------------------------------
qc <- olink_qc_plot(npx_data1, color_g = "QC_Warning", IQR_outlierDef = 3, median_outlierDef = 3)
qc$data %>% filter(Outlier == 1) %>% select(SampleID, Panel, IQR, sample_median, Outlier)

## ----message=FALSE, fig.width=8, fig.height=4---------------------------------
# perform t-test
ttest_results <- olink_ttest(df = npx_data1,
                             variable = 'Treatment')
# select names of proteins to show
top_10_name <- ttest_results %>%
  slice_head(n = 10) %>%
  pull(OlinkID)
# volcano plot
olink_volcano_plot(p.val_tbl = ttest_results,
                   x_lab = 'Treatment',
                   olinkid_list = top_10_name)

## ----message=FALSE, fig.width=8, fig.height=4---------------------------------
npx_data1 %>% 
  filter(OlinkID == 'OID01216') %>% 
  ggplot(aes(x = Treatment, y = NPX, fill = Treatment)) +
  geom_boxplot() +
  set_plot_theme()

## ----message=FALSE, fig.width=8, fig.height=4---------------------------------
npx_data1 %>% 
  filter(OlinkID == 'OID01216') %>% 
  ggplot(aes(x = Treatment, y = NPX, fill = Treatment)) +
  geom_boxplot() +
  set_plot_theme() +
  olink_fill_discrete()

