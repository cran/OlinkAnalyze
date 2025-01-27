## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 6,
  fig.height = 3,
  fig.align = "center",
  collapse = TRUE,
  comment = "#>")

options(tibble.print_min = 4L, tibble.print_max = 4L)

## ----eval=FALSE---------------------------------------------------------------
# install.packages("OlinkAnalyze")

## ----message=FALSE, warning=FALSE---------------------------------------------
# Load OlinkAnalyze
library(OlinkAnalyze)

# Load other libraries used in Vignette
library(dplyr)
library(ggplot2)
library(stringr)

## ----message=FALSE, eval=FALSE------------------------------------------------
# data <- read_NPX("~/NPX_file_location.xlsx")

## ----message=FALSE, eval=FALSE------------------------------------------------
# # Read in multiple NPX files in .csv format
# data <-  list.files(path = "path/to/dir/with/NPX/files",
#                     pattern = "csv$",
#                     full.names = TRUE) |>
#          lapply(FUN = function(x){
#     OlinkAnalyze::read_NPX(x) |>
#       # Optionally add additional columns to add file identifiers
#       dplyr::mutate(File = x)
#       }  |>
#          # optional to return a single data frame of all files instead of a list of data frames
#          dplyr::bind_rows()
# 
# # Read in multiple NPX files in .parquet format
# data <-  list.files(path = "path/to/dir/with/NPX/files",
#                     pattern = "parquet$",
#                     full.names = TRUE) |>
#          lapply(OlinkAnalyze::read_NPX)  |>
#          dplyr::bind_rows()
# 
# # Read in multiple NPX files in either format
# data <-  list.files(path = "path/to/dir/with/NPX/files",
#                     pattern = "parquet$|csv$",
#                     full.names = TRUE) |>
#          lapply(OlinkAnalyze::read_NPX)  |>
#          dplyr::bind_rows()

## ----message=FALSE, eval=FALSE------------------------------------------------
# olink_ttest(df = npx_data1,
#             variable = 'Treatment')

## ----message=FALSE, eval=FALSE------------------------------------------------
# olink_wilcox(df = npx_data1,
#              variable = 'Treatment')

## ----message=FALSE, eval=FALSE------------------------------------------------
# # Remove control samples and assays
# npx_data1_no_controls <- npx_data1 |>
#   dplyr::filter(!stringr::str_detect(SampleID,
#                                      regex("control|ctrl",
#                                            ignore_case = TRUE))) |>
#   dplyr::filter(!stringr::str_detect(Assay,
#                                      regex("control|ctrl",
#                                            ignore_case = TRUE)))
# 
# # One-way ANOVA, no covariates
# anova_results_oneway <- olink_anova(df = npx_data1_no_controls,
#                                     variable = 'Site')
# # Two-way ANOVA, no covariates
# anova_results_twoway <- olink_anova(df = npx_data1_no_controls,
#                                     variable = c('Site', 'Time'))
# # One-way ANOVA, Treatment as covariates
# anova_results_oneway <- olink_anova(df = npx_data1_no_controls,
#                                     variable = 'Site',
#                                     covariates = 'Treatment')

## ----message=FALSE, eval=FALSE------------------------------------------------
# # calculate the p-value for the ANOVA
# anova_results_oneway <- olink_anova(df = npx_data1_no_controls,
#                                     variable = 'Site')
# 
# # extracting the significant proteins
# anova_results_oneway_significant <- anova_results_oneway %>%
#   dplyr::filter(Threshold == 'Significant') %>%
#   dplyr::pull(OlinkID)
# 
# anova_posthoc_oneway_results <- olink_anova_posthoc(df = npx_data1_no_controls,
#                                                     olinkid_list = anova_results_oneway_significant,
#                                                     variable = 'Site',
#                                                     effect = 'Site')

## ----message=FALSE, eval=FALSE------------------------------------------------
# if (requireNamespace("lme4", quietly = TRUE) & requireNamespace("lmerTest", quietly = TRUE)){
# # Linear mixed model with one variable.
#   lmer_results_oneway <- olink_lmer(df = npx_data1,
#                                     variable = 'Site',
#                                     random = 'Subject')
# # Linear mixed model with two variables.
#   lmer_results_twoway <- olink_lmer(df = npx_data1,
#                                     variable = c('Site', 'Treatment'),
#                                     random = 'Subject')
# }

## ----message=FALSE, eval=FALSE------------------------------------------------
# if (requireNamespace("lme4", quietly = TRUE) & requireNamespace("lmerTest", quietly = TRUE)){
#   # Linear mixed model with two variables.
#   lmer_results_twoway <- olink_lmer(df = npx_data1,
#                                     variable = c('Site', 'Treatment'),
#                                     random = 'Subject')
#   # extracting the significant proteins
#   lmer_results_twoway_significant <- lmer_results_twoway %>%
#     dplyr::filter(Threshold == 'Significant', term == 'Treatment') %>%
#     dplyr::pull(OlinkID)
#   # performing post-hoc analysis
#   lmer_posthoc_twoway_results <- olink_lmer_posthoc(df = npx_data1,
#                                                     olinkid_list = lmer_results_twoway_significant,
#                                                     variable = c('Site', 'Treatment'),
#                                                     random = 'Subject',
#                                                     effect = 'Treatment')
# }

## ----message=FALSE------------------------------------------------------------
npx_df <- npx_data1 |>
  dplyr::filter(!grepl("control", SampleID, ignore.case = TRUE))

ttest_results <- olink_ttest(
  df = npx_df,
  variable = "Treatment",
  alternative = "two.sided")

try({ # This expression might fail if dependencies are not installed
gsea_results <- olink_pathway_enrichment(data = npx_data1, test_results = ttest_results)
ora_results <- olink_pathway_enrichment(
  data = npx_data1,
  test_results = ttest_results, method = "ORA")
}, silent = TRUE)

## ----message=FALSE, eval = FALSE----------------------------------------------
# if (requireNamespace("umap", quietly = TRUE) ){
# npx_data1 %>%
#   dplyr::filter(!stringr::str_detect(SampleID, 'CONTROL_SAMPLE')) %>%
#   olink_umap_plot(df = .,
#                  color_g = "QC_Warning", byPanel = TRUE)
#   }

## ----message=FALSE, echo = FALSE----------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/olink_umap_plot.png"),
                        error = FALSE)

## ----message=FALSE, eval =FALSE-----------------------------------------------
# npx_data1 %>%
#   dplyr::filter(!stringr::str_detect(SampleID, 'CONTROL_SAMPLE')) %>%
#   olink_umap_plot(df = .,
#                  color_g = "QC_Warning", byPanel = TRUE)

## ----message=FALSE, eval = FALSE----------------------------------------------
# # Remove control samples and assays
# npx_data1_no_controls <- npx_data1 |>
#   dplyr::filter(!stringr::str_detect(SampleID,
#                                      regex("control|ctrl",
#                                            ignore_case = TRUE))) |>
#   dplyr::filter(!stringr::str_detect(Assay,
#                                      regex("control|ctrl",
#                                            ignore_case = TRUE)))
# 
# plot <- npx_data1_no_controls |>
#   dplyr::filter(!is.na(Site)) |> # removing missing values which exist for Site
#   olink_boxplot(variable = "Site",
#                 olinkid_list = c("OID00488", "OID01276"),
#                 number_of_proteins_per_plot  = 2)
# 
# plot[[1]]

## ----message=FALSE, echo = FALSE----------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/olink_boxplot.png"),
                        error = FALSE)

## ----message=FALSE, eval = FALSE----------------------------------------------
# anova_posthoc_results<-npx_data1_no_controls %>%
#   olink_anova_posthoc(olinkid_list = c("OID00488", "OID01276"),
#                       variable = 'Site',
#                       effect = 'Site')
# 
# plot2 <- npx_data1_no_controls %>%
#   na.omit() %>% # removing missing values which exists for Site
#   olink_boxplot(variable = "Site",
#                 olinkid_list = c("OID00488", "OID01276"),
#                 number_of_proteins_per_plot  = 2,
#                 posthoc_results = anova_posthoc_results)
# 
# plot2[[1]]
# 

## ----message=FALSE, echo = FALSE----------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/olink_boxplot_anova_posthoc.png"),
                        error = FALSE)

## ----message=FALSE, eval = FALSE----------------------------------------------
# if (requireNamespace("lme4", quietly = TRUE) & requireNamespace("lmerTest", quietly = TRUE)){
#   plot <- olink_lmer_plot(df = npx_data1,
#                           olinkid_list = c("OID01216", "OID01217"),
#                           variable = c('Site', 'Treatment'),
#                           x_axis_variable =  'Site',
#                           col_variable = 'Treatment',
#                           random = 'Subject')
#   plot[[1]]
# }

## ----message=FALSE, fig.width = 8, echo = FALSE-------------------------------
knitr::include_graphics(normalizePath("../man/figures/olink_lmer_plot.png"),
                        error = FALSE)

## ----message=FALSE, eval = FALSE----------------------------------------------
# # GSEA Heatmap from t-test results
# try({ # This expression might fail if dependencies are not installed
# olink_pathway_heatmap(enrich_results = gsea_results,
#                       test_results = ttest_results)
# })

## ----message=FALSE, echo = FALSE----------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/olink_pathway_heatmap_gsea.png"), error = FALSE)

## ----message = FALSE, fig.height=4, fig.width=8, eval = FALSE-----------------
# # ORA Heatmap from t-test results with cell keyword
# try({ # This expression might fail if dependencies are not installed
# olink_pathway_heatmap(enrich_results = ora_results,
#                       test_results = ttest_results,
#                       method = "ORA", keyword = "immune")
# })

## ----message=FALSE, echo = FALSE----------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/olink_pathway_heatmap_ora.png"), error = FALSE)

## ----message=FALSE,  eval = FALSE---------------------------------------------
# first10 <- npx_data1 %>%
#   dplyr::pull(OlinkID) %>%
#   unique() %>%
#   head(10)
# 
# first15samples <- npx_data1$SampleID %>%
#   unique() %>%
#   head(15)
# 
# npx_data_small <- npx_data1 %>%
#   dplyr::filter(!str_detect(SampleID, 'CONT')) %>%
#   dplyr::filter(OlinkID %in% first10) %>%
#   dplyr::filter(SampleID %in% first15samples)
# 
# try({ #This statement might fail if dependencies are not installed
#   olink_heatmap_plot(npx_data_small, variable_row_list =  'Treatment')
#   })

## ----message=FALSE, fig.height=4, echo = FALSE--------------------------------
knitr::include_graphics(normalizePath("../man/figures/olink_heatmap_plot.png"), error = FALSE)

## ----message=FALSE, eval = FALSE----------------------------------------------
# # perform t-test
# ttest_results <- olink_ttest(df = npx_data1,
#                              variable = 'Treatment')
# # select names of proteins to show
# top_10_name <- ttest_results %>%
#   dplyr::slice_head(n = 10) %>%
#   dplyr::pull(OlinkID)
# # volcano plot
# olink_volcano_plot(p.val_tbl = ttest_results,
#                    x_lab = 'Treatment',
#                    olinkid_list = top_10_name)

## ----message=FALSE, echo = FALSE----------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/olink_volcano_plot.png"), error = FALSE)

## ----message=FALSE, eval = FALSE----------------------------------------------
# npx_data1 %>%
#   dplyr::filter(!is.na(Treatment)) |>
#   dplyr::filter(OlinkID == 'OID01216') %>%
#   ggplot(aes(x = Treatment, y = NPX, fill = Treatment)) +
#   geom_boxplot() +
#   set_plot_theme()

## ----message=FALSE, echo = FALSE----------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/set_plot_theme_boxplot.png"), error = FALSE)

## ----message=FALSE, eval = FALSE----------------------------------------------
# npx_data1 %>%
#   dplyr::filter(!is.na(Treatment)) |>
#   dplyr::filter(OlinkID == 'OID01216') %>%
#   ggplot(aes(x = Treatment, y = NPX, fill = Treatment)) +
#   geom_boxplot() +
#   set_plot_theme() +
#   olink_fill_discrete()

## ----message=FALSE, echo = FALSE----------------------------------------------
knitr::include_graphics(normalizePath("../man/figures/olink_fill_discrete_boxplot.png"), error = FALSE)

## ----message=FALSE, eval = FALSE----------------------------------------------
# 
# npx_ht <- data_exploreht |>
#   dplyr::filter(SampleType == "SAMPLE") |>
#   dplyr::mutate(Project = "data1")
# 
# npx_3072 <- data_explore3072 |>
#   dplyr::filter(SampleType == "SAMPLE") |>
#   dplyr::mutate(Project = "data2")
# 
# overlapping_samples <- unique(intersect(npx_ht |>
#                                           dplyr::distinct(SampleID) |>
#                                           dplyr::pull(),
#                                         npx_3072 |>
#                                           dplyr::distinct(SampleID) |>
#                                           dplyr::pull()))
# 
# npx_br_data <- OlinkAnalyze::olink_normalization(
#   df1 = npx_ht,
#   df2 = npx_3072,
#   overlapping_samples_df1 = overlapping_samples,
#   df1_project_nr = "Explore HT",
#   df2_project_nr = "Explore 3072",
#   reference_project = "Explore HT"
#   )
# 
# 
# npx_br_data_bridgeable_plt <- olink_bridgeability_plot(
#   data = npx_br_data,
#   median_counts_threshold = 150,
#   min_count = 10
#   )
# 
# npx_br_data_bridgeable_plt[[1]]

## ----message=FALSE, echo = FALSE, out.width = "600px"-------------------------
knitr::include_graphics(normalizePath("../man/figures/bridgeable_plt_MedianCenter.png"), error = FALSE)

