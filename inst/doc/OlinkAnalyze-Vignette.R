## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 6,
  fig.height = 3,
  fig.align = "center",
  collapse = TRUE,
  comment = "#>"
)

options(tibble.print_min = 4L, tibble.print_max = 4L)

## ----eval=FALSE---------------------------------------------------------------
# install.packages("OlinkAnalyze")

## ----message=FALSE, eval=FALSE------------------------------------------------
# data <- OlinkAnalyze::read_NPX("~/NPX_file_location.xlsx")

## ----message=FALSE, eval=FALSE------------------------------------------------
# # Read in multiple NPX files in .csv format
# data <- list.files(
#   path = "path/to/dir/with/NPX/files",
#   pattern = "csv$",
#   full.names = TRUE
# ) |>
#   lapply(FUN = function(x) {
#     df_tmp <- OlinkAnalyze::read_NPX(x) |>
#       # Optionally add additional columns to add file identifiers
#       dplyr::mutate(File = x)
#     return(df_tmp)
#   })  |>
#   # optional to return a single data frame of all files instead of a list of dfs
#   dplyr::bind_rows()
# 
# # Read in multiple NPX files in .parquet format
# data <- list.files(
#   path = "path/to/dir/with/NPX/files",
#   pattern = "parquet$",
#   full.names = TRUE
# ) |>
#   lapply(
#     OlinkAnalyze::read_NPX
#   )  |>
#   dplyr::bind_rows()
# 
# # Read in multiple NPX files in either format
# data <- list.files(
#   path = "path/to/dir/with/NPX/files",
#   pattern = "parquet$|csv$",
#   full.names = TRUE
# ) |>
#   lapply(
#     OlinkAnalyze::read_NPX
#   )  |>
#   dplyr::bind_rows()

## ----message=FALSE, eval=FALSE------------------------------------------------
# # Check NPX data quality and format
# check_npx_result <- OlinkAnalyze::check_npx(
#   df = OlinkAnalyze::npx_data1
# )

## ----message=FALSE, eval=FALSE------------------------------------------------
# # Clean the NPX data using the check_npx output
# npx_clean <- OlinkAnalyze::clean_npx(
#   df = OlinkAnalyze::npx_data1,
#   check_log = check_npx_result
# )

## ----message=FALSE, eval=FALSE------------------------------------------------
# # Check NPX data quality and format
# check_npx_clean <- OlinkAnalyze::check_npx(
#   df = npx_clean
# )

## ----message=FALSE, eval=FALSE------------------------------------------------
# # Run check_npx and clean_npx before analysis
# OlinkAnalyze::olink_ttest(
#   df = npx_clean,
#   variable = "Treatment",
#   check_log = check_npx_clean
# )

## ----message=FALSE, eval=FALSE------------------------------------------------
# OlinkAnalyze::olink_wilcox(
#   df = npx_clean,
#   variable = "Treatment",
#   check_log = check_npx_clean
# )

## ----message=FALSE, eval=FALSE------------------------------------------------
# # One-way ANOVA, no covariates
# anova_results_oneway <- OlinkAnalyze::olink_anova(
#   df = npx_clean,
#   variable = "Site",
#   check_log = check_npx_clean
# )
# # Two-way ANOVA, no covariates
# anova_results_twoway <- OlinkAnalyze::olink_anova(
#   df = npx_clean,
#   variable = c("Site", "Time"),
#   check_log = check_npx_clean
# )
# # One-way ANOVA, Treatment as covariates
# anova_results_oneway <- OlinkAnalyze::olink_anova(
#   df = npx_clean,
#   variable = "Site",
#   covariates = "Treatment",
#   check_log = check_npx_clean
# )

## ----message=FALSE, eval=FALSE------------------------------------------------
# # calculate the p-value for the ANOVA
# anova_results_oneway <- OlinkAnalyze::olink_anova(
#   df = npx_clean,
#   variable = "Site",
#   check_log = check_npx_clean
# )
# 
# # extracting the significant proteins
# anova_results_oneway_sign <- anova_results_oneway |>
#   dplyr::filter(
#     .data[["Threshold"]] == "Significant"
#   ) |>
#   dplyr::pull(
#     .data[["OlinkID"]]
#   )
# 
# anova_posthoc_oneway_results <- OlinkAnalyze::olink_anova_posthoc(
#   df = npx_clean,
#   olinkid_list = anova_results_oneway_sign,
#   variable = "Site",
#   effect = "Site",
#   check_log = check_npx_clean
# )

## ----message=FALSE, eval=FALSE------------------------------------------------
# # Linear mixed model with one variable.
# lmer_results_oneway <- OlinkAnalyze::olink_lmer(
#   df = npx_clean,
#   variable = "Site",
#   random = "Subject",
#   check_log = check_npx_clean
# )
# 
# # Linear mixed model with two variables.
# lmer_results_twoway <- OlinkAnalyze::olink_lmer(
#   df = npx_clean,
#   variable = c("Site", "Treatment"),
#   random = "Subject",
#   check_log = check_npx_clean
# )

## ----message=FALSE, eval=FALSE------------------------------------------------
# # Linear mixed model with two variables.
# lmer_results_twoway <- OlinkAnalyze::olink_lmer(
#   df = npx_clean,
#   variable = c("Site", "Treatment"),
#   random = "Subject",
#   check_log = check_npx_clean
# )
# 
# # extracting the significant proteins
# lmer_results_twoway_sign <- lmer_results_twoway |>
#   dplyr::filter(
#     .data[["Threshold"]] == "Significant" &
#       .data[["term"]] == "Treatment"
#   ) |>
#   dplyr::pull(
#     .data[["OlinkID"]]
#   )
# 
# # performing post-hoc analysis
# lmer_posthoc_twoway_results <- OlinkAnalyze::olink_lmer_posthoc(
#   df = npx_clean,
#   olinkid_list = lmer_results_twoway_sign,
#   variable = c("Site", "Treatment"),
#   random = "Subject",
#   effect = "Treatment",
#   check_log = check_npx_clean
# )

## ----message=FALSE, eval=FALSE------------------------------------------------
# ttest_results <- OlinkAnalyze::olink_ttest(
#   df = npx_df_clean,
#   variable = "Treatment",
#   alternative = "two.sided",
#   check_log = check_npx_clean
# )
# 
# gsea_results <- OlinkAnalyze::olink_pathway_enrichment(
#   df = npx_df_clean,
#   test_results = ttest_results
# )
# 
# ora_results <- OlinkAnalyze::olink_pathway_enrichment(
#   df = npx_df_clean,
#   test_results = ttest_results,
#   method = "ORA"
# )

## ----message=FALSE, eval=FALSE------------------------------------------------
# OlinkAnalyze::olink_umap_plot(
#   df = npx_clean,
#   color_g = "QC_Warning",
#   byPanel = TRUE,
#   check_log = check_npx_clean
# )

## ----message=FALSE, echo = FALSE----------------------------------------------
knitr::include_graphics(
  path = normalizePath("../man/figures/olink_umap_plot.png"),
  error = FALSE
)

## ----message=FALSE, eval=FALSE------------------------------------------------
# plot <- npx_clean |>
#   dplyr::filter(
#     !is.na(.data[["Site"]])
#   ) |> # removing missing values which exist for Site
#   OlinkAnalyze::olink_boxplot(
#     variable = "Site",
#     olinkid_list = c("OID00488", "OID01276"),
#     number_of_proteins_per_plot = 2L,
#     check_log = check_npx_clean
#   )
# 
# plot[[1L]]

## ----message=FALSE, echo = FALSE----------------------------------------------
knitr::include_graphics(
  path = normalizePath("../man/figures/olink_boxplot.png"),
  error = FALSE
)

## ----message=FALSE, eval=FALSE------------------------------------------------
# anova_posthoc_results <- npx_clean |>
#   OlinkAnalyze::olink_anova_posthoc(
#     olinkid_list = c("OID00488", "OID01276"),
#     variable = "Site",
#     effect = "Site",
#     check_log = check_npx_clean
#   )
# 
# plot2 <- npx_clean |>
#   stats::na.omit() |> # removing missing values which exists for Site
#   OlinkAnalyze::olink_boxplot(
#     variable = "Site",
#     olinkid_list = c("OID00488", "OID01276"),
#     number_of_proteins_per_plot = 2L,
#     posthoc_results = anova_posthoc_results,
#     check_log = check_npx_clean
#   )
# 
# plot2[[1L]]

## ----message=FALSE, echo=FALSE------------------------------------------------
knitr::include_graphics(
  path = normalizePath("../man/figures/olink_boxplot_anova_posthoc.png"),
  error = FALSE
)

## ----message=FALSE, eval=FALSE------------------------------------------------
# plot <- OlinkAnalyze::olink_lmer_plot(
#   df = npx_clean,
#   olinkid_list = c("OID01216", "OID01217"),
#   variable = c("Site", "Treatment"),
#   x_axis_variable = "Site",
#   col_variable = "Treatment",
#   random = "Subject",
#   check_log = check_npx_clean
# )
# plot[[1L]]

## ----message=FALSE, fig.width=8, echo=FALSE-----------------------------------
knitr::include_graphics(
  path = normalizePath("../man/figures/olink_lmer_plot.png"),
  error = FALSE
)

## ----message=FALSE, eval=FALSE------------------------------------------------
# OlinkAnalyze::olink_pathway_heatmap(
#   enrich_results = gsea_results,
#   test_results = ttest_results
# )

## ----message=FALSE, echo=FALSE------------------------------------------------
knitr::include_graphics(
  path = normalizePath("../man/figures/olink_pathway_heatmap_gsea.png"),
  error = FALSE
)

## ----message=FALSE, fig.height=4, fig.width=8, eval=FALSE---------------------
# OlinkAnalyze::olink_pathway_heatmap(
#   enrich_results = ora_results,
#   test_results = ttest_results,
#   method = "ORA",
#   keyword = "immune"
# )

## ----message=FALSE, echo=FALSE------------------------------------------------
knitr::include_graphics(
  path = normalizePath("../man/figures/olink_pathway_heatmap_ora.png"),
  error = FALSE
)

## ----message=FALSE, eval=FALSE------------------------------------------------
# first10 <- npx_clean |>
#   dplyr::pull(
#     .data[["OlinkID"]]
#   ) |>
#   unique() |>
#   utils::head(n = 10L)
# 
# first15samples <- npx_clean |>
#   dplyr::pull(
#     .data[["SampleID"]]
#   ) |>
#   unique() |>
#   utils::head(n = 15L)
# 
# npx_data_small <- npx_clean |>
#   dplyr::filter(
#     .data[["OlinkID"]] %in% .env[["first10"]]
#   ) |>
#   dplyr::filter(
#     .data[["SampleID"]] %in% .env[["first15samples"]]
#   )
# 
# OlinkAnalyze::olink_heatmap_plot(
#   df = npx_data_small,
#   variable_row_list = "Treatment",
#   check_log = check_npx_clean
# )

## ----message=FALSE, fig.height=4, echo=FALSE----------------------------------
knitr::include_graphics(
  path = normalizePath("../man/figures/olink_heatmap_plot.png"),
  error = FALSE
)

## ----message=FALSE, eval=FALSE------------------------------------------------
# # perform t-test
# ttest_results <- OlinkAnalyze::olink_ttest(
#   df = npx_clean,
#   variable = "Treatment",
#   check_log = check_npx_clean
# )
# 
# # select names of proteins to show
# top_10_name <- ttest_results |>
#   dplyr::slice_head(
#     n = 10L
#   ) |>
#   dplyr::pull(
#     .data[["OlinkID"]]
#   )
# 
# # volcano plot
# OlinkAnalyze::olink_volcano_plot(
#   p.val_tbl = ttest_results,
#   x_lab = "Treatment",
#   olinkid_list = top_10_name
# )

## ----message=FALSE, echo=FALSE------------------------------------------------
knitr::include_graphics(
  path = normalizePath("../man/figures/olink_volcano_plot.png"),
  error = FALSE
)

## ----message=FALSE, eval=FALSE------------------------------------------------
# OlinkAnalyze::npx_data1 |>
#   dplyr::filter(
#     !is.na(.data[["Treatment"]])
#   ) |>
#   dplyr::filter(
#     .data[["OlinkID"]] == "OID01216"
#   ) |>
#   ggplot2::ggplot(
#     ggplot2::aes(
#       x = .data[["Treatment"]],
#       y = .data[["NPX"]],
#       fill = .data[["Treatment"]]
#     )
#   ) +
#   ggplot2::geom_boxplot() +
#   OlinkAnalyze::set_plot_theme()

## ----message=FALSE, echo=FALSE------------------------------------------------
knitr::include_graphics(
  path = normalizePath("../man/figures/set_plot_theme_boxplot.png"),
  error = FALSE
)

## ----message=FALSE, eval=FALSE------------------------------------------------
# OlinkAnalyze::npx_data1 |>
#   dplyr::filter(
#     !is.na(.data[["Treatment"]])
#   ) |>
#   dplyr::filter(
#     .data[["OlinkID"]] == "OID01216"
#   ) |>
#   ggplot2::ggplot(
#     mapping = ggplot2::aes(
#       x = .data[["Treatment"]],
#       y = .data[["NPX"]],
#       fill = .data[["Treatment"]]
#     )
#   ) +
#   ggplot2::geom_boxplot() +
#   OlinkAnalyze::set_plot_theme() +
#   OlinkAnalyze::olink_fill_discrete()

## ----message=FALSE, echo=FALSE------------------------------------------------
knitr::include_graphics(
  path = normalizePath("../man/figures/olink_fill_discrete_boxplot.png"),
  error = FALSE
)

## ----message=FALSE, eval=FALSE------------------------------------------------
# npx_ht <- data_exploreht |>
#   dplyr::filter(
#     .data[["SampleType"]] == "SAMPLE"
#   ) |>
#   dplyr::mutate(
#     Project = "data1"
#   )
# 
# check_npx_ht <- OlinkAnalyze::check_npx(
#   df = npx_ht
# )
# 
# npx_3072 <- data_explore3072 |>
#   dplyr::filter(
#     .data[["SampleType"]] == "SAMPLE"
#   ) |>
#   dplyr::mutate(
#     Project = "data2"
#   )
# 
# check_npx_3072 <- OlinkAnalyze::check_npx(
#   df = npx_3072
# )
# 
# overlapping_samples <- unique(
#   intersect(
#     x = npx_ht |> dplyr::distinct(.data[["SampleID"]]) |> dplyr::pull(),
#     y = npx_3072 |> dplyr::distinct(.data[["SampleID"]]) |> dplyr::pull()
#   )
# )
# 
# npx_br_data <- OlinkAnalyze::olink_normalization(
#   df1 = npx_ht,
#   df2 = npx_3072,
#   overlapping_samples_df1 = overlapping_samples,
#   df1_project_nr = "Explore HT",
#   df2_project_nr = "Explore 3072",
#   reference_project = "Explore HT",
#   format = FALSE,
#   df1_check_log = check_npx_ht,
#   df2_check_log = check_npx_3072
# )
# 
# check_npx_br_data <- OlinkAnalyze::check_npx(
#   df = npx_br_data
# )
# 
# npx_br_data_bridgeable_plt <- OlinkAnalyze::olink_bridgeability_plot(
#   df = npx_br_data,
#   median_counts_threshold = 150L,
#   min_count = 10L,
#   check_log = check_npx_br_data
# )
# 
# npx_br_data_bridgeable_plt[[1L]]

## ----message=FALSE, echo=FALSE, out.width="600px"-----------------------------
knitr::include_graphics(
  path = normalizePath("../man/figures/bridgeable_plt_MedianCenter.png"),
  error = FALSE
)

