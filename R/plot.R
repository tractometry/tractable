# Use colorblind-friendly palette from http://jfly.iam.u-tokyo.ac.jp/color/
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#' Retrieve tract name from abbreviation
#'
#' @param tract AFQ tract abbreviation
#'
#' @return Formatted tract name
#' @export
#'
#' @examples
#' name <- tract_name("OR")
#' name <- tract_name("CST_L")
tract_name <- function(tract) {
  name <- switch(
    tract,
    "OR" = "Optic Radiation",
    "CST_R" = "Right Corticospinal",
    "CST_L" = "Left Corticospinal",
    "UNC_R" = "Right Uncinate",
    "UNC_L" = "Left Uncinate",
    "IFO_L" = "Left IFOF",
    "IFO_R" = "Right IFOF",
    "ARC_R" = "Right Arcuate",
    "ARC_L" = "Left Arcuate",
    "ATR_R" = "Right Thalamic Radiation",
    "ATR_L" = "Left Thalamic Radiation",
    "CGC_R" = "Right Cingulum Cingulate",
    "CGC_L" = "Left Cingulum Cingulate",
    "HCC_R" = "Right Cingulum Hippocampus",
    "HCC_L" = "Left Cingulum Hippocampus",
    "FP" = "Callosum Forceps Major",
    "FA" = "Callosum Forceps Minor",
    "ILF_R" = "RightILF",
    "ILF_L" = "LeftILF",
    "SLF_R" = "RightSLF",
    "SLF_L" = "LeftSLF",
    tract
  )

  return(name)
}


utils::globalVariables(c("metric", "tracts", "nodeID", "value"))

#' Plot GAM splines for each group
#'
#' @param gam_model GAM object, produced by gam/bam
#' @param tract AFQ tract name
#' @param df_tract A dataframe of AFQ nodes for certain tract
#' @param dwi_metric Diffusion MRI metric (e.g. FA, MD)
#' @param covariates List of strings of GAM covariates,
#'     not including the smoothing terms over nodes and the random effect due
#'     to subjectID.
#' @param participant_id The name of the column that encodes participant ID
#' @param group_by The grouping variable used to group nodeID smoothing terms
#' @param out_dir directory in which to save plots
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' gam_fit <- fit_gam(df_afq,
#'                    target = "dti_fa",
#'                    covariates = list("group", "sex"),
#'                    family = "gamma",
#'                    k = 40)
#' plot_gam_splines(gam_model = gam_fit,
#'                  tract = "OR",
#'                  df_tract = df_afq,
#'                  dwi_metric = "dti_fa",
#'                  covariates = c("group", "sex"),
#'                  out_dir = ".")
#' }
plot_gam_splines <- function(
    gam_model,
    tract,
    df_tract,
    dwi_metric,
    covariates,
    group_by       = "group",
    participant_id = "subjectID",
    out_dir
) {
  # generate predictions
  values <- vector(mode = "list", length = length(covariates))
  names(values) <- covariates
  df_pred <- mgcv::predict.bam(
    gam_model,
    exclude_terms = c(covariates, "subjectID"),
    values = values,
    se.fit = T,
    type = "response"
  )

  # convert predictions to dataframe
  df_pred <- data.frame(
    nodeID = df_tract$nodeID,
    fit = df_pred$fit,
    se.fit = df_pred$se.fit
  )

  for (covar in covariates) {
    df_pred[[covar]] <- df_tract[[covar]]
  }

  if (!is.null(group_by)) {
    df_pred[[group_by]] <- df_tract[[group_by]]
  }

  df_pred[[participant_id]] <- df_tract[[participant_id]]

  # set up for plot
  h_tract <- tract_name(tract)
  h_title <- paste0("GAM fit of ", h_tract, " ", dwi_metric, " values")

  # draw plot
  options(warn = -1)
  p <- ggplot2::ggplot(data = df_pred) +
    ggplot2::geom_smooth(mapping = ggplot2::aes_string(
      x = "nodeID", y = "fit", color = group_by
    )) +
    ggplot2::ggtitle(h_title) +
    ggplot2::ylab(paste0("Fit ", dwi_metric)) +
    ggplot2::xlab("Tract Node") +
    ggplot2::theme(text = ggplot2::element_text(
      family = "Times New Roman", face = "bold", size = 14
    ))
  options(warn = 0)

  # Use colorblind palette for fills and lines
  p + ggplot2::scale_color_manual(
    values = cbbPalette,
  )

  plot_filename <- file.path(out_dir,
                             paste0("plot_gam_", sub(" ", "_", tract), ".png"))

  ggplot2::ggsave(
    plot_filename,
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}

#' Calculate and plot difference between two splines
#'
#' This function both
#'   1) draws a spline-difference plot for 2 splines,
#'   2) and returns a dataframe of differences at each node
#'
#' @param gam_model GAM object, produced by gam/bam
#' @param tract AFQ tract name
#' @param group_by The grouping variable used to group nodeID smoothing terms
#' @param factor_a First group factor, string
#' @param factor_b Second group factor, string
#' @param save_output Boolean. If TRUE, save plot output
#' @param sim.ci Logical: Using simultaneous confidence intervals or not
#'   (default set to FALSE). The implementation of simultaneous CIs follows
#'   Gavin Simpson's blog of December 15, 2016:
#'   http://www.fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/.
#'   This interval is calculated from simulations based. Please specify a seed
#'   (e.g., set.seed(123)) for reproducable results. Note: in contrast with
#'   Gavin Simpson's code, here the Bayesian posterior covariance matrix of
#'   the parameters is uncertainty corrected (unconditional=TRUE) to reflect
#'   the uncertainty on the estimation of smoothness parameters.
#' @param out_dir Directory in which to save plots
#'
#' @return A dataframe of spline differences at each node
#' @export
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' gam_fit <- fit_gam(df_afq,
#'                    target = "dti_fa",
#'                    covariates = list("group", "sex"),
#'                    family = "gamma",
#'                    k = 40)
#' df_diff = spline_diff(gam_model = gam_fit,
#'                       tract = "OR",
#'                       factor_a = "0",
#'                       factor_b = "1",
#'                       out_dir = ".")
#' }
spline_diff <- function(gam_model,
                        tract,
                        group_by = "group",
                        factor_a,
                        factor_b,
                        save_output = TRUE,
                        sim.ci = FALSE,
                        out_dir) {
  # determine bottom of plot
  comp <- list(c(factor_a, factor_b))
  names(comp) <- c(group_by)

  df_pair <- itsadug::plot_diff(
    gam_model,
    view = "nodeID",
    comp = comp,
    rm.ranef = TRUE,
    plot = FALSE,
    print.summary = FALSE,
    sim.ci = sim.ci,
    n.grid = 100,
  )

  h_min <- min(df_pair$est)
  h_ci <- df_pair[which(df_pair$est == h_min), ]$CI
  min_val <- h_min - h_ci

  # add comparison column to df
  df_pair$comp <- paste0(factor_a, "/", factor_b)

  if (save_output) {
    # set output
    grDevices::png(
      filename = file.path(out_dir, paste0(
        "plot_diff_", sub(" ", "_", tract), "_pair.png"
      )),
      width = 600, height = 600
    )

    # draw plot
    graphics::par(mar = c(5, 5, 4, 2), family = "Times New Roman")

    p_summary <- utils::capture.output(itsadug::plot_diff(
      gam_model,
      view = "nodeID",
      comp = comp,
      # sim.ci = TRUE,
      n.grid = 100,
      rm.ranef = TRUE,
      print.summary = TRUE,
      main = paste0(tract, "Difference Scores, ", factor_a, "-", factor_b),
      ylab = "Est. difference",
      xlab = "Tract Node",
      cex.lab = 2,
      cex.axis = 2,
      cex.main = 2,
      cex.sub = 1.5,
      col.diff = "red"
    ))

    # determine sig nodes
    if (p_summary[14] != "Difference is not significant.") {
      sig_regions <- p_summary[15:length(p_summary)]
      sig_regions <- gsub("\\t", "", sig_regions)

      # make list of start and end nodes, for shading
      sig_list <- as.list(strsplit(sig_regions, " - "))
      start_list <- as.numeric(sapply(sig_list, "[[", 1))
      end_list <- as.numeric(sapply(sig_list, "[[", 2))

      # shade significant regions
      for (h_ind in 1:length(start_list)) {
        graphics::polygon(
          x = c(rep(start_list[h_ind],2), rep(end_list[h_ind], 2)),
          y = c(0, min_val, min_val, 0),
          col = grDevices::rgb(1, 0, 0, 0.2),
          border = NA
        )
      }
    }

    graphics::par(mar = c(5, 4, 4, 2))
    grDevices::dev.off()
  }

  return(df_pair)
}


#' Plot tract profiles for each bundle as a facet and each metric as a figure.
#'
#' @param df Data frame.
#' @param metrics Name(s) of the metrics to plot per figure, character vector.
#'          By default, will be all diffusion metrics in the provided data frame.
#' @param bundles Name(s) of the tract bundles to plot per facet, character
#'          vector. By default, will be all tract bundles in the provided data
#'          frame.
#' @param bundles_col Name of the column in the provided data frame with the
#'          tract bundles.
#' @param group_col Name of the column in the data frame to group by as a color,
#'          character. By default, no grouping variable is provided.
#' @param line_func Line function that provides the line positioning. See
#'          \link[ggplot2]{stat_summary} for more information.
#' @param linewidth Line thickness of the tract profile line.
#' @param ribbon_func Ribbon function that provides the range for the ribbon.
#'          See \link[ggplot2]{stat_summary} for more information.
#' @param ribbon_alpha Ribbon alpha level.
#' @param n_groups Number of groups to split a numeric grouping variable.
#' @param pal_name Grouping color palette name, character. Default is colorblind.
#' @param out_dir Output directory of saved plots.
#' @param figsize Figure size. A numeric vector of (width, height) in inches.
#'
#' @return List of plot handles corresponding to the specified metrics.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df <- read_afq_sarica()
#'
#'plot_tract_profiles(
#'   df,
#'   metrics = c("fa"),
#'   bundles = c("Left Corticospinal", "Right Corticospinal"),
#'   group_col = "class",
#'   figsize   = c(8, 6)
#')
#'
#'plot_tract_profiles(
#'   df,
#'   metrics = c("fa"),
#'   bundles = c("Left Corticospinal", "Right Corticospinal"),
#'   group_col = "age",
#'   n_groups  = 3,
#'   pal_name  = "Spectral",
#'   figsize   = c(8, 6)
#')
#'}
plot_tract_profiles <- function (
    df,
    metrics      = NULL,
    bundles      = NULL,
    bundles_col  = "tractID",
    group_col    = NULL,
    line_func    = "mean",
    linewidth    = 1,
    ribbon_func  = "mean_cl_boot",
    ribbon_alpha = 0.25,
    n_groups     = 3,
    pal_name     = "colorblind",
    out_dir      = getwd(),
    figsize      = c(8, 11.5)
) {

  # argument preparation
  if (is.null(metrics)) {
    metrics <- names(df) # extract all column names from given data frame
    indx <- sapply(metrics, function(x) any(startsWith(x, c("csd", "dki", "dti", "fwdti"))))
    metrics <- metrics[indx]
  }

  if (is.null(bundles)) {
    bundles <- unique(df[[bundles_col]])
  }

  if (is.null(group_col)) {
    group_col <- "_group"
    df[[group_col]] <- "placeholder"
  }

  ribbon_func <- switch(
    ribbon_func,
    "mean_cl_boot"   = ggplot2::mean_cl_boot,
    "mean_cl_normal" = ggplot2::mean_cl_normal,
    "mean_sdl"       = ggplot2::mean_sdl,
    "median_hilow"   = ggplot2::median_hilow
  )

  # prepare data.frame for plotting
  plot_df <- df %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(metrics), names_to = "metric") %>%
    dplyr::rename(tracts = tidyselect::all_of(bundles_col)) %>%
    dplyr::filter(tracts %in% bundles, metric %in% metrics)

  # factorized grouping variable, split into groups if numeric
  if (is.numeric(plot_df[[group_col]])) {
    plot_df[[group_col]] <- Hmisc::cut2(plot_df[[group_col]], g = n_groups)
  } else if (is.character(plot_df[[group_col]])) {
    plot_df[[group_col]] <- forcats::fct(plot_df[[group_col]])
  }

  # declare color palette for plotting
  if (pal_name == "colorblind") {
    color_palette <- cbPalette
  } else if (pal_name == "colorblindblack") {
    color_palette = cbbPalette
  } else {
    n <- RColorBrewer::brewer.pal.info[pal_name, "maxcolors"]
    color_palette <- RColorBrewer::brewer.pal(n, pal_name)
  }

  plot_handles <- list() # initialize
  for (curr_metric in metrics) {
    # create current metric figure handle
    plot_handle <- plot_df %>%
      dplyr::filter(metric == curr_metric) %>%
      ggplot2::ggplot(ggplot2::aes_string(x = "nodeID",
                                          y = "value",
                                          group = group_col,
                                          color = group_col,
                                          fill = group_col)) +
                      ggplot2::stat_summary(
                          color = NA,
                          geom = "ribbon",
                          fun.data = ribbon_func,
                          alpha = ribbon_alpha) +
                      ggplot2::stat_summary(
                          geom = "line",
                          fun = line_func,
                          linewidth = linewidth) +
                      ggplot2::scale_x_continuous(name = "") +
                      ggplot2::scale_y_continuous(name = curr_metric) +
                      ggplot2::scale_color_manual(values = color_palette) +
                      ggplot2::scale_fill_manual(values = color_palette) +
                      ggplot2::facet_wrap(~ tracts) +
                      ggplot2::theme_bw()

    # remove legend if no group
    # if (group_col == "_group") {
    #   plot_handle <- plot_handle + ggplot::theme(legend.position = "none")
    # }

    # save tract profile figure
    plot_fname <- paste0("tract-profile_by-", group_col, "_",
                         stringr::str_replace_all(curr_metric, "_", "-"), ".png")
    ggplot2::ggsave(
      filename = file.path(out_dir, plot_fname),
      plot     = plot_handle,
      width    = figsize[1],
      height   = figsize[2],
      units    = "in",
      device   = "png"
    )

    # collect plot handles by metric
    plot_handles <- c(plot_handles, list(plot_handle))
  }

  # assign names to plot handles and return list
  names(plot_handles) <- metrics
  return(plot_handles)
}
