# Use colorblind-friendly palette from http://jfly.iam.u-tokyo.ac.jp/color/
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Global Variables
PALETTE_NAMES    <- c("cb", "cbb", "colorblind", "colorblind_black", 
                      rownames(RColorBrewer::brewer.pal.info))
RIBBON_FUNCTIONS <- c("mean_cl_boot", "mean_cl_normal", "mean_sdl", 
                      "median_hilow")

# Set Dummy Global Variables to appease R CMD check
utils::globalVariables(c("x", "ymin", "ymax", "group"))

#' Retrieve tract name from abbreviation
#'
#' @param tract_abbr AFQ tract abbreviation
#'
#' @return Formatted tract name
#'
#' @examples
#' full_name <- tract_name("OR")
#' full_name <- tract_name("CST_L")
#' @export
tract_name <- function(tract_abbr) {
  name <- switch(
    tract_abbr,
    "OR"    = "Optic Radiation",
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
    "FP"    = "Callosum Forceps Major",
    "FA"    = "Callosum Forceps Minor",
    "ILF_R" = "RightILF",
    "ILF_L" = "LeftILF",
    "SLF_R" = "RightSLF",
    "SLF_L" = "LeftSLF",
    tract_abbr
  )
  return(name)
}


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
#' @param output_dir directory in which to save plots
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' 
#' gam_fit <- fit_gam(
#'   df_afq,
#'   target = "dti_fa",
#'   covariates = list("group", "sex"),
#'   family = "gamma",
#'   k = 40
#' )
#' 
#' plot_gam_splines(
#'   gam_model = gam_fit,
#'   tract = "OR",
#'   df_tract = df_afq,
#'   dwi_metric = "dti_fa",
#'   covariates = c("group", "sex"),
#'   output_dir = ".")}
#' @export
plot_gam_splines <- function(
    gam_model,
    tract,
    df_tract,
    dwi_metric,
    covariates,
    group_by       = "group",
    participant_id = "subjectID",
    output_dir
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

  plot_filename <- file.path(output_dir,
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
#'   (e.g., set.seed(123)) for reproducible results. Note: in contrast with
#'   Gavin Simpson's code, here the Bayesian posterior covariance matrix of
#'   the parameters is uncertainty corrected (unconditional=TRUE) to reflect
#'   the uncertainty on the estimation of smoothness parameters.
#' @param output_dir Directory in which to save plots
#'
#' @return A dataframe of spline differences at each node
#' @export
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' 
#' gam_fit <- fit_gam(
#'   df_afq,
#'   target = "dti_fa",
#'   covariates = list("group", "sex"),
#'   family = "gamma",
#'   k = 40
#' )
#' 
#' df_diff <- spline_diff(
#'   gam_model = gam_fit,
#'   tract = "OR",
#'   factor_b = "1",
#'   output_dir = "."
#' )}
spline_diff <- function(
  gam_model,
  tract,
  group_by = "group",
  factor_a,
  factor_b,
  save_output = TRUE,
  sim.ci = FALSE,
  output_dir = getwd()
) {
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
      filename = file.path(output_dir, paste0(
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


#' Plot tract profiles
#' 
#' @description
#' Create tract profile figures for each tract as a facet and for each `y`
#' value as a figure.
#'
#' @param df              The input dataframe. 
#' @param y               Column name(s) of the variables to plot on the y-axis.
#' @param tracts          Name(s) of the tract tracts to plot per facet. 
#'                        Default: NULL \cr \cr
#'                        If NULL, will be all tracts in the data frame.
#' @param tract_col       The column name that encodes the tract information.
#'                        Default: tractID
#' @param node_col        The column name that encodes tract node positions.
#'                        Default: "nodeID".
#' @param group_col       Column name that encodes group information. Will be 
#'                        drawn by color. Default: NULL. 
#'                        \itemize{
#'                          \item  If `group_col` is a factor, will use unique
#'                                 values as groups.
#'                          \item  If `group_col` is a numeric, will use
#'                                 [ggplot2::cut_interval] to create equal range
#'                                 `n_groups` as groups.
#'                        }
#' @param n_groups        Number of groups to split a numeric grouping variable.
#'                        Only used if `group_col` is numeric. 
#'                        Default: NULL.
#' @param group_pal       Grouping color palette name. Valid options include: 
#'                        "cb" (colorblind), "cbb" (colorblind_black), and all 
#'                        named [RColorBrewer] palettes.
#'                        Default: "cb" (colorblind)
#' @param participant_col The column name that encodes participant ID.
#'                        Default: "subjectID". 
#' @param ribbon_func     Ribbon summarizing function that provides the y, ymin,
#'                        and ymax for the ribbon. See `fun.data` argument from
#'                        [ggplot2::stat_summary] for more information.
#'                        Default: "mean_cl_boot"
#' @param ribbon_alpha    Ribbon alpha level. Default: 0.25
#' @param linewidth       Line thickness of the tract profile line. Default: 1.
#' @param save_figure     Boolean flag. If TRUE, save tract profiles. If FALSE, 
#'                        do not save tract profiles. Default: FALSE
#' @param output_dir      Output directory for the figure image(s).
#'                        Default: [getwd] (current working directory). 
#' @param ...             Keyword arguments to be passed to [ggplot2::ggsave].
#'                      
#' @return Named list of plot handles corresponding to the specified y values.
#' 
#' @details
#' If `save_figure` is TRUE, the naming convention is as follows:
#' 
#' - If `group_col`, "tract_by-(group_col)_param-(y)_profiles.png".
#' - If `group_col == NULL`, "tract_param-(y)_profiles.png"
#' 
#'
#' @examples
#' \dontrun{
#' df_sarica <- read_afq_sarica(na_omit = TRUE)
#'
#' plot_handle <- plot_tract_profiles(
#'   df        = df,
#'   y         = c("fa", "md"),
#'   tracts    = c("Left Corticospinal", "Right Corticospinal"),
#'   width     = 12, 
#'   height    = 6,
#'   units     = "in"
#' )
#' 
#' plot_handle <- plot_tract_profiles(
#'   df        = df,
#'   y         = c("fa", "md"),
#'   tracts    = c("Left Corticospinal", "Right Corticospinal"),
#'   group_col = "group", 
#'   width     = 12, 
#'   height    = 6,
#'   units     = "in"
#' )
#'
#' plot_handle <- plot_tract_profiles(
#'   df        = df,
#'   y         = "fa",
#'   tracts    = c("Left Corticospinal", "Right Corticospinal"),
#'   group_col = "age", 
#'   n_groups  = 5, 
#'   group_pal = "Spectral", 
#'   width     = 12, 
#'   height    = 6,
#'   units     = "in"
#' )}
#' @export
plot_tract_profiles <- function (
    df,
    y, 
    tracts          = NULL,
    tract_col       = "tractID",
    node_col        = "nodeID", 
    group_col       = NULL,
    n_groups        = NULL,
    group_pal       = "cb",
    participant_col = "subjectID", 
    ribbon_func     = "mean_cl_boot",
    ribbon_alpha    = 0.25,
    linewidth       = 1,
    save_figure     = TRUE, 
    output_dir      = getwd(),
    ... 
) {
  # argument input control
  stopifnot("`df` must be a class data.frame or tibble" = 
    any(class(df) %in% c("data.frame", "tbl_df")))
  stopifnot("`y` must be a character" = is.character(y))
  if (!is.null(tracts)) {
    stopifnot("`tracts` must be a character" = is.character(tracts))
  }
  stopifnot("`tract_col` must be a character" = is.character(tract_col))
  stopifnot("There can be only one `tract_col`" = length(tract_col) == 1)
  stopifnot("`node_col` must be a character" = is.character(node_col))
  stopifnot("There can be only one `node_col`" = length(node_col) == 1)
  if (!is.null(group_col)) { # group_col is NOT null 
    stopifnot("`group_col` must be a character" = is.character(group_col))
    stopifnot("There can be only one `group_col`" = length(group_col) == 1)
    stopifnot("`group_col` must be a column in the dataframe `df`" =
      group_col %in% colnames(df))
    if (is.numeric(df[[group_col]])) {
      stopifnot("`n_groups` must be a integer" = round(n_groups) == n_groups)
    }
    rlang::arg_match(group_pal, values = PALETTE_NAMES)
  }
  stopifnot("`linewidth` must be a numeric" = is.numeric(linewidth))
  rlang::arg_match(ribbon_func, values = RIBBON_FUNCTIONS) 
  stopifnot("`ribbon_alpha` must be a numeric" = is.numeric(ribbon_alpha))
  stopifnot("`save_figure` must be a logical." = is.logical(save_figure))
  if (save_figure) {
    stopifnot("`output_dir` must exist." = fs::is_dir(output_dir))
  }

  # prepare tracts, if NULL all unique tracts
  if (is.null(tracts)) {
    tracts <- unique(df[[tract_col]])
  } 
   
  # prepare grouping column all to factors
  if (!is.null(group_col)) { # if not NULL
    group_values  <- df[[group_col]] # grouping values
    color_palette <- `_get_palette`(group_pal) # get group color palette
    if (is.numeric(group_values)) {
      df[[group_col]] <- ggplot2::cut_interval(group_values, n = n_groups)
    } else if (is.character(group_values)) {
      df[[group_col]] <- forcats::fct(group_values)
    }  
  }

  # prepare data frame for plotting
  keep_cols <- c(participant_col, tract_col, node_col, group_col, y)
  df_plot <- df[df[[tract_col]] %in% tracts, ] %>%
    dplyr::select(tidyselect::all_of(keep_cols)) %>% 
    dplyr::rename(
      x = tidyselect::all_of(node_col), 
      tracts = tidyselect::all_of(tract_col)
    )

  # prepare summarizing function 
  ribbon_func <- `_get_ribbon_func`(ribbon_func) # get ribbon function

  plot_handles <- list() # initialize
  for (y_curr in y) { # for each y-axis variable
    if (is.null(group_col)) {
      # prepare current y-axis values to plot
      df_curr <- df_plot %>% 
        dplyr::rename(y = tidyselect::all_of(y_curr)) %>% 
        dplyr::group_by(x, tracts) %>% 
        dplyr::summarize(ribbon_func(y), .groups = "drop")

      # create current metric figure handle
      plot_handle <- df_curr %>% 
        ggplot2::ggplot(ggplot2::aes(x = x, y = y, ymin = ymin, ymax = ymax)) +
        ggplot2::geom_ribbon(color = NA, alpha = ribbon_alpha) +
        ggplot2::geom_line(linewidth = linewidth) + 
        ggplot2::scale_x_continuous(name = "Node Position") + 
        ggplot2::scale_y_continuous(name = stringr::str_to_upper(y_curr)) + 
        ggplot2::facet_wrap(~ tracts) + 
        ggplot2::theme_bw()

      # prepare the saved figure file name
      output_fname <- sprintf("tracts_param-%s_profile.png", y_curr)
    } else {
      # prepare current y-axis values to plot
      df_curr <- df_plot %>% 
        dplyr::rename(
          y = tidyselect::all_of(y_curr), 
          group = tidyselect::all_of(group_col)
        ) %>% 
        dplyr::group_by(x, group, tracts) %>% 
        dplyr::summarize(ribbon_func(y), .groups = "drop")
  
      # create current metric figure handle
      plot_handle <- df_curr %>% 
        ggplot2::ggplot(ggplot2::aes(x = x, y = y, ymin = ymin, ymax = ymax, 
          group = group, color = group, fill = group)) +
        ggplot2::geom_ribbon(color = NA, alpha = ribbon_alpha) +
        ggplot2::geom_line(linewidth = linewidth) + 
        ggplot2::scale_x_continuous(name = "Node Position") + 
        ggplot2::scale_y_continuous(name = stringr::str_to_upper(y_curr)) + 
        ggplot2::scale_color_manual(name = group_col, values = color_palette) +
        ggplot2::scale_fill_manual(name = group_col, values = color_palette) +
        ggplot2::facet_wrap(~ tracts) + 
        ggplot2::theme_bw()
  
      # prepare the saved figure file name
      output_fname <- sprintf("tracts_by-%s_param-%s_profile.png", group_col, y_curr)
    }
      
    # save the figure 
    if (save_figure) {
      ggplot2::ggsave(
        filename = file.path(output_dir, output_fname),
        plot     = plot_handle,
        ...      = ...
      )  
    }

    # collect plot handles by metric
    plot_handles <- c(plot_handles, list(plot_handle))
  }

  # assign names to plot handles and return
  names(plot_handles) <- y

  return(plot_handles)
}


# Helper Functions that are NOT exported
`_get_ribbon_func` <- function(func_name) {
  rlang::arg_match(func_name, values = RIBBON_FUNCTIONS)
  ribbon_func <- switch(
    func_name,
    "mean_cl_boot"   = ggplot2::mean_cl_boot,
    "mean_cl_normal" = ggplot2::mean_cl_normal,
    "mean_sdl"       = ggplot2::mean_sdl,
    "median_hilow"   = ggplot2::median_hilow
  )
  return(ribbon_func)
}


`_get_palette` <- function(palette_name) {
  rlang::arg_match(palette_name, values = PALETTE_NAMES)
  color_palette <- switch(
    palette_name, 
    "cb"               = cbPalette, 
    "cbb"              = cbbPalette, 
    "colorblind"       = cbPalette, 
    "colorblind_black" = cbbPalette, 
    { # default RColorBrewer palettes
      n <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
      RColorBrewer::brewer.pal(n, palette_name)
    }
  )
  return(color_palette)
}