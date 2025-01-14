#' Shuffle an AFQ dataframe
#'
#' @description
#' This function shuffles participants' demographic information (i.e., age, 
#' group, sex), thereby destroying correlations between participants' tract 
#' profiles and phenotypic data.
#'
#' @param df              The input dataframe. 
#' @param target          The column name that encodes the metric to model.
#' @param shuffle_cols    Column names that should be shuffled.
#' @param node_col        The column name that encodes tract node positions.
#'                        Default: "nodeID"
#' @param node_group      The column name to group the tract node smooth by.
#'                        Default: NULL.
#' @param tract_col       The column name that encodes tract names.
#'                        Default: "tractID"
#' @param participant_col The column name that encodes participant ID.
#'                        Default: "subjectID".
#' @param sample_uniform  Boolean flag. If TRUE, shuffling should sample 
#'                        uniformly from the unique values in the columns. If 
#'                        FALSE, shuffling will shuffle without replacement.
#'
#' @return A shuffled AFQ dataframe
#'
#' @examples
#' \dontrun{
#' df_afq <- read_csv("/path/to/afq/output.csv")
#' df_shuffled <- shuffle_df(df_afq, target = "dti_fa")}
#' @export
shuffle_df <- function(
  df,
  target,
  shuffle_cols    = NULL,
  node_col        = "nodeID", 
  node_group      = NULL,
  tract_col       = "tractID", 
  participant_col = "subjectID",
  sample_uniform  = FALSE
) {
  # argument input control 
  stopifnot("`df` must be a class data.frame or tibble" = 
    any(class(df) %in% c("data.frame", "tbl_df")))
  stopifnot("`target` must be a character" = is.character(target))
  if (!is.null(shuffle_cols)) {
    stopifnot("`shuffle_cols` must be a character" = is.character(shuffle_cols))
  }
  stopifnot("`node_col` must be a character" = is.character(node_col))
  if (!is.null(node_group)) {
    stopifnot("`node_group` must be a character" = is.character(node_group))
  }
  stopifnot("`tract_col` must be a character" = is.character(tract_col))
  stopifnot("`participant_col` must be a character" = is.character(participant_col))
  stopifnot("`sample_uniform` must be a logical" = is.logical(sample_uniform))

  # pivot data frame to one row per participant
  df_wide <- tidyr::pivot_wider(
    data        = df,
    names_from  = tidyselect::all_of(node_col),
    values_from = tidyselect::all_of(target)
  )

  # if not given, determine shuffle columns
  original_colnames <- colnames(df)
  if (is.null(shuffle_cols)) {
    static_cols <- c(participant_col, node_col, tract_col, node_group, target)
    shuffle_cols <- original_colnames[!original_colnames %in% static_cols]
  }

  # shuffle participants' shuffle_cols and the grouping variable
  for (svar in unique(c(shuffle_cols, node_group))) {
    x <- df_wide[[svar]] # current values to shuffle
    if (sample_uniform) {
      # sample uniformly from the unique values (with replacement)
      df_wide[[svar]] <- sample(unique(x), length(x), replace = TRUE)
    } else {
      # sample and shuffle the existing values
      df_wide[[svar]] <- sample(x, length(x))
    }
  }

  # return to long format (one row per node)
  df_shuffled <- tidyr::pivot_longer(
    data      = df_wide,
    cols      = tidyselect::all_of(as.character(unique(df[[node_col]]))),
    names_to  = node_col, 
    values_to = target
  ) %>% 
    dplyr::select(tidyselect::all_of(original_colnames)) 

  # format column class to match original
  for (var in original_colnames) {
    class(df_shuffled[[var]]) <- class(df[[var]])
  }

  return(df_shuffled)
}


#' Bootstrap an AFQ dataframe
#'
#' @description
#' This function bootstrap samples an AFQ dataframe by participant. That is, it 
#' first pivots to wide format with one row per participant, bootstrap samples, 
#' and finally pivots back to long format.
#'
#' @param df              The input dataframe. 
#' @param target          The column name that encodes the metric to model.
#' @param node_col        The column name that encodes tract node positions.
#'                        Default: "nodeID"
#' @param node_group      The column name to group the tract node smooth by.
#'                        Default: NULL.
#' @param participant_col The column name that encodes participant ID.
#'                        Default: "subjectID".
#' 
#' @return A shuffled AFQ dataframe
#' @export
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' df_bootstrap <- bootstrap_df(df_afq, "dti_fa")}
bootstrap_df <- function(
  df,
  target,
  node_col        = "nodeID", 
  node_group      = "group",
  participant_col = "subjectID"
) {
  # argument input control 
  stopifnot("`df` must be a class data.frame or tibble" = 
    any(class(df) %in% c("data.frame", "tbl_df")))
  stopifnot("`target` must be a character" = is.character(target))
  stopifnot("`node_col` must be a character" = is.character(node_col))
  if (!is.null(node_group)) {
    stopifnot("`node_group` must be a character" = is.character(node_group))
  }
  stopifnot("`participant_col` must be a character" = is.character(participant_col))

  # pivot data frame to one row per participant
  df_wide <- tidyr::pivot_wider(
    data        = df,
    names_from  = tidyselect::all_of(node_col), 
    values_from = tidyselect::all_of(target)
  ) %>% 
    dplyr::slice_sample(prop = 1, replace = TRUE)

  # determine columns not used for pivoting
  original_colnames <- colnames(df)
  static_cols <- original_colnames[!original_colnames %in% c(node_col, target)]

  # return to long format (one row per node)
  df_bootstrap <- tidyr::pivot_longer(
    data      = df_wide,
    cols      = -tidyselect::all_of(static_cols),
    names_to  = node_col,
    values_to = target
  ) %>% 
    dplyr::select(tidyselect::all_of(original_colnames)) 

  # format column class to match original
  for (var in original_colnames) {
    class(df_bootstrap[[var]]) <- class(df[[var]])
  }

  return(df_bootstrap)
}


#' Perform repeated sampling tests on an AFQ dataframe.
#'
#' When permutation_test == FALSE (the default), this function bootstrap
#' samples from an AFQ dataframe and returns pairwise differences at each node
#' for each bootstrap sample. These results can then be used to construct
#' bootstrap confidence intervals for the node-wise differences.
#'
#' When permutation_test == TRUE, this function simulates the null distribution
#' using permutation testing. That is, it shuffles to destroy any relationships
#' between covariates and dwi_metrics. It then computes node-wise differences
#' for each shuffled sample.
#'
#' @param df AFQ Dataframe of node metric values for single tract
#' @param n_samples Number of sample tests to perform
#' @param target The diffusion metric to model (e.g. "FA", "MD")
#' @param tract AFQ tract name
#' @param group_by The grouping variable used to group nodeID smoothing terms
#' @param participant_col The name of the column that encodes participant ID
#' @param sample_uniform Boolean flag. If TRUE, shuffling should sample
#'     uniformly from the unique values in the columns. If FALSE, shuffling
#'     will shuffle without replacement.
#' @param covariates List of strings of GAM covariates, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#' @param smooth_terms Smoothing terms, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#' @param k Dimension of the basis used to represent the node smoothing term
#' @param family Distribution to use for the gam. Must be 'gamma' or 'beta'
#' @param formula Optional explicit formula to use for the GAM. If provided,
#'     this will override the dynamically generated formula build from the
#'     target, covariate, and k inputs. Default = NULL.
#' @param factor_a First group factor, string
#' @param factor_b Second group factor, string
#' @param permute Boolean flag. If TRUE, perform a permutation test.
#'     Otherwise, do a bootstrap simulation.
#'
#' @return Dataframe with bootstrap or permutation test coefficients
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' df_tract <- df_afq[which(df_afq$tractID == tract), ]
#' bootstrap_coefs <- sampling_test(df_afq,
#'   dwi_metric = "dti_fa",
#'   covariates = list("group", "sex"),
#'   family = "gamma",
#'   k = 40,
#'   n_samples = 1000
#' )
#' }
#' #@export
# sampling_test <- function(
#   df_tract,
#   n_samples,
#   target,
#   tract,
#   node_col = "nodeID", 
#   node_group = "group",
#   participant_col = "subjectID",
#   sample_uniform = FALSE,
#   covariates = NULL,
#   smooth_terms = NULL,
#   k = NULL,
#   family = NULL,
#   formula = NULL,
#   factor_a = NULL,
#   factor_b = NULL,
#   permute = FALSE
# ) {
#   if (group_by %in% covariates) {
#     coefs <- vector(mode = "list", length = n_samples)
#     pvalues <- vector(mode = "list", length = n_samples)
#   }
#   node_diffs <- data.frame(
#     nodeID = numeric(0),
#     est = numeric(0),
#     permIdx = numeric(0)
#   )

#   pb <- progress::progress_bar$new(total = n_samples)
#   for (idx in 1:n_samples) {
#     pb$tick()

#     if (permute) {
#       df_shuffle <- shuffle_df(
#         df              = df,
#         target          = target,
#         shuffle_cols    = covariates,
#         node_col        = node_col, 
#         node_group      = node_group,
#         participant_col = participant_col,
#         sample_uniform  = sample_uniform
#       )
#     } else {
#       df_shuffle <- bootstrap_df(
#         df              =  df,
#         target          = target,
#         node_group      = node_group,
#         participant_col = participant_col
#       )
#     }

#     gam_shuffle <- fit_gam(
#       df_tract = df_shuffle,
#       target = dwi_metric,
#       covariates = covariates,
#       smooth_terms = smooth_terms,
#       group_by = group_by,
#       participant_id = participant_id,
#       formula = formula,
#       k = k,
#       family = family
#     )
#     ff <- summary(gam_shuffle)
#     pvalues[[idx]] <- ff$p.table[, "Pr(>|t|)"][[paste0(group_by, factor_b)]]
#     if (group_by %in% covariates) {
#       coef_name <- grep(paste0("^", group_by),
#         names(gam_shuffle$coefficients),
#         value = TRUE
#       )

#       coefs[[idx]] <- gam_shuffle$coefficients[[coef_name]]
#     }

#     if (!is.null(factor_a) & !is.null(factor_b)) {
#       df_pair <- spline_diff(
#         gam_model = gam_shuffle,
#         tract = tract,
#         group_by = group_by,
#         factor_a,
#         factor_b,
#         save_output = FALSE,
#         out_dir = NULL
#       )

#       df_pair$permIdx <- idx
#       df_pair <- dplyr::select(df_pair, c("nodeID", "est", "permIdx"))
#       node_diffs <- rbind(node_diffs, df_pair)
#     }
#   }

#   df_sampling_test <- tidyr::pivot_wider(
#     node_diffs,
#     names_from = "nodeID",
#     names_prefix = "node_",
#     values_from = "est",
#     id_cols = "permIdx"
#   )

#   if (group_by %in% covariates) {
#     df_sampling_test$group_coefs <- unlist(coefs)
#   }

#   df_sampling_test$pvalue <- unlist(pvalues)
#   return(df_sampling_test)
# }
