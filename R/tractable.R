#' Analyze group differences in a single dMRI tract profile using GAMs
#'
#' @param df_afq Input AFQ dataframe. If NULL, this function will load data
#'     using read_afq_files and the additional arguments in ...
#' @param tract Abbreviated tract name, e.g., "CST_L" or "OR"
#' @param dwi_metric The diffusion metric to model (e.g. "FA", "MD")
#' @param participant_id The name of the column that encodes participant ID
#' @param group_by The grouping variable used to group nodeID smoothing terms
#' @param covariates List of strings of GAM covariates,
#'     not including the smoothing terms over nodes and the random effect due
#'     to subjectID.
#' @param smooth_terms Smoothing terms, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#' @param k Dimension of the basis used to represent the node smoothing term,
#'     If k = 'auto', this function will attempt to find the best value
#' @param family Distribution to use for the gam. Must be either 'gamma',
#'     'beta', or 'auto'. If 'auto', this function will select the best fit
#'     between beta and gamma distributions.
#' @param ... Arguments to pass to fit_gam
#'
#' @examples
#' df_sarica <- read_afq_sarica()
#'
#' tractable_single_tract(
#'   df = df_sarica,
#'   tract = "Right Corticospinal",
#'   participant_id = "subjectID",
#'   group_by = "group",
#'   covariates = c("age","group"),
#'   dwi_metric = "fa"
#' )
#' @export
tractable_single_tract <- function(
  tract, 
  df, 
  tract_col = "tractID",
  ...
) {
  # argument input control
  stopifnot("`tract` must be a character" = is.character(tract))
  stopifnot("`tract_col` must be a character" = is.character(tract_col))
  stopifnot("`tract_col` must be a column in the dataframe `df`" = 
    tract_col %in% df)
  stopifnot("`tract` must be a valid option inside `tract_col`" = 
    any(df[[tract_col]] == tract))

  # subset the full dataset by the given tract name
  df_tract <- df[df[[tract_col]] == tract,]

  # call fit_gam 
  tract_fit <- fit_gam(
    df = df_tract, 
    ...
  )

  return(gam_fit)
}
