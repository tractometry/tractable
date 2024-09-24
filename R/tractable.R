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
#' @export
#'
#' @examples
#' sarica <- read_afq_sarica()
#' sarica$group <- factor(sarica$class)
#' sarica$subjectID <- unclass(factor(sarica$subjectID))
#' tractable_single_bundle(df_afq = sarica,
#'                      tract = "Right Corticospinal",
#'                      participant_id = "subjectID",
#'                      group_by = "group",
#'                      covariates = c("age","group"),
#'                      dwi_metric = "fa")
tractable_single_bundle <- function(df_afq,
                                 tract,
                                 dwi_metric,
                                 participant_id = "subjectID",
                                 group_by = "group",
                                 covariates = c(group_by),
                                 smooth_terms = NULL,
                                 k = "auto",
                                 family = "auto",
                                 ...) {
  selected <- select_bundle(
    df_afq = df_afq,
    tract = tract,
    dwi_metric = dwi_metric,
    covariates = covariates,
    participant_id = participant_id,
    group_by = group_by)

  df_tract <- selected$df_tract

  gam_fit <- fit_gam(df_tract = df_tract,
                     target = dwi_metric,
                     covariates = covariates,
                     smooth_terms = smooth_terms,
                     group_by = group_by,
                     participant_id = participant_id,
                     k = k,
                     family = family,
                     ... = ...)

  return(gam_fit)
}
