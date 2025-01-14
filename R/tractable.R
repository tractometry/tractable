#' Analyze group differences in a single dMRI tract profile using GAMs
#'
#' @param df         Input data frame.
#' @param tract      Tract name to subset from the dataframe.
#' @param tract_col  The column name that encodes the tract information.
#'                   Default: tractID
#' @inheritDotParams fit_gam 
#' 
#' @return GAM model fit object.
#' 
#' @examples
#' df_sarica <- read_afq_sarica(na_omit = TRUE)
#'
#' tractable_single_tract(
#'   df         = df_sarica, 
#'   tract      = "Right Corticospinal",
#'   target     = "fa",
#'   regressors = c("age", "group"),
#'   node_group = "group"
#' )
#' @export
tractable_single_tract <- function(
  df, 
  tract, 
  tract_col = "tractID",
  ...
) {
  # argument input control
  stopifnot("`df` must be a class data.frame or tibble" = 
    any(class(df) %in% c("data.frame", "tbl_df")))
  stopifnot("`tract` must be a character" = is.character(tract))
  stopifnot("`tract_col` must be a character" = is.character(tract_col))
  stopifnot("`tract_col` must be a column in the dataframe `df`" = 
    tract_col %in% names(df))
  stopifnot("`tract` must be a valid option inside `tract_col`" = 
    any(df[[tract_col]] == tract))

  # subset the full dataset by the given tract name
  df_tract <- df[df[[tract_col]] == tract,]

  # extract kwargs for fit_gam versus other kwargs
  vargs <- list(...) # listify kwargs
  fit_gam_func <- ifelse(
    "target" %in% names(vargs), fit_gam.default, fit_gam.formula)
  fit_gam_kwargs <- `_get_function_kwargs`(fit_gam_func, vargs)
  other_kwargs <- vargs[setdiff(names(vargs), names(fit_gam_kwargs))]

  # define fit_gam arguments
  fit_gam_kwargs <- c(list(df = df_tract), fit_gam_kwargs, other_kwargs)

  # fit the model and return model fit
  model_fit <- do.call(fit_gam_func, fit_gam_kwargs)

  return(model_fit)
}