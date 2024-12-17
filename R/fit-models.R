#' Build a GAM formula.
#' 
#' @description
#' Function used to build a generic GAM formula string from the given arguments. 
#' The function automatically includes a tract node smoothing term (with or 
#' without a grouping term) and a participant random effects (random intercept). 
#' Any regressor terms are included as additive effects.
#' 
#' @param target          The column name that encodes the metric to model.
#' @param regressors      Column name or list of column names to use as 
#'                        regressors, not including nodes smoothing terms and 
#'                        the participant random effect. This list can also 
#'                        include smoothing terms. Default: NULL.
#' @param node_col        The column name that encodes tract node positions.
#'                        Default: "nodeID".
#' @param node_k          The basis dimensions used to represent the node 
#'                        smoother. If `node_group`, the basis value is applied
#'                        to the group as well. Default: 10. 
#' @param node_group      The column name to group the tract node smooth by.
#'                        Default: NULL.
#' @param participant_col The column name that encodes participant ID.
#'                        Default: "subjectID".
#'
#' @return A GAM formula string.
#' @export
#'
#' @examples
#' formula <- build_formula(target = "dti_fa", node_k = 40)
#' 
#' formula <- build_formula(target = "dki_md",
#'                          regressors = c("group", "sex"), 
#'                          node_k = 32, 
#'                          node_group = "group")
build_formula <- function(
  target, 
  regressors = NULL, 
  node_col = "nodeID", 
  node_k = 10, 
  node_group = NULL, 
  participant_col = "subjectID"
) {
  # argument input control
  stopifnot("`target` must be a character" = is.character(target))
  if (!is.null(regressors)) {
    stopifnot("`regressors` must be a character or character vector" = 
    is.character(regressors))
  }
  stopifnot("`node_col` must be a character" = is.character(node_col))
  stopifnot("`node_k` must be a numeric" = is.numeric(node_k))
  stopifnot("`node_k` must be a integer value" = (node_k %% 1) == 0)
  if (!is.null(node_group)) {
    stopifnot("`node_group` must be a character" = is.character(node_group))
    stopifnot("There can be only one `node_group`" = length(node_group) == 1)
  }
  stopifnot("`participant_col` must be a character" = 
    is.character(participant_col))

  # define node smooth term (with or without group)
  if (!is.null(node_group)) {
    node_smoother <- sprintf("s(%s, by = %s, k = %d)", 
      node_col, node_group, node_k)
  } else {
    node_smoother <- sprintf("s(%s, k = %d)", node_col, node_k)
  }

  # define random effects (intercept) of participant
  participant_random_effect <- sprintf("s(%s, bs = 're')", participant_col)

  # define formula depending on regressors
  if (!is.null(regressors)) { 
    # remove node and participant columns from regressor list
    regressors <- regressors[! regressors %in% c(node_col, participant_col)]

    # define regressor effect(s) as additive effects ONLY
    regressor_effects <- stringr::str_flatten(regressors, collapse = " + ")

    # define formula string with regressor effect(s)
    formula_string <- sprintf("%s ~ %s + %s + %s", target, regressor_effects, 
      node_smoother, participant_random_effect)
  } else {
    # define formula string without regressor effect(s)
    formula_string <- sprintf("%s ~ %s + %s", target, node_smoother, 
      participant_random_effect)
  }
  
  return(formula_string)
}


#' Fit a GAM for tract node metrics (e.g., FA, MD)
#'
#' @param formula         Explicit formula to use for the GAM.
#' @param df              The data frame contains GAM metrics.  
#' @param node_col        The column name that encodes tract node positions.
#'                        Default: "nodeID".
#' 
#' @param target          The column name that encodes the metric to model.
#' @param regressors      Column name or list of column names to use as 
#'                        regressors, not including nodes smoothing terms and 
#'                        the participant random effect. This list can also 
#'                        include smoothing terms. Default: NULL.
#' @param node_k          The basis dimensions used to represent the node 
#'                        smoother. If `node_group`, the basis value is applied
#'                        to the group as well. Default: 10. 
#' @param node_group      The column name to group the tract node smooth by
#'                        (i.e., `s(node_col, by = node_group, k = node_k)``).
#'                        Default: NULL.
#' @param participant_col The column name that encodes participant ID.
#'                        Default: "subjectID".
#' @param autocor         Whether to account for autocorrelation in the tract 
#'                        profile with the AR(1) autocorrelation model. 
#'                        Default: TRUE.
#' @param family          Name or function of the distribution to use for the 
#'                        GAM dependent variable.
#'                        \cr \cr
#'                        If name, possible values: {'auto', 'beta', 'gamma', 
#'                        'norm'}. If 'auto', will automatically determine the 
#'                        distribution of best fit between mgcv::betar ('beta'), 
#'                        stats::Gamma ('gamma'), or stats::gaussian ('norm').
#'                        Default: 'auto'
#'                        \cr \cr
#'                        See [stats::family()] or [mgcv::family.mgcv()] for 
#'                        more \code{family} or \code{extended.family} functions.
#' @param method          GAM fitting method, string; passed to [mgcv::bam()]. 
#'                        Default: "fREML"
#' @param discrete        With method = "fREML" it is possible to discretize 
#'                        covariates for storage and efficiency reasons. See 
#'                        [mgcv::bam()] for more information. 
#'                        Default: TRUE
#' @param ...             Further keyword arguments passed to [mgcv::bam()]
#' @return Fitted GAM model
#'
#' @details
#' This function has a series of steps:
#' \enumerate{
#'   \item If family == "auto", choose the distribution (either 'beta', 'gamma',
#'         or 'norm') that has the lowest AIC when fitting to the GAM dependent
#'         variable data using [fitdistrplus::fitdist()].
#'   \item If node_k == "auto", build an initial GAM model with k = 2 and
#'      continue to double the k value until gam.check shows that k is
#'      large enough.
#'   \item Fit a GAM model according to the formula given or build a formula 
#'         that follow the below pattern: 
#'      \cr \cr
#'      target ~ regressor1 (+ regressor2 + ...) +
#'               s(node_col, by = node_group, k = node_k) +
#'               s(participant_col, bs = "re")
#'      \cr
#'   \item Optionally, save the output of gam.check and summary to files.
#' }
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' 
#' df_tract <- df_afq %>% 
#'   filter(tractID == "CST_L")
#' 
#' gam_fit <- fit_gam(
#'   target = "dti_fa", 
#'   regressors = c("sex", "group"), 
#'   df = df_tract, 
#'   family = "gamma", 
#'   node_k = "auto"
#' )
#' }
#' @export
fit_gam <- function(...) {
  UseMethod("fit_gam")
}

#' @export
fit_gam.default <- function(
  formula, 
  df, 
  node_col = "nodeID", 
  autocor  = TRUE, 
  family   = "auto", 
  method   = "fREML",
  discrete = TRUE, 
  ...
) {
  # argument input control
  stopifnot("`formula` must be a class formulas" = class(formula) == "formula")
  stopifnot("`df` must be a class data.frame or tibble" = 
    any(class(df) %in% c("data.frame", "tbl_df")))
  stopifnot("`node_col` must be a character" = is.character(node_col))
  stopifnot("`autocor` must be a logical" = is.logical(autocor))
  if (is.character(family)) { family <- stringr::str_to_lower(family) }
  stopifnot("`family` must be a recognized family (auto, gamma, beta) or a family class function" = 
    any(family %in% c("auto", "gamma", "beta", "norm")) || 
    any(class(family) %in% c("family", "extended.family")))
  stopifnot("`method` must be a recognized method" = 
    any(method %in% c("GCV.Cp", "GACV.Cp", "REML", "P-REML", "fREML")))
  stopifnot("`discrete` must be a logical" = is.logical(discrete))
  
  # define GAM linkfamily function
  if (is.character(family)) {
    # if automatically determining linkfamily function 
    if (family == "auto") {
      values <- df[[all.vars(formula)[attr(stats::terms(formula), "response")]]]
      distr_aic <- sapply(c("beta", "gamma", "norm"), 
        function(x) fitdistrplus::fitdist(values, x)$aic)
      family <- names(distr_aic)[which.min(distr_aic)]
    }
    # set linkfamily based on family character
    if (family == "beta") {
      linkfamily <- mgcv::betar(link = "logit")
    } else if (family == "gamma") {
      linkfamily <- stats::Gamma(link = "logit")
    } else if (family == "gaussian" || family == "norm" ) {
      linkfamily <- stats::gaussian(link = "identity")
    }
  } else { # family function provided
    linkfamily <- family
  }

  if (autocor) {
    # define AR1.start as first nodeID position
    df$ar_start <- df[[node_col]] == 0

    # fit GAM model without autocorrelation 
    gam_fit_wo_rho <- mgcv::bam(
      formula  = formula, 
      family   = linkfamily, 
      data     = df, 
      method   = method, 
      rho      = 0, 
      AR.start = ar_start, 
      discrete = discrete, 
      ...      = ...
    )

    # determine autocorrelation parameter, rho
    rho <- itsadug::start_value_rho(gam_fit_wo_rho)

    # fit GAM model with autocorrelation 
    gam_fit <- mgcv::bam(
      formula  = formula, 
      family   = linkfamily, 
      data     = df, 
      method   = method, 
      rho      = 0, 
      AR.start = ar_start, 
      discrete = discrete, 
      ...      = ...
    )
  } else {
    # fit GAM model 
    gam_fit <- mgcv::bam(
      formula  = formula,
      family   = linkfamily, 
      data     = data, 
      method   = method, 
      discrete = discrete, 
      ...      = ...
    )
  }

  return(gam_fit)
}

#' @rdname fit_gam
#' @export
fit_gam.formula <- function(
  formula, 
  df, 
  node_col = "nodeID", 
  autocor  = TRUE, 
  family   = "auto", 
  method   = "fREML",
  discrete = TRUE, 
  ...
) {
  fit_gam.default(
    formula  = formula, 
    df       = df, 
    autocor  = autocor, 
    family   = family, 
    method   = method, 
    discrete = discrete, 
    ...      = ...
  )
}

#' @rdname fit_gam
#' @export
fit_gam.character <- function(
  target, 
  df, 
  regressors = NULL, 
  node_k = 10, 
  node_col = "nodeID", 
  node_group = NULL, 
  participant_col = "subjectID", 
  autocor = TRUE, 
  family = "auto", 
  method = "fREML",
  discrete = TRUE, 
  ...
) {
  # use build formula to create formula string
  formula <- tractable::build_formula(
    target = target, 
    regressors = regressors, 
    node_k = node_k, 
    node_col = node_col, 
    node_group = node_group, 
    participant_col = participant_col
  )

  # call default fit_gam() method
  fit_gam.default(
    formula  = as.formula(formula, env = .GlobalEnv),
    df       = df, 
    autocor  = autocor, 
    family   = family, 
    method   = method,
    discrete = discrete, 
    ...      = ...
  )
}

#     if (k == "auto") {
#       # Initial k value. This will be multiplied by 2 in the while loop
#       k_model <- 4

#       # Initialize k check results to enter the while loop
#       k_indices <- list(0.0, 0.0)
#       k_pvals <- list(0.0, 0.0)

#       while (any(k_indices < 1.0) | any(k_pvals <= 0.05)) {
#         k_model <- k_model * 2
#         formula <- build_formula(target = target,
#                                  covariates = covariates,
#                                  smooth_terms = smooth_terms,
#                                  group_by = group_by,
#                                  participant_id = participant_id,
#                                  k = k_model)

#         # Fit the gam
#         gam_fit <- mgcv::bam(
#           stats::as.formula(formula),
#           data = df_tract,
#           family = linkfamily,
#           method = method,
#           ... = ...
#         )

#         k_check <- mgcv::k.check(gam_fit)
#         k_indices <- stats::na.omit(k_check[, "k-index"])
#         k_pvals <- stats::na.omit(k_check[, "p-value"])
#        }
#     } else {
#       k_model <- k
#     }
#     formula <- build_formula(target = target,
#                              covariates = covariates,
#                              smooth_terms = smooth_terms,
#                              group_by = group_by,
#                              participant_id = participant_id,
#                              k = k_model)
#   }




# save_gam_outputs <- function(gam_fit, out_dir, tract_name){
#     utils::capture.output(
#       mgcv::gam.check(gam_fit, rep = 500),
#       file = file.path(out_dir, paste0(
#         "k_check_gam_", family, "_", sub(" ", "_", tract_name), ".txt"
#       ))
#     )
#     gam_summary <- summary(gam_fit)
#     utils::capture.output(
#       gam_summary,
#       file = file.path(out_dir, paste0(
#         "fit_summary_gam_", family, "_", sub(" ", "_", tract_name), ".txt"
#       ))
#     )

#     utils::write.csv(
#       t(gam_summary$p.table[, "Pr(>|t|)"]),
#       file.path(out_dir, paste0(
#         "gam_pvals_", family, "_", sub(" ", "_", tract_name), ".csv"
#         )))

#   }
  