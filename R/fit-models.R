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

#' Estimate smoothing parameter (k) for smoothing terms
#' 
#' This function
#' @param name description
#' @param name description
#' 

#' Fit a GAM for tract node metrics (e.g., FA, MD)
#'
#' This function has a series of steps:
#' \itemize{
#'   \item If family == "auto", choose the distribution (either beta or gamma)
#'      that has the lowest AIC when fitting to the dMRI metric data
#'   \item If k == "auto", build an initial GAM model with k = 16 and
#'      continue to double the k value until gam.check shows that k is
#'      large enough
#'   \item Fit a GAM model such that: \cr \cr
#'      target ~ covariates +
#'               s(nodeID, by=group, k = k_value) +
#'               s(subjectID, bs = "re")
#'      \cr
#'   \item Optionally save the output of gam.check and summary to files.
#' }
#'
#' @param target The diffusion metric to model (e.g. "FA", "MD"). If this is
#'  set, `formula` must NOT be set.
#' @param regressors List of strings of GAM regressors, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#'     If this is set, `formula` must NOT be set.
#' @param smooth_terms Smoothing terms, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#'     If this is set, `formula` must NOT be set.
#' @param node_group The grouping variable used to group nodeID smoothing terms
#'     If this is set, `formula` must NOT be set.
#' @param df_tract AFQ Dataframe of node metric values for single tract
#' @param participant_col The name of the column that encodes participant ID
#'     If this is set, `formula` must NOT be set.
#' @param formula Optional explicit formula to use for the GAM. If provided,
#'     this will override the dynamically generated formula build from the
#'     target and covariate inputs. Default = NULL. If this is set, all other
#'     inputs that determine the formula must be set to NULL.
#' @param node_k Dimension of the basis used to represent the node smoothing term,
#'     If k = 'auto' (default), this function will attempt to find the best
#'     value.
#' @param autocor Whether to account for autocorrelation in the tract profile.
#'     Default = TRUE.
#' @param family Distribution to use for the gam. Must be either 'gamma',
#'     'beta', or 'auto'. If 'auto', this function will select the best fit
#'     between beta and gamma distributions.
#' @param method String, fitting method passed to mgcv::bam
#' @param ... Further keyword arguments passed to mgcv::bam
#'
#' @return Fitted GAM model
#' @export
#'
#' @examples
#' \dontrun{
#' df_afq <- read.csv("/path/to/afq/output.csv")
#' tract <- "CST_L"
#' df_tract <- df_afq[which(df_afq$tractID == tract), ]
#' gam_fit <- fit_gam(df_tract,
#'                    target = "dti_fa",
#'                    covariates = list("group", "sex"),
#'                    family = "gamma",
#'                    k = "auto")
#' }

# define the s3 method 
fit_gam <- function(...) {
  UseMethod("fit_gam")
}

fit_gam.default <- function(
  formula, 
  df, 
  autocor = TRUE, 
  family = "auto", 
  method = "fREML",
  discrete = TRUE, 
  ...
) {
  # argument input control
  stopifnot("`autocor` must be a logical" = is.logical(autocor))
  stopifnot("`df` must be a class data.frame or tibble" = 
    any(class(df) %in% c("data.frame", "tbl_df")))
  if (is.character(family)) { family <- stringr::str_to_lower(family) }
  stopifnot("`family` must be a recognized family (auto, gamma, beta) or a family class function" = 
    any(family %in% c("auto", "gamma", "beta")) || 
    any(class(family) %in% c("family", "extended.family")))

  # if determining linkfamily function 
  if (family == "auto") {
    values <- df[[all.vars(formula)[attr(stats::terms(formula), "response")]]]
    distr_aic <- sapply(c("beta", "gamma", "norm"), 
      function(x) fitdistrplus::fitdist(values, x)$aic)
    family <- names(distr_aic)[which.min(distr_aic)]
  } 
  
  # define GAM linkfamily function
  if (family == "beta") {
    linkfamily <- mgcv::betar(link = "logit")
  } else if (family == "gamma") {
    linkfamily <- stats::Gamma(link = "logit")
  } else if (family == "gaussian" || family == "norm" ) {
    linkfamily <- stats::gaussian(link = "identity")
  } else { # family function provided
    linkfamily <- family
  }

  if (autocor) {
    # define AR1.start as first nodeID position
    df <- df %>% dplyr::mutate(ar_start = nodeID == 0)

    # fit GAM model without autocorrelation 
    gam_fit_wo_rho <- mgcv::bam(
      formula = formula,
      family = linkfamily, 
      data = df, 
      method = method, 
      rho = 0, 
      AR.start = ar_start, 
      discrete = discrete, 
      ... = ...
    )

    # determine autocorrelation parameter, rho
    rho <- itsadug::start_value_rho(gam_fit_wo_rho)

    # fit GAM model with autocorrelation 
    gam_fit <- mgcv::bam(
      formula = formula,
      family = linkfamily, 
      data = df, 
      method = method, 
      rho = rho, 
      AR.start = ar_start, 
      discrete = discrete, 
      ... = ...
    )
  } else {
    # fit GAM model 
    gam_fit <- mgcv::bam(
      formula = formula,
      family = linkfamily, 
      data = data, 
      method = method, 
      discrete = discrete, 
      ... = ...
    )
  }

  return(gam_fit)
}

fit_gam.formula <- function(...) {
  fit_gam.default(...)
}

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
    formula  = as.formula(formula), 
    df       = df, 
    autocor  = autocor, 
    family   = family, 
    method   = method,
    discrete = discrete, 
  )
}


# fit_gam <- function(
#   df_tract,
#   target = NULL,
#   regressors = NULL,
#   smooth_terms = NULL,
#   participant_id = NULL,
#   formula = NULL,
#   node_k = NULL,
#   autocor = TRUE,
#   family = "auto",
#   method = "fREML",
#   ...
# ) {
#   # Check that if formula is non-NULL, all the other formula-setting inputs
#   # are null
#   if (!is.null(formula)) {
#     if (!(is.null(target) &&
#           is.null(covariates) &&
#           is.null(smooth_terms) &&
#           is.null(group_by) &&
#           is.null(k))) {
#     stop(
#     "If `formula` is provided no other formula-setting input may be provided")
#     } else{
#       # If it's a string input, we'll cast it into a formula
#     if (is.character(formula) & length(formula) == 1) {
#       formula <- as.formula(formula)
#     }
#     # Get the target from the formula for use below
#     target <- terms(formula)[[2]]
#   }}

#   # Set other defaults
#   if (is.null(group_by)) {
#     group_by <- "group"
#   }
#   if (is.null(participant_id)) {
#     participant_id <- "subjectID"
#   }

#   if (is.null(k)) {
#     k <- "auto"
#   }


#   if (is.null(formula)) {
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
  