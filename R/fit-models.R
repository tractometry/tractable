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

estimate_distribution <- function(
  x, distr = c("beta", "gamma", "norm"), eval_metric = "aic"
) {
  distr_eval <- sapply(distr, 
    function(d) fitdistrplus::fitdist(x, d)[[eval_metric]])
  distr_name <- names(distr_eval)[which.min(distr_eval)]
  return(distr_name)
}

#' Estimate smoothing basis dimensions for GAM smoothers.
#' 
#' @description
#' 
#' @param target          The column name that encodes the metric to model.
#' @param regressors      Column name or list of column names to use as 
#'                        regressors, not including nodes smoothing terms and 
#'                        the participant random effect. This list can also 
#'                        include smoothing terms. Default: NULL.
#' @param family          The family
#' @return description
#' @export 
estimate_smooth_basis <- function(
  target, 
  df, 
  regressors = NULL, 
  family = "auto", 
  ...
) {
  # input argument control 
  smooth_terms <- list(...)

  # determine family function 
  if (family == "auto") {
    family <- estimate_distribution(df[[target]], 
      distr = c("beta", "gamma", "norm"))
  }

  smooth_terms <- list(
    age = list(k_end = 10), 
    nodeID = list(k_start = 4, k_end = 50, bs = "tp")
  )
  term_names <- names(smooth_terms)
  all(term_names %in% names(df))

  default_args <- list(
    k_start = 1, 
    k_end = 50, 
    k_step = 1, 
    bs = "tp", 
    k_index = 0.5, 
    p_value_thr = 0.05
  )
  
  df_k <- tibble::tibble() # initialize results
  for (i in 1:length(smooth_terms)) { # for each smooth term
    curr_term <- term_names[i] # current smoothing term name

    curr_args <- default_args # copy default arguments
    curr_args[names(smooth_terms[[i]])] <- smooth_terms[[i]]
    
    k_index <- 0.0
    p_value <- 1.0
    curr_k <- curr_args$k_start
    while (threshold_has_not_been_reached) {
      curr_smooth <- sprintf("s(%s, k = %d)", curr_term, curr_k)

      formula <- sprintf("%s ~ %s + %s", target, 
        stringr::str_c(regressors, sep = " + "), curr_smooth)
      
      model_fit <- mgcv::bam(
        formula = as.formula(formula, env = .GlobalEnv), 
        df = df, 
        family = linkfamily, 
        method = method,
        discrete = discrete
      )

      k_check <- mgcv::k.check(model_fit)
      k_index <- stats::na.omit(kcheck[, "k-index"])
      p_value <- stats::na.omit(k_check[, "p-value"])
    }

    # collect term, k, k-index, and p-value
    curr_row <- tibble::tibble(
      term = curr_term, 
      k = curr_k, 
      k_index = k_index, 
      p_value = p_value, 
      regressor = if_else(curr_k == 1, curr_term, 
        sprintf("s(%s, bs = %s, k = %d)", curr_term, curr_args$bs, curr_k))
    )
    
    if (nrow(df) < 1) {
      df_k <- curr_row
    } else {
      df_k <- df_k %>% tibble::add_row(curr_row)
    }
    
  }
}


#' Fit a GAM model.
#' 
#' @description 
#' Function used to fit a GAM formula from an explicit formula or builds a 
#' formula from the given arguments. 
#'
#' @param formula         Explicit formula to use for the GAM.
#' @param df              The data frame contains GAM metrics.  
#' @param node_col        The column name that encodes tract node positions.
#'                        Default: "nodeID".
#' @param target          The column name that encodes the metric to model.
#' @param regressors      Column name or list of column names to use as 
#'                        regressors, not including nodes smoothing terms and 
#'                        the participant random effect. This list can also 
#'                        include smoothing terms. Default: NULL.
#' @param node_k          The basis dimensions used to represent the node 
#'                        smoother. If `node_group`, the basis value is applied
#'                        to the group as well. Default: 'auto'. 
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
  # if automatically determining linkfamily function 
  if (family == "auto") {
    values <- df[[all.vars(formula)[attr(stats::terms(formula), "response")]]]
    family <- estimate_distribution(values, 
      distr = c("beta", "gamma", "norm"))
  }

  # call fit_gam.default
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
  node_k = "auto", 
  node_col = "nodeID", 
  node_group = NULL, 
  participant_col = "subjectID", 
  autocor = TRUE, 
  family = "auto", 
  method = "fREML",
  discrete = TRUE, 
  ... = ...
) {
  # if automatically determining linkfamily function 
  if (family == "auto") {
    family <- estimate_distribution(df[[target]], 
      distr = c("beta", "gamma", "norm"))
  }

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

#' DESCRIPTION
#' 
#' @description
#' A short description...
#' 
#' @param gam_model
#' @param file
#' @param model_rdata
#' @param model_summary
#' @param model_check 
#' @return None. 
#' 
#' @example
#' 
#' @export
save_gam <- function(
  gam_model, 
  file, 
  model_rdata   = TRUE, 
  model_summary = TRUE,
  model_check   = FALSE 
) {
  # argument input control
  stopifnot("`gam_model` must be a class gam or bam" = 
    any(class(gam_model) %in% c("bam", "gam")))
  stopifnot("`model_rdata` must be a logical" = is.logical(model_rdata))
  stopifnot("`model_summary` must be a logical" = is.logical(model_summary))
  stopifnot("`model_check` must be a logical" = is.logical(model_check))

  # save model as Rdata file
  if (model_rdata) {
    saveRDS(gam_model, file = file)
  }

  # save model summary as text file
  if (model_summary) {
    save_name <- stringr::str_replace(
      file, ".(RData|Rdata|rdata|rda)$", "_Summary.txt")
    utils::capture.output(summary(gam_model), file = save_name)
  }

  # save gam.check output (text and figures)
  if (model_check) {
    save_name <- stringr::str_replace(
      file, ".(RData|Rdata|rdata|rda)$", "_GamCheck.txt")
    utils::capture.output(
      mgcv::gam.check(gam_model, rep = 500, plot = FALSE), 
      file = save_name,
    )
  }
}
