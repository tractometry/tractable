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
#' @return A GAM formula. 
#'
#' @examples
#' formula <- build_formula(target = "dti_fa", node_k = 40)
#' 
#' formula <- build_formula(target = "dki_md",
#'                          regressors = c("group", "sex"), 
#'                          node_k = 32, 
#'                          node_group = "group")
#' @export
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

  # define random effects (intercept and node smooth) of participant
  participant_random_effect <- sprintf("s(%s, bs = 're')", participant_col)
  participant_node_smooth <- sprintf("s(%s, by = %s, bs = 'fs')", 
    node_col, participant_col)

  # define additional regressors (additive effects ONLY)
  if (!is.null(regressors)) { 
    # remove node and participant columns from regressor list
    regressors <- regressors[! regressors %in% c(node_col, participant_col)]
  } 

  # define formula (in global environment)
  formula <- stats::reformulate(
    termlabels = c(regressors, node_smoother, participant_random_effect,
                   participant_node_smooth), 
    response = target, 
    env = .GlobalEnv
  )
  
  return(formula)
}

#' Estimate distribution function. 
#' 
#' @param x           A numeric vector that will be evaluated. 
#' @param distr_names A vector of distribution names to evaluate the values. \cr
#'                    Possible options: ('beta', 'gamma', 'norm') \cr
#'                    Default: c("beta", "gamma", "norm")
#' @param eval_metric The distribution evaluation metric names. \cr
#'                    Possible options: ('aic', 'bic', 'loglik') \cr
#'                    Default: "aic"
#' @details
#' See [fitdistrplus::fitdist()] for additional information.
#' 
#' @return The best fitting family function to the values. 
#' 
#' @examples
#' x <- rnorm(1000)
#' estimate_distribution(x)
#' 
#' x <- rgamma(1000, shape = 12)
#' estimate_distribution(x)
#' @export
estimate_distribution <- function(
  x, 
  distr_names = c("beta", "gamma", "norm"),
  eval_metric = "aic"
) {
  # argument input control 
  stopifnot("`x` must be a numeric" = is.numeric(x))
  stopifnot("`distr_names` must be a recognized distribution name" = 
    is.character(distr_names) && 
    all(distr_names %in% c("beta", "gamma", "norm")))
  stopifnot("`eval_metric` must be a recognized evauluation metric" = 
    is.character(eval_metric) && length(eval_metric) == 1 && 
    any(eval_metric %in% c("aic", "bic", "loglik")))

  # estimate distribution function based on evaluation metric
  distr_eval <- sapply(distr_names, function(d) {
    tryCatch({invisible(utils::capture.output(
      y <- fitdistrplus::fitdist(x, d)[[eval_metric]]))
      y
    }, error = function(conds) {
      NA_real_
    })
  })
  
  # determine best distribution based on evaluation metric 
  eval_func  <- if (eval_metric == "loglik") which.max else which.min
  best_distr <- names(distr_eval)[eval_func(distr_eval)]

  # based on evaluated distribution, return family function
  if (best_distr == "beta") {
    linkfamily <- mgcv::betar(link = "logit")
  } else if (best_distr == "gamma") {
    linkfamily <- stats::Gamma(link = "logit")
  } else if (best_distr == "norm") {
    linkfamily <- stats::gaussian(link = "identity")
  }
  return(linkfamily)
}


#' Estimate smoothing basis dimensions for GAM smoothers.
#' 
#' @description 
#' asdfasdfa
#' 
#' @param target          The column name that encodes the metric to model.
#' @param regressors      Column name or list of column names to use as 
#'                        regressors, not including nodes smoothing terms and 
#'                        the participant random effect. This list can also 
#'                        include smoothing terms. Default: NULL.
#' @param family          The family
#' @return description
#' @export 
estimate_smooth_basis <- function(...) {
  UseMethod("estimate_smooth_basis")
}

#' @export
estimate_smooth_basis.default <- function() {

}


#' @export
estimate_smooth_basis.formula <- function(
  formula, 
  df, 
  kindex_thr = 0.95, 
  pvalue_thr = 0.05, 
  k_values   = 1:50,
  bs         = "ts",
  family     = "auto", 
  method     = "fREML",
  discrete   = TRUE, 
  ...
) {

  formula <- as.formula("fa ~ sex + 
    s(age, k = c(1, 3, 5)) + 
    s(nodeID, bs = 'fs') + 
    s(participantID, bs = 're')")

  #TODO: consider interaction smoothers, and ti, te
  formula <- as.formula("fa ~ sex + 
    s(nodeID, age, bs = 'tp') + 
    s(participantID, bs = 're')")

  # define default arguments for estimate smooth basis
  default_args <- list(
    k_values = k_values,
    bs = bs, 
    kindex_thr = kindex_thr, 
    pvalue_thr = pvalue_thr
  )

  # extract formula terms
  formula_terms <- labels(terms(formula))
    
  # determine which terms have smoothers from the formula 
  f_terms <- stringr::str_replace_all(formula_terms, "\\s", "")
  smooth_index <- (stringr::str_detect(f_terms, "^s\\(.+\\)$") &
    stringr::str_detect(f_terms, "bs=.re.", negate = TRUE))
  
  # define the target, regressor, and smoothing terms
  target_term     <- all.vars(formula)[attr(stats::terms(formula), "response")]
  regressor_terms <- f_terms[!smooth_index]
  smooth_terms    <- f_terms[smooth_index]

  # for each smooth term
  for (i in 1:length(smooth_terms)) {
    # extract information from the current smooth term
    # curr_term <- stringr::str_replace_all(smooth_terms[i], "^s\\(|\\)$", "")
    # curr_term <- stringr::str_split(curr_term, ",")[[1]]
    # curr_term

    # curr_variable <- curr_term[1]

    # curr <- smooth_terms[1]
    # call <- substitute(expr(curr))
    

    # curr_args <- default_args # initialize
    # if (length(curr_term) > 1) { # if additional smoothing terms
    #   # for each additional argument
    #   for (term in curr_term[2:length(curr_term)]) {
    #     term <- stringr::str_split(term, "=")[[1]] # (key, value)
    #     if (term[1] == "k") {
    #       curr_args$k_values <- eval(parse(text = term[2]))
    #     } else {
    #       curr_args[[term[1]]] <- term[2]
    #     }
    #   }
    # } 

    # declare initial values of k_index and p_value
    k_index <- 0.0 # desire towards 1.0
    p_value <- 0.0 # desire towards 1.0
    while (k_index < curr_args$kindex_thr && p_value < curr_args$pvalue_thr) {


      # fit current model with regressors and current smooth term
      curr_model <- mgcv::bam(
        formula  = as.formula(formula, env = .GlobalEnv), 
        df       = df, 
        family   = linkfamily, 
        method   = method,
        discrete = discrete
      )

      # call k.check on model fit
      k_check <- mgcv::k.check(curr_model)

      # redefine k.check metrics for all smooth coefficients
      smooth_rowname <- sprintf("s(%s)", curr_variable)
      k_index <- k_check[smooth_rowname, "k-index"]
      p_value <- k_check[smooth_rowname, "p-value"]
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
#'                        (i.e., `s(node_col, by = node_group, k = node_k)`).
#'                        Default: NULL.
#' @param participant_col The column name that encodes participant ID.
#'                        Default: "subjectID".
#' @param autocor         Whether to account for autocorrelation in the tract 
#'                        profile with the AR(1) autocorrelation model. 
#'                        Default: TRUE.
#' @param family          Name or function of the distribution to use for the 
#'                        GAM dependent variable.
#'                        \cr \cr
#'                        If name, possible values: ('auto', 'beta', 'gamma', 
#'                        'norm'). If 'auto', will automatically determine the 
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
      distr_names = c("beta", "gamma", "norm"))
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
  ...
) {
  # if automatically determining linkfamily function 
  if (family == "auto") {
    family <- tractable::estimate_distribution(df[[target]], 
      distr_names = c("beta", "gamma", "norm"))
  }

  # use build formula to create formula string
  formula <- build_formula(
    target = target, 
    regressors = regressors, 
    node_k = node_k, 
    node_col = node_col, 
    node_group = node_group, 
    participant_col = participant_col
  )

  # call default fit_gam() method
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

#' Save GAM model outputs.
#' 
#' @description
#' This function saves the GAM model as RData (.rda) file, the model summary as
#' a text (.txt) file, and the [mgcv::gam.check()] figures as a PNG (.png) file.
#' 
#' @param gam_model asdfasdf
#' @param file asdfasdf
#' @param model_rdata asdfasdf
#' @param model_summary asdfasdf
#' @param model_check  asdfasdf
#' @return None. 
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
