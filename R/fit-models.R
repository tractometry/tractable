# Global Variables
FAMILY_FUNCTION_NAMES <- c("beta", "gamma" , "gaussian", "norm")
MGCV_METHODS <- c("GCV.Cp", "GACV.Cp", "REML", "P-REML", "fREML")

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
#'                        regressors, not including node smoothing terms and 
#'                        the participant effects. This list can also include 
#'                        smoothing terms. Default: NULL.
#' @param node_col        The column name that encodes tract node positions.
#'                        Default: "nodeID".
#' @param node_k          The basis dimensions used to represent the node 
#'                        smoother. If \code{node_group}, the basis value is 
#'                        applied to the group as well. Default: NULL
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
#' formula <- build_formula(
#'   target     = "dki_md",
#'   regressors = c("group", "sex"), 
#'   node_k     = 32, 
#'   node_group = "group"
#' )
#' @export
build_formula <- function(
  target, 
  regressors      = NULL, 
  node_col        = "nodeID", 
  node_k          = NULL, 
  node_group      = NULL, 
  participant_col = "subjectID"
) {
  # argument input control
  stopifnot("`target` must be a character" = is.character(target))
  if (!is.null(regressors)) {
    stopifnot("`regressors` must be a character or character vector" = 
    is.character(regressors))
  }
  stopifnot("`node_col` must be a character" = is.character(node_col))
  if (!is.null(node_k)) {
    stopifnot("`node_k` must be a numeric" = is.numeric(node_k))
    stopifnot("`node_k` must be a integer value" = (node_k %% 1) == 0)
  } 
  if (!is.null(node_group)) {
    stopifnot("`node_group` must be a character" = is.character(node_group))
    stopifnot("There can be only one `node_group`" = length(node_group) == 1)
  }
  stopifnot("`participant_col` must be a character" = 
    is.character(participant_col))

  # define node smooth term (with or without group)
  if(!is.null(node_group) && !is.null(node_k)) {
    node_smoother <- sprintf("s(%s, by = %s, bs = 'fs', k = %d)", 
      node_col, node_group, node_k)
  } else if (!is.null(node_group)) {
    node_smoother <- sprintf("s(%s, by = %s, bs = 'fs')", node_col, node_group)
  } else if (!is.null(node_k)) {
    node_smoother <- sprintf("s(%s, k = %d)", node_col, node_k)
  } else {
    node_smoother <- sprintf("s(%s)", node_col)
  }

  # define random effects (intercept) of participant
  participant_random_effect <- sprintf("s(%s, bs = 're')", participant_col)

  # define additional regressors (additive effects ONLY)
  if (!is.null(regressors)) { 
    # remove node and participant columns from regressor list
    regressors <- regressors[! regressors %in% c(node_col, participant_col)]
  } 

  # define formula (in global environment)
  formula <- stats::reformulate(
    termlabels = c(regressors, node_smoother, participant_random_effect), 
    response = target, 
    env = .GlobalEnv
  )
  
  return(formula)
}


#' Estimate distribution function from values.
#' 
#' @param x           A numeric vector that will be evaluated. 
#' @param distr_names A vector of distribution names to evaluate the values. \cr
#'                    Possible options: ("beta", "gamma", "gaussian") \cr
#'                    Default: c("beta", "gamma", "gaussian")
#' @param eval_metric The distribution evaluation metric names. \cr
#'                    Possible options: ("aic", "bic", "loglik") \cr
#'                    Default: "aic"
#' 
#' @return The estimated \code{family} or \code{extended.family} function.
#' 
#' @seealso [fitdistrplus::fitdist]
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
  distr_options = c("beta", "gamma", "gaussian"),
  eval_metric   = "aic"
) {
  # argument input control 
  stopifnot("`x` must be a numeric" = is.numeric(x))
  stopifnot("`distr_options` must be a character class" = 
    is.character(distr_options))
  rlang::arg_match(distr_options, values = FAMILY_FUNCTION_NAMES, multiple = TRUE)
  stopifnot("`eval_metric` must be a recognized evauluation metric" = 
    is.character(eval_metric) && length(eval_metric) == 1)
  rlang::arg_match(eval_metric, values = c("aic", "bic", "loglik"))

  # edit "gaussian" distribution name to "norm"
  if (any(distr_options == "gaussian")) {
    distr_options[distr_options == "gaussian"] <- "norm"
  }

  # estimate distribution function based on evaluation metric
  distr_eval <- sapply(distr_options, function(d) {
    tryCatch({invisible(utils::capture.output(
      y <- fitdistrplus::fitdist(x, d)[[eval_metric]]))
      y
    }, error = function(conds) {
      NA_real_
    })
  })

  # error condition: all NAs for considered distributions
  stopifnot("All distributions considered evaluated with errors or NAs. 
    Please provide different distriubutions to `distr_options`." = 
    !all(is.na(distr_eval)))

  # determine best distribution based on evaluation metric 
  eval_func  <- ifelse(eval_metric == "loglik", which.max, which.min)
  best_distr <- names(distr_eval)[eval_func(distr_eval)]
  
  # based on evaluated distribution, return family function
  linkfamily <- `_get_family_function`(best_distr)

  return(linkfamily)
}


#' Estimate smoothing basis dimensions for GAM smoothers.
#' 
#' @description 
#' Function used to estimate smoothing basis, \code{k}, for each smooth term. 
#'
#' @param target         The column name that encodes the metric to model.
#' @param smooth_terms   List of smooth terms to estimate smoothing basis. See 
#'                       details for examples of smoothing terms.
#' @param df             The data frame that contains the GAM metrics.  
#' @param regressors     Column name or list of column names to use as 
#'                       regressors. This list can also include smoothing terms.
#'                       Default: NULL.
#' @param k_values       A list of k values to consider. Default: 1:50
#' @param bs             The name of the default smoothing basis. Default: "tp"
#' @param kindex_thr     The k-index threshold. Default: 0.95
#' @param pvalue_thr     The p-value threshold. Default: 0.05
#' @param family         Name or family function of the distribution to use for
#'                       modeling the GAM dependent variable. 
#'                       \itemize{
#'                         \item If name, the possible values: ("auto", "beta", 
#'                               "gamma", "gaussian").
#'                         \item If "auto", will automatically determine the 
#'                               distribution of best fit between 
#'                               [mgcv::betar] ("beta"), [stats::Gamma] 
#'                               ("gamma"), or [stats::gaussian] ("gaussian").
#'                         \item If function, see [family][stats::family]
#'                               or [family.mgcv][mgcv::family.mgcv] for 
#'                               more \code{family} or \code{extended.family} 
#'                               class functions.
#'                       } 
#' @param method         GAM fitting method passed to [bam][mgcv::bam()]. 
#'                       Default: "fREML"
#' @param discrete       With \code{method} is "fREML" it is possible to 
#'                       discretize covariates for storage and efficiency 
#'                       reasons. See [bam][mgcv::bam] for more information. 
#'                       Default: TRUE
#' @param ...            Further keyword arguments to be passed to 
#'                       [bam][mgcv::bam]
#' 
#' @details 
#' ## Smooth terms specification
#' Smooth terms can be specified as: \cr
#' 
#' \code{s(x)} \cr
#' The smooth term will be estimated with the defaults from \code{k_values} and
#' \code{bs}. \cr \cr
#' \code{s(x, k = 1:10, bs = 'cp')} \cr
#' The smooth term will be estimated with \code{k} values 1 to 10 and a basis 
#' set of 'cp'. \cr \cr
#' \code{s(x, by = group, k = c(2, 7))} \cr 
#' The smooth term will estimate with \code{k} value of 2 and 7 and using the
#' \code{by} variable \code{group}. \cr \cr
#' \code{s(x, y, bs = 'fs', m = 3)} \cr
#' The smooth term over two variables, \code{x} and \code{y}, will be estimated 
#' with the default \code{k_values} with additional arguments of basis set of 
#' 'fs' and \code{m} of 3. \cr
#' 
#' Not shown are other mgcv smoothers, such as [te][mgcv::te], [ti][mgcv::ti], 
#' and [t2][mgcv::t2], which are also available.
#' 
#' ## Estimation process
#' For each \code{smooth_term}, the function will iteratively fit a GAM model
#' following the formula while incrementing through \code{k_values}:
#' 
#' \tabular{lll}{ \tab \tab \code{target ~ regressor_terms + smooth_term} }
#' 
#' where \code{target} is the dependent variable, \code{regressor_terms} are the
#' additive effects that should be accounted for while estimating the smoothing
#' term, and \code{smooth_term} is the smoothing term that is currently being 
#' estimated.
#' 
#' ## Stopping criterion
#' The estimation process has two stopping criterion: 
#'
#' \itemize{
#'   \item The procedure will stop once the k-index value exceeds 
#'         \code{kindex_thr} \strong{AND} the p-value exceeds the 
#'         \code{pvalue_thr}. \cr
#'   \item If the thresholds are not met, the procedure will stop once it runs 
#'         through all of the \code{k_values}. \cr
#' }
#' 
#' @return A list containing two items: \tabular{llll}{
#'   \tab \code{est_terms} \tab \tab List of smoothing terms with "best" estimated 
#'                         smoothing basis. \cr
#'   \tab \code{k_estimates} \tab \tab A data frame with contains all estimated 
#'                         smoothing terms and corresponding k-index and p-values. \cr
#' }
#' @seealso [choose.k][mgcv::choose.k]
#' @export
estimate_smooth_basis <- function(...) {
  UseMethod("estimate_smooth_basis")
}

#' @rdname estimate_smooth_basis
#' @export
estimate_smooth_basis.default <- function(
  target, 
  smooth_terms, 
  df, 
  regressors = NULL, 
  k_values   = 1:50,
  bs         = "tp",
  kindex_thr = 0.95, 
  pvalue_thr = 0.05, 
  family     = "auto", 
  method     = "fREML",
  discrete   = TRUE, 
  ...
) {
  # argument input control 
  stopifnot("`target` must be a character" = is.character(target))
  stopifnot("`smooth_terms` must be a character vector" = is.character(smooth_terms))
  stopifnot("`df` must be a class data.frame or tibble" = 
    any(class(df) %in% c("data.frame", "tbl_df")))
  if (!is.null(regressors)) {
    stopifnot("`regressors` must be a character class" = is.character(regressors))
  }
  stopifnot("`k_values` must be numeric values" = is.numeric(k_values))
  stopifnot("`kindex_thr` must be a numeric value" = 
    length(kindex_thr) == 1 && is.numeric(kindex_thr))
  stopifnot("`pvalue_thr` must be a numeric value" = 
    length(pvalue_thr) == 1 && is.numeric(pvalue_thr))
  rlang::arg_match(bs, values = c("tp", "ts", "ds", "cr", "cs", "cc", "sos", 
    "bs", "ps", "cp", "mrf", "gp", "so", "sw", "sf", "ad", "sz", "fs"))
  if (is.character(family)) { 
    family <- stringr::str_to_lower(family) 
    rlang::arg_match(family, values = c("auto", FAMILY_FUNCTION_NAMES))
  } else {
    stopifnot("`family` must be a family class function" = 
      class(family) %in% c("family", "extended.family"))
  }
  rlang::arg_match(method, values = MGCV_METHODS)
  stopifnot("`discrete` must be a logical" = is.logical(discrete))

  # prepare regressor terms 
  regressor_terms <- regressors
  if (!is.null(regressor_terms)) {
    regressor_terms <- stringr::str_c(regressor_terms, collapse = " + ")
  }

  # prepare default arguments for estimate smooth basis
  default_args <- list(
    k_values   = k_values,
    kindex_thr = kindex_thr, 
    pvalue_thr = pvalue_thr, 
    kwargs     = list("bs" = bs)
  )

  # prepare link family function
  if (is.character(family) && family == "auto") {
    linkfamily <- estimate_distribution(df[[target]])
  } else if (is.character(family)) {
    linkfamily <- `_get_family_function`(family)
  } else {
    linkfamily <- family
  }

  # estimate smoothing basis for each smooth term
  k_estimates <- tibble::tibble() # create empty tibble
  est_smooth_terms <- rep("", length(smooth_terms)) # initialize
  for (i in 1:length(smooth_terms)) { # for each smooth term
    # extract information from the current smooth term (using a defused call)
    defused_call     <- rlang::parse_expr(smooth_terms[i])
    smooth_func      <- rlang::call_name(defused_call)
    smooth_args      <- rlang::call_args(defused_call)
    smooth_arg_names <- names(smooth_args)

    # define current smooth variable(s)
    smooth_vars <- as.character(smooth_args[smooth_arg_names == ""]) 
    by_index <- smooth_arg_names == "by"
    if (any(by_index)) { # if by argument, concatename variable to smooth_vars
      by_var <- as.character(smooth_args[by_index])
      smooth_vars <- c(smooth_vars, stringr::str_c("by = ", by_var))
    }
    smooth_vars <- stringr::str_flatten_comma(smooth_vars)

    # define smooth function additional arguments
    curr_args <- default_args # initialize with default arguments
    var_index <- smooth_arg_names %in% c("", "by") # variable args to skip
    if (length(smooth_args) > sum(var_index)) { # if additional smoothing arguments
      for (j in (sum(var_index) + 1):length(smooth_args)) {
        curr_key <- smooth_arg_names[j] # extract argument name
        curr_values <- eval(smooth_args[[j]]) # evaluate argument value
        if (curr_key == "k") { # if current value is `k`
          curr_args$k_values <- sort(curr_values)
        } else { # else assign argument to kwargs list
          curr_args$kwargs[[curr_key]] <- eval(smooth_args[[j]]) 
        }
      }
    } 

    # define k_check smoothing term rowname(s)
    smooth_vars_strip <- stringr::str_replace_all(smooth_vars, "\\s", "")
    if (stringr::str_detect(smooth_vars, "by = ")) {
      s_vars <- stringr::str_replace(smooth_vars_strip, sprintf(",by=%s", by_var), "")
      by_levels <- levels(df_sarica[[by_var]]) # extract group by levels
      smooth_rowname <- sprintf("%s(%s):%s%s", smooth_func, s_vars, by_var, by_levels)
    } else { # else, no group by levels to consider
      smooth_rowname <- sprintf("%s(%s)", smooth_func, smooth_vars_strip)
    }

    # start estimate smoothing basis paramter, k
    k_values <- curr_args$k_values # current k values
    for (j in 1:length(k_values)) { # for each k value
      # define smooth term k value and additional kwargs
      curr_kwargs   <- curr_args$kwargs # initialize
      curr_kwargs$k <- curr_args$k_values[j] # increment k value
      curr_kwargs   <- lapply(curr_kwargs, function(x) { 
        ifelse(is.character(x), sprintf("'%s'", x), x) })

      # define current smoothing variable
      curr_smooth_term <- sprintf("%s(%s, %s)",
        smooth_func, smooth_vars, 
        stringr::str_c(names(curr_kwargs), " = ", curr_kwargs, collapse = ", ")
      )

      # define formula with regressor terms and ONLY current smoothing term
      curr_formula <- stats::reformulate(
        termlabels = c(regressor_terms, curr_smooth_term), 
        response = target, 
        env = .GlobalEnv
      )

      # fit current model with regressors and current smooth term
      suppressWarnings(
        curr_model <- mgcv::bam(
          formula  = curr_formula, 
          data     = df, 
          family   = linkfamily, 
          method   = method,
          discrete = discrete
        )
      )

      # call k.check on model fit
      k_check <- mgcv::k.check(curr_model)

      # redefine k.check metrics for all smooth coefficients, averaged by levels
      k_index <- mean(k_check[smooth_rowname, "k-index"], na.rm = TRUE)
      p_value <- mean(k_check[smooth_rowname, "p-value"], na.rm = TRUE)

      # store k estimates into data frame
      k_estimates <- k_estimates %>% 
        dplyr::bind_rows(tibble::tibble_row(
          term = curr_smooth_term, k_index = k_index, p_value = p_value))

      # stopping condition based on k-index and p-value
      if (k_index > curr_args$kindex_thr && 
        p_value > curr_args$pvalue_thr) break;
    }

    # collect "best" estimated smoothing term
    est_smooth_terms[i] <- curr_smooth_term
  }

  return(list(est_terms = est_smooth_terms, k_estimates = k_estimates))
}


#' @rdname estimate_smooth_basis
#' @export
estimate_smooth_basis.formula <- function(
  formula, 
  df, 
  k_values   = 1:50,
  bs         = "tp",
  kindex_thr = 0.95, 
  pvalue_thr = 0.05, 
  family     = "auto", 
  method     = "fREML",
  discrete   = TRUE, 
  ...
) {
  # argument input control
  stopifnot("`formula` must be a formula class" = class(formula) == "formula")

  # extract formula terms
  formula_terms <- `_get_formula_terms`(formula)
    
  # determine which terms have smoothers from the formula 
  formula_terms <- stringr::str_replace_all(formula_terms, "\\s", "")
  smooth_index  <- (stringr::str_detect(formula_terms, 
    "^(s|te|ti|t2)\\(.+\\)$") &
    stringr::str_detect(formula_terms, "bs=['\"]re['\"]", negate = TRUE))
  
  # define the target, regressor, and smoothing terms
  target       <- `_get_formula_target`(formula)
  regressors   <- formula_terms[!smooth_index]
  smooth_terms <- formula_terms[smooth_index]

  # prepare regressors if empty)
  if (length(regressors) == 0) { regressors <- NULL }

  estimate_smooth_basis.default(
    target       = target, 
    smooth_terms = smooth_terms, 
    regressors   = regressors, 
    df           = df, 
    k_values     = k_values, 
    kindex_thr   = kindex_thr, 
    pvalue_thr   = pvalue_thr, 
    bs           = bs,
    family       = family,
    method       = method,
    ...          = ...
  )
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
#' @return Fit GAM model
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

#' @rdname fit_gam
#' @export
fit_gam.default <- function(
  target, 
  df, 
  regressors      = NULL, 
  node_col        = "nodeID", 
  node_k          = "auto", 
  node_group      = NULL, 
  participant_col = "subjectID", 
  autocor         = TRUE, 
  family          = "auto", 
  method          = "fREML",
  discrete        = TRUE, 
  ...
) {
  # argument input control
  stopifnot("`target` must be a character" = is.character(target))
  stopifnot("There can be only one `target`" = length(target) == 1)
  stopifnot("`df` must be a class data.frame or tibble" = 
    any(class(df) %in% c("data.frame", "tbl_df")))
  if (!is.null(regressors)) {
    stopifnot("`regressors` must be a character" = is.character(regressors))
  }
  stopifnot("`node_col` must be a character" = is.character(node_col))
  stopifnot("`node_col` must be a column in the dataframe." = 
    node_col %in% names(df))
  if (is.character(node_k)) {
    stopifnot("`node_k` must be 'auto' if character" = node_k == "auto")
  } else {
    stopifnot("`node_k` must be a numeric" = is.numeric(node_k))
    stopifnot("`node_k` must be a integer value" = (node_k %% 1) == 0)
  }
  if (!is.null(node_group)) {
    stopifnot("`node_group` must be a character" = is.character(node_group))
    stopifnot("There can be only one `node_group`" = length(node_group) == 1)
  }
  stopifnot("`participant_col` must be a character" = is.character(participant_col))
  stopifnot("`autocor` must be a logical" = is.logical(autocor))
  if (is.character(family)) { 
    family <- stringr::str_to_lower(family) 
    rlang::arg_match(family, values = c("auto", FAMILY_FUNCTION_NAMES))
  } else {
    stopifnot("`family` must be a family class function" = 
      class(family) %in% c("family", "extended.family"))
  }
  rlang::arg_match(method, values = MGCV_METHODS)
  stopifnot("`discrete` must be a logical" = is.logical(discrete))

  # prepare variables arguments
  vargs <- list(...) # listify variable argument
  mgcv_kwargs <- `_get_function_kwargs`(mgcv::bam, vargs)

  # define GAM linkfamily function
  if (is.character(family) && family == "auto") {
    linkfamily <- estimate_distribution(df[[target]])
  } else if (is.character(family)) {
    linkfamily <- `_get_family_function`(family)
  } else { # family function provided
    linkfamily <- family
  }

  # if estimating node smoother
  if (node_k == "auto") {
    # define default formula with given variables
    base_formula <- build_formula(
      target          = target, 
      regressors      = regressors, 
      node_col        = node_col, 
      node_group      = node_group, 
      participant_col = participant_col
    )

    # extract node term (smoother(node_col)) and regressor terms (everything else)
    formula_terms <- `_get_formula_terms`(base_formula)
    node_index <- stringr::str_starts(formula_terms, sprintf(".+\\(%s", node_col))
    node_terms <- formula_terms[node_index]
    regressors <- formula_terms[!node_index]
    if (length(regressors) == 0) { regressors <- NULL }

    # prepare estimate_smooth_basis arguments
    est_kwargs <- `_get_function_kwargs`(estimate_smooth_basis.default, vargs)
    est_kwargs <- c(list(
      target       = target, 
      smooth_terms = node_terms, 
      df           = df, 
      regressors   = regressors, 
      family       = linkfamily, 
      method       = method, 
      discrete     = discrete
    ), est_kwargs)

    # estimate node term smoothing basis
    node_terms <- do.call(estimate_smooth_basis.default, est_kwargs)$est_terms

    # reformulate the formula with estimated node smoother
    formula_terms[node_index] <- node_terms
    formula <- stats::reformulate(
      termlabels = formula_terms, 
      response = target, 
      env = .GlobalEnv
    )

  } else {
    # define GAM formula with the provided inputs
    formula <- build_formula(
      target     = target, 
      regressors = regressors, 
      node_col   = node_col, 
      node_k     = node_k, 
      node_group = node_group, 
      participant_col = participant_col
    )
  }

  # start model fitting process (with or without autocorrelations)
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
      ...      = mgcv_kwargs
    )

    # determine autocorrelation parameter, rho
    rho <- itsadug::start_value_rho(gam_fit_wo_rho)

    # fit GAM model with autocorrelation 
    gam_fit <- mgcv::bam(
      formula  = formula, 
      family   = linkfamily, 
      data     = df, 
      method   = method, 
      rho      = rho, 
      AR.start = ar_start, 
      discrete = discrete, 
      ...      = mgcv_kwargs
    )
  } else {
    # fit GAM model without autocorrelations
    gam_fit <- mgcv::bam(
      formula  = formula,
      family   = linkfamily, 
      data     = df, 
      method   = method, 
      discrete = discrete, 
      ...      = mgcv_kwargs
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
  # argument input control
  stopifnot("`formula` must be a formula class" = class(formula) == "formula")
  stopifnot("`df` must be a class data.frame or tibble" = 
    any(class(df) %in% c("data.frame", "tbl_df")))
  stopifnot("`node_col` must be a character" = is.character(node_col))
  stopifnot("`node_col` must be a column in the dataframe." = 
    node_col %in% names(df))
  stopifnot("`autocor` must be a logical" = is.logical(autocor))
  if (is.character(family)) { 
    family <- stringr::str_to_lower(family) 
    rlang::arg_match(family, values = c("auto", FAMILY_FUNCTION_NAMES))
  } else {
    stopifnot("`family` must be a family class function" = 
      class(family) %in% c("family", "extended.family"))
  }
  rlang::arg_match(method, values = MGCV_METHODS)
  stopifnot("`discrete` must be a logical" = is.logical(discrete))

  # prepare variables arguments
  vargs <- list(...) # listify variable argument
  mgcv_kwargs <- `_get_function_kwargs`(mgcv::bam, vargs)

  # define GAM linkfamily function
  if (is.character(family) && family == "auto") {
    target <- `_get_formula_target`(formula)
    linkfamily <- estimate_distribution(df[[target]])
  } else if (is.character(family)) {
    linkfamily <- `_get_family_function`(family)
  } else { # family function provided
    linkfamily <- family
  }

  # extract node term (smoother(node_col)) and regressor terms (everything else)
  formula_target <- `_get_formula_target`(formula)
  formula_terms  <- `_get_formula_terms`(formula)
  node_index <- stringr::str_starts(formula_terms, sprintf(".+\\(%s", node_col))
  node_terms <- formula_terms[node_index]
  regressors <- formula_terms[!node_index]
  if (length(regressors) == 0) { regressors <- NULL }

  # stop if more than one node smoother given in formula
  stopifnot("Cannot specify more than one node smoother in formula." = 
    length(node_terms) <= 1)

  # determine if estimating node smoother (k does not already have a value)
  defused_call <- rlang::parse_expr(node_terms)
  node_args    <- rlang::call_args(defused_call)
  node_k       <- eval(node_args$k)

  if (length(node_terms) == 1 && (is.null(node_k) || length(node_k) > 1)) {
    # prepare estimate_smooth_basis arguments
    est_kwargs <- `_get_function_kwargs`(estimate_smooth_basis.default, vargs)
    est_kwargs <- c(list(
      target       = formula_target, 
      smooth_terms = node_terms, 
      df           = df, 
      regressors   = regressors, 
      family       = linkfamily, 
      method       = method, 
      discrete     = discrete
    ), est_kwargs)
      
    # estimate node term smoothing basis
    node_terms <- do.call(estimate_smooth_basis.default, est_kwargs)$est_terms

    # reformulate the formula with estimated node smoother
    formula_terms[node_index] <- node_terms
    formula <- stats::reformulate(
      termlabels = formula_terms, 
      response = formula_target, 
      env = .GlobalEnv
    )
  }
  
  # start model fitting process (with or without autocorrelations)
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
      ...      = mgcv_kwargs
    )

    # determine autocorrelation parameter, rho
    rho <- itsadug::start_value_rho(gam_fit_wo_rho)

    # fit GAM model with autocorrelation 
    gam_fit <- mgcv::bam(
      formula  = formula, 
      family   = linkfamily, 
      data     = df, 
      method   = method, 
      rho      = rho, 
      AR.start = ar_start, 
      discrete = discrete, 
      ...      = mgcv_kwargs
    )
  } else {
    # fit GAM model without autocorrelations
    gam_fit <- mgcv::bam(
      formula  = formula,
      family   = linkfamily, 
      data     = df, 
      method   = method, 
      discrete = discrete, 
      ...      = mgcv_kwargs
    )
  }

  return(gam_fit)
}


#' Save GAM model outputs.
#' 
#' @description
#' This function saves the GAM model as RData (.rda) file, the model summary as
#' a text (.txt) file, and the [mgcv::gam.check()] figures as a PNG (.png) file.
#' 
#' @param gam_model asdfasdf
#' @param output_file asdfasdf
#' @param model_summary asdfasdf
#' @param model_check  asdfasdf
#' @return None. 
#'  
#' @export
save_gam <- function(
  gam_model, 
  output_file, 
  model_summary = TRUE,
  model_check   = FALSE 
) {
  # argument input control
  stopifnot("`gam_model` must be a class gam or bam" = 
    any(class(gam_model) %in% c("bam", "gam")))
  stopifnot("`output_file` must be a character" = is.character(output_file))
  if (fs::path_ext(output_file) == "") { # if missing file extension
    output_file <- fs::path_ext_set(output_file, "rda") # append
  }
  stopifnot("`output_file` file extension should be Rdata, Rdata, rdata, or rda" = 
    fs::path_ext(output_file) %in% c("RData", "Rdata", "rdata", "rda"))
  stopifnot("`model_summary` must be a logical" = is.logical(model_summary))
  stopifnot("`model_check` must be a logical" = is.logical(model_check))

  # save model as Rdata file
  saveRDS(gam_model, file = output_file)

  # save model summary as text file
  if (model_summary) {
    save_name <- stringr::str_replace(
      output_file, ".(RData|Rdata|rdata|rda)$", "_Summary.txt")
    utils::capture.output(summary(gam_model), file = save_name)
  }

  # save gam.check output (figure and text)
  if (model_check) {
    # save gam.check figures (via gratia package)
    save_name <- stringr::str_replace(
      output_file, ".(RData|Rdata|rdata|rda)$", "_GamCheck.png")
    gam_check <- gratia::appraise(gam_model)
    ggplot2::ggsave(save_name, plot = gam_check, scale = 2.5)

    # save k.check component
    save_name <- stringr::str_replace(
      output_file, ".(RData|Rdata|rdata|rda)$", "_GamCheck.txt")
    utils::capture.output(
      mgcv::k.check(gam_model, n.rep = 500),
      file = save_name
    )
  }
}


# Helper Functions that are NOT exported
`_get_formula_target` <- function(formula) {
  return(all.vars(formula)[attr(stats::terms(formula), "response")])
}

`_get_formula_terms` <- function(formula) {
  return(labels(terms(formula)))
}
  
`_get_family_function` <- function(distr_name) {
  rlang::arg_match(distr_name, values = FAMILY_FUNCTION_NAMES)
  family_func <- switch(distr_name, 
    "beta"     = mgcv::betar(link = "logit"), 
    "gamma"    = stats::Gamma(link = "logit"), 
    "gaussian" = stats::gaussian(link = "identity"),
    "norm"     = stats::gaussian(link = "identity")
  )
  return(family_func)
}