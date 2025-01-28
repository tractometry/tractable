
`_get_function_kwargs` <- function(func, vargs) {
  func_kwargs <- rlang::fn_fmls_names(func)
  func_kwargs <- sapply(func_kwargs, function(x) vargs[[x]])
  func_kwargs <- func_kwargs[!sapply(func_kwargs, is.null)]
  if (length(func_kwargs) == 0) { return(NULL) }
  return(func_kwargs)
}

.collap_slice_resample <- function(df){
  if(inherits(df, c("GRP_df", "grouped_df"))) {
    indices <- unlist(lapply(collapse::GRPN(df, expand=FALSE), sample.int, replace = TRUE), recursive = FALSE, use.names = FALSE)
  } else {
    indices <- sample.int(nrow(df), replace = TRUE)
  }
  collapse::ss(df, i = indices)
}


.cast_as_thing.factor <- function(from, to=NULL, toClass=NULL, conv_fun_namer = function(cls) Filter(\(x)nchar(x)>3,paste0(c("as.", "as_"),Filter(\(x)nchar(x)>0,c(cls[[1]],paste0(toupper(substr(cls[[1]],1,1)),substring(cls[[1]],2)))))), ..., .use_collapse=FALSE) {
  stopifnot(inherits(from, "factor"))
  if (is.null(from)) {
    return(from)
  }
  if (rlang::obj_address(from) == rlang::obj_address(to)) {
    return(from)
  }
  if (!is.null(to)) {
    if (all(class(from) == class(to))) {
      return(from)
    }
  }
  res <- NULL
  cls <- if (is.null(to)) toClass else class(to)
  if (inherits(from, cls[[1]])) {
    return(from)
  }
  if (!is.null(res <- tryCatch(methods::as(from, cls[[1]]), error = \(e)NULL))) {
    return(res)
  }
  pkg_ns <- tryCatch(asNamespace(attr(methods::getClass(cls), "package", exact = TRUE)), error = \(...)NULL)
  conv_fn_names <- conv_fun_namer(cls[[1]])
  for (conv_fn_name in conv_fn_names) {
    conv_fn <- get0(conv_fn_name)
    if (is.null(conv_fn) && !is.null(pkg_ns)) conv_fn <- get0(conv_fn_name, envir = pkg_ns)
    if (!is.null(conv_fn)) break
  }

  if (is.null(conv_fn)) {
    for (conv_fn_name in conv_fn_names) {
      found <- utils::getAnywhere(conv_fn_name)
      if (length(found$objs) > 0) {
        conv_fn <- conv_fns$objs[[1]]
        break
      }
    }
  }
  if (!is.null(conv_fn)) res <- conv_fn(from)

  if (!is.null(to)) {
    if (!.use_collapse) {
      unique_vals <- unique(to)
      names(unique_vals) <- as.character(unique_vals)
      stopifnot(!anyNA(names(unique_vals)))
      res <- unique_vals[match(as.character(from), table = names(unique_vals))]
    } else {
      unique_vals <- collapse::funique(to)
      names(unique_vals) <- as.character(unique_vals)
      stopifnot(length(collapse::whichNA(names(unique_vals))) == 0)
      res <- unique_vals[collapse::fmatch(collapse::as_character_factor(from), table = names(unique_vals))]
    }
  }
  stopifnot(!is.null(res))
  dn <- attr(res, "dimnames", exact = TRUE)
  d <- attr(res, "dim", exact = TRUE)
  attributes(res) <- attributes(to)
  attr(res, "dimnames") <- dn
  attr(res, "dimnames") <- d
  res
}

