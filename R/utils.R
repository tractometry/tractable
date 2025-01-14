
`_get_function_kwargs` <- function(func, vargs) {
  func_kwargs <- rlang::fn_fmls_names(func)
  func_kwargs <- sapply(func_kwargs, function(x) vargs[[x]])
  func_kwargs <- func_kwargs[!sapply(func_kwargs, is.null)]
  if (length(func_kwargs) == 0) { return(NULL) }
  return(func_kwargs)
}
