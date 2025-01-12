#' Create a merged AFQ and phenotype dataframe
#' 
#' @description 
#' There should a description of what an "AFQ" dataset is
#'
#' @param nodes_file    Path to a nodes file.
#' @param pheno_file    Path to a phenotypic file. Leave NULL to return an 
#'                      "unsupervised" AFQ dataset
#' @param by            Column used for merging the nodes and phenotypic files. 
#'                      If variable names differ between the nodes and 
#'                      phenotypic datasets, use a named character vector like 
#'                      by = c("x_a" = "y_a", "x_b" = "y_b").
#'                      Default: "subjectID"
#' @param tract_col   description
#' @param node_col    description
#' @param other_cols  description
#' @param na_omit     description
#' @param ...
#'
#' @examples
#' \dontrun{
#' df_afq <- read_afq_files(
#'   nodes_file = "path/to/nodes.csv",
#'   pheno_file = "path/to/pheno.csv"
#' )
#' }
#' @export
read_afq_files <- function(
  nodes_file, 
  pheno_file  = NULL, 
  by          = "subjectID", 
  tract_col   = "tractID",
  node_col    = "nodeID", 
  other_cols  = NULL, 
  na_omit     = FALSE,
  ...
) {
  # argument input control
  stopifnot("`nodes_file` must be a character" = is.character(nodes_file))
  stopifnot("`nodes_file` must be a recognized extension (csv, tsv)" = 
    fs::path_ext(nodes_file) %in% c("csv", "tsv"))
  if (!is.null(pheno_file)) {
    stopifnot("`pheno_file` must be a character" = is.character(pheno_file))
    stopifnot("`pheno_file` must be a recognized extension (csv, tsv)" = 
      fs::path_ext(pheno_file) %in% c("csv", "tsv"))
  }
  stopifnot("`by` must be a character" = is.character(by))
  stopifnot("`node_col` must be a character" = is.character(node_col))
  if (!is.null(other_cols)) {
    stopifnot("`other_cols` must be a character" = is.character(other_cols))
  }
  stopifnot("`na_omit` must be a logical" = is.logical(na_omit))

  # read nodes file
  df_nodes <- `_read_afq_file`(nodes_file, ... = ...)

  if (!is.null(pheno_file)) {
    # read the phenotypic file
    df_pheno <- `_read_afq_file`(pheno_file, ... = ...)

    # merge the nodes and phenotypic datasets
    df_merged <- dplyr::full_join(x = df_nodes, y = df_pheno, by = by)
  } else {
    # place index first in "merged" dataframe (nodes dataset)
    df_merged <- df_nodes %>% dplyr::relocate(all_of(by), 1)
  }

  # subset columns by requested values
  if (!is.null(names(by))) { by <- names(by) }
  keep_cols <- c(by, tract_col, node_col, other_cols)
  df <- df_merged %>% dplyr::select(tidyselect::all_of(keep_cols))

  # remove any row with NAs 
  if (na_omit) {
    df <- df %>% tidyr::drop_na()
  }

  return(df)
}


#' Load tract profiles from Sarica et al.
#'
#' @param ... arguments to be passed to read.csv
#'
#' @return A merged dataframe with data from Sarica et al.
#'
#' @examples
#' df_sarica <- read_afq_sarica()
#' @export
read_afq_sarica  <- function(na_omit = TRUE, ...) {
  # define node and pheno S3 URLs
  nodes_file <- "https://github.com/yeatmanlab/Sarica_2017/raw/gh-pages/data/nodes.csv"
  pheno_file <- "https://github.com/yeatmanlab/Sarica_2017/raw/gh-pages/data/subjects.csv"
    
  # read and prepare sarica dataset
  suppressMessages(
    df <- read_afq_files(
      nodes_file   = nodes_file,
      pheno_file   = pheno_file,
      other_cols   = c("fa", "md", "age", "class", "gender"), 
      na_omit      = na_omit, 
      ...          = ...
    ) %>% 
      dplyr::rename(group = class) %>% 
      dplyr::mutate(
        dplyr::across(tidyselect::all_of(c("subjectID", "group", "gender")), 
        ~ factor(.x))
      )
  )

  return(df)
}


#' Load tract profiles from Yeatman et al.
#'
#' @param ... arguments to be passed to read.csv
#'
#' @return A merged dataframe with data from Yeatman et al.
#' @examples
#' df_weston_havens <- read_afq_weston_havens()
#' @export
read_afq_weston_havens <- function(na_omit = TRUE, ...) {
  # define node and pheno S3 URLs
  nodes_file <- "https://yeatmanlab.github.io/AFQBrowser-demo/data/nodes.csv"
  pheno_file <- "https://yeatmanlab.github.io/AFQBrowser-demo/data/subjects.csv"
  
  # read and prepare yeatman et al. dataset
  suppressMessages(
    df <- read_afq_files(
      nodes_file = nodes_file,
      pheno_file = pheno_file,
      other_cols = c("fa", "md", "Age", "Gender", "IQ"), 
      na_omit    = na_omit,
      ...        = ...
    ) %>% 
      dplyr::mutate(
        dplyr::across(tidyselect::all_of(c("subjectID", "Gender")), 
        ~ factor(.x))
      )
  )

  return(df)
}


#' Load tract profiles from the Healthy Brain Network dataset
#'
#' @param truncate   Truncate the data to 49 rows, boolean.
#'                   Default: FALSE.
#' @param ...        Arguments to be passed to read.csv
#'
#' @return A merged dataframe with data from HBN
#'
#' @examples
#' \dontrun{
#'   df_hbn <- read_afq_hbn()
#' }
#' @export
read_afq_hbn <- function(truncate = FALSE, na_omit = TRUE, ...) {
  # define node and pheno S3 URLs
  if (truncate) { # if truncating HBN download
    nodes_file <- "s3://fcp-indi/data/Projects/HBN/BIDS_curated/derivatives/afq/.truncated_tract_profiles.csv"
    pheno_file <- "s3://fcp-indi/data/Projects/HBN/BIDS_curated/derivatives/afq/.truncated_participants.tsv"
  } else { 
    nodes_file <- "s3://fcp-indi/data/Projects/HBN/BIDS_curated/derivatives/afq/combined_tract_profiles.csv"
    pheno_file <- "s3://fcp-indi/data/Projects/HBN/BIDS_curated/derivatives/qsiprep/participants.tsv"
  } 

  # read and prepare HBN dataset
  suppressMessages({
    # read nodes file to edit the subject id column
    df_nodes <- `_read_afq_file`(nodes_file) %>% 
      dplyr::mutate(subjectID = stringr::str_c("sub-", subjectID))

    # save the edit nodes file to a temporary file
    temp_file <- stringr::str_c(tempfile(), ".csv")
    readr::write_csv(df_nodes, file = temp_file)

    # read afq files with the temporary and phenotypic files
    df <- read_afq_files(
      nodes_file   = temp_file,
      pheno_file   = pheno_file,
      by           = c("subjectID" = "subject_id"), 
      other_cols   = c("dki_fa", "dki_md", "age", "sex"),
      na_omit      = na_omit, 
      ...          = ...
    ) %>% 
      dplyr::mutate(
        dplyr::across(tidyselect::all_of(c("subjectID", "sex")), 
        ~ factor(.x))
      )
  })

  return(df)
}


# Helper Functions that are NOT exported
`_parse_s3_uri` <- function(s3_uri) {
  s3_uri_regex <- stringr::regex("
    ^s3://  # initial S3 prefix
    (?<bucket>(?!(xn--|-s3alias$))[a-z0-9][a-z0-9-]{1,61}[a-z0-9]) # bucket name
    /  # separating /
    (?<object>[a-zA-Z0-9!-_.*'()][a-zA-Z0-9!-_.*'()/]+)$  # object name
  ", 
    comments = TRUE
  )

  m <- stringr::str_match(s3_uri, s3_uri_regex)[1, ]
  bucket <- getElement(m, "bucket")
  object <- getElement(m, "object")

  s3_info <- list("bucket" = bucket, "object" = object)
  return(s3_info)
}


`_read_s3_file` <- function(csv, func, ...) {
  parsed_uri <- `_parse_s3_uri`(csv)
  df <- aws.s3::s3read_using(
      FUN = func,  
      bucket = parsed_uri[["bucket"]],
      object = parsed_uri[["object"]],
      show_col_types = FALSE, 
      ... = ...
  )
  return(df) 
}


`_get_file_function` <- function(ext) {
  rlang::arg_match(ext, values = c("csv", "tsv"))
  file_func <- switch(ext, 
    "csv" = readr::read_csv, 
    "tsv" = readr::read_tsv
  )
  return(file_func)
}


`_read_afq_file` <- function(file, ...) {
  # determine file extension and reader function
  file_ext  <- fs::path_ext(file)
  file_func <- `_get_file_function`(file_ext)

  # load file from s3, local, or url
  if (stringr::str_detect(file, "^s3://")) {
    df <- `_read_s3_file`(file, file_func, ... = ...)
  } else { # read from local or URL
    df <- file_func(file, show_col_types = FALSE, ... = ...)
  }

  return(df)
}
