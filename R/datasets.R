#' Parse an S3 URI
#'
#' @param uri An AWS S3 URI
#'
#' @return A named list with bucket and object elements
#' @export
#'
#' @examples
#' parse_s3_uri("s3://bucket-name/a/path/to/an/object.txt")
parse_s3_uri <- function(uri) {
  s3_uri_regex <- stringr::regex("
    ^s3://  # initial S3 prefix
    (?<bucket>(?!(xn--|-s3alias$))[a-z0-9][a-z0-9-]{1,61}[a-z0-9]) # bucket name
    /  # separating /
    (?<object>[a-zA-Z0-9!-_.*'()][a-zA-Z0-9!-_.*'()/]+)$  # object name
    ", comments = TRUE)

  m <- stringr::str_match(uri, s3_uri_regex)[1, ]
  bucket <- getElement(m, "bucket")
  object <- getElement(m, "object")

  return(c("bucket" = bucket, "object" = object))
}

#' Create a merged AFQ/phenotype dataframe
#'
#' @param nodes_csv path to a nodes file
#' @param pheno_csv path to a phenotypic file, leave NULL to return an "unsupervised" AFQ dataset
#' @param index specification of the column used for merging
#' @param index_nodes specification of the column used for merging
#' @param index_pheno specification of the column used for merging
#' @param dwi_metrics which diffusion metrics should be retained
#' @param factor_cols which columns should be treated as factors
#' @param pheno_cols which columns to include from pheno file
#' @param ... arguments to be passed to read.csv
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   df_afq <- read_afq_files(nodes_csv = "a/path/to/nodes.csv",
#'                            pheno_csv = "a/path/to/pheno.csv")
#' }
read_afq_files <- function(nodes_csv, pheno_csv = NULL, index = "subjectID", 
                           index_nodes = index, index_pheno = index, 
                           dwi_metrics = NULL, factor_cols = NULL, 
                           pheno_cols = NULL,  ...) {
  if (grepl("^s3://", nodes_csv)) {
    # Read from S3
    uri <- parse_s3_uri(nodes_csv)
    nodes_df <- aws.s3::s3read_using(FUN = utils::read.csv,
                                     bucket = uri[["bucket"]],
                                     object = uri[["object"]],
                                     ... = ...)
  } else {
    # Read from local or URL
    nodes_df <- utils::read.csv(nodes_csv, ... = ...)
  }

  # Select only the user's desired DWI metrics
  if (!is.null(dwi_metrics)) {
    node_cols <- unique(c(index_nodes,
                          "siteID",
                          "tractID",
                          "nodeID",
                          dwi_metrics))
    nodes_df <- dplyr::select(nodes_df, dplyr::any_of(node_cols))
  }

  # Add "sub-" prefix to participant IDs if not already there
  nodes_df[[index_nodes]] <- trimws(as.character(nodes_df[[index_nodes]]))
  nodes_df[[index_nodes]] <- ifelse(grepl("^sub-", nodes_df[[index_nodes]]),
                                    nodes_df[[index_nodes]],
                                    paste0("sub-", nodes_df[[index_nodes]]))

  if (!is.null(pheno_csv)) {
    # Treat some pheno columns as factors
    if (is.null(factor_cols)) {
      colClasses <- NA
    } else {
      colClasses <- rep("factor", length(factor_cols))
      names(colClasses) <- factor_cols
    }

    # Determine separator from file suffix
    if (grepl("\\.tsv$", pheno_csv)) {
      sep <- "\t"
    } else {
      sep <- ","
    }

    if (grepl("^s3://", pheno_csv)) {
      # Read from S3
      uri <- parse_s3_uri(pheno_csv)
      pheno_df <- aws.s3::s3read_using(FUN = utils::read.csv,
                                      bucket = uri[["bucket"]],
                                      object = uri[["object"]],
                                      sep = sep,
                                      ... = ...)
    } else {
      # Read from local or URL
      pheno_df <- utils::read.csv(pheno_csv,
                          check.names = FALSE,
                          colClasses = colClasses,
                          sep = sep,
                          ... = ...)
    }

    # Select only the user supplied pheno columns
    if (!is.null(pheno_cols)) {
      pheno_df <- dplyr::select(pheno_df, unique(c(index_pheno, pheno_cols)))
    }

    # Add "sub-" prefix to participant IDs if not already there
    pheno_df[[index_pheno]] <- trimws(as.character(pheno_df[[index_pheno]]))
    pheno_df[[index_pheno]] <- ifelse(grepl("^sub-", pheno_df[[index_pheno]]),
                                      pheno_df[[index_pheno]],
                                      paste0("sub-", pheno_df[[index_pheno]]))

    # Merge pheno is provided
    merged_df <- merge(nodes_df, pheno_df, by.x = index_nodes, 
                       by.y = index_pheno, incomparables = NA)
    return(merged_df)
  }

  # If we get here, then we are in the unsupervised case, return just the nodes_df.
  return(nodes_df)
}

#' Load tract profiles from Sarica et al.
#'
#' @param ... arguments to be passed to read.csv
#'
#' @return A merged dataframe with data from Sarica et al.
#' @export
#'
#' @examples
#' df_sarica <- read_afq_sarica()
read_afq_sarica  <- function(...) {
  nodes_url <- "https://github.com/yeatmanlab/Sarica_2017/raw/gh-pages/data/nodes.csv"
  pheno_url <- "https://github.com/yeatmanlab/Sarica_2017/raw/gh-pages/data/subjects.csv"
  df <- read_afq_files(nodes_csv   = nodes_url,
                       pheno_csv   = pheno_url,
                       dwi_metrics = c("fa", "md"),
                       factor_cols = c("class", "gender"),
                       pheno_cols  = c("age", "class", "gender"),
                       ... = ...)
  return(df)
}

#' Load tract profiles from Yeatman et al.
#'
#' @param ... arguments to be passed to read.csv
#'
#' @return A merged dataframe with data from Yeatman et al.
#' @export
#'
#' @examples
#' df_weston_havens <- read_afq_weston_havens()
read_afq_weston_havens <- function(...) {
  nodes_url <- "https://yeatmanlab.github.io/AFQBrowser-demo/data/nodes.csv"
  pheno_url <- "https://yeatmanlab.github.io/AFQBrowser-demo/data/subjects.csv"
  df <- read_afq_files(nodes_csv   = nodes_url,
                       pheno_csv   = pheno_url,
                       dwi_metrics = c("fa", "md"),
                       factor_cols = c("Gender"),
                       pheno_cols  = c("Age", "Gender", "IQ"),
                       ... = ...)
  return(df)
}

#' Load tract profiles from the Healthy Brain Network dataset
#'
#' @param truncate if TRUE, truncate the data to 49 rows. default = FALSE
#' @param ... arguments to be passed to read.csv
#'
#' @return A merged dataframe with data from HBN
#' @export
#'
#' @examples
#' \dontrun{
#'   df_hbn <- read_afq_hbn()
#' }
read_afq_hbn <- function(truncate = FALSE, ...) {
  if (truncate) {
    nodes_url <- "s3://fcp-indi/data/Projects/HBN/BIDS_curated/derivatives/afq/.truncated_tract_profiles.csv"
    pheno_url <- "s3://fcp-indi/data/Projects/HBN/BIDS_curated/derivatives/afq/.truncated_participants.tsv"
  } else { # nocov start
    nodes_url <- "s3://fcp-indi/data/Projects/HBN/BIDS_curated/derivatives/afq/combined_tract_profiles.csv"
    pheno_url <- "s3://fcp-indi/data/Projects/HBN/BIDS_curated/derivatives/qsiprep/participants.tsv"
  } # nocov end

  df <- read_afq_files(nodes_csv   = nodes_url,
                       pheno_csv   = pheno_url,
                       index_pheno = "subject_id",
                       dwi_metrics = c("dki_fa", "dki_md"),
                       factor_cols = c("sex"),
                       pheno_cols  = c("age", "sex"),
                       ... = ...)
  return(df)
}
