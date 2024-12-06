#' Select a single tract from an AFQ dataframe
#'
#' @param df_afq The input AFQ dataframe
#' @param tract Abbreviated tract name, e.g., "CST_L" or "OR"
#' @param dwi_metric The diffusion MRI metric (e.g. "FA", "MD")
#' @param covariates List of strings of GAM covariates, not including
#'     the smoothing terms over nodes and the random effect due to subjectID.
#'     This list can also include smoothing terms.
#' @param participant_id The name of the column that encodes participant ID
#' @param group_by The grouping variable used to group nodeID smoothing terms
#'
#' @return An AFQ dataframe containing only the selected tract
#' @export
select_tract <- function(df_afq,
                          tract,
                          dwi_metric,
                          covariates,
                          participant_id = "subjectID",
                          group_by = "group") {

  if (!is.null(group_by)) {
    df_afq[[group_by]] <- factor(df_afq[[group_by]])
  }

  df_afq[[participant_id]] <- factor(df_afq[[participant_id]])

  cols <- unique(c(participant_id,
                   "nodeID",
                   "tractID",
                   dwi_metric,
                   group_by,
                   covariates))
  df_afq <- dplyr::select(df_afq, dplyr::any_of(cols))

  if (tract != "all") {
    df_tract <- df_afq[which(df_afq$tractID == tract), ]
    df_tract <- stats::na.omit(df_tract)
    tract_names <- c(tract)
  } else {
    df_tract <- df_afq
    old_group_by <- group_by
    group_by <- paste0("tractID", group_by)
    df_tract[[group_by]] <- as.factor(paste(df_tract$tractID,
                                            df_tract[[old_group_by]]))
    df_tract <- stats::na.omit(df_tract)
    tract_names <- unique(df_tract$tractID)
  }
  return(list(df_tract = df_tract,
              tract_names = tract_names))
}