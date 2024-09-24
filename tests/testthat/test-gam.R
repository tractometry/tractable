test_that("tractable_single_bundle runs as expected", {
 sarica <- read_afq_sarica()
 sarica$group <- factor(sarica$class)
 sarica$subjectID <- unclass(factor(sarica$subjectID))
 gam_fit <- expect_no_error(
    tractable_single_bundle(df_afq = sarica,
                         tract = "Right Corticospinal",
                         participant_id = "subjectID",
                         group_by = "group",
                         covariates = c("age", "group"),
                         dwi_metric = "fa"))

expect_length(gam_fit["terms"], 1)

                      })