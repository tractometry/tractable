test_that("select_tract runs as expected", {
 sarica <- read_afq_sarica()
 sarica$group <- factor(sarica$class)
 sarica$subjectID <- unclass(factor(sarica$subjectID))

 selected <- expect_no_error(
    select_tract(df_afq = sarica,
                  tract = "Right Corticospinal",
                  dwi_metric = "fa",
                  covariates = c("age", "group"),
                  participant_id = "subjectID",
                  group_by = "group")
                )
df_tract <- selected$df_tract
tract_names <- selected$tract_names

expect_identical(unique(df_tract$tractID), "Right Corticospinal")
expect_identical(tract_names, c("Right Corticospinal"))
expect_identical(length(colnames(df_tract)), as.integer(6))

 selected <- expect_no_error(
    select_tract(df_afq = sarica,
                  tract = "all",
                  dwi_metric = "fa",
                  covariates = c("age", "group"),
                  participant_id = "subjectID",
                  group_by = "group")
                )
df_tract <- selected$df_tract
tract_names <- selected$tract_names

expect_identical(length(tract_names), as.integer(20))
    })
