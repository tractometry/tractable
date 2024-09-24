test_that("fit_gam runs as expected", {

  sarica <- read_afq_sarica()
  sarica$group <- factor(sarica$class)
  sarica$subjectID <- unclass(factor(sarica$subjectID))

  selected <- select_bundle(
    df_afq = sarica,
    tract = "Right Corticospinal",
    dwi_metric = "fa",
    covariates = c("age", "group"),
    participant_id = "subjectID",
    group_by = "group"
  )

  df_tract <- selected$df_tract
  tract_names <- selected$tract_names

  formula <- tractable::build_formula(
    "fa",
    covariates = c("age", "group"),
    group_by = "group",
    participant_id = "subjectID",
    k = 40
  )

  # One and only one of "target" and "formula" should be set to non-NULL
  gam_fit <- expect_error(tractable::fit_gam(df_tract = df_tract,
                                             target = NULL,
                                             formula = NULL))

  gam_fit <- expect_error(tractable::fit_gam(df_tract = df_tract,
                                             target = "fa",
                                             formula = formula))

  gam_fit <- expect_no_error(tractable::fit_gam(df_tract = df_tract,
                                                formula = formula))

  # Check that formula can be passed as a string
  string_formula = 'fa ~ age + group + s(nodeID, by = group, k = 40) + 
    s(subjectID, bs = "re")'
  gam_fit <- expect_no_error(tractable::fit_gam(df_tract = df_tract,
                                            formula = string_formula))

})