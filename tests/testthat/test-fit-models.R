test_that("build_formula runs as expected", {

  expect_identical(
    tractable::build_formula(target = "fa"),
    "fa ~ s(nodeID, k = 10) + s(subjectID, bs = 're')"
  )

  expect_identical(
    tractable::build_formula(target = "fa", node_k = 15, node_col = "x", 
                             participant_col = "participant_id"),
    "fa ~ s(x, k = 15) + s(participant_id, bs = 're')"
  )

  expect_identical(
    tractable::build_formula(target = "fa", regressors = "group", 
                             node_k = 15, node_col = "x", 
                             participant_col = "participant_id"),
    "fa ~ group + s(x, k = 15) + s(participant_id, bs = 're')"
  )

  expect_identical(
    tractable::build_formula(target = "fa", regressors = "group", 
                             node_k = 15, node_col = "x", node_group = "group",
                             participant_col = "participant_id"),
    "fa ~ group + s(x, by = group, k = 15) + s(participant_id, bs = 're')"
  )

  expect_identical(
    tractable::build_formula(target = "fa", 
                             regressors = c("group", "sex"), 
                             node_k = 15, node_col = "x", 
                             participant_col = "participant_id"),
    "fa ~ group + sex + s(x, k = 15) + s(participant_id, bs = 're')"
  )

  expect_identical(
    tractable::build_formula(target = "fa", 
                             regressors = c("group", "sex", "nodeID"), 
                             node_k = 15, node_group = "group"), 
    "fa ~ group + sex + s(nodeID, by = group, k = 15) + s(subjectID, bs = 're')"
  )

  expect_identical(
    tractable::build_formula(
      target = "fa", regressors = c("age", "group", "s(age, by = group, k = 2)")), 
    "fa ~ age + group + s(age, by = group, k = 2) + s(nodeID, k = 10) + s(subjectID, bs = 're')"
  )
})

# test_that("fit_gam runs as expected", {

#   sarica <- read_afq_sarica()
#   sarica$group <- factor(sarica$class)
#   sarica$subjectID <- unclass(factor(sarica$subjectID))

#   selected <- select_bundle(
#     df_afq = sarica,
#     tract = "Right Corticospinal",
#     dwi_metric = "fa",
#     covariates = c("age", "group"),
#     participant_id = "subjectID",
#     group_by = "group"
#   )

#   df_tract <- selected$df_tract
#   tract_names <- selected$tract_names

#   formula <- tractable::build_formula(
#     target = "fa",
#     regressors = c("age", "group"),
#     node_k = 40
#   )

#   # One and only one of "target" and "formula" should be set to non-NULL
#   gam_fit <- expect_error(tractable::fit_gam(df_tract = df_tract,
#                                              target = NULL,
#                                              formula = NULL))

#   gam_fit <- expect_error(tractable::fit_gam(df_tract = df_tract,
#                                              target = "fa",
#                                              formula = formula))

#   gam_fit <- expect_no_error(tractable::fit_gam(df_tract = df_tract,
#                                                 formula = formula))

#   gam_fit <- expect_no_error(tractable::fit_gam(df_tract = df_tract,
#                                                 autocor = FALSE,
#                                                 formula = formula))

#   # Check that formula can be passed as a string
#   string_formula = "fa ~ age + group + s(nodeID, k = 40) + s(subjectID, bs = 're')"
#   gam_fit <- expect_no_error(tractable::fit_gam(df_tract = df_tract,
#                                             formula = string_formula))

# })