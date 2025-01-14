test_that("tractable_single_bundle run as expected", {
  df_sarica <- tractable::read_afq_sarica(na_omit = TRUE)

  expect_no_error(
    model_fit <- tractable_single_tract(
      df = df_sarica,
      tract = "Right Corticospinal",
      target = "fa",
      regressors = c("age", "group"),
      node_group = "group"
    )
  )

  formula <- "fa ~ age + group + s(nodeID, by = group, bs = 'fs', k = 9) + s(subjectID, bs = 're')"
  expect_identical(model_fit$formula, as.formula(formula, env = .GlobalEnv))

  expect_no_error(
    model_fit <- tractable_single_tract(
      formula = fa ~ age + group + s(nodeID, by = group, bs = "fs") + s(subjectID, bs = "re"),
      df = df_sarica,
      tract = "Right Corticospinal"
    )
  )

  formula <- "fa ~ age + group + s(nodeID, by = group, bs = 'fs', k = 9) + s(subjectID, bs = 're')"
  expect_identical(model_fit$formula, as.formula(formula, env = .GlobalEnv))
})
