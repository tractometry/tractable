test_that("build_formula runs as expected", {
  expect_identical(
    build_formula(target = "fa"),
    as.formula("fa ~ s(nodeID) + s(subjectID, bs = 're')", env = .GlobalEnv)
  )

  expect_identical(
    build_formula(target = "fa", node_k = 15, node_col = "x", 
                  participant_col = "participant_id"),
    as.formula("fa ~ s(x, k = 15) + s(participant_id, bs = 're')", 
               env = .GlobalEnv)
  )

  expect_identical(
    build_formula(target = "fa", regressors = "group", 
                  node_k = 15, node_col = "x", 
                  participant_col = "participant_id"),
    as.formula("fa ~ group + s(x, k = 15) + s(participant_id, bs = 're')", 
               env = .GlobalEnv)
  )

  expect_identical(
    build_formula(target = "fa", regressors = "group", 
                  node_k = 15, node_col = "x", node_group = "group",
                  participant_col = "participant_id"),
    as.formula("fa ~ group + s(x, by = group, bs = 'fs', k = 15) + s(participant_id, bs = 're')", 
               env = .GlobalEnv)
  )

  expect_identical(
    build_formula(target = "fa", 
                  regressors = c("group", "sex"), 
                  node_k = 15, node_col = "x", 
                  participant_col = "participant_id"),
   as.formula("fa ~ group + sex + s(x, k = 15) + s(participant_id, bs = 're')", 
              env = .GlobalEnv)
  )

  expect_identical(
    build_formula(target = "fa", 
                  regressors = c("group", "sex", "nodeID"), 
                  node_k = 15, node_group = "group"), 
    as.formula("fa ~ group + sex + s(nodeID, by = group, bs = 'fs', k = 15) + s(subjectID, bs = 're')", 
               env = .GlobalEnv)
  )

  expect_identical(
    build_formula(
      target = "fa", regressors = c("age", "group", "s(age, by = group, k = 2)")), 
    as.formula("fa ~ age + group + s(age, by = group, k = 2) + s(nodeID) + s(subjectID, bs = 're')", 
               env = .GlobalEnv)
  )
})


test_that("estimate_distribution runs as expected", {
  x <- stats::rnorm(1000)
  expect_identical(estimate_distribution(x), 
                   stats::gaussian(link = "identity"))
  
  expect_error(estimate_distribution("gaussian"))
  expect_error(estimate_distribution(
    x, distr_options = c("logit", "gamma")))
  expect_error(estimate_distribution(x, eval_metric = "ML"))

  x <- stats::rgamma(1000, shape = 12)
  expect_identical(estimate_distribution(x), 
                   stats::Gamma(link = "logit"))
})


test_that("estimate_smooth_basis.default runs as expected", {
  df_sarica <- tractable::read_afq_sarica(na_omit = TRUE) %>% 
    dplyr::filter(tractID == "Right Corticospinal")

  estimated_information <- estimate_smooth_basis(
    target       = "md", 
    smooth_terms = "s(nodeID, bs = 'ts', k = c(12, 20))", 
    df           = df_sarica, 
  )
  expected_term <- "s(nodeID, bs = 'ts', k = 20)"
  expect_identical(estimated_information$est_terms, expected_term)
  expected_values <- estimated_information$k_estimates %>% 
    dplyr::filter(term == expected_term)
  expect_true(expected_values$k_index > 0.95 && expected_values$p_value > 0.05)
  
  estimated_information <- estimate_smooth_basis(
    target       = "fa", 
    smooth_terms = "s(nodeID, by = group, bs = 'fs')", 
    df           = df_sarica, 
    regressors   = c("age", "group"), 
  )
  expected_term <- "s(nodeID, by = group, bs = 'fs', k = 9)"
  expect_identical(estimated_information$est_terms, expected_term)
  expected_values <- estimated_information$k_estimates %>% 
    dplyr::filter(term == expected_term)
  expect_true(expected_values$k_index > 0.95 && expected_values$p_value > 0.05)
  
  estimated_information <- estimate_smooth_basis(
    target       = "fa", 
    smooth_terms = "s(nodeID, by = group, bs = 'fs', k = c(4, 11))", 
    df           = df_sarica, 
    regressors   = c("age", "group"), 
    kindex_thr   = 0.98,
    pvalue_thr   = 0.1
  )
  expected_term <- "s(nodeID, by = group, bs = 'fs', k = 11)"
  expect_identical(estimated_information$est_terms, expected_term)
  expected_values <- estimated_information$k_estimates %>% 
    dplyr::filter(term == expected_term)
  expect_true(expected_values$k_index > 0.98 && expected_values$p_value > 0.1)

  estimated_information <- estimate_smooth_basis(
    target       = "fa", 
    smooth_terms = c("s(age, k = 1:10)", "s(nodeID, by = group, bs = 'fs', k = c(4,9), m = 3)"),
    df           = df_sarica, 
    regressors   = "group"
  )
  expected_term <- c("s(age, bs = 'tp', k = 10)", "s(nodeID, by = group, bs = 'fs', m = 3, k = 9)")
  expect_identical(estimated_information$est_terms, expected_term)
  expected_values <- estimated_information$k_estimates %>% 
    dplyr::filter(term == expected_term[1]) # age term
  expect_false(expected_values$k_index > 0.95 && expected_values$p_value > 0.05)
  expected_values <- estimated_information$k_estimates %>% 
    dplyr::filter(term == expected_term[2]) # node term
  expect_true(expected_values$k_index > 0.95 && expected_values$p_value > 0.05)
})


test_that("estimate_smooth_basis.formula runs as expected", {
  df_sarica <- tractable::read_afq_sarica(na_omit = TRUE) %>% 
    dplyr::filter(tractID == "Right Corticospinal")

  estimated_information <- estimate_smooth_basis(
    formula = md ~ s(nodeID, bs = "ts", k = k = c(12, 20)), 
    df      = df_sarica
  )
  expected_term <- "s(nodeID, bs = 'ts', k = 20)"
  expect_identical(estimated_information$est_terms, expected_term)
  expected_values <- estimated_information$k_estimates %>% 
    dplyr::filter(term == expected_term)
  expect_true(expected_values$k_index > 0.95 && expected_values$p_value > 0.05)
  
  estimated_information <- estimate_smooth_basis(
    formula = fa ~ age + group + s(nodeID, by = group, bs = "fs"), 
    df      = df_sarica, 
  )
  expected_term <- "s(nodeID, by = group, bs = 'fs', k = 9)"
  expect_identical(estimated_information$est_terms, expected_term)
  expected_values <- estimated_information$k_estimates %>% 
    dplyr::filter(term == expected_term)
  expect_true(expected_values$k_index > 0.95 && expected_values$p_value > 0.05)
  
  estimated_information <- estimate_smooth_basis(
    formula    = fa ~ age + group + s(nodeID, by = group, bs = "fs", k = c(4, 11)),
    df         = df_sarica, 
    kindex_thr = 0.98,
    pvalue_thr = 0.1
  )
  expected_term <- "s(nodeID, by = group, bs = 'fs', k = 11)"
  expect_identical(estimated_information$est_terms, expected_term)
  expected_values <- estimated_information$k_estimates %>% 
    dplyr::filter(term == expected_term)
  expect_true(expected_values$k_index > 0.98 && expected_values$p_value > 0.1)

  estimated_information <- estimate_smooth_basis(
    formula = fa ~ group + s(age, k = 1:10) 
      + s(nodeID, by = group, bs = "fs", k = c(4,9), m = 3), 
    df      = df_sarica, 
  )
  expected_term <- c("s(age, bs = 'tp', k = 10)", "s(nodeID, by = group, bs = 'fs', m = 3, k = 9)")
  expect_identical(estimated_information$est_terms, expected_term)
  expected_values <- estimated_information$k_estimates %>% 
    dplyr::filter(term == expected_term[1]) # age term
  expect_false(expected_values$k_index > 0.95 && expected_values$p_value > 0.05)
  expected_values <- estimated_information$k_estimates %>% 
    dplyr::filter(term == expected_term[2]) # node term
  expect_true(expected_values$k_index > 0.95 && expected_values$p_value > 0.05)
})


test_that("fit_gam.default runs as expected", {
  df_sarica <- tractable::read_afq_sarica(na_omit = TRUE) %>% 
    dplyr::filter(tractID == "Right Corticospinal")

  expect_no_error({
    model_fit <- fit_gam(
      target     = "fa", 
      df         = df_sarica, 
      regressors = c("age", "group"), 
      node_group = "group", 
    )
    model_summary <- summary(model_fit)
  })
  expect_identical(model_summary$formula, as.formula(
    "fa ~ age + group + s(nodeID, by = group, bs = 'fs', k = 9) + s(subjectID, bs = 're')",
    env = .GlobalEnv))

  expect_no_error({
    model_fit <- fit_gam(
      target     = "fa", 
      df         = df_sarica, 
      regressors = c("age", "group"), 
      k_values   = c(2, 9, 13)
    )
    model_summary <- summary(model_fit)
  })
  expect_identical(model_summary$formula, as.formula(
    "fa ~ age + group + s(nodeID, bs = 'tp', k = 13) + s(subjectID, bs = 're')",
    env = .GlobalEnv))

  expect_no_error({
    model_fit <- fit_gam(
      target     = "fa", 
      df         = df_sarica, 
      regressors = c("age", "group"), 
      autocor    = FALSE, 
    )
    model_summary <- summary(model_fit)
  })
  expect_identical(model_fit$AR1.rho, 0)
  
  expect_no_error({
    model_fit <- fit_gam(
      target     = "fa", 
      df         = df_sarica, 
      regressors = c("age", "group"),
      node_k     = 10, 
    )
    model_summary <- summary(model_fit)
  })
  expect_identical(model_summary$formula, as.formula(
    "fa ~ age + group + s(nodeID, k = 10) + s(subjectID, bs = 're')",
    env = .GlobalEnv))

})


test_that("fit_gam.formula runs as expected", {
  df_sarica <- tractable::read_afq_sarica(na_omit = TRUE) %>% 
    dplyr::filter(tractID == "Right Corticospinal")

  expect_no_error({
    model_fit <- fit_gam(
      formula = fa ~ age + group 
        + s(nodeID, by = group, bs = "fs") 
        + s(subjectID, bs = "re"),
      df = df_sarica, 
    )
    model_summary <- summary(model_fit)
  })
  expect_identical(model_summary$formula, as.formula(
    "fa ~ age + group + s(nodeID, by = group, bs = 'fs', k = 9) + s(subjectID, bs = 're')",
    env = .GlobalEnv))

  expect_no_error({
    model_fit <- fit_gam(
      formula = fa ~ age + group 
        + s(nodeID) 
        + s(subjectID, bs = "re"),
      df = df_sarica, 
      k_values = c(1, 9, 13)
    )
    model_summary <- summary(model_fit)
  })
  expect_identical(model_summary$formula, as.formula(
    "fa ~ age + group + s(nodeID, bs = 'tp', k = 13) + s(subjectID, bs = 're')",
    env = .GlobalEnv))

  expect_no_error({
    model_fit <- fit_gam(
      formula = fa ~ age + group +
        + s(nodeID, bs = "ts") 
        + s(subjectID, bs = "re"),
      df = df_sarica, 
      autocor = FALSE
    )
    model_summary <- summary(model_fit)
  })
  expect_identical(model_fit$AR1.rho, 0)
  
  expect_no_error({
    model_fit <- fit_gam(
      formula = fa ~ age + group
        + s(nodeID, k = 10) 
        + s(subjectID, bs = "re"),
      df = df_sarica
    )
    model_summary <- summary(model_fit)
  })
  expect_identical(model_summary$formula, as.formula(
    "fa ~ age + group + s(nodeID, k = 10) + s(subjectID, bs = 're')",
    env = .GlobalEnv))

})


test_that("save_gam runs as expected", {
  df_sarica <- tractable::read_afq_sarica(na_omit = TRUE) %>% 
    dplyr::filter(tractID == "Right Corticospinal")

  model_fit <- fit_gam(
    target     = "fa", 
    df         = df_sarica, 
    regressors = c("age", "group"), 
    node_group = "group", 
  )

  expect_error(save_gam(model_fit, output_file = "model.txt"))
  expect_no_error(save_gam(model_fit, output_file = tempfile()))
  expect_no_error(save_gam(model_fit, output_file = tempfile(), model_summary = FALSE))
  expect_no_error(save_gam(model_fit, output_file = tempfile(), model_check = TRUE))
})