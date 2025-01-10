# test_that("build_formula runs as expected", {

#   expect_identical(
#     tractable::build_formula(target = "fa"),
#     as.formula("fa ~ s(nodeID) + s(subjectID, bs = 're')", env = .GlobalEnv)
#   )

#   expect_identical(
#     tractable::build_formula(target = "fa", node_k = 15, node_col = "x", 
#                              participant_col = "participant_id"),
#     as.formula("fa ~ s(x, k = 15) + s(participant_id, bs = 're')", 
#                env = .GlobalEnv)
#   )

#   expect_identical(
#     tractable::build_formula(target = "fa", regressors = "group", 
#                              node_k = 15, node_col = "x", 
#                              participant_col = "participant_id"),
#     as.formula("fa ~ group + s(x, k = 15) + s(participant_id, bs = 're')", 
#                env = .GlobalEnv)
#   )

#   expect_identical(
#     tractable::build_formula(target = "fa", regressors = "group", 
#                              node_k = 15, node_col = "x", node_group = "group",
#                              participant_col = "participant_id"),
#     as.formula("fa ~ group + s(x, by = group, bs = 'fs', k = 15) + s(participant_id, bs = 're')", 
#                env = .GlobalEnv)
#   )

#   expect_identical(
#     tractable::build_formula(target = "fa", 
#                              regressors = c("group", "sex"), 
#                              node_k = 15, node_col = "x", 
#                              participant_col = "participant_id"),
#    as.formula("fa ~ group + sex + s(x, k = 15) + s(participant_id, bs = 're')", 
#               env = .GlobalEnv)
#   )

#   expect_identical(
#     tractable::build_formula(target = "fa", 
#                              regressors = c("group", "sex", "nodeID"), 
#                              node_k = 15, node_group = "group"), 
#     as.formula("fa ~ group + sex + s(nodeID, by = group, bs = 'fs', k = 15) + s(subjectID, bs = 're')", 
#                env = .GlobalEnv)
#   )

#   expect_identical(
#     tractable::build_formula(
#       target = "fa", regressors = c("age", "group", "s(age, by = group, k = 2)")), 
#     as.formula("fa ~ age + group + s(age, by = group, k = 2) + s(nodeID) + s(subjectID, bs = 're')", 
#                env = .GlobalEnv)
#   )
# })


# test_that("estimate_distribution runs as expected", {
#   x <- stats::rnorm(1000)
#   expect_identical(tractable::estimate_distribution(x), 
#                    stats::gaussian(link = "identity"))
  
#   expect_error(tractable::estimate_distribution("gaussian"))
#   expect_error(tractable::estimate_distribution(
#     x, distr_options = c("logit", "gamma")))
#   expect_error(tractable::estimate_distribution(x, eval_metric = "ML"))

#   x <- stats::rgamma(1000, shape = 12)
#   expect_identical(tractable::estimate_distribution(x), 
#                    stats::Gamma(link = "logit"))
# })


# test_that("estimate_smooth_basis.default runs as expected", {
#   df_sarica <- tractable::read_afq_sarica() %>% 
#     dplyr::filter(tractID == "Right Corticospinal") %>% 
#     dplyr::mutate(subjectID = factor(subjectID), 
#                   gender = factor(gender), 
#                   group = factor(class)) %>% 
#     tidyr::drop_na() 

#   estimated_information <- tractable::estimate_smooth_basis(
#     target       = "md", 
#     smooth_terms = "s(nodeID, bs = 'ts', k = seq(2, 20, 2))", 
#     df           = df_sarica, 
#   )
#   expected_term <- "s(nodeID, bs = 'ts', k = 18)"
#   expect_identical(estimated_information$est_terms, expected_term)
#   expected_values <- estimated_information$k_estimates %>% 
#     dplyr::filter(term == expected_term)
#   expect_identical(expected_values$k_index > 0.95 && expected_values$p_value > 0.05, TRUE)
  
#   estimated_information <- tractable::estimate_smooth_basis(
#     target       = "fa", 
#     smooth_terms = "s(nodeID, by = group, bs = 'fs')", 
#     df           = df_sarica, 
#     regressors   = c("age", "group"), 
#   )
#   expected_term <- "s(nodeID, by = group, bs = 'fs', k = 9)"
#   expect_identical(estimated_information$est_terms, expected_term)
#   expected_values <- estimated_information$k_estimates %>% 
#     dplyr::filter(term == expected_term)
#   expect_identical(expected_values$k_index > 0.95 && expected_values$p_value > 0.05, TRUE)
  
#   estimated_information <- tractable::estimate_smooth_basis(
#     target       = "fa", 
#     smooth_terms = "s(nodeID, by = group, bs = 'fs')", 
#     df           = df_sarica, 
#     regressors   = c("age", "group"), 
#     kindex_thr   = 0.98,
#     pvalue_thr   = 0.1
#   )
#   expected_term <- "s(nodeID, by = group, bs = 'fs', k = 11)"
#   expect_identical(estimated_information$est_terms, expected_term)
#   expected_values <- estimated_information$k_estimates %>% 
#     dplyr::filter(term == expected_term)
#   expect_identical(expected_values$k_index > 0.98 && expected_values$p_value > 0.1, TRUE)

#   estimated_information <- tractable::estimate_smooth_basis(
#     target       = "fa", 
#     smooth_terms = c("s(age, k = 1:10)", "s(nodeID, by = group, bs = 'fs', m = 3)"),
#     df           = df_sarica, 
#     regressors   = "group"
#   )
#   expected_term <- c("s(age, bs = 'tp', k = 10)", "s(nodeID, by = group, bs = 'fs', m = 3, k = 9)")
#   expect_identical(estimated_information$est_terms, expected_term)
#   expected_values <- estimated_information$k_estimates %>% 
#     dplyr::filter(term == expected_term[1]) # age term
#   expect_identical(expected_values$k_index > 0.95 && expected_values$p_value > 0.05, FALSE)
#   expected_values <- estimated_information$k_estimates %>% 
#     dplyr::filter(term == expected_term[2]) # node term
#   expect_identical(expected_values$k_index > 0.95 && expected_values$p_value > 0.05, TRUE)
# })


# test_that("estimate_smooth_basis.formula runs as expected", {
#   df_sarica <- tractable::read_afq_sarica() %>% 
#     dplyr::filter(tractID == "Right Corticospinal") %>% 
#     dplyr::mutate(subjectID = factor(subjectID), 
#                   gender = factor(gender), 
#                   group = factor(class)) %>% 
#     tidyr::drop_na() 

#   estimated_information <- tractable::estimate_smooth_basis(
#     formula = md ~ s(nodeID, bs = "ts", k = seq(2, 20, 2)), 
#     df      = df_sarica
#   )
#   expected_term <- "s(nodeID, bs = 'ts', k = 18)"
#   expect_identical(estimated_information$est_terms, expected_term)
#   expected_values <- estimated_information$k_estimates %>% 
#     dplyr::filter(term == expected_term)
#   expect_identical(expected_values$k_index > 0.95 && expected_values$p_value > 0.05, TRUE)
  
#   estimated_information <- tractable::estimate_smooth_basis(
#     formula = fa ~ age + group + s(nodeID, by = group, bs = "fs"), 
#     df      = df_sarica, 
#   )
#   expected_term <- "s(nodeID, by = group, bs = 'fs', k = 9)"
#   expect_identical(estimated_information$est_terms, expected_term)
#   expected_values <- estimated_information$k_estimates %>% 
#     dplyr::filter(term == expected_term)
#   expect_identical(expected_values$k_index > 0.95 && expected_values$p_value > 0.05, TRUE)
  
#   estimated_information <- tractable::estimate_smooth_basis(
#     formula    = fa ~ age + group + s(nodeID, by = group, bs = "fs"),
#     df         = df_sarica, 
#     kindex_thr = 0.98,
#     pvalue_thr = 0.1
#   )
#   expected_term <- "s(nodeID, by = group, bs = 'fs', k = 11)"
#   expect_identical(estimated_information$est_terms, expected_term)
#   expected_values <- estimated_information$k_estimates %>% 
#     dplyr::filter(term == expected_term)
#   expect_identical(expected_values$k_index > 0.98 && expected_values$p_value > 0.1, TRUE)

#   estimated_information <- tractable::estimate_smooth_basis(
#     formula = fa ~ group + s(age, k = 1:10) + s(nodeID, by = group, bs = "fs", m = 3), 
#     df      = df_sarica, 
#   )
#   expected_term <- c("s(age, bs = 'tp', k = 10)", "s(nodeID, by = group, bs = 'fs', m = 3, k = 9)")
#   expect_identical(estimated_information$est_terms, expected_term)
#   expected_values <- estimated_information$k_estimates %>% 
#     dplyr::filter(term == expected_term[1]) # age term
#   expect_identical(expected_values$k_index > 0.95 && expected_values$p_value > 0.05, FALSE)
#   expected_values <- estimated_information$k_estimates %>% 
#     dplyr::filter(term == expected_term[2]) # node term
#   expect_identical(expected_values$k_index > 0.95 && expected_values$p_value > 0.05, TRUE)
# })




# test_that("fit_gam run as expected", {

#   df_sarica <- tractable::read_afq_sarica() %>% 
#     dplyr::filter(tractID == "Right Corticospinal") %>% 
#     dplyr::mutate(subjectID = factor(subjectID), 
#                   gender = factor(gender), 
#                   group = factor(class)) %>% 
#     tidyr::drop_na() 

# model_fit <- tractable::fit_gam(
#   formula = fa ~ s(age, k = 3), 
#   df = df_sarica, 
#   family = "auto"
# )

# model_fit <- tractable::fit_gam(
#   formula = fa ~ s(nodeID, k = 3) + s(subjectID, bs = "re"), 
#   df = df_sarica, 
# )

# model_fit <- tractable::fit_gam(
#   target = "fa",
#   df = df_sarica, 
#   node_k = 14, 
#   family = "auto"
# )

# estimate_smooth_basis(
#   target = "fa", 
#   df = df_sarica, 
#   regressors = "sex",
  
#   age    = list(k_start = 2, k_end = 5), 
#   nodeID = list(k_start = 2, k_end = 50, bs = "tp")
# )

#   model0 <- tractable::fit_gam(
#     formula = fa ~ group + 
#       s(nodeID, by = group, k = 16) + 
#       s(subjectID, bs = "re") + 
#       s(nodeID, subjectID, bs = "fs", m = 1),
#     df = df_sarica, 
#     family = "norm"
#   )
#   summary(model0)

#   model1 <- tractable::fit_gam(
#     formula = fa ~ group + 
#       s(nodeID, by = group, k = 16) + 
#       s(subjectID, bs = "re") + 
#       s(nodeID, subjectID, bs = "fs", m = 3),
#     df = df_sarica, 
#     family = "norm"
#   )
#   summary(model1)

#   itsadug::compareML(model0, model1)

# })

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