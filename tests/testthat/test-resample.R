test_that("shuffle_df shuffles an afq dataframe as expected", {
  df_afq <- readr::read_csv("afq_12_rows.csv", show_col_types = FALSE)

  set.seed(0) # set seed for shuffling
  df_shuffled <- shuffle_df(
    df_afq, target = "dti_fa", node_group = "group", sample_uniform = FALSE)
  df_reference <- readr::read_csv("shuffle_12.csv", show_col_types = FALSE)
  expect_equal(df_shuffled, df_reference)
})


test_that("shuffle_df with sample_uniform samples from an afq dataframe", {
  df_afq <- readr::read_csv("afq_12_rows.csv", show_col_types = FALSE)

  set.seed(10) # set seed for shuffling
  df_shuffled <- shuffle_df(
    df_afq, target = "dti_fa", node_group = "group", sample_uniform = TRUE)
  df_reference <- readr::read_csv("sample_12.csv", show_col_types = FALSE)
  expect_equal(df_shuffled, df_reference)
})


test_that("bootstrap_df resamples an afq dataframe", {
  df_afq <- readr::read_csv("afq_12_rows.csv", show_col_types = FALSE)

  set.seed(0) # set seed for shuffling
  df_bootstrap <- bootstrap_df(df_afq, target = "dti_fa", node_group = "group")
  df_reference <- readr::read_csv("boot_12.csv", show_col_types = FALSE)
  expect_equal(df_bootstrap, df_reference)
})
