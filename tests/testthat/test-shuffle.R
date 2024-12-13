library(testthat)

test_that("shuffle_df shuffles an afq dataframe", {
  df_afq <- read.csv("afq_12_rows.csv")

  set.seed(0)
  df_shuffle <- shuffle_df(df_afq, "dti_fa")
  df_ref <- read.csv("shuffle_12.csv")

  expect_true(dplyr::all_equal(df_shuffle, df_ref))
})

test_that("shuffle_df with sample_uniform samples from an afq dataframe", {
  df_afq <- read.csv("afq_12_rows.csv")

  set.seed(0)
  df_shuffle <- shuffle_df(df_afq, "dti_fa")

  df_ref <- read.csv("sample_12.csv")
  expect_true(dplyr::all_equal(df_shuffle, df_ref))
})

test_that("bootstrap_df resamples an afq dataframe", {
  df_afq <- read.csv("afq_12_rows.csv")

  set.seed(0)
  df_boot <- bootstrap_df(df_afq, "dti_fa")
  df_ref <- read.csv("boot_12.csv")

  expect_true(dplyr::all_equal(df_boot, df_ref))
})
