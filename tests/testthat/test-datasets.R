test_that("read_afq_files returns an unsupervised dataset when pheno_csv is NULL", {
  nodes_file <- "https://github.com/yeatmanlab/Sarica_2017/raw/gh-pages/data/nodes.csv"
  
  # cannot request columns that do not exist in the dataset
  expect_error(
    df <- read_afq_files(
      nodes_file   = nodes_file,
      other_cols = c("fa", "md", "age", "class", "gender")
    )
  )

  df <- read_afq_files(
    nodes_file   = nodes_file,
    other_cols = c("fa", "md")
  )
  expect_equal(nrow(df), 96000)
  expect_equal(ncol(df), 5)
  expect_equal(length(unique(df$subjectID)), 48)
})


test_that("read_afq_sarica load the sarica dataset", {
  df_sarica <- read_afq_sarica(na_omit = FALSE)
  expect_equal(nrow(df_sarica), 96000)
  expect_equal(ncol(df_sarica), 8)
  expect_equal(length(unique(df_sarica$subjectID)), 48)

  df_sarica <- read_afq_sarica(na_omit = TRUE)
  expect_equal(nrow(df_sarica), 93377)
  expect_equal(ncol(df_sarica), 8)
  expect_equal(length(unique(df_sarica$subjectID)), 48)
})


test_that("read_afq_weston_havens loads the WH dataset", {
  df_wh <- read_afq_weston_havens(na_omit = FALSE)
  expect_equal(nrow(df_wh), 154000)
  expect_equal(ncol(df_wh), 8)
  expect_equal(length(unique(df_wh$subjectID)), 77)

  df_wh <- read_afq_weston_havens(na_omit = TRUE)
  expect_equal(nrow(df_wh), 124389)
  expect_equal(ncol(df_wh), 8)
  expect_equal(length(unique(df_wh$subjectID)), 63)
})


test_that("read_afq_hbn loads the HBN dataset", {
  df_hbn <- read_afq_hbn(truncate = TRUE, na_omit = TRUE)
  expect_equal(nrow(df_hbn), 49)
  expect_equal(ncol(df_hbn), 7)
  expect_equal(length(unique(df_hbn$subjectID)), 1)
})

