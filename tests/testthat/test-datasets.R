test_that("read_afq_sarica loads the sarica dataset", {
  df_sarica <- read_afq_sarica()

  expect_equal(nrow(df_sarica), 96000)
  expect_equal(ncol(df_sarica), 8)
  expect_equal(length(unique(df_sarica$subjectID)), 48)
})

test_that("read_afq_weston_havens loads the WH dataset", {
  df_hbn <- read_afq_weston_havens()

  expect_equal(nrow(df_hbn), 154000)
  expect_equal(ncol(df_hbn), 8)
  expect_equal(length(unique(df_hbn$subjectID)), 77)
})

test_that("read_afq_hbn loads the HBN dataset", {
  df_hbn <- read_afq_hbn(truncate = TRUE)

  expect_equal(nrow(df_hbn), 49)
  expect_equal(ncol(df_hbn), 7)
  expect_equal(length(unique(df_hbn$subjectID)), 1)
})

test_that("read_afq_files returns an unsupervised dataset when pheno_csv is NULL", {
  nodes_url <- "https://github.com/yeatmanlab/Sarica_2017/raw/gh-pages/data/nodes.csv"
  df <- read_afq_files(nodes_csv   = nodes_url,
                       dwi_metrics = c("fa", "md"),
                       pheno_cols  = c("age", "class", "gender"))

  expect_equal(! "age" %in% colnames(df), TRUE)
  expect_equal(! "class" %in% colnames(df), TRUE)
  expect_equal(! "gender" %in% colnames(df), TRUE)
})