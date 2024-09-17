test_that("read.afq.sarica loads the sarica dataset", {
  df_sarica <- read.afq.sarica()

  expect_equal(nrow(df_sarica), 96000)
  expect_equal(ncol(df_sarica), 8)
  expect_equal(length(unique(df_sarica$subjectID)), 48)
})

test_that("read.afq.weston.havens loads the WH dataset", {
  df_hbn <- read.afq.weston.havens()

  expect_equal(nrow(df_hbn), 154000)
  expect_equal(ncol(df_hbn), 8)
  expect_equal(length(unique(df_hbn$subjectID)), 77)
})

test_that("read.afq.hbn loads the HBN dataset", {
  df_hbn <- read.afq.hbn(truncate = TRUE)

  expect_equal(nrow(df_hbn), 49)
  expect_equal(ncol(df_hbn), 7)
  expect_equal(length(unique(df_hbn$subjectID)), 1)
})

test_that("read.afq.files returns an unsupervised dataset when pheno_csv is NULL", {
  url.nodes <- "https://github.com/yeatmanlab/Sarica_2017/raw/gh-pages/data/nodes.csv"
  df <- read.afq.files(nodes_csv = url.nodes,
                       dwi_metrics = c("fa", "md"),
                       pheno_cols = c("age", "class", "gender"))

  expect(! "age" %in% colnames(df), "Age should not be in the dataset")
  expect(! "class" %in% colnames(df), "Class should not be in the dataset")
  expect(! "gender" %in% colnames(df), "gender should be in the dataset")
})