test_that("SOMA object can be sliced by dimension", {
  uri <- withr::local_tempdir("soma-dim-slice1")

  pbmc_small_rna <- pbmc_small[["RNA"]]
  var_ids <- c("PPBP", "VDAC3")
  obs_ids <- c("GAGTTGTGGTAGCT", "ATTACCTGCCTTAT", "CTTGATTGATCTTC")

  soma <- SOMA$new(uri = uri, verbose = TRUE)
  soma$from_seurat_assay(pbmc_small_rna, obs = pbmc_small[[]])

  # slice by obs
  soma$set_query(obs_ids = obs_ids[1])

  obs <- soma$obs$to_dataframe()
  expect_is(obs, "data.frame")
  expect_equal(nrow(obs), 1)

  mat_counts <- soma$X$members$counts$to_matrix(transpose = TRUE)
  expect_true(is_matrix(mat_counts))
  expect_equal(ncol(mat_counts), 1)

  # slice by obs and var
  soma$set_query(obs_ids = obs_ids, var_ids = var_ids)

  obs <- soma$obs$to_dataframe()
  expect_equal(nrow(obs), 3)
  var <- soma$var$to_dataframe()
  expect_equal(nrow(var), 2)

  mat_counts <- soma$X$members$counts$to_matrix(transpose = TRUE)
  expect_equal(dim(mat_counts), c(2, 3))

  pbmc_small_rna2 <- soma$to_seurat_assay()
  expect_equal(dim(pbmc_small_rna2), c(2, 3))


  # var attribute filter (obs attributes are at Seurat object level)
  soma <- SOMA$new(uri = uri)
  soma$set_query(var_attr_filter = vst.mean > 5)

  obs <- soma$obs$to_dataframe()
  expect_equal(nrow(obs), 80)
  var <- soma$var$to_dataframe()
  expect_equal(nrow(var), 12)

  mat_counts <- soma$X$members$counts$to_matrix(transpose = TRUE)
  expect_equal(dim(mat_counts), c(12, 80))

  # var attribute filter + obs/var dimension slicing
  soma$set_query(
    var_attr_filter = vst.mean > 5,
    obs_ids = obs_ids,
    # NCOA4 doesn't meet the attribute filter
    var_ids = c("COTL1", "CST3", "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "NCOA4")
  )

  obs <- soma$obs$to_dataframe()
  expect_equal(nrow(obs), 3)
  var <- soma$var$to_dataframe()
  expect_equal(nrow(var), 5)

  mat_counts <- soma$X$members$counts$to_matrix(transpose = TRUE)
  expect_equal(dim(mat_counts), c(5, 3))

  pbmc_small_rna2 <- soma$to_seurat_assay()
  expect_equal(dim(pbmc_small_rna2), c(5, 3))

  # reset the query
  soma$reset_query()
  pbmc_small_rna3 <- soma$to_seurat_assay()
  expect_equal(dim(pbmc_small_rna3), dim(pbmc_small_rna))
})
