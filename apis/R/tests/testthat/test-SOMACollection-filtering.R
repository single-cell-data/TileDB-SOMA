test_that("SOMACollection can be sliced by dimension", {
  uri <- withr::local_tempdir("soco")

  var_ids <- c("PPBP", "VDAC3")
  obs_ids <- c("GAGTTGTGGTAGCT", "ATTACCTGCCTTAT", "CTTGATTGATCTTC")

  soco <- SOMACollection$new(uri = uri, verbose = TRUE)
  expect_error(
    soco$set_query(obs_ids = "foo"),
    "SOMACollection must contain a SOMA to query"
  )

  soco$from_seurat(pbmc_small)

  # selected ranges are passed to the SOMA object
  soco$set_query(obs_ids = obs_ids, var_ids = var_ids)

  obs <- soco$somas$RNA$obs$to_dataframe()
  expect_equal(nrow(obs), 3)
  var <- soco$somas$RNA$var$to_dataframe()
  expect_equal(nrow(var), 2)

  pbmc_small2 <- soco$to_seurat()
  expect_equal(dim(pbmc_small2), c(2, 3))

  # clear the query
  soco$reset_query()
  pbmc_small2 <- soco$to_seurat()
  expect_equal(dim(pbmc_small2), dim(pbmc_small))

  # attribute filters
  soco$set_query(
    var_attr_filter = vst.mean > 5,
    obs_attr_filter = groups == "g1"
  )

  pbmc_small2 <- soco$to_seurat()

  expect_equal(
    ncol(pbmc_small2),
    sum(pbmc_small[[]]$groups == "g1")
  )

  expect_equal(
    nrow(pbmc_small2),
    sum(pbmc_small[["RNA"]][[]]$vst.mean > 5)
  )

  # reset the query
  soco$reset_query()
  pbmc_small2 <- soco$to_seurat()
  expect_equal(dim(pbmc_small2), dim(pbmc_small))
})
