test_that("Write SummarizedExperiment mechanics", {
  skip_if(!extended_tests() || covr_tests())
  suppressWarnings(suppressMessages(skip_if_not_installed(
    "SummarizedExperiment",
    "1.28.0"
  )))
  skip_if_not_installed("pbmc3k")

  se <- pbmc3k_sce()
  var_df <- SummarizedExperiment::rowData(se)
  features <- rownames(se)

  se <- as(se, "SummarizedExperiment")
  SummarizedExperiment::rowData(se) <- var_df
  rownames(se) <- features

  uri <- withr::local_tempdir("summarized-experiment")

  expect_no_condition(uri <- suppressMessages(write_soma(se, uri, "RNA")))

  expect_type(uri, "character")
  expect_true(grepl("^summarized-experiment", basename(uri)))

  expect_no_condition(experiment <- SOMAExperimentOpen(uri))
  expect_s3_class(experiment, "SOMAExperiment")
  on.exit(experiment$close())

  expect_no_error(experiment$ms)

  expect_equal(experiment$mode(), "READ")
  expect_s3_class(experiment, "SOMAExperiment")
  expect_true(grepl("^summarized-experiment", basename(experiment$uri)))

  expect_s3_class(experiment$obs, "SOMADataFrame")

  expect_identical(experiment$ms$names(), "RNA")
  expect_s3_class(ms <- experiment$ms$get("RNA"), "SOMAMeasurement")

  expect_s3_class(ms$var, "SOMADataFrame")

  expect_identical(
    sort(ms$X$names()),
    sort(SummarizedExperiment::assayNames(se))
  )

  expect_error(ms$obsm)
  expect_error(ms$varm)
  expect_error(ms$obsp)
  expect_error(ms$varp)

  expect_identical(
    setdiff(experiment$obs$attrnames(), "obs_id"),
    names(SummarizedExperiment::colData(se))
  )

  expect_identical(
    setdiff(ms$var$attrnames(), "var_id"),
    names(SummarizedExperiment::rowData(se))
  )

  # Verify X data values round-trip correctly
  for (assay in SummarizedExperiment::assayNames(se)) {
    original <- SummarizedExperiment::assay(se, assay)
    stored <- ms$X$get(assay)$read()$sparse_matrix()$concat()
    expect_equal(
      sum(stored != 0),
      sum(original != 0),
      label = sprintf("non-zero count for assay '%s'", assay)
    )
  }

  # Test ms_name assertions
  expect_error(write_soma(se, uri))
  expect_error(write_soma(se, uri, ""))
  expect_error(write_soma(se, uri, NA_character_))
  expect_error(write_soma(se, uri, c("a", "b")))
  expect_error(write_soma(se, uri, 1))
  expect_error(write_soma(se, uri, TRUE))
})

test_that("Resume-mode adds a second measurement to an existing experiment", {
  suppressWarnings(suppressMessages(skip_if_not_installed(
    "SummarizedExperiment",
    "1.28.0"
  )))

  n_obs <- 5L
  obs <- S4Vectors::DataFrame(id = paste0("s", seq_len(n_obs)))
  rownames(obs) <- obs$id

  make_se <- function(n_var, prefix) {
    x <- matrix(seq_len(n_obs * n_var), nrow = n_var, ncol = n_obs)
    rownames(x) <- paste0(prefix, seq_len(n_var))
    colnames(x) <- obs$id
    var <- S4Vectors::DataFrame(id = rownames(x))
    rownames(var) <- var$id
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = x),
      colData = obs,
      rowData = var
    )
  }

  se1 <- make_se(10L, "a")
  se2 <- make_se(20L, "b")

  uri <- withr::local_tempdir("multi-ms-se")
  write_soma(se1, uri, ms_name = "ms1", ingest_mode = "write")
  write_soma(se2, uri, ms_name = "ms2", ingest_mode = "resume")

  exp <- SOMAExperimentOpen(uri)
  on.exit(exp$close(), add = TRUE, after = FALSE)

  expect_setequal(exp$ms$names(), c("ms1", "ms2"))

  mat2 <- exp$ms$get("ms2")$X$get("counts")$read()$sparse_matrix()$concat()
  # Verify that the second measurement's X data is not all zeros (CX-279)
  expect_true(sum(mat2 != 0) > 0)
})

test_that("Write SummarizedExperiment relatively (SOMA-906)", {
  skip_if(!extended_tests())
  suppressWarnings(suppressMessages(skip_if_not_installed(
    "SummarizedExperiment",
    "1.28.0"
  )))
  skip_if_not_installed("pbmc3k")

  se <- pbmc3k_sce()
  var_df <- SummarizedExperiment::rowData(se)
  features <- rownames(se)

  se <- as(se, "SummarizedExperiment")
  SummarizedExperiment::rowData(se) <- var_df
  rownames(se) <- features

  uri <- tempfile(pattern = "summarizedexperiment-relative-")

  for (i in c(TRUE, FALSE)) {
    expect_error(
      write_soma(se, uri, relative = i, ms_name = "RNA"),
      regexp = "^The dots '\\.\\.\\.' must be empty when",
      label = sprintf("write_soma.SummarizedExperiment(relative = %s)", i)
    )
  }
})
