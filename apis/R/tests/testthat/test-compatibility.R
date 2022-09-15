setup({
  zipfile <- testthat::test_path("test-data/scdataset-pbmc-small_0-1-2.zip")
  unzip(zipfile, exdir = tempdir())
  scdataset_uri <<- file.path(
    tempdir(),
    "tests/testthat/test-data",
    sub("\\.zip$", "", basename(zipfile))
  )
})

teardown({
  unlink(scdataset_uri, recursive = TRUE)
})

test_that("compatibility is maintained with v0.1.2 SCDataset", {

  expect_warning(
    scdataset <- SCDataset$new(scdataset_uri, verbose = FALSE)
  )
  expect_true(inherits(scdataset, "SCDataset"))

  # verify uns wasn't created in the presence of an existing misc group
  expect_false(dir.exists(file.path(scdataset$uri, "uns")))

  # verify misc group contains only a single array
  expect_identical(
    suppressWarnings(scdataset$misc$count_members()),
    1
  )

  # Members are SCGroups, not SOMAs
  expect_true(inherits(scdataset$members$RNA, "SCGroup"))

  # SCGroup's misc group is accessible via misc
  expect_identical(
    suppressWarnings(scdataset$members$RNA$misc),
    scdataset$members$RNA$uns
  )

  # misc is an alias for uns
  expect_identical(
    suppressWarnings(scdataset$misc$uri),
    scdataset$uns$uri
  )

  expect_silent(pbmc_small2 <- scdataset$to_seurat())
})

test_that("compatibility is maintained with v0.1.2 SCGroup", {
  scgroup_uri <- file.path(scdataset_uri, "scgroup_RNA")
  expect_warning(
    scgroup <- SCGroup$new(scgroup_uri, verbose = FALSE)
  )
  expect_true(inherits(scgroup, "SCGroup"))

  # verify uns wasn't created in the presence of an existing misc group
  expect_false(dir.exists(file.path(scgroup$uri, "uns")))

  # obs dataframe manually added to the SCGroup misc group is present
  expect_identical(
    suppressWarnings(scgroup$misc$count_members()),
    1
  )

  # misc is an alias for uns
  expect_identical(
    suppressWarnings(scgroup$misc$uri),
    scgroup$uns$uri
  )

  expect_silent(
    pbmc_small_assay2 <- scgroup$to_seurat_assay()
  )
})
