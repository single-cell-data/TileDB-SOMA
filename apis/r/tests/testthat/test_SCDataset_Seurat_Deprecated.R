test_that("SCDataset can be created from a Seurat object", {
  tdb_uri <- withr::local_tempdir("test-scdataset")

  expect_warning(
    scdataset <- SCDataset$new(uri = tdb_uri, verbose = TRUE),
    "'SCDataset' is deprecated."
  )
  expect_true(inherits(scdataset, "SCDataset"))

  # uns is created by default
  expect_equal(scdataset$count_members(), 1L)
  expect_equal(scdataset$list_members()$NAME, "uns")

  scdataset$from_seurat(pbmc_small)
  expect_true(inherits(scdataset$members$RNA, "SOMA"))
  expect_true(inherits(scdataset$members$uns, "TileDBGroup"))
  expect_warning(scdataset$scgroups, "'scgroups' is deprecated")
  expect_mapequal(
    suppressWarnings(scdataset$scgroups),
    scdataset$members["RNA"]
  )

  # check for dimensionality reduction results
  expect_identical(
    names(scdataset$members$RNA$obsm$members),
    c("dimreduction_pca", "dimreduction_tsne")
  )
  expect_identical(
    names(scdataset$members$RNA$varm$members),
    "dimreduction_pca"
  )

  # check for graph results
  expect_identical(
    names(scdataset$members$RNA$obsp$members),
    c("graph_snn")
  )

  # create a new SCDataset from an existing TileDB group
  expect_warning(
    scdataset2 <- SCDataset$new(uri = tdb_uri, verbose = TRUE),
    "'SCDataset' is deprecated."
  )
  expect_true(inherits(scdataset2, "SCDataset"))
  expect_true(inherits(scdataset2$members$RNA, "SOMA"))
  expect_false(inherits(scdataset2$members$uns, "SOMA"))
  expect_true(inherits(scdataset2$members$uns, "TileDBGroup"))

  # check for auxiliary arrays
  scgroup <- scdataset2$get_member("RNA")

  expect_equal(scgroup$obsm$count_members(), 2)
  expect_equal(scgroup$varm$count_members(), 1)
  expect_equal(scgroup$obsp$count_members(), 1)

  # validate restored aux data
  pbmc_small2 <- scdataset2$to_seurat()

  reductions <- SeuratObject::Reductions(pbmc_small)
  for (r in reductions) {
    reduc1 <- SeuratObject::Reductions(pbmc_small, r)
    reduc2 <- SeuratObject::Reductions(pbmc_small2, r)

    load1 <- SeuratObject::Loadings(reduc1)
    load2 <- SeuratObject::Loadings(reduc2)
    expect_identical(load2[rownames(load1), ], load1)

    embed1 <- SeuratObject::Embeddings(reduc1)
    embed2 <- SeuratObject::Embeddings(reduc2)
    expect_identical(embed2[rownames(embed1), ], embed1)
  }

  # check for cell identities
  # factors are stored in tiledb as character vectors
  obs_ids <- SeuratObject::Cells(pbmc_small)
  expect_equal(
    as.character(SeuratObject::Idents(pbmc_small2)[obs_ids]),
    as.character(SeuratObject::Idents(pbmc_small)[obs_ids])
  )

  expect_identical(
    SeuratObject::Graphs(pbmc_small2),
    SeuratObject::Graphs(pbmc_small)
  )

  # check for commands
  command_names <- SeuratObject::Command(object=pbmc_small)
  command_names2 <- SeuratObject::Command(object=pbmc_small2)
  expect_identical(command_names, command_names2)

  # calling from_seurat again will update the existing data
  scdataset2$from_seurat(pbmc_small)
  expect_identical(
    suppressWarnings(scdataset2$scgroups$RNA$obs$fragment_count()),
    2
  )
})

test_that("a dataset containing an assay with empty cells is fully retrieved", {
  uri <- withr::local_tempdir("assay-with-empty-cells")

  cell_ids <- SeuratObject::Cells(pbmc_small)

  # remove all counts for a subset of cells
  counts2 <- SeuratObject::GetAssayData(pbmc_small[["RNA"]], "counts")
  counts2[, 1:10] <- 0
  pbmc_small[["RNA2"]] <- SeuratObject::CreateAssayObject(counts = counts2)

  expect_warning(
    scdataset <- SCDataset$new(uri = uri, verbose = FALSE),
    "'SCDataset' is deprecated."
  )
  scdataset$from_seurat(pbmc_small)

  # Should not trigger error:
  # Cannot add a different number of cells than already present
  expect_silent(scdataset$to_seurat())
})

test_that("a dataset with empty cell identities is retrieved", {
  uri <- withr::local_tempdir("empty-cell-identities")
  assay_counts <- SeuratObject::CreateSeuratObject(
    counts = GetAssayData(pbmc_small[["RNA"]])
  )

  # verify identities is unset
  testthat::expect_length(nlevels(SeuratObject::Idents(assay_counts)), 1L)

  expect_warning(
    scdataset <- SCDataset$new(uri, verbose = FALSE),
    "'SCDataset' is deprecated."
  )
  scdataset$from_seurat(assay_counts)

  # Should not trigger error
  expect_silent(scdataset$to_seurat())
})
