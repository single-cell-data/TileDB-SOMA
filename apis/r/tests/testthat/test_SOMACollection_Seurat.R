test_that("SOMACollection can be created from a Seurat object", {
  tdb_uri <- withr::local_tempdir("test-soco")
  soco <- SOMACollection$new(uri = tdb_uri, verbose = TRUE)
  expect_true(inherits(soco, "SOMACollection"))

  # uns is created by default
  expect_equal(soco$count_members(), 1L)
  expect_equal(soco$list_members()$NAME, "uns")

  soco$from_seurat(pbmc_small)
  expect_true(inherits(soco$members$RNA, "SOMA"))
  expect_true(inherits(soco$members$uns, "TileDBGroup"))
  expect_mapequal(soco$somas, soco$members["RNA"])

  # check for dimensionality reduction results
  expect_identical(
    names(soco$members$RNA$obsm$members),
    c("dimreduction_pca", "dimreduction_tsne")
  )
  expect_identical(
    names(soco$members$RNA$varm$members),
    "dimreduction_pca"
  )

  # check for graph results
  expect_identical(
    names(soco$members$RNA$obsp$members),
    c("graph_snn")
  )

  # create a new SOMACollection from an existing TileDB group
  soco2 <- SOMACollection$new(uri = tdb_uri, verbose = TRUE)
  expect_true(inherits(soco2, "SOMACollection"))
  expect_true(inherits(soco2$members$RNA, "SOMA"))
  expect_false(inherits(soco2$members$uns, "SOMA"))
  expect_true(inherits(soco2$members$uns, "TileDBGroup"))

  # check for auxiliary arrays
  soma <- soco2$get_member("RNA")

  expect_equal(soma$obsm$count_members(), 2)
  expect_equal(soma$varm$count_members(), 1)
  expect_equal(soma$obsp$count_members(), 1)

  # validate restored aux data
  pbmc_small2 <- soco2$to_seurat()

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
  soco2$from_seurat(pbmc_small)
  expect_identical(
    soco2$somas$RNA$obs$fragment_count(),
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

  soco <- SOMACollection$new(uri, verbose = FALSE)
  soco$from_seurat(pbmc_small)

  # Should not trigger error:
  # Cannot add a different number of cells than already present
  expect_silent(soco$to_seurat())
})

test_that("a dataset with empty cell identities is retrieved", {
  uri <- withr::local_tempdir("empty-cell-identities")
  assay_counts <- SeuratObject::CreateSeuratObject(
    counts = GetAssayData(pbmc_small[["RNA"]])
  )

  # verify identities is unset
  testthat::expect_length(nlevels(SeuratObject::Idents(assay_counts)), 1L)

  soco <- SOMACollection$new(uri, verbose = FALSE)
  soco$from_seurat(assay_counts)

  # Should not trigger error
  expect_silent(soco$to_seurat())
})
