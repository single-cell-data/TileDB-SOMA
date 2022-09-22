setup({
  tdb_uri <<- file.path(tempdir(), "test-soma-seurat")
  assay1 <<- pbmc_small[["RNA"]]
})

teardown({
  tiledb::tiledb_vfs_remove_dir(tdb_uri)
})


test_that("SOMA object can be created from a Seurat assay", {

  soma <<- SOMA$new(uri = tdb_uri, verbose = TRUE)

  # should fail if not provided a Seurat assay
  expect_error(
    soma$from_seurat_assay(pbmc_small),
    "A SOMA must be created from a Seurat Assay"
  )

  soma$from_seurat_assay(assay1, obs = pbmc_small[[]])

  # SOMA slot classes
  expect_true(inherits(soma$X, "AssayMatrixGroup"))
  expect_true(inherits(soma$obs, "AnnotationDataframe"))
  expect_true(inherits(soma$var, "AnnotationDataframe"))
  expect_true(inherits(soma$obsm, "AnnotationMatrixGroup"))
  expect_true(inherits(soma$varm, "AnnotationMatrixGroup"))
  expect_true(inherits(soma$obsp, "AnnotationPairwiseMatrixGroup"))
  expect_true(inherits(soma$varp, "AnnotationPairwiseMatrixGroup"))
  expect_true(inherits(soma$uns, "TileDBGroup"))

  # AnnotationGroup dimensions
  expect_equal(soma$X$dimension_name, c("obs_id", "var_id"))
  expect_equal(soma$obsm$dimension_name, "obs_id")
  expect_equal(soma$varm$dimension_name, "var_id")
  expect_equal(soma$obsp$dimension_name, "obs_id")
  expect_equal(soma$varp$dimension_name, "var_id")
})

test_that("Seurat Assay can be recreated from an existing SOMA", {
  soma <- SOMA$new(uri = tdb_uri, verbose = TRUE)

  # SOMA slot classes are restored
  expect_true(inherits(soma$X, "AssayMatrixGroup"))
  expect_true(inherits(soma$obs, "AnnotationDataframe"))
  expect_true(inherits(soma$var, "AnnotationDataframe"))
  expect_true(inherits(soma$obsm, "AnnotationMatrixGroup"))
  expect_true(inherits(soma$varm, "AnnotationMatrixGroup"))
  expect_true(inherits(soma$obsp, "AnnotationPairwiseMatrixGroup"))
  expect_true(inherits(soma$varp, "AnnotationPairwiseMatrixGroup"))
  expect_true(inherits(soma$uns, "TileDBGroup"))

  # AnnotationGroup dimensions are restored
  expect_equal(soma$X$dimension_name, c("obs_id", "var_id"))
  expect_equal(soma$obsm$dimension_name, "obs_id")
  expect_equal(soma$varm$dimension_name, "var_id")
  expect_equal(soma$obsp$dimension_name, "obs_id")
  expect_equal(soma$varp$dimension_name, "var_id")

  # Seurat assay conversion
  assay2 <- soma$to_seurat_assay()

  expect_s4_class(assay2, "Assay")
  expect_equivalent(slot(assay2, "key"), slot(assay1, "key"))

  # use feature/sample names to ensure objects being compared are sorted
  var_ids <- rownames(assay1)
  obs_ids <- colnames(assay1)

  # validate sample metadata, which isn't part of the Seurat Assay so we grab
  # it from the Seurat object
  obs <- pbmc_small[[]]

  # factors are stored in tiledb as character vectors
  expect_equal(
    fac2char(obs)[obs_ids, ],
    soma$obs$to_dataframe()[obs_ids, ]
  )

  # validate feature metadata
  # (manually remove vst.variable column because logicals are returned as ints)
  expect_equal(
    assay2[[]][var_ids, -5],
    assay1[[]][var_ids, -5]
  )

  # validate variable features
  expect_identical(
    SeuratObject::VariableFeatures(assay2),
    sort(SeuratObject::VariableFeatures(assay1))
  )

  # validate raw counts matrix
  expect_identical(
    SeuratObject::GetAssayData(assay2, "counts")[var_ids, obs_ids],
    SeuratObject::GetAssayData(assay1, "counts")[var_ids, obs_ids]
  )

  # validate normalized data matrix
  expect_identical(
    SeuratObject::GetAssayData(assay2, "data")[var_ids, obs_ids],
    SeuratObject::GetAssayData(assay1, "data")[var_ids, obs_ids]
  )

  # Seurat assay conversion with batch mode on
  with_allocation_size_preference(1e4)
  expect_message(
    assay3 <- soma$to_seurat_assay(batch_mode = TRUE, layers = "counts"),
    "...reading in batches"
  )

  expect_identical(
    SeuratObject::GetAssayData(assay3, "counts")[var_ids, obs_ids],
    SeuratObject::GetAssayData(assay1, "counts")[var_ids, obs_ids]
  )
})

test_that("Individual layers can be retrieved from an existing SOMA", {
  soma <<- SOMA$new(uri = tdb_uri, verbose = TRUE)
  expect_s4_class(soma$to_seurat_assay(layers = "counts"), "Assay")
  expect_s4_class(soma$to_seurat_assay(layers = "data"), "Assay")
  expect_error(
    soma$to_seurat_assay(layers = "scale.data"),
    "Creation of a Seurat Assay requires either 'counts' or 'data'"
  )
})

test_that("obs and var are created when even no annotations are present", {
  uri <- withr::local_tempdir("assay-with-no-annotations")

  assay <- SeuratObject::CreateAssayObject(
    counts = SeuratObject::GetAssayData(pbmc_small[["RNA"]], "counts")
  )
  expect_equal(ncol(assay[[]]), 0L)

  soma <- SOMA$new(uri = uri)
  soma$from_seurat_assay(assay)

  obs <- soma$obs$to_dataframe()
  expect_length(obs, 0)
  expect_setequal(rownames(obs), colnames(assay))

  var <- soma$var$to_dataframe()
  expect_length(var, 0)
  expect_setequal(rownames(var), rownames(assay))
})

test_that("dimensional reduction data can be stored and retrieved", {
  soma <- SOMA$new(uri = tdb_uri)

  # obsm/varm are empty
  expect_length(soma$obsm$members, 0L)
  expect_length(soma$varm$members, 0L)

  user_md <- list(foo = "bar")
  pca1 <- SeuratObject::Reductions(pbmc_small, slot = "pca")
  soma$add_seurat_dimreduction(pca1, technique = "pca", metadata = user_md)

  # obsm/varm are discovered
  soma <- SOMA$new(uri = tdb_uri)
  expect_length(soma$obsm$members, 1L)
  expect_length(soma$varm$members, 1L)

  # check dimreduction metadata
  expect_identical(
    soma$obsm$members[[1]]$get_metadata(key = "dimreduction_technique"),
    "pca"
  )
  expect_identical(
    soma$obsm$members[[1]]$get_metadata(key = "dimreduction_key"),
    "PC_"
  )
  expect_identical(
    soma$obsm$members[[1]]$get_metadata(key = "foo"),
    "bar"
  )

  # validate recreated dimreduction data
  pca2 <- soma$get_seurat_dimreduction()

  var_ids <- rownames(SeuratObject::Loadings(pca1))
  expect_identical(
    SeuratObject::Loadings(pca2)[var_ids, ],
    SeuratObject::Loadings(pca1)[var_ids, ]
  )

  obs_ids <- SeuratObject::Cells(pca1)
  expect_identical(
    SeuratObject::Embeddings(pca2)[obs_ids, ],
    SeuratObject::Embeddings(pca1)[obs_ids, ]
  )

  # tsne results only include cell-aligned Embeddings
  tsne1 <- SeuratObject::Reductions(pbmc_small, slot = "tsne")
  soma$add_seurat_dimreduction(tsne1, technique = "tsne")
  tsne2 <- soma$get_seurat_dimreduction(technique = "tsne")

  expect_identical(
    SeuratObject::Embeddings(tsne2)[obs_ids, ],
    SeuratObject::Embeddings(tsne1)[obs_ids, ]
  )
})

test_that("creation from a Seurat Assay without scale.data", {
  uri <- withr::local_tempdir()

  assay1 <- SeuratObject::SetAssayData(
    assay1,
    slot = "scale.data",
    new.data = new(Class = "matrix")
  )

  soma <- SOMA$new(uri = uri, verbose = FALSE)
  testthat::expect_silent(soma$from_seurat_assay(assay1))

  assay2 <- soma$to_seurat_assay()
  testthat::expect_equal(
    SeuratObject::GetAssayData(assay2, "scale.data"),
    SeuratObject::GetAssayData(assay1, "scale.data")
  )
})

test_that("an assay with scale.data containing all features", {
  uri <- withr::local_tempdir()

  # create a seurat object containing on
  object1 <- pbmc_small[SeuratObject::VariableFeatures(pbmc_small),]
  SeuratObject::VariableFeatures(object1) <- character()
  assay1 <- object1[["RNA"]]

  expect_identical(
    dim(SeuratObject::GetAssayData(assay1, "counts")),
    dim(SeuratObject::GetAssayData(assay1, "scale.data"))
  )

  soma <- SOMA$new(uri = uri, verbose = FALSE)
  testthat::expect_silent(soma$from_seurat_assay(assay1))

  var_ids <- rownames(assay1)
  obs_ids <- colnames(assay1)
  assay2 <- soma$to_seurat_assay()
  testthat::expect_equal(
    SeuratObject::GetAssayData(assay2, "scale.data")[var_ids, obs_ids],
    SeuratObject::GetAssayData(assay1, "scale.data")[var_ids, obs_ids]
  )
})


test_that("an assay with empty counts slot can be converted", {
  uri <- withr::local_tempdir("assay-with-empty-counts")
  assay <- SeuratObject::SetAssayData(
    pbmc_small[["RNA"]],
    slot = "counts",
    new.data = new(Class = "matrix")
  )

  expect_true(is_empty(SeuratObject::GetAssayData(assay, "counts")))

  soma <- SOMA$new(uri, verbose = FALSE)
  expect_silent(soma$from_seurat_assay(assay))
  expect_match(tiledb::tiledb_object_type(soma$X$uri), "GROUP")

  assay2 <- soma$to_seurat_assay()

  rlabs <- rownames(assay)
  clabs <- colnames(assay)

  expect_identical(
    SeuratObject::GetAssayData(assay2, "data")[rlabs, clabs],
    SeuratObject::GetAssayData(assay, "data")[rlabs, clabs]
  )

  expect_identical(
    SeuratObject::GetAssayData(assay2, "counts"),
    SeuratObject::GetAssayData(assay, "counts")
  )
})

test_that("an assay with empty feature metdata can be converted", {
  uri <- withr::local_tempdir("assay-without-feature-metadata")

  assay <- SeuratObject::CreateAssayObject(
    counts = SeuratObject::GetAssayData(pbmc_small[["RNA"]], "counts")
  )
  expect_length(assay[[]], 0L)

  soma <- SOMA$new(uri, verbose = FALSE)
  expect_silent(soma$from_seurat_assay(assay))
  assay2 <- soma$to_seurat_assay()
  expect_identical(assay2[[]][rownames(assay),], assay[[]])
})

test_that("individual layers can be added or updated", {
  uri <- withr::local_tempdir("assay-with-individual-layers")

  assay <- SeuratObject::CreateAssayObject(
    counts = SeuratObject::GetAssayData(pbmc_small[["RNA"]], "counts")
  )

  soma <- SOMA$new(uri, verbose = TRUE)
  soma$from_seurat_assay(assay)

  # only counts was created
  expect_equal(names(soma$members$X$members), "counts")

  # add data layer
  soma$from_seurat_assay(assay, var = FALSE, layers = "data")
  expect_equal(names(soma$members$X$members), c("counts", "data"))
  # X$counts, obs and var were not updated
  expect_identical(soma$obs$fragment_count(), 1)
  expect_identical(soma$var$fragment_count(), 1)
  expect_identical(soma$X$members$counts$fragment_count(), 1)

  # update data layer
  soma$from_seurat_assay(assay, var = FALSE, layers = "data")
  expect_equal(names(soma$members$X$members), c("counts", "data"))
  expect_identical(soma$X$members$counts$fragment_count(), 1)
  expect_identical(soma$X$members$data$fragment_count(), 2)
})
