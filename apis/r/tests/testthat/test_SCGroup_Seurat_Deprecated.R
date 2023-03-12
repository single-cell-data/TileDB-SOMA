setup({
  tdb_uri <<- file.path(tempdir(), "test-scgroup")
  assay1 <<- pbmc_small[["RNA"]]
})

teardown({
  unlink(tdb_uri, recursive = TRUE)
})


test_that("SCGroup object can be created from a Seurat assay", {

  expect_warning(
    scgroup <- SCGroup$new(uri = tdb_uri, verbose = TRUE),
    "'SCGroup' is deprecated."
  )
  expect_true(inherits(scgroup, "SCGroup"))

  # default groups were created
  expect_equal(scgroup$count_members(), 6L)
  expect_setequal(
    scgroup$list_members()$NAME,
    c("X", "obsm", "varm", "obsp", "varp", "uns")
  )

  expect_error(
    scgroup$from_seurat_assay(pbmc_small),
    "A SOMA must be created from a Seurat Assay"
  )

  scgroup$from_seurat_assay(assay1, obs = pbmc_small[[]])

  # SCGroup slot classes
  expect_true(inherits(scgroup$X, "AssayMatrixGroup"))
  expect_true(inherits(scgroup$obs, "AnnotationDataframe"))
  expect_true(inherits(scgroup$var, "AnnotationDataframe"))
  expect_true(inherits(scgroup$obsm, "AnnotationMatrixGroup"))
  expect_true(inherits(scgroup$varm, "AnnotationMatrixGroup"))
  expect_true(inherits(scgroup$obsp, "AnnotationPairwiseMatrixGroup"))
  expect_true(inherits(scgroup$varp, "AnnotationPairwiseMatrixGroup"))
  expect_true(inherits(scgroup$uns, "TileDBGroup"))

  # AnnotationGroup dimensions
  expect_equal(scgroup$X$dimension_name, c("obs_id", "var_id"))
  expect_equal(scgroup$obsm$dimension_name, "obs_id")
  expect_equal(scgroup$varm$dimension_name, "var_id")
  expect_equal(scgroup$obsp$dimension_name, "obs_id")
  expect_equal(scgroup$varp$dimension_name, "var_id")
})

test_that("Seurat Assay can be recreated from an existing SCGroup", {
  expect_warning(
    scgroup <- SCGroup$new(uri = tdb_uri, verbose = TRUE),
    "'SCGroup' is deprecated."
  )


  # SCGroup slot classes are restored
  expect_true(inherits(scgroup$X, "AssayMatrixGroup"))
  expect_true(inherits(scgroup$obs, "AnnotationDataframe"))
  expect_true(inherits(scgroup$var, "AnnotationDataframe"))
  expect_true(inherits(scgroup$obsm, "AnnotationMatrixGroup"))
  expect_true(inherits(scgroup$varm, "AnnotationMatrixGroup"))
  expect_true(inherits(scgroup$obsp, "AnnotationPairwiseMatrixGroup"))
  expect_true(inherits(scgroup$varp, "AnnotationPairwiseMatrixGroup"))
  expect_true(inherits(scgroup$uns, "TileDBGroup"))

  # AnnotationGroup dimensions are restored
  expect_equal(scgroup$X$dimension_name, c("obs_id", "var_id"))
  expect_equal(scgroup$obsm$dimension_name, "obs_id")
  expect_equal(scgroup$varm$dimension_name, "var_id")
  expect_equal(scgroup$obsp$dimension_name, "obs_id")
  expect_equal(scgroup$varp$dimension_name, "var_id")

  # Seurat assay conversion
  assay2 <- scgroup$to_seurat_assay()

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
    scgroup$obs$to_dataframe()[obs_ids, ]
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
})

test_that("Individual layers can be retrieved from an existing SCGroup", {
    expect_warning(
    scgroup <- SCGroup$new(uri = tdb_uri, verbose = TRUE),
    "'SCGroup' is deprecated."
  )
  expect_s4_class(scgroup$to_seurat_assay(layers = "counts"), "Assay")
  expect_s4_class(scgroup$to_seurat_assay(layers = "data"), "Assay")
  expect_error(
    scgroup$to_seurat_assay(layers = "scale.data"),
    "Creation of a Seurat Assay requires either 'counts' or 'data'"
  )
})

test_that("obs and var are created when even no annotations are present", {
  uri <- withr::local_tempdir("assay-with-no-annotations")

  assay <- SeuratObject::CreateAssayObject(
    counts = SeuratObject::GetAssayData(pbmc_small[["RNA"]], "counts")
  )
  expect_equal(ncol(assay[[]]), 0L)

  expect_warning(
    scgroup <- SCGroup$new(uri = uri),
    "'SCGroup' is deprecated."
  )
  scgroup$from_seurat_assay(assay)

  obs <- scgroup$obs$to_dataframe()
  expect_length(obs, 0)
  expect_setequal(rownames(obs), colnames(assay))

  var <- scgroup$var$to_dataframe()
  expect_length(var, 0)
  expect_setequal(rownames(var), rownames(assay))
})

test_that("dimensional reduction data can be stored and retrieved", {
  expect_warning(
    scgroup <- SCGroup$new(uri = tdb_uri),
    "'SCGroup' is deprecated."
  )

  # obsm/varm are empty
  expect_length(scgroup$obsm$members, 0L)
  expect_length(scgroup$varm$members, 0L)

  user_md <- list(foo = "bar")
  pca1 <- SeuratObject::Reductions(pbmc_small, slot = "pca")
  scgroup$add_seurat_dimreduction(pca1, technique = "pca", metadata = user_md)

  # obsm/varm are discovered
  expect_warning(
    scgroup <- SCGroup$new(uri = tdb_uri),
    "'SCGroup' is deprecated."
  )
  expect_length(scgroup$obsm$members, 1L)
  expect_length(scgroup$varm$members, 1L)

  # check dimreduction metadata
  expect_identical(
    scgroup$obsm$members[[1]]$get_metadata(key = "dimreduction_technique"),
    "pca"
  )
  expect_identical(
    scgroup$obsm$members[[1]]$get_metadata(key = "dimreduction_key"),
    "PC_"
  )
  expect_identical(
    scgroup$obsm$members[[1]]$get_metadata(key = "foo"),
    "bar"
  )

  # validate recreated dimreduction data
  pca2 <- scgroup$get_seurat_dimreduction()

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
  scgroup$add_seurat_dimreduction(tsne1, technique = "tsne")
  tsne2 <- scgroup$get_seurat_dimreduction(technique = "tsne")

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

  expect_warning(
    scgroup <- SCGroup$new(uri = uri, verbose = FALSE),
    "'SCGroup' is deprecated."
  )
  testthat::expect_silent(scgroup$from_seurat_assay(assay1))

  assay2 <- scgroup$to_seurat_assay()
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

  expect_warning(
    scgroup <- SCGroup$new(uri = uri, verbose = FALSE),
    "'SCGroup' is deprecated."
  )
  testthat::expect_silent(scgroup$from_seurat_assay(assay1))

  var_ids <- rownames(assay1)
  obs_ids <- colnames(assay1)
  assay2 <- scgroup$to_seurat_assay()
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

  expect_warning(
    scgroup <- SCGroup$new(uri = uri, verbose = FALSE),
    "'SCGroup' is deprecated."
  )
  expect_silent(scgroup$from_seurat_assay(assay))
  expect_match(tiledb::tiledb_object_type(scgroup$X$uri), "GROUP")

  assay2 <- scgroup$to_seurat_assay()

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

  expect_warning(
    scgroup <- SCGroup$new(uri = uri, verbose = FALSE),
    "'SCGroup' is deprecated."
  )
  expect_silent(scgroup$from_seurat_assay(assay))
  assay2 <- scgroup$to_seurat_assay()
  expect_identical(assay2[[]][rownames(assay),], assay[[]])
})

test_that("individual layers can be added or updated", {
  uri <- withr::local_tempdir("assay-with-individual-layers")

  assay <- SeuratObject::CreateAssayObject(
    counts = SeuratObject::GetAssayData(pbmc_small[["RNA"]], "counts")
  )

  expect_warning(
    scgroup <- SCGroup$new(uri = uri, verbose = FALSE),
    "'SCGroup' is deprecated."
  )
  scgroup$from_seurat_assay(assay)

  # only counts was created
  expect_equal(names(scgroup$members$X$members), "counts")

  # add data layer
  scgroup$from_seurat_assay(assay, var = FALSE, layers = "data")
  expect_equal(names(scgroup$members$X$members), c("counts", "data"))
  # X$counts, obs and var were not updated
  expect_identical(scgroup$obs$fragment_count(), 1)
  expect_identical(scgroup$var$fragment_count(), 1)
  expect_identical(scgroup$X$members$counts$fragment_count(), 1)

  # update data layer
  scgroup$from_seurat_assay(assay, var = FALSE, layers = "data")
  expect_equal(names(scgroup$members$X$members), c("counts", "data"))
  expect_identical(scgroup$X$members$counts$fragment_count(), 1)
  expect_identical(scgroup$X$members$data$fragment_count(), 2)
})
