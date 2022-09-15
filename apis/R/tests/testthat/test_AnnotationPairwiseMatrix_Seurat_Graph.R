test_that("graph data can be stored and retrieved", {
  uri <- withr::local_tempdir("test-soma-graph")
  soma <- SOMA$new(uri = uri)

  soma$from_seurat_assay(pbmc_small[["RNA"]], obs = pbmc_small[[]])

  # obsp/varp are empty
  expect_length(soma$obsp$arrays, 0L)
  expect_length(soma$varp$arrays, 0L)

  graph1 <- SeuratObject::Graphs(pbmc_small, slot = "RNA_snn")
  soma$obsp$add_seurat_graph(graph1, technique = "snn")

  # obsp/varp are discovered
  soma2 <- SOMA$new(uri = uri)
  expect_length(soma2$obsp$members, 1L)
  expect_length(soma2$varp$members, 0L)

  # validate recreated graph
  graph2 <- soma2$obsp$members$graph_snn$to_seurat_graph()
  expect_identical(
    SeuratObject::DefaultAssay(graph2),
    SeuratObject::DefaultAssay(graph1)
  )
  labs <- rownames(graph1)
  expect_identical(
    graph2[labs, labs],
    graph1[labs, labs]
  )
})
