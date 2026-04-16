test_that("SOMATileDBContext is defunct", {
  lifecycle::expect_defunct(SOMATileDBContext$new())
})

test_that("soma_context() is defunct", {
  lifecycle::expect_defunct(soma_context())
})
