.match <- function(x, table) match(x, table, nomatch = 0L) - 1L

test_that("IntIndexer mechanics", {
  keys <- 1L
  lookups <- rep_len(1L, length.out = 4L)
  expect_s3_class(indexer <- IntIndexer(keys), "IntIndexer")
  expect_equal(
    indexer$get_indexer(lookups),
    .match(lookups, keys)
  )

  keys <- c(-1, 1, 2, 3, 4, 5)
  lookups <- unlist(replicate(n = 4L, c(-1L, 1:5), simplify = FALSE))
  expect_no_condition(indexer <- IntIndexer$new(keys))
  expect_equal(
    indexer$get_indexer(lookups),
    .match(lookups, keys)
  )

  keys <- c(-10000, -100000, 200000, 5, 1, 7)
  expect_no_condition(indexer <- IntIndexer$new(keys))
  expect_equal(
    indexer$get_indexer(lookups),
    .match(lookups, keys)
  )

  keys <- c(-10000, -200000, 1000, 3000, 1, 2)
  expect_no_condition(indexer <- IntIndexer$new(keys))
  expect_equal(
    indexer$get_indexer(lookups),
    .match(lookups, keys)
  )
  keys <- bit64::as.integer64(c(-10000, -200000, 1000, 3000, 1, 2))
  lookups <- bit64::as.integer64(unlist(replicate(n = 4L, c(-1L, 1:5), simplify = FALSE)))
  expect_no_condition(indexer <- IntIndexer$new(keys))
  expect_equal(
    indexer$get_indexer(lookups),
    .match(lookups, keys)
  )

  keys <- seq.int(1L, 10000L)
  lookups <- seq.int(1L, 10L)
  expect_no_condition(indexer <- IntIndexer$new(keys))
  expect_equal(
    indexer$get_indexer(lookups),
    bit64::match(lookups, table = keys, nomatch = 0L) - 1L
  )

  keys <- arrow::Array$create(seq.int(1L, 10000L - 1L))
  lookups <- arrow::Array$create(seq.int(1L, 10000L - 1L))
  expect_no_condition(indexer <- IntIndexer$new(keys))
  expect_equal(
    indexer$get_indexer(lookups),
    arrow::Array$create(.match(lookups$as_vector(), keys$as_vector()))
  )

  keys <- arrow::chunked_array(list(
    seq.int(1L, 10000L - 1L),
    seq.int(10000L, 20000L - 1L),
    seq.int(30000L, 40000L - 1L)
  ))
  lookups <- arrow::chunked_array(list(
    seq.int(1L, 10000L - 1L),
    seq.int(10000L, 20000L - 1L),
    seq.int(30000L, 40000L - 1L)
  ))
  expect_no_condition(indexer <- IntIndexer$new(keys))
  expect_equal(
    indexer$get_indexer(lookups),
    arrow::Array$create(.match(lookups$as_vector(), keys$as_vector()))
  )

  # Test assertions
  expect_error(IntIndexer$new(TRUE))
  expect_error(IntIndexer$new(1.1))
  expect_error(IntIndexer$new(list(1L)))
})
