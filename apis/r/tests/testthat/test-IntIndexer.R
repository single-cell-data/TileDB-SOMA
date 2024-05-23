.match <- function(x, table) bit64::match(x, table, nomatch = 0L) - 1L

test_that("IntIndexer mechanics", {
  keys <- 1L
  lookups <- rep_len(1L, length.out = 4L)
  expect_s3_class(indexer <- IntIndexer$new(keys), "IntIndexer")
  expect_s3_class(val <- indexer$get_indexer(lookups), 'integer64')
  expect_equal(val, .match(lookups, keys))

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
    .match(lookups, keys)
  )

  keys <- seq.int(1L, 10000L - 1L)
  lookups <- arrow::Array$create(seq.int(1L, 10000L - 1L))
  expect_no_condition(indexer <- IntIndexer$new(keys))
  expect_equal(
   indexer$get_indexer(lookups),
   .match(lookups$as_vector(), keys)
  )

  keys <- c(
    seq.int(1L, 10000L - 1L),
    seq.int(10000L, 20000L - 1L),
    seq.int(30000L, 40000L - 1L)
  )
  lookups <- arrow::chunked_array(list(
    seq.int(1L, 10000L - 1L),
    seq.int(10000L, 20000L - 1L),
    seq.int(30000L, 40000L - 1L)
  ))
  expect_no_condition(indexer <- IntIndexer$new(keys))
  expect_equal(
   indexer$get_indexer(lookups),
   .match(unlist(lookups$as_vector()), keys)
  )

  # Test assertions
  expect_error(IntIndexer$new(arrow::Array$create(seq.int(1L, 10000L - 1L))))
  expect_error(IntIndexer$new(TRUE))
  expect_error(IntIndexer$new(1.1))
  expect_error(IntIndexer$new(list(1L)))
})

test_that("IndIndexer nomatch", {
  keys <- c(-10000, -200000, 1000, 3000)
  lookups <- unlist(replicate(n = 4L, c(-1L, 1:5), simplify = FALSE))
  expect_s3_class(indexer <- IntIndexer$new(keys), "IntIndexer")
  expect_s3_class(vals <- indexer$get_indexer(lookups), "integer64")
  expect_equal(vals, rep_len(-1L, length.out = length(lookups)))

  expect_s3_class(
    vals <- indexer$get_indexer(lookups, nomatch_na = TRUE),
    "integer64"
  )
  expect_true(all(is.na(vals)))

  keys <- c(-10000, -200000, 1000, 3000, 1L, 2L)
  lookups <- 1:5
  expect_s3_class(indexer <- IntIndexer$new(keys), "IntIndexer")
  expect_s3_class(vals <- indexer$get_indexer(lookups), "integer64")
  expect_equal(sum(vals == -1L), 3L)

  expect_s3_class(
    vals <- indexer$get_indexer(lookups, nomatch_na = TRUE),
    "integer64"
  )
  expect_equal(sum(is.na(vals)), 3L)

  # Test assertions
  expect_error(indexer$get_indexer(lookups, nomatch_na = NA))
  expect_error(indexer$get_indexer(lookups, nomatch_na = c(TRUE, TRUE)))
  expect_error(indexer$get_indexer(lookups, nomatch_na = 1.1))
  expect_error(indexer$get_indexer(lookups, nomatch_na = list(1L)))
})
