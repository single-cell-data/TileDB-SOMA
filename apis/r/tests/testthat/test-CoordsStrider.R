test_that("CoordsStrider start/end mechanics", {
  skip_if_not_installed('iterators')
  skip_if_not_installed('itertools')
  start <- 1L
  end <- 200L
  # Test no stride
  expect_s3_class(strider <- CoordsStrider$new(start = start, end = end), "CoordsStrider")
  expect_equal(strider$stride, end)
  expect_null(strider$coords)
  expect_type(coords <- as.list(strider), "list")
  expect_length(coords, 1L)
  expect_length(coords[[1L]], end)
  expect_equal(coords[[1L]], seq.int(start, end))
  # Test with even stride
  expect_s3_class(
    strider <- CoordsStrider$new(start = start, end = end, stride = 10L),
    "CoordsStrider"
  )
  expect_true(strider$has_next())
  expect_true(itertools::hasNext(strider))
  expect_type(coords <- as.list(strider), "list")
  expect_length(coords, end / strider$stride)
  for (i in seq_along(coords)) {
    expect_length(coords[[i]], as.integer(strider$stride))
  }
  expect_equal(unlist64(coords), seq.int(start, end))
  expect_error(strider$next_element(), "StopIteration", class = "stopIteration")
  expect_error(iterators::nextElem(strider), "StopIteration", class = "stopIteration")
  expect_false(strider$has_next())
  expect_false(itertools::hasNext(strider))
  # Test with uneven stride
  expect_s3_class(
    strider <- CoordsStrider$new(start = start, end = end, stride = 15L),
    "CoordsStrider"
  )
  expect_type(coords <- as.list(strider), "list")
  expect_length(coords, ceiling(end / strider$stride))
  for (i in seq_along(coords)) {
    target <- ifelse(
      i == length(coords),
      yes = end %% as.integer(strider$stride),
      no = as.integer(strider$stride)
    )
    expect_length(coords[[i]], target)
  }
  expect_equal(unlist64(coords), seq.int(start, end))
  # Test with larger stride
  expect_s3_class(
    strider <- CoordsStrider$new(start = start, end = end, stride = 300),
    "CoordsStrider"
  )
  expect_type(coords <- as.list(strider), "list")
  expect_length(coords, 1L)
  expect_length(coords[[1L]], end)
  expect_equal(coords[[1L]], seq.int(start, end))
})

test_that("CoordsStrider coodinate mechanics", {
  skip_if_not_installed('iterators')
  skip_if_not_installed('itertools')
  init <- seq.int(1L, 205L, 3L)
  # Test no stride
  expect_s3_class(strider <- CoordsStrider$new(init), "CoordsStrider")
  expect_equal(strider$stride, length(init))
  expect_equal(strider$start, min(init))
  expect_equal(strider$end, max(init))
  expect_type(coords <- as.list(strider), "list")
  expect_length(coords, 1L)
  expect_length(coords[[1L]], length(init))
  expect_equal(coords[[1L]], init)
  # Test with even stride
  expect_s3_class(
    strider <- CoordsStrider$new(init, stride = 3L),
    "CoordsStrider"
  )
  expect_true(strider$has_next())
  expect_true(itertools::hasNext(strider))
  expect_type(coords <- as.list(strider), "list")
  expect_length(coords, length(init) / as.integer(strider$stride))
  for (i in seq_along(coords)) {
    expect_length(coords[[i]], as.integer(strider$stride))
  }
  expect_equal(unlist64(coords), init)
  expect_error(strider$next_element(), "StopIteration", class = "stopIteration")
  expect_error(iterators::nextElem(strider), "StopIteration", class = "stopIteration")
  expect_false(strider$has_next())
  expect_false(itertools::hasNext(strider))
  # Test with uneven stride
  expect_s3_class(
    strider <- CoordsStrider$new(init, stride = 15L),
    "CoordsStrider"
  )
  expect_type(coords <- as.list(strider), "list")
  expect_length(coords, ceiling(length(init) / as.integer(strider$stride)))
  for (i in seq_along(coords)) {
    target <- ifelse(
      i == length(coords),
      yes = length(init) %% as.integer(strider$stride),
      no = as.integer(strider$stride)
    )
    expect_length(coords[[i]], target)
  }
  expect_equal(unlist64(coords), init)
  # Test with larger stride
  expect_s3_class(
    strider <- CoordsStrider$new(init, stride = 300),
    "CoordsStrider"
  )
  expect_type(coords <- as.list(strider), "list")
  expect_length(coords, 1L)
  expect_length(coords[[1L]], length(init))
  expect_equal(coords[[1L]], init)
})
