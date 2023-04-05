test_that("MappingBase virtual class", {
  expect_error(MappingBase$new())
})

test_that("ScalarMap mechanics", {
  map <- ScalarMap$new()
  expect_equal(map$type, 'any')
  expect_error(map$type <- 'value')
  expect_output(print(map))
  # Ensure empty get
  expect_error(map$get('a'))
  expect_null(map$get('a', default = NULL))
  # Check set
  expect_no_condition(map$set('a', 1L))
  expect_no_condition(map$get('a'))
  expect_no_condition(map$set('b', c(2L, 3L)))
  # Check properties
  expect_equal(map$keys(), c('a', 'b'))
  expect_identical(names(map), map$keys())
  expect_equal(map$get('a'), 1L)
  expect_equal(map$get('b'), c(2L, 3L))
  expect_length(map$items(), 2L)
  expect_equal(map$length(), 2L)
  expect_identical(length(map), map$length())
  # Check removing a value
  expect_no_condition(map$set('a', NULL))
  expect_length(map$items(), 1L)
  expect_identical(map$items(), list(b=c(2L, 3L)))
  # Check [[ and [[<-
  expect_no_condition(map[['a']] <- 1L)
  expect_equal(map[['a']], 1L)
  # Check setv
  expect_no_condition(map$setv(b = 2L, c = 3L))
  expect_length(map$items(), 3L)
  expect_mapequal(map$items(), list(a = 1L, b = 2L, c = 3L))
  expect_equal(sort(map$keys()), c('a', 'b', 'c'))
  # Check update
  nm <- ScalarMap$new()
  nm$setv(x = TRUE, y = FALSE)
  expect_no_condition(map$update(nm))
  expect_mapequal(
    map$items(),
    list(a = 1L, b = 2L, c = 3L, x = TRUE, y = FALSE)
  )
  mm <- ScalarMap$new()
  mm$setv(a = 'a', z = 'z')
  expect_no_condition(map$update(mm))
  expect_mapequal(
    map$items(),
    list(a = 'a', b = 2L, c = 3L, x = TRUE, y = FALSE, z = 'z')
  )
  expect_no_condition(map$remove("b"))
  expect_equal(map$keys(), c('a', 'c', 'x', 'y', 'z'))
  expect_no_condition(map$remove('x')$remove('z'))
  expect_mapequal(
    map$items(),
    list(a = 'a', c = 3L, y = FALSE)
  )
})

test_that("Scalar Map types", {
  atomics <- c('numeric', 'integer', 'character', 'logical')
  for (i in atomics) {
    expect_no_condition(map <- ScalarMap$new(type = i))
    expect_no_condition(map$set('a', vector(mode = i, length = 1L)))
    expect_equal(map$type, i)
    expect_type(map$get('a'), type = ifelse(i == 'numeric', 'double', i))
    expect_no_condition(map$set('a', NULL))
    for (j in setdiff(atomics, i)) {
      expect_error(map$set('b', vector(mode = j, length = 1L)))
    }
  }
  non_atomics <- c('list', 'data.frame', 'factor', 'matrix', 'array')
  for (i in non_atomics) {
    expect_error(map <- ScalarMap$new(type = i))
  }
})
