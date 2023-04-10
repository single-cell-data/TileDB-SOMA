test_that("ConfigList mechanics", {
  cfg <- ConfigList$new()
  expect_output(print(cfg))
  # Check `set`
  expect_no_condition(cfg$set('op1', 'a', 1L))
  expect_equal(cfg$keys(), 'op1')
  expect_length(cfg, 1L)
  expect_s3_class(map <- cfg$get('op1'), 'ScalarMap')
  expect_length(map, 1L)
  expect_mapequal(map$items(), list(a = 1L))
  # Check `set` with map
  map <- ScalarMap$new()
  map$setv(a = 'a', b = 'b')
  expect_no_condition(cfg$set('op2', value = map))
  expect_length(cfg, 2L)
  expect_equal(cfg$keys(), c('op1', 'op2'))
  expect_s3_class(map2 <- cfg$get('op2'), 'ScalarMap')
  expect_length(map2, 2L)
  # Check `set` errors
  expect_error(cfg$set(c('op1', 'op2'), 'a', 1L))
  expect_error(cfg$set('op1', c('a', 'b'), c(1L, 2L)))
})
