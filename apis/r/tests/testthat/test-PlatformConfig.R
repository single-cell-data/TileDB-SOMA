test_that("PlatformConfig mechanics", {
  cfg <- PlatformConfig$new()
  expect_output(print(cfg))
  # Check `set`
  expect_no_condition(cfg$set('plat1', 'op1', 'a', 1L))
  expect_length(cfg, 1L)
  expect_equal(cfg$keys(), 'plat1')
  expect_identical(names(cfg), cfg$keys())
  expect_error(cfg$get('platform1'))
  expect_null(cfg$get('platform1', default = NULL))
  expect_equal(cfg$get('platform1', default = 3L), 3L)
  # Check `get`
  expect_s3_class(opcfg <- cfg$get('plat1'), 'ConfigList')
  expect_length(opcfg, 1L)
  expect_s3_class(map <- opcfg$get('op1'), 'ScalarMap')
  expect_length(map, 1L)
  expect_identical(cfg$get('plat1', 'op1'), map)
  expect_equal(cfg$get('plat1', 'op1', 'a'), 1L)
  expect_error(cfg$get('plat1', 'op1', 'b'))
  expect_null(cfg$get('plat1', 'op1', 'b', default = NULL))
  # Check `set` with map
  map <- ScalarMap$new()
  map$setv(a = TRUE, b = FALSE)
  expect_no_condition(cfg$set('plat1', 'op2', value = map))
  expect_length(cfg$get('plat1'), 2L)
  expect_s3_class(cfg$get('plat1', 'op2'), 'ScalarMap')
  # Check PlatformConfig information
  cfg$set('plat2', 'op1', 'a', 1L)
  expect_length(cfg, 2L)
  expect_equal(cfg$platforms(), c('plat1', 'plat2'))
  expect_equal(cfg$params(), cfg$get('plat1')$keys())
  expect_equal(cfg$params('plat2'), cfg$get('plat2')$keys())
  expect_equal(
    cfg$params(TRUE),
    union(cfg$get('plat1')$keys(), cfg$get('plat2')$keys())
  )
  expect_s3_class(cfg$get_params('plat1'), 'ConfigList')
  expect_error(cfg$get_params('platform1'))
})
