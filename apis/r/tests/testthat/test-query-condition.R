test_that("DataFrame Factory", {
    skip_if(!extended_tests())
    uri <- tempfile()
    if (dir.exists(uri)) unlink(uri, recursive=TRUE)

    ctx <- soma_context()

    sch <- arrow::schema(
      arrow::field("soma_joinid", arrow::int64()),
      arrow::field("int8",   arrow::int8()),
      arrow::field("int16",  arrow::int16()),
      arrow::field("int32",  arrow::int32()),
      arrow::field("int64",  arrow::int64()),
      arrow::field("uint8",  arrow::uint8()),
      arrow::field("uint16", arrow::uint16()),
      arrow::field("uint32", arrow::uint32()),
      arrow::field("uint64", arrow::uint64()),
      arrow::field("string", arrow::string()),
      arrow::field("large_utf8", arrow::large_utf8()),
      arrow::field("enum",
          arrow::dictionary(
              index_type = arrow::int8(),
              value_type = arrow::utf8(),
              ordered = TRUE)),
      arrow::field("float32", arrow::float32()),
      arrow::field("float64", arrow::float64())
      # arrow::field("float_column", arrow::float64(), nullable = bl),
      # arrow::field("string_column", arrow::large_utf8(), nullable = bl)
      # XXX MORE TO DO
    )

    sdf <- SOMADataFrameCreate(uri, sch, index_column_names = "soma_joinid")
    expect_true(sdf$exists())
    expect_true(dir.exists(uri))

    tbl <- arrow::arrow_table(
        soma_joinid = 1L:10L,
        int32 = 101L:110L,
        string = c("apple", "ball", "cat", "dog", "egg", "fig", "goose", "hay", "ice", "jam"),
        large_utf8 = c("APPLE", "BALL", "CAT", "DOG", "EGG", "FIG", "GOOSE", "HAY", "ICE", "JAM"),
        enum = factor(
            c("red", "yellow", "green", "red", "red", "red", "yellow", "green", "red", "green"),
            levels = c("red", "yellow", "green")),
        float32 = 1.5:10.5,
        float64 = 11.5:20.5,
        schema = sch)
    sdf$write(tbl)
    sdf$close()

    sdf$reopen("READ")

    cases <- list(
      'soma_joinid > 5' = function(df) {
          expect_true(all(df$soma_joinid == 6:10))
          expect_true(all(df$int32 == 106:110))
      },
      'soma_joinid == 10' = function(df) {
          expect_true(all(df$soma_joinid == 10))
          expect_true(all(df$int32 == 110))
          expect_true(all(df$enum == "green"))
      },
      'soma_joinid > 4 && soma_joinid < 8' = function(df) {
          expect_true(all(df$soma_joinid == 5:7))
          expect_true(all(df$string == c("egg", "fig", "goose")))
          expect_true(all(df$large_utr8 == c("EGG", "FIG", "GOOSE")))
      },
      'soma_joinid < 4 || soma_joinid > 8' = function(df) {
          expect_true(all(df$soma_joinid == c(1:3, 9:10)))
      },

      'string == "dog"' = function(df) {
          expect_true(df$soma_joinid == c(4))
      },
      'string %in% c("fig", "dog")' = function(df) {
          expect_true(all(df$soma_joinid == c(4, 6)))
      },
      'string %nin% c("fig", "dog")' = function(df) {
          expect_true(all(df$soma_joinid == c(1, 2, 3, 5, 7, 8, 9, 10)))
      },

      'enum == "red"' = function(df) {
          expect_true(all(df$soma_joinid == c(1, 4, 5, 6, 9)))
      },
      'enum != "red"' = function(df) {
          expect_true(all(df$soma_joinid == c(2, 3, 7, 8, 10)))
      },
      'enum == "orange"' = function(df) {
          expect_true(all(df$soma_joinid == c(4, 6)))
      },
      'enum != "orange"' = function(df) {
          expect_true(all(df$soma_joinid == 1:10))
      },
      'enum %in% c("red", "green")' = function(df) {
          expect_true(all(df$soma_joinid == c(1, 3, 4, 5, 6, 8, 9, 10)))
      }, 
      'enum %nin% c("red", "green")' = function(df) {
          expect_true(all(df$soma_joinid == c(2, 7)))
      },
      'enum %in% c("orange", "green")' = function(df) {
          expect_true(all(df$soma_joinid == c(3, 8, 10)))
      },
      'enum %nin% c("orange", "green")' = function(df) {
          expect_true(all(df$soma_joinid == c(1, 2, 4, 5, 6, 7, 9)))
      },
      'enum %in% c("orange", "yellow")' = function(df) {
          expect_true(all(df$soma_joinid == c()))
      },
      'enum %nin% c("orange", "purple")' = function(df) {
          expect_true(all(df$soma_joinid == 1:10))
      }
    )

    for (query_string in names(cases)) {
        parsed <- do.call(
            what = tiledbsoma:::parse_query_condition_new,
            args = list(expr=str2lang(query_string), schema=sch, somactx=ctx))
        clib_value_filter <- parsed@ptr

        sr <- sr_setup(uri = sdf$uri, ctx, qc=clib_value_filter)
        iter <- TableReadIter$new(sr)
        tbl <- iter$read_next()
        expect_true(iter$read_complete())
        df <- as.data.frame(tbl)
        cases[[query_string]](df)
    }

    sdf$close()
})
