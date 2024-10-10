test_that("DataFrame Factory", {
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
      # TODO: for a follow-up PR
      # arrow::field("timestamp_s", arrow::timestamp(unit="s")),
      # arrow::field("timestamp_ms", arrow::timestamp(unit="ms")),
      # arrow::field("timestamp_us", arrow::timestamp(unit="us")),
      # arrow::field("timestamp_ns", arrow::timestamp(unit="ns"))
      # Not supported in libtiledbsoma
      # arrow::field("datetime_day", arrow::date32())
    )

    sdf <- SOMADataFrameCreate(uri, sch, index_column_names = "soma_joinid")
    expect_true(sdf$exists())
    expect_true(dir.exists(uri))

    tbl <- arrow::arrow_table(
        soma_joinid = 1L:10L,
        int8   =  -11L:-20L,
        int16  = -201L:-210L,
        int32  = -301L:-310L,
        int64  = -401L:-410L,
        uint8  =   11L:20L,
        uint16 =  201L:210L,
        uint32 =  301L:310L,
        uint64 =  401L:410L,
        string = c("apple", "ball", "cat", "dog", "egg", "fig", "goose", "hay", "ice", "jam"),
        large_utf8 = c("APPLE", "BALL", "CAT", "DOG", "EGG", "FIG", "GOOSE", "HAY", "ICE", "JAM"),
        enum = factor(
            c("red", "yellow", "green", "red", "red", "red", "yellow", "green", "red", "green"),
            levels = c("red", "yellow", "green")),
        float32 = 1.5:10.5,
        float64 = 11.5:20.5,
        # TODO: for a follow-up PR
        # timestamp_s  = as.POSIXct(as.numeric(3600 + 1:10), tz="GMT"),
        # timestamp_ms = as.POSIXct(as.numeric(3600*1000 + 1:10), tz="GMT"),
        # timestamp_us = as.POSIXct(as.numeric(3600*1000*1000 + 1:10), tz="GMT"),
        # timestamp_ns = as.POSIXct(as.numeric(3600*1000*1000*1000 + 1:10), tz="GMT"),
        schema = sch)
    sdf$write(tbl)
    sdf$close()

    sdf$reopen("READ")

    good_cases <- list(
      'soma_joinid > 5' = function(df) {
          expect_equal(df$soma_joinid, 6:10)
          expect_equal(df$int32, -306:-310)
      },
      'soma_joinid == 10' = function(df) {
          expect_equal(df$soma_joinid, 10)
          expect_equal(df$int32, -310)
          expect_equal(as.character(df$enum), c("green"))
      },
      'soma_joinid > 4 && soma_joinid < 8' = function(df) {
          expect_equal(df$soma_joinid, 5:7)
          expect_equal(df$string, c("egg", "fig", "goose"))
          expect_equal(df$large_utf8, c("EGG", "FIG", "GOOSE"))
      },
      'soma_joinid < 4 || soma_joinid > 8' = function(df) {
          expect_equal(df$soma_joinid, c(1:3, 9:10))
      },

      'int8 == 8' = function(df) {
          expect_equal(length(df$soma_joinid), 0)
      },
      'int8 == -12' = function(df) {
          expect_equal(df$soma_joinid, c(2))
      },
      'int16 > -203' = function(df) {
          expect_equal(df$soma_joinid, c(1, 2))
      },
      'uint16 < 204' = function(df) {
          expect_equal(df$soma_joinid, c(1, 2, 3))
      },
      'int32 > -303' = function(df) {
          expect_equal(df$soma_joinid, c(1, 2))
      },
      'uint32 < 304' = function(df) {
          expect_equal(df$soma_joinid, c(1, 2, 3))
      },
      'int64 > -403' = function(df) {
          expect_equal(df$soma_joinid, c(1, 2))
      },
      'uint64 < 404' = function(df) {
          expect_equal(df$soma_joinid, c(1, 2, 3))
      },

      'float32 < 4.5' = function(df) {
          expect_equal(df$soma_joinid, c(1, 2, 3))
      },
      'float64 < 14.5' = function(df) {
          expect_equal(df$soma_joinid, c(1, 2, 3))
      },

      'string == "dog"' = function(df) {
          expect_equal(df$soma_joinid, c(4))
      },
      'string %in% c("fig", "dog")' = function(df) {
          expect_equal(df$soma_joinid, c(4, 6))
      },
      'string %nin% c("fig", "dog")' = function(df) {
          expect_equal(df$soma_joinid, c(1, 2, 3, 5, 7, 8, 9, 10))
      },

      'enum == "red"' = function(df) {
          expect_equal(df$soma_joinid, c(1, 4, 5, 6, 9))
      },
      'enum != "red"' = function(df) {
          expect_equal(df$soma_joinid, c(2, 3, 7, 8, 10))
      },
      'enum == "orange"' = function(df) {
          expect_equal(length(df$soma_joinid), 0)
      },
      'enum != "orange"' = function(df) {
          expect_equal(df$soma_joinid, 1:10)
      },
      'enum %in% c("red", "green")' = function(df) {
          expect_equal(df$soma_joinid, c(1, 3, 4, 5, 6, 8, 9, 10))
      }, 
      'enum %nin% c("red", "green")' = function(df) {
          expect_equal(df$soma_joinid, c(2, 7))
      },
      'enum %in% c("orange", "green")' = function(df) {
          expect_equal(df$soma_joinid, c(3, 8, 10))
      },
      'enum %nin% c("orange", "green")' = function(df) {
          expect_equal(df$soma_joinid, c(1, 2, 4, 5, 6, 7, 9))
      },
      'enum %in% c("orange", "purple")' = function(df) {
          expect_equal(length(df$soma_joinid), 0)
      },
      'enum %nin% c("orange", "purple")' = function(df) {
          expect_equal(df$soma_joinid, 1:10)
      }

      # TODO: for a follow-up PR
      # 'timestamp_s < "1969-12-31 20:01:04 EST"' = function(df) {
      #     expect_equal(df$soma_joinid, 1:3)
      # },
      # 'timestamp_ms != "1970-02-11 11:00:05 EST"' = function(df) {
      #     expect_equal(df$soma_joinid, 1:10)
      # },
      # 'timestamp_us > "1970-01-01 00:00:01 GMT"' = function(df) {
      #     expect_equal(df$soma_joinid, 1:10)
      # },
      # 'timestamp_ns > "1970-01-01 00:00:01 GMT"' = function(df) {
      #     expect_equal(df$soma_joinid, 1:10)
      # }
    )

    for (query_string in names(good_cases)) {
        parsed <- do.call(
            what = tiledbsoma:::parse_query_condition_new,
            args = list(expr=str2lang(query_string), schema=sch, somactx=ctx))
        clib_value_filter <- parsed@ptr

        sr <- sr_setup(uri = sdf$uri, ctx, qc=clib_value_filter)
        iter <- TableReadIter$new(sr)
        tbl <- iter$read_next()
        expect_true(iter$read_complete())
        df <- as.data.frame(tbl)
        # Call the validator
        good_cases[[query_string]](df)
    }

    bad_cases <- list(
      '',
      ' ',
      'nonesuch < 10',
      'soma_joinid << 10',
      'soma_joinid',
      'soma_joinid < 4 or soma_joinid > 8'
    )

    for (query_string in names(bad_cases)) {
        expect_error(
            do.call(
                what = tiledbsoma:::parse_query_condition_new,
                args = list(expr=str2lang(query_string), schema=sch, somactx=ctx)))
    }

    sdf$close()
})
