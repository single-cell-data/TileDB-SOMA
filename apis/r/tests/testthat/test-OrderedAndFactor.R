test_that("SOMADataFrame round-trip with factor and ordered", {
    skip_if(!extended_tests())

    uri <- tempfile()

    ## borrowed from tiledb-r test file test_ordered.R
    ## A data.frame with an ordered column, taken from package `earth` and its `etitanic` cleaned

    ## dataset of Titanic survivors (with NAs removed).
    ##
    ## et <- earth::etitanic
    ## et$pclass <- as.ordered(et$pclass)
    ## set.seed(42)
    ## et <- et[sort(sample(nrow(et), 100)), ]
    ## dput(et)
    ##
    ## Slightly edited (for code alignment) `dput(et)` output below
    et <- structure(list(pclass = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                              1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                              1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                                              2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L,
                                              3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L,
                                              3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L,
                                              3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L),
                                            levels = c("1st", "2nd", "3rd"), class = c("ordered", "factor")),
                     survived = c(0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L,
                                  1L, 0L, 0L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 1L, 1L,
                                  0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L,
                                  1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                  0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                  0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                  0L, 0L, 0L),
                     sex = structure(c(1L, 2L, 1L, 1L, 1L, 2L, 1L, 2L,
                                       2L, 1L, 1L, 1L, 2L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 2L, 2L, 1L,
                                       2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 1L,
                                       2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 1L,
                                       2L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 1L,
                                       2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L,
                                       2L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L),
                                     levels = c("female", "male"), class = "factor"),
                     age = c(2, 24, 29, 58, 59, 28, 36,
                             27, 39, 27, 48, 24, 19, 22, 48, 35, 38, 16, 65, 28.5, 35, 34,
                             32, 43, 49, 31, 30, 18, 28, 32, 19, 40, 0.833299994, 19, 37,
                             32, 34, 54, 8, 27, 34, 16, 21, 62, 21, 23, 36, 29, 41, 33, 25,
                             25, 18.5, 13, 20, 6, 32, 21, 18, 26, 32, 29, 18.5, 21, 17, 37,
                             35, 30, 22, 47, 26, 21, 28, 25, 28, 43, 22, 30, 20.5, 51, 35,
                             28, 19, 28, 29, 41, 19, 28, 8, 39, 2, 45, 30, 33, 21, 24, 11.5,
                             18, 36, 45.5),
                     sibsp = c(1L, 0L, 0L, 0L, 2L, 0L, 1L, 1L, 1L,
                               1L, 1L, 3L, 3L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L,
                               0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 1L,
                               0L, 1L, 0L, 2L, 2L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 4L, 1L,
                               0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 2L, 2L,
                               0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 4L,
                               0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L),
                     parch = c(2L, 1L, 0L, 0L, 0L, 0L, 2L, 0L, 0L, 2L, 0L, 2L, 2L, 2L, 0L, 0L, 0L, 1L,
                               0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 2L, 0L,
                               0L, 0L, 1L, 0L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 2L,
                               0L, 0L, 0L, 2L, 0L, 2L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L,
                               0L, 0L, 0L, 0L, 0L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                               1L, 0L, 4L, 5L, 0L, 0L, 1L, 5L, 1L, 4L, 0L, 0L, 0L, 0L, 1L, 0L,
                               0L, 0L)),
                row.names = c("3", "17", "25", "34", "43", "53", "58",
                              "65", "85", "91", "100", "112", "115", "123", "146", "165", "169",
                              "188", "206", "223", "258", "260", "279", "282", "295", "299",
                              "324", "327", "335", "337", "338", "353", "360", "365", "369",
                              "390", "397", "398", "399", "402", "415", "417", "420", "433",
                              "445", "448", "449", "453", "533", "543", "556", "568", "569",
                              "602", "616", "624", "656", "676", "677", "678", "685", "689",
                              "693", "697", "701", "711", "730", "761", "786", "794", "804",
                              "807", "839", "854", "864", "869", "953", "975", "978", "980",
                              "996", "1022", "1051", "1084", "1101", "1107", "1109", "1127",
                              "1146", "1147", "1157", "1212", "1219", "1223", "1225", "1238",
                              "1264", "1289", "1299", "1302"),
                class = "data.frame")
    expect_true(is.data.frame(et))

    ett <- data.frame(soma_joinid=bit64::as.integer64(seq(1, nrow(et))), et)
    ## quick write with tiledb-r so that we get a schema from the manifested array
    ## there should possibly be a helper function to create the schema from a data.frame
    turi <- tempfile()
    expect_silent(tiledb::fromDataFrame(ett, turi, col_index="soma_joinid"))

    tsch <- tiledb::schema(turi)
    expect_true(inherits(tsch, "tiledb_array_schema"))

    sch <- tiledbsoma:::arrow_schema_from_tiledb_schema(tsch)
    expect_true(inherits(sch, "Schema"))

    att <- arrow::as_arrow_table(ett)
    expect_true(inherits(att, "Table"))

    lvls <- tiledbsoma:::extract_levels(att)
    expect_true(is.list(lvls))
    expect_equal(length(lvls), ncol(et))  # et, not ett or tsch or sch as no soma_joinid
    expect_equal(names(lvls), colnames(et))

    sdf <- SOMADataFrameCreate(uri, sch, levels=lvls)
    expect_true(inherits(sdf, "SOMADataFrame"))

    sdf$write(att)

    op <- getOption("arrow.int64_downcast")
    options("arrow.int64_downcast"=FALSE) # else it becomes int
    ndf <- SOMADataFrameOpen(uri)$read()$concat()
    expect_true(inherits(ndf, "Table"))

    expect_equivalent(tibble::as_tibble(ndf), tibble::as_tibble(att))

    options("arrow.int64_downcast"=op)

})
