library(tiledbsoma)

create_arrow_schema <- function() {
  arrow::schema(
    arrow::field("foo", arrow::int32(), nullable = FALSE),
    arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
    arrow::field("bar", arrow::float64(), nullable = FALSE),
    arrow::field("baz", arrow::large_utf8(), nullable = FALSE)
  )
}


uri <- "sdf"
asch <- create_arrow_schema()

sdf <- SOMADataFrame$new(uri)
sdf$create(asch, index_column_names = "soma_joinid")

tbl0 <- arrow::arrow_table(foo = 1L:36L,
                           soma_joinid = 1L:36L,
                           bar = 1.1:36.1,
                           baz = c("á", "ą", "ã", "à", "å", "ä", "æ", "ç", "ć", "Ç", "í",
                                   "ë", "é", "è", "ê", "ł", "Ł", "ñ", "ń", "ó", "ô", "ò",
                                   "ö", "ø", "Ø", "ř", "š", "ś", "ş", "Š", "ú", "ü", "ý",
                                   "ź", "Ž", "Ż"),
                           schema = asch)

sdf$write(tbl0)
