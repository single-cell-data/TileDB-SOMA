#include <Rcpp/Lighter>
#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"

// [[Rcpp::export]]
void soma_collection_create(
    const std::string& uri,
    Rcpp::XPtr<somactx_wrap_t> context,
    Rcpp::Nullable<Rcpp::DatetimeVector> timestamp = R_NilValue) {
    auto timestamp_range = makeTimestampRange(timestamp);
    tiledbsoma::SOMACollection::create(uri, context->ctxptr, timestamp_range);
}
