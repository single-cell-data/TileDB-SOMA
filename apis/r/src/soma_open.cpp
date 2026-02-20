#include <Rcpp/Lighter>

#include <optional>

#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"
#include "xptr-utils.h"  // xptr taggging utilities

// [[Rcpp::export]]
Rcpp::XPtr<tiledbsoma::SOMADataFrame> open_dataframe_handle(
    const std::string& uri,
    const std::string& mode,
    Rcpp::XPtr<somactx_wrap_t> soma_context,
    Rcpp::Nullable<Rcpp::DatetimeVector> tiledb_timestamp) {
    auto open_mode = get_open_mode(mode);
    auto timestamp_range = makeTimestampRange(tiledb_timestamp);
    return make_xptr<tiledbsoma::SOMADataFrame>(
        new tiledbsoma::SOMADataFrame(open_mode, uri, soma_context->ctxptr, timestamp_range));
}

// [[Rcpp::export]]
Rcpp::XPtr<tiledbsoma::SOMASparseNDArray> open_sparse_ndarray_handle(
    const std::string& uri,
    const std::string& mode,
    Rcpp::XPtr<somactx_wrap_t> soma_context,
    Rcpp::Nullable<Rcpp::DatetimeVector> tiledb_timestamp) {
    auto open_mode = get_open_mode(mode);
    auto timestamp_range = makeTimestampRange(tiledb_timestamp);
    return make_xptr<tiledbsoma::SOMASparseNDArray>(
        new tiledbsoma::SOMASparseNDArray(open_mode, uri, soma_context->ctxptr, timestamp_range));
}

// [[Rcpp::export]]
Rcpp::XPtr<tiledbsoma::SOMADenseNDArray> open_dense_ndarray_handle(
    const std::string& uri,
    const std::string& mode,
    Rcpp::XPtr<somactx_wrap_t> soma_context,
    Rcpp::Nullable<Rcpp::DatetimeVector> tiledb_timestamp) {
    auto open_mode = get_open_mode(mode);
    auto timestamp_range = makeTimestampRange(tiledb_timestamp);
    return make_xptr<tiledbsoma::SOMADenseNDArray>(
        new tiledbsoma::SOMADenseNDArray(open_mode, uri, soma_context->ctxptr, timestamp_range));
}

// [[Rcpp::export]]
Rcpp::XPtr<tiledbsoma::SOMACollection> open_collection_handle(
    const std::string& uri,
    const std::string& mode,
    Rcpp::XPtr<somactx_wrap_t> soma_context,
    Rcpp::Nullable<Rcpp::DatetimeVector> tiledb_timestamp) {
    auto open_mode = get_open_mode(mode);
    auto timestamp_range = makeTimestampRange(tiledb_timestamp);
    return make_xptr<tiledbsoma::SOMACollection>(
        new tiledbsoma::SOMACollection(open_mode, uri, soma_context->ctxptr, timestamp_range));
}

// [[Rcpp::export]]
Rcpp::XPtr<tiledbsoma::SOMAMeasurement> open_measurement_handle(
    const std::string& uri,
    const std::string& mode,
    Rcpp::XPtr<somactx_wrap_t> soma_context,
    Rcpp::Nullable<Rcpp::DatetimeVector> tiledb_timestamp) {
    auto open_mode = get_open_mode(mode);
    auto timestamp_range = makeTimestampRange(tiledb_timestamp);
    return make_xptr<tiledbsoma::SOMAMeasurement>(
        new tiledbsoma::SOMAMeasurement(open_mode, uri, soma_context->ctxptr, timestamp_range));
}

// [[Rcpp::export]]
Rcpp::XPtr<tiledbsoma::SOMAExperiment> open_experiment_handle(
    const std::string& uri,
    const std::string& mode,
    Rcpp::XPtr<somactx_wrap_t> soma_context,
    Rcpp::Nullable<Rcpp::DatetimeVector> tiledb_timestamp) {
    auto open_mode = get_open_mode(mode);
    auto timestamp_range = makeTimestampRange(tiledb_timestamp);
    return make_xptr<tiledbsoma::SOMAExperiment>(
        new tiledbsoma::SOMAExperiment(open_mode, uri, soma_context->ctxptr, timestamp_range));
}
