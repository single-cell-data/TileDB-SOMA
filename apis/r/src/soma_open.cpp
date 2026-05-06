#include <Rcpp/Lighter>

#include <optional>

#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"
#include "xptr-utils.h"  // xptr taggging utilities

// [[Rcpp::export]]
Rcpp::XPtr<somaobj_wrap_t> open_soma_handle(
    const std::string& uri,
    const std::string& mode,
    Rcpp::XPtr<somactx_wrap_t> soma_context,
    Rcpp::Nullable<Rcpp::DatetimeVector> tiledb_timestamp) {
    auto open_mode = get_open_mode(mode);
    auto timestamp_range = makeTimestampRange(tiledb_timestamp);

    std::shared_ptr<tiledbsoma::SOMAObject> handle = tiledbsoma::SOMAObject::open(
        uri, open_mode, soma_context->ctxptr, timestamp_range);
    return make_xptr<somaobj_wrap_t>(new SOMAWrapper(handle));
}

// [[Rcpp::export]]
Rcpp::XPtr<somaobj_wrap_t> open_dataframe_handle(
    const std::string& uri,
    const std::string& mode,
    Rcpp::XPtr<somactx_wrap_t> soma_context,
    Rcpp::Nullable<Rcpp::DatetimeVector> tiledb_timestamp) {
    auto open_mode = get_open_mode(mode);
    auto timestamp_range = makeTimestampRange(tiledb_timestamp);

    std::shared_ptr<tiledbsoma::SOMAObject> handle = std::make_shared<tiledbsoma::SOMADataFrame>(
        open_mode, uri, soma_context->ctxptr, timestamp_range);
    return make_xptr<somaobj_wrap_t>(new SOMAWrapper(handle));
}

// [[Rcpp::export]]
Rcpp::XPtr<somaobj_wrap_t> open_sparse_ndarray_handle(
    const std::string& uri,
    const std::string& mode,
    Rcpp::XPtr<somactx_wrap_t> soma_context,
    Rcpp::Nullable<Rcpp::DatetimeVector> tiledb_timestamp) {
    auto open_mode = get_open_mode(mode);
    auto timestamp_range = makeTimestampRange(tiledb_timestamp);

    std::shared_ptr<tiledbsoma::SOMAObject> handle = std::make_shared<tiledbsoma::SOMASparseNDArray>(
        open_mode, uri, soma_context->ctxptr, timestamp_range);
    return make_xptr<somaobj_wrap_t>(new SOMAWrapper(handle));
}

// [[Rcpp::export]]
Rcpp::XPtr<somaobj_wrap_t> open_dense_ndarray_handle(
    const std::string& uri,
    const std::string& mode,
    Rcpp::XPtr<somactx_wrap_t> soma_context,
    Rcpp::Nullable<Rcpp::DatetimeVector> tiledb_timestamp) {
    auto open_mode = get_open_mode(mode);
    auto timestamp_range = makeTimestampRange(tiledb_timestamp);

    std::shared_ptr<tiledbsoma::SOMAObject> handle = std::make_shared<tiledbsoma::SOMADenseNDArray>(
        open_mode, uri, soma_context->ctxptr, timestamp_range);
    return make_xptr<somaobj_wrap_t>(new SOMAWrapper(handle));
}

// [[Rcpp::export]]
Rcpp::XPtr<somaobj_wrap_t> open_collection_handle(
    const std::string& uri,
    const std::string& mode,
    Rcpp::XPtr<somactx_wrap_t> soma_context,
    Rcpp::Nullable<Rcpp::DatetimeVector> tiledb_timestamp) {
    auto open_mode = get_open_mode(mode);
    auto timestamp_range = makeTimestampRange(tiledb_timestamp);

    std::shared_ptr<tiledbsoma::SOMAObject> handle = std::make_shared<tiledbsoma::SOMACollection>(
        open_mode, uri, soma_context->ctxptr, timestamp_range);
    return make_xptr<somaobj_wrap_t>(new SOMAWrapper(handle));
}

// [[Rcpp::export]]
Rcpp::XPtr<somaobj_wrap_t> open_measurement_handle(
    const std::string& uri,
    const std::string& mode,
    Rcpp::XPtr<somactx_wrap_t> soma_context,
    Rcpp::Nullable<Rcpp::DatetimeVector> tiledb_timestamp) {
    auto open_mode = get_open_mode(mode);
    auto timestamp_range = makeTimestampRange(tiledb_timestamp);

    std::shared_ptr<tiledbsoma::SOMAObject> handle = std::make_shared<tiledbsoma::SOMAMeasurement>(
        open_mode, uri, soma_context->ctxptr, timestamp_range);
    return make_xptr<somaobj_wrap_t>(new SOMAWrapper(handle));
}

// [[Rcpp::export]]
Rcpp::XPtr<somaobj_wrap_t> open_experiment_handle(
    const std::string& uri,
    const std::string& mode,
    Rcpp::XPtr<somactx_wrap_t> soma_context,
    Rcpp::Nullable<Rcpp::DatetimeVector> tiledb_timestamp) {
    auto open_mode = get_open_mode(mode);
    auto timestamp_range = makeTimestampRange(tiledb_timestamp);

    std::shared_ptr<tiledbsoma::SOMAObject> handle = std::make_shared<tiledbsoma::SOMAExperiment>(
        open_mode, uri, soma_context->ctxptr, timestamp_range);
    return make_xptr<somaobj_wrap_t>(new SOMAWrapper(handle));
}
