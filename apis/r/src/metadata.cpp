#include <Rcpp/Lighter>  // for R interface to C++
#include <sstream>

#include <nanoarrow/r.h>            // for C/C++ interface to Arrow (via header exported from the R package)
#include <RcppInt64>                // for fromInteger64
#include <nanoarrow/nanoarrow.hpp>  // for C/C++ interface to Arrow (vendored)

#include <tiledbsoma/reindexer/reindexer.h>
#include <tiledbsoma/tiledbsoma>

#include "metadata.h"
#include "rutilities.h"  // local declarations
#include "xptr-utils.h"  // xptr taggging utilities

namespace tdbs = tiledbsoma;

std::unique_ptr<tdbs::SOMAObject> getObjectUniquePointer(
    bool is_array,
    OpenMode mode,
    std::string& uri,
    std::shared_ptr<tdbs::SOMAContext> ctx,
    Rcpp::Nullable<Rcpp::DatetimeVector> tsvec = R_NilValue) {
    // optional timestamp range
    std::optional<tdbs::TimestampRange> tsrng = makeTimestampRange(tsvec);

    if (is_array) {
        return tdbs::SOMAArray::open(mode, uri, ctx, tsrng);
    } else {
        return tdbs::SOMAGroup::open(mode, uri, ctx, "unnamed", tsrng);
    }
}

// Get number of metadata items
//
// @param uri The array URI
// @param is_array A boolean to indicate array or group
// @param ctxxp An external pointer to the SOMAContext wrapper
//
// @export
//
// [[Rcpp::export]]
int32_t get_metadata_num(std::string& uri, bool is_array, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // SOMA Object unique pointer (aka soup)
    auto soup = getObjectUniquePointer(is_array, OpenMode::soma_read, uri, sctx);
    int32_t nb = soup->metadata_num();
    return nb;
}

// Read all metadata (as named list)
//
// This function currently supports metadata as either a string or an 'int64'
// (or 'int32'). It will error if a different datatype is encountered.
//
// @param uri The array URI
// @param is_array A boolean to indicate array or group
// @param ctxxp An external pointer to the SOMAContext wrapper
//
// @export
//
// [[Rcpp::export]]
Rcpp::List get_all_metadata(std::string& uri, bool is_array, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;

    std::shared_ptr<tiledbsoma::SOMAObject> soup = getObjectUniquePointer(is_array, OpenMode::soma_read, uri, sctx);
    return soma_get_all_metadata(make_xptr<somaobj_wrap_t>(new SOMAWrapper(soup)));
}

// Read metadata (as a string)
//
// @param uri The array URI
// @param key The array metadata key
// @param is_array A boolean to indicate array or group
// @param ctxxp An external pointer to the SOMAContext wrapper
//
// @export
//
// [[Rcpp::export]]
std::string get_metadata(std::string& uri, std::string& key, bool is_array, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;

    std::shared_ptr<tiledbsoma::SOMAObject> soup = getObjectUniquePointer(is_array, OpenMode::soma_read, uri, sctx);
    return Rcpp::as<std::string>(soma_get_metadata(make_xptr<somaobj_wrap_t>(new SOMAWrapper(soup)), key));
}

// Check for metadata given key
//
// @param uri The array URI
// @param key The array metadata key
// @param is_array A boolean to indicate array or group
// @param ctxxp An external pointer to the SOMAContext wrapper
//
// @export
//
// [[Rcpp::export]]
bool has_metadata(std::string& uri, std::string& key, bool is_array, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;

    std::shared_ptr<tiledbsoma::SOMAObject> soup = getObjectUniquePointer(is_array, OpenMode::soma_read, uri, sctx);
    return soma_has_metadata(make_xptr<somaobj_wrap_t>(new SOMAWrapper(soup)), key);
}

// Delete metadata for given key
//
// @param uri The array URI
// @param key The array metadata key
// @param is_array A boolean to indicate array or group
// @param ctxxp An external pointer to the SOMAContext wrapper
//
// @export
//
// [[Rcpp::export]]
void delete_metadata(std::string& uri, std::string& key, bool is_array, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;

    std::shared_ptr<tiledbsoma::SOMAObject> soup = getObjectUniquePointer(is_array, OpenMode::soma_write, uri, sctx);
    return soma_delete_metadata(make_xptr<somaobj_wrap_t>(new SOMAWrapper(soup)), key);
}

// Set metadata (as a string)
//
// @param uri The array URI
// @param key The array metadata key
// @param valuesxp The metadata value
// @param type The datatype
// @param is_array A boolean to indicate array or group
// @param ctxxp An external pointer to the SOMAContext wrapper
// @param tsvec An optional two-element datetime vector
//
// @export
//
// [[Rcpp::export]]
void set_metadata(
    std::string& uri,
    std::string& key,
    SEXP valuesxp,
    std::string& type,
    bool is_array,
    Rcpp::XPtr<somactx_wrap_t> ctxxp,
    Rcpp::Nullable<Rcpp::DatetimeVector> tsvec = R_NilValue) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;

    std::shared_ptr<tiledbsoma::SOMAObject> soup = getObjectUniquePointer(
        is_array, OpenMode::soma_read, uri, sctx, tsvec);
    return soma_set_metadata(make_xptr<somaobj_wrap_t>(new SOMAWrapper(soup)), key, valuesxp);
}

// [[Rcpp::export]]
void soma_delete_metadata(Rcpp::XPtr<somaobj_wrap_t> objxp, const std::string& key) {
    objxp->ptr()->delete_metadata(key);
}

// [[Rcpp::export]]
bool soma_has_metadata(Rcpp::XPtr<somaobj_wrap_t> objxp, const std::string& key) {
    return objxp->ptr()->has_metadata(key);
}

// [[Rcpp::export]]
Rcpp::List soma_get_all_metadata(Rcpp::XPtr<somaobj_wrap_t> objxp) {
    auto metadata = objxp->ptr()->get_metadata();

    std::vector<std::string> namvec;
    Rcpp::List lst;

    for (const auto& [key, value] : metadata) {
        namvec.push_back(key);
        lst.push_back(_metadata_to_sexp(value));
    }

    lst.attr("names") = Rcpp::CharacterVector(namvec.begin(), namvec.end());
    return lst;
}

// [[Rcpp::export]]
SEXP soma_get_metadata(Rcpp::XPtr<somaobj_wrap_t> objxp, const std::string& key) {
    auto result = objxp->ptr()->get_metadata(key);
    if (!result.has_value()) {
        // Rcpp::stop("No value for '%s'", key.c_str());
        return R_NilValue;
    }

    return _metadata_to_sexp(result.value());
}

// [[Rcpp::export]]
void soma_set_metadata(Rcpp::XPtr<somaobj_wrap_t> objxp, const std::string& key, SEXP value) {
    switch (TYPEOF(value)) {
        case VECSXP: {
            Rcpp::stop("List objects are not supported.");
            break;  // not reached
        }
        case REALSXP: {
            Rcpp::NumericVector v(value);
            if (Rcpp::isInteger64(value)) {
                objxp->ptr()->set_metadata(key, tdbs::common::DataType::int64, v.size(), v.begin());
            } else {
                objxp->ptr()->set_metadata(key, tdbs::common::DataType::float64, v.size(), v.begin());
            }
            break;
        }
        case INTSXP: {
            Rcpp::IntegerVector v(value);
            objxp->ptr()->set_metadata(key, tdbs::common::DataType::int32, v.size(), v.begin());
            break;
        }
        case STRSXP: {
            Rcpp::CharacterVector v(value);
            std::string s(v[0]);
            objxp->ptr()->set_metadata(key, tdbs::common::DataType::string_utf8, s.length(), s.c_str());
            break;
        }
        case LGLSXP: {  // experimental: map R logical (ie TRUE, FALSE, NA) to
                        // int8
            Rcpp::stop("Logical vector objects are not supported.");
            break;  // not reached
        }
        default: {
            Rcpp::stop("No support (yet) for type '%s'.", Rf_type2char(TYPEOF(value)));
            break;  // not reached
        }
    }
}

// helper function to convert_metadata
SEXP _metadata_to_sexp(const tdbs::common::MetadataValue& value) {
    // This supports a limited set of basic types as the metadata
    // annotation is not meant to support complete serialization

    return std::visit(
        [&](auto&& arg) -> SEXP {
            using T = std::decay_t<decltype(arg)>;

            if constexpr (std::is_same_v<T, std::string>) {
                return Rcpp::wrap(arg);
            } else if constexpr (std::is_same_v<T, int8_t> || std::is_same_v<T, bool>) {
                return Rcpp::LogicalVector({{static_cast<bool>(arg)}});
            } else if constexpr (
                std::is_same_v<T, uint8_t> || std::is_same_v<T, uint16_t> || std::is_same_v<T, int16_t> ||
                std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>) {
                return copy_int_vector<T>({{arg}});
            } else if constexpr (std::is_same_v<T, int32_t>) {
                return Rcpp::IntegerVector({{arg}});
            } else if constexpr (std::is_same_v<T, int64_t>) {
                return Rcpp::toInteger64({{arg}});
            } else if constexpr (std::is_same_v<T, float>) {
                return Rcpp::NumericVector({{static_cast<double>(arg)}});
            } else if constexpr (std::is_same_v<T, double>) {
                return Rcpp::NumericVector({{arg}});
            } else if constexpr (std::is_same_v<T, std::vector<int8_t>> || std::is_same_v<T, std::vector<bool>>) {
                Rcpp::LogicalVector vec(arg.size());
                for (size_t i = 0; i < arg.size(); i++)
                    vec[i] = static_cast<bool>(arg[i]);
                return vec;
            } else if constexpr (
                std::is_same_v<T, std::vector<uint8_t>> || std::is_same_v<T, std::vector<uint16_t>> ||
                std::is_same_v<T, std::vector<int16_t>> || std::is_same_v<T, std::vector<uint32_t>> ||
                std::is_same_v<T, std::vector<uint64_t>>) {
                return copy_int_vector<typename T::value_type>(arg);
            } else if constexpr (std::is_same_v<T, std::vector<int32_t>>) {
                Rcpp::IntegerVector vec(arg.size());
                std::copy(arg.cbegin(), arg.cend(), vec.begin());
                return vec;
            } else if constexpr (std::is_same_v<T, std::vector<int64_t>>) {
                return Rcpp::toInteger64(arg);
            } else if constexpr (std::is_same_v<T, std::vector<float>>) {
                Rcpp::NumericVector vec(arg.size());
                for (size_t i = 0; i < arg.size(); i++)
                    vec[i] = static_cast<double>(arg[i]);
                return vec;
            } else if constexpr (std::is_same_v<T, std::vector<double>>) {
                Rcpp::NumericVector vec(arg.size());
                std::copy(arg.cbegin(), arg.cend(), vec.begin());
                return vec;
            } else {
                Rcpp::stop("No support yet for TileDB data type %s", tdbs::common::demangle_name(typeid(T).name()));
            }
        },
        value);
}