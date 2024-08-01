#include <Rcpp/Lighter>                         // for R interface to C++
#include <nanoarrow/r.h>                        // for C/C++ interface to Arrow (via header exported from the R package)
#include <nanoarrow/nanoarrow.hpp>              // for C/C++ interface to Arrow (vendored)
#include <RcppInt64>                            // for fromInteger64

#include <tiledbsoma/tiledbsoma>
#include <tiledbsoma/reindexer/reindexer.h>

#include "rutilities.h"         				// local declarations
#include "xptr-utils.h"         				// xptr taggging utilities

namespace tdbs = tiledbsoma;

std::unique_ptr<tdbs::SOMAObject> getObjectUniquePointer(bool is_array, OpenMode mode,
                                                         std::string& uri,
                                                         std::shared_ptr<tdbs::SOMAContext> ctx) {
    if (is_array) {
        return tdbs::SOMAArray::open(mode, uri, ctx);
    } else {
        return tdbs::SOMAGroup::open(mode, uri, ctx);
    }
}

//' Get nnumber of metadata items
//' @param uri The array URI
//' @param ctxxp An external pointer to the SOMAContext wrapper
//' @export
// [[Rcpp::export]]
int32_t get_metadata_num(std::string& uri, bool is_array, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // SOMA Object unique pointer (aka soup)
    auto soup = getObjectUniquePointer(is_array, OpenMode::read, uri, sctx);
    int32_t nb = soup->metadata_num();
    return nb;
}

//' Read all metadata (as named list)
//'
//' This function currently supports metadata as either a string or an 'int64' (or 'int32').
//' It will error if a different datatype is encountered.
//' @param uri The array URI
//' @param ctxxp An external pointer to the SOMAContext wrapper
//' @export
// [[Rcpp::export]]
Rcpp::List get_all_metadata(std::string& uri, bool is_array, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;

    // SOMA Object unique pointer (aka soup)
    auto soup = getObjectUniquePointer(is_array, OpenMode::read, uri, sctx);
    auto mvmap = soup->get_metadata();

    std::vector<std::string> namvec;
    Rcpp::List lst;
    for (auto it = mvmap.begin(); it != mvmap.end(); it++) {
        std::string key = it->first;
        namvec.push_back(key);
        tdbs::MetadataValue val = it->second;
        auto dtype = std::get<0>(val);
        auto len = std::get<1>(val);
        const void* ptr = std::get<2>(val);
        if (dtype == TILEDB_STRING_UTF8 || dtype == TILEDB_STRING_ASCII) {
            auto str = std::string((char*) ptr, len);
            lst.push_back(str);
        } else if (dtype == TILEDB_INT64) {
            std::vector<int64_t> v(len);
            std::memcpy(&(v[0]), ptr, len*sizeof(int64_t));
            lst.push_back(Rcpp::toInteger64(v));
        } else if (dtype == TILEDB_INT32) {
            Rcpp::IntegerVector v(len);
            std::memcpy(v.begin(), ptr, len*sizeof(int32_t));
        } else {
            auto txt = tiledb::impl::type_to_str(dtype);
            Rcpp::stop("Currently unsupported type '%s'", txt.c_str());
        }
    }
    lst.attr("names") = Rcpp::CharacterVector(namvec.begin(), namvec.end());
    return lst;
}

//' Read metadata (as a string)
//'
//' @param uri The array URI
//' @param key The array metadata key
//' @param ctxxp An external pointer to the SOMAContext wrapper
//' @export
// [[Rcpp::export]]
std::string get_metadata(std::string& uri, std::string& key, bool is_array,
                         Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;

    // SOMA Object unique pointer (aka soup)
    auto soup = getObjectUniquePointer(is_array, OpenMode::read, uri, sctx);
    auto mv = soup->get_metadata(key);
    if (!mv.has_value()) {
        Rcpp::stop("No value for '%s'", key.c_str());
    }
    tdbs::MetadataValue val = *mv;
    auto dtype = std::get<0>(val);
    auto txt = tiledb::impl::type_to_str(dtype);
    if (txt != "STRING_UTF8" && txt != "STRING_ASCII") {
        Rcpp::stop("Currently unsupported type '%s'", txt.c_str());
    }
    auto len = std::get<1>(val);
    const void* ptr = std::get<2>(val);
    auto str = std::string((char*) ptr, len);
    return str;
}


//' Check for metadata given key
//'
//' @param uri The array URI
//' @param key The array metadata key
//' @param ctxxp An external pointer to the SOMAContext wrapper
//' @export
// [[Rcpp::export]]
bool has_metadata(std::string& uri, std::string& key, bool is_array,
                  Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // SOMA Object unique pointer (aka soup)
    auto soup = getObjectUniquePointer(is_array, OpenMode::read, uri, sctx);
    return soup->has_metadata(key);
}


//' Delete metadata for given key
//'
//' @param uri The array URI
//' @param key The array metadata key
//' @param ctxxp An external pointer to the SOMAContext wrapper
//' @export
// [[Rcpp::export]]
void delete_metadata(std::string& uri, std::string& key, bool is_array,
                     Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // SOMA Object unique pointer (aka soup)
    auto soup = getObjectUniquePointer(is_array, OpenMode::write, uri, sctx);
    soup->delete_metadata(key);
    soup->close();
}


//' Set metadata (as a string)
//'
//' @param uri The array URI
//' @param key The array metadata key
//' @param value The metadata value
//' @
//' @param ctxxp An external pointer to the SOMAContext wrapper
//' @export
// [[Rcpp::export]]
void set_metadata(std::string& uri, std::string& key, SEXP valuesxp, std::string& type,
                  bool is_array, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // SOMA Object unique pointer (aka soup)
    auto soup = getObjectUniquePointer(is_array, OpenMode::write, uri, sctx);

    if (type == "character") {
        const tiledb_datatype_t value_type = TILEDB_STRING_UTF8;
        std::string value = Rcpp::as<std::string>(valuesxp);
        spdl::debug("[set_metadata] key {} value {} is_array {} type {}", key, value, is_array, type);
        soup->set_metadata(key, value_type, value.length(), (void*) value.c_str(), true);
    } else if (type == "integer64") {
        const tiledb_datatype_t value_type = TILEDB_INT64;
        double dv = Rcpp::as<double>(valuesxp);
        int64_t value = Rcpp::fromInteger64(dv);
        spdl::debug("[set_metadata] key {} value {} is_array {} type {}", key, value, is_array, type);
        soup->set_metadata(key, value_type, 1, (void*) &value, true);
    } else {
        Rcpp::stop("Unsupported type '%s'", type);
    }
    soup->close();
}
