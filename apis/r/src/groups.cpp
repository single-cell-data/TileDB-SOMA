#include <Rcpp/Lighter>  // for R interface to C++
#include <sstream>

#include <nanoarrow/r.h>            // for C/C++ interface to Arrow (via header exported from the R package)
#include <RcppInt64>                // for fromInteger64
#include <nanoarrow/nanoarrow.hpp>  // for C/C++ interface to Arrow (vendored)

#include <tiledbsoma/reindexer/reindexer.h>
#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"  // local declarations
#include "xptr-utils.h"  // xptr taggging utilities

namespace tdbs = tiledbsoma;

// [[Rcpp::export]]
Rcpp::List soma_group_get_members(Rcpp::XPtr<tiledbsoma::SOMAGroup> group) {
    auto mm = group->members_map();
    auto n = mm.size();
    Rcpp::List lst(n);
    Rcpp::CharacterVector names(n);
    int i = 0;
    for (auto const& key : mm) {
        names[i] = key.first;
        Rcpp::List row = Rcpp::List::create(
            Rcpp::Named("type") = std::string(key.second.second),
            Rcpp::Named("uri") = std::string(key.second.first),
            Rcpp::Named("name") = key.first);
        lst[i] = row;
        i++;
    }
    lst.attr("names") = names;
    return lst;
}

std::map<int, URIType> uritypemap = {{0, URIType::automatic}, {1, URIType::absolute}, {2, URIType::relative}};

// [[Rcpp::export]]
void soma_group_set(
    Rcpp::XPtr<tiledbsoma::SOMAGroup> group,
    const std::string& uri,
    int uri_type_int,  // "automatic", "absolute", "relative"
    const std::string& name,
    const std::string& soma_type) {
    if (!group) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    group->set(uri, uritypemap[uri_type_int], name, soma_type);
}

// [[Rcpp::export]]
void soma_group_remove_member(Rcpp::XPtr<tiledbsoma::SOMAGroup> group, const std::string& name) {
    if (!group) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    group->del(name);
}

// [[Rcpp::export]]
void soma_group_put_metadata(Rcpp::XPtr<tiledbsoma::SOMAGroup> group, std::string key, SEXP obj) {
    if (!group) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    // we implement a simpler interface here as the 'type' is given from the
    // supplied SEXP, as is the extent
    switch (TYPEOF(obj)) {
        case VECSXP: {
            Rcpp::stop("List objects are not supported.");
            break;  // not reached
        }
        case REALSXP: {
            Rcpp::NumericVector v(obj);
            if (Rcpp::isInteger64(obj)) {
                group->set_metadata(key, tiledbsoma::common::DataType::int64, v.size(), v.begin());
            } else {
                group->set_metadata(key, tiledbsoma::common::DataType::float64, v.size(), v.begin());
            }
            break;
        }
        case INTSXP: {
            Rcpp::IntegerVector v(obj);
            group->set_metadata(key, tiledbsoma::common::DataType::int32, v.size(), v.begin());
            break;
        }
        case STRSXP: {
            Rcpp::CharacterVector v(obj);
            std::string s(v[0]);
            // We use TILEDB_CHAR interchangeably with TILEDB_STRING_ASCII is
            // this best string type?
            // Use TILEDB_STRING_UTF8 for compatibility with Python API
            group->set_metadata(key, tiledbsoma::common::DataType::string_utf8, s.length(), s.c_str());
            break;
        }
        case LGLSXP: {  // experimental: map R logical (ie TRUE, FALSE, NA) to
                        // int8
            Rcpp::stop("Logical vector objects are not supported.");
            break;  // not reached
        }
        default: {
            Rcpp::stop("No support (yet) for type '%s'.", Rf_type2char(TYPEOF(obj)));
            break;  // not reached
        }
    }
}
