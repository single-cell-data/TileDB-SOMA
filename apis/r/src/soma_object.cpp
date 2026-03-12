#include <Rcpp/Lighter>
#include <RcppInt64>

#include <optional>

#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"
#include "xptr-utils.h"  // xptr taggging utilities

// [[Rcpp::export]]
Rcpp::List soma_object_get_metadata(Rcpp::XPtr<tiledbsoma::SOMAObject> soma_object) {
    if (!soma_object) {
        throw Rcpp::exception("Internal error: SOMA object handle is not initialized.");
    }
    auto metadata = soma_object->get_metadata();
    return metadata_as_rlist(metadata);
}

// [[Rcpp::export]]
void soma_object_set_metadata(Rcpp::XPtr<tiledbsoma::SOMAObject> soma_object, const std::string& key, SEXP value) {
    if (!soma_object) {
        throw Rcpp::exception("Internal error: SOMA object handle is not initialized.");
    }
    switch (TYPEOF(value)) {
        case VECSXP: {
            Rcpp::stop("Writing List values to metadata is not supported.");
            break;  // not reached
        }
        case REALSXP: {
            Rcpp::NumericVector v(value);
            if (Rcpp::isInteger64(value)) {
                soma_object->set_metadata(key, tiledbsoma::common::DataType::int64, v.size(), v.begin());
            } else {
                soma_object->set_metadata(key, tiledbsoma::common::DataType::float64, v.size(), v.begin());
            }
            break;
        }
        case INTSXP: {
            Rcpp::IntegerVector v(value);
            soma_object->set_metadata(key, tiledbsoma::common::DataType::int32, v.size(), v.begin());
            break;
        }
        case STRSXP: {
            Rcpp::CharacterVector v(value);
            std::string s(v[0]);
            soma_object->set_metadata(key, tiledbsoma::common::DataType::string_utf8, s.length(), s.c_str());
            break;
        }
        case LGLSXP: {  // experimental: map R logical (ie TRUE, FALSE, NA) toint8
            Rcpp::stop("Writing logical vectors to metadata is not supported.");
            break;
        }
        default: {
            Rcpp::stop("Support for writing '%s' to metadata is not yet implemented.", Rf_type2char(TYPEOF(value)));
            break;  // not reached
        }
    }
}

// [[Rcpp::export]]
void soma_object_close(Rcpp::XPtr<tiledbsoma::SOMAObject> soma_object) {
    if (soma_object) {
        soma_object->close();
    }
}

// [[Rcpp::export]]
bool soma_object_is_open(Rcpp::XPtr<tiledbsoma::SOMAObject> soma_object) {
    if (soma_object) {
        return soma_object->is_open();
    }
    return false;
}

// [[Rcpp::export]]
std::string soma_object_open_mode(Rcpp::XPtr<tiledbsoma::SOMAObject> soma_object) {
    if (!soma_object || !soma_object->is_open()) {
        return "CLOSED";
    }
    return get_open_mode_string(soma_object->mode());
}
