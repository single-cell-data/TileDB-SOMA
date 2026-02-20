#include <Rcpp/Lighter>

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
