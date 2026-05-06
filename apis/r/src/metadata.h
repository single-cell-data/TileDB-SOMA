#pragma once

#include <Rcpp/Lighter>  // for R interface to C++
#include <sstream>

#include <nanoarrow/r.h>            // for C/C++ interface to Arrow (via header exported from the R package)
#include <RcppInt64>                // for fromInteger64
#include <nanoarrow/nanoarrow.hpp>  // for C/C++ interface to Arrow (vendored)

#include <tiledbsoma/reindexer/reindexer.h>
#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"  // local declarations
#include "xptr-utils.h"  // xptr taggging utilities

void soma_delete_metadata(Rcpp::XPtr<somaobj_wrap_t> objxp, const std::string& key);

bool soma_has_metadata(Rcpp::XPtr<somaobj_wrap_t> objxp, const std::string& key);

Rcpp::List soma_get_all_metadata(Rcpp::XPtr<somaobj_wrap_t> objxp);

SEXP soma_get_metadata(Rcpp::XPtr<somaobj_wrap_t> objxp, const std::string& key);

void soma_set_metadata(Rcpp::XPtr<somaobj_wrap_t> objxp, const std::string& key, SEXP value);

// helper function to convert_metadata
SEXP _metadata_to_sexp(const tdbs::common::MetadataValue& value);