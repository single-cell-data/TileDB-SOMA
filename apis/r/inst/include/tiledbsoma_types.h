#pragma once

#include <Rcpp.h> 

// defining this prevents spdlog to use stderr -- see bottom of spdlog/logger-inl.h
#define USING_R
#define R_R_H

// it also needs these R headers to define REprintf and ::R_FlushConsole
#include <R.h>
#include <Rinterface.h>
#include <R_ext/Print.h>

// we currently get deprecation warnings by default which are noisy
// this turns them off for RcppExports.cpp
#ifndef TILEDB_NO_API_DEPRECATION_WARNINGS
#define TILEDB_NO_API_DEPRECATION_WARNINGS
#endif

#include <nanoarrow/nanoarrow.h>            // for C interface to Arrow
#include <tiledb/tiledb>					// for QueryCondition etc
#define ARROW_SCHEMA_AND_ARRAY_DEFINED 1
#include <tiledbsoma/tiledbsoma>
#include <tiledbsoma/reindexer/reindexer.h>
#include "rutilities.h"

#include <concepts>

namespace tdbs = tiledbsoma;

// We create a struct containing a shared pointer to a Context, this will have standard C++
// semantics on its reference count. To not have the R external pointer reference count get
// in the way of the Context object, the external pointer we use to interface is on the struct
struct ContextWrapper {
    //ContextWrapper(std::shared_ptr<tiledb::Context> ctx_ptr_) : ctxptr(std::move(ctx_ptr_)) {}
    ContextWrapper(std::shared_ptr<tiledb::Context> ctx_ptr_) : ctxptr(ctx_ptr_) {}
    std::shared_ptr<tiledb::Context> ctxptr;
};
typedef struct ContextWrapper ctx_wrap_t;

struct SOMAContextWrapper {
    SOMAContextWrapper(std::shared_ptr<tdbs::SOMAContext> ctx_ptr_) : ctxptr(ctx_ptr_) {}
    std::shared_ptr<tdbs::SOMAContext> ctxptr;
};
typedef struct SOMAContextWrapper somactx_wrap_t;

struct SOMAWrapper {
    SOMAWrapper(std::shared_ptr<tdbs::SOMAObject> objptr) : objptr_(objptr) {}

    template <class T = tdbs::SOMAObject>
    requires std::derived_from<T, tdbs::SOMAObject>
    std::shared_ptr<T> ptr() {
        auto result = std::dynamic_pointer_cast<T>(objptr_);

        if (objptr_.get() == nullptr) {
            throw std::runtime_error("Unexpected SOMA handle null pointer");
        }

        if (result == nullptr) {
            throw std::runtime_error("SOMAObject cast failed from " + tdbs::common::demangle_name(typeid(objptr_.get()).name()) + " to " + tdbs::common::demangle_name(typeid(T*).name()));
        }

        return result;
    }
   private:
    std::shared_ptr<tdbs::SOMAObject> objptr_;
};
typedef struct SOMAWrapper somaobj_wrap_t;


// make the function signature nicer as using an uppercase SEXP 'screams'
// we can not tag these as we do in xptrUtils.h they pass through to the
// nanoarrow helpers and are expected to be plain EXTPTR SEXP types.
typedef SEXP naxpArray;
typedef SEXP naxpSchema;
