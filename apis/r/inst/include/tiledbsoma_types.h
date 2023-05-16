
// this file can make headers available for the generated file RcppExports.cpp

// defining this prevents spdlog to use stderr -- see bottom of spdlog/logger-inl.h
#define USING_R
#define R_R_H
// it also needs these R headers to define REprintf and ::R_FlushConsole
#include <R.h>
#include <Rinterface.h>
#include <R_ext/Print.h>

#include <nanoarrow.h>          			// for C interface to Arrow
#include <tiledb/tiledb>					// for QueryCondition etc
#define ARROW_SCHEMA_AND_ARRAY_DEFINED 1
#include <tiledbsoma/tiledbsoma>

namespace tdbs = tiledbsoma;
