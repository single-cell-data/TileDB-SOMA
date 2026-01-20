#include <Rcpp.h>                   // for R interface to C++
#include <nanoarrow/r.h>            // for C interface to Arrow (via R package)
#include <RcppInt64>                // for fromInteger64
#include <nanoarrow/nanoarrow.hpp>  // for C/C++ interface to Arrow
#include <sstream>

// we currently get deprecation warnings by default which are noisy
#ifndef TILEDB_NO_API_DEPRECATION_WARNINGS
#define TILEDB_NO_API_DEPRECATION_WARNINGS
#endif

// We get these via nanoarrow and must cannot include carrow.h again
#define ARROW_SCHEMA_AND_ARRAY_DEFINED 1
#include <tiledb/tiledb_experimental>
#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"  // local declarations
#include "xptr-utils.h"  // xptr taggging utilities

// Helper
tiledb_query_condition_combination_op_t _tiledb_query_string_to_condition_combination_op(const std::string& opstr) {
    if (opstr == "AND") {
        return TILEDB_AND;
    } else if (opstr == "OR") {
        return TILEDB_OR;
    } else if (opstr == "NOT") {
        return TILEDB_NOT;
    } else {
        Rcpp::stop("Unknown TileDB combination op string '%s'", opstr.c_str());
    }
}

// Helper
tiledb_query_condition_op_t _op_name_to_tdb_op(const std::string& opstr) {
    if (opstr == "LT") {
        return TILEDB_LT;
    } else if (opstr == "LE") {
        return TILEDB_LE;
    } else if (opstr == "GT") {
        return TILEDB_GT;
    } else if (opstr == "GE") {
        return TILEDB_GE;
    } else if (opstr == "EQ") {
        return TILEDB_EQ;
    } else if (opstr == "NE") {
        return TILEDB_NE;
    } else if (opstr == "IN") {
        return TILEDB_IN;
    } else if (opstr == "NOT_IN") {
        return TILEDB_NOT_IN;
    } else {
        Rcpp::stop("Unknown TileDB op string '%s'", opstr.c_str());
    }
}

// [[Rcpp::export]]
Rcpp::XPtr<tdbs::QueryCondition> libtiledbsoma_empty_query_condition(Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // Shared pointer to SOMAContext from external pointer wrapper:
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // Shared pointer to TileDB Context from SOMAContext:
    std::shared_ptr<tiledb::Context> ctx = sctx->tiledb_ctx();
    // Core constructor
    return make_xptr<tdbs::QueryCondition>(new tdbs::QueryCondition(*ctx.get()));
}

// [[Rcpp::export]]
void libtiledbsoma_query_condition_from_triple(
    Rcpp::XPtr<tdbs::QueryCondition> query_cond,
    const std::string& attr_name,
    SEXP condition_value,
    const std::string& arrow_type_name,
    const std::string& cond_op_string) {
    // No such:
    // print(arrow::large_string()$name)
    // print(arrow::double()$name)

    // print(arrow::int64()$name)         [1] "int64"
    // print(arrow::uint64()$name)        [1] "uint64"
    // print(arrow::int32()$name)         [1] "int32"
    // print(arrow::uint32()$name)        [1] "uint32"
    // print(arrow::int16()$name)         [1] "int16"
    // print(arrow::uint16()$name)        [1] "uint16"
    // print(arrow::int8()$name)          [1] "int8"
    // print(arrow::uint8()$name)         [1] "uint8"
    // print(arrow::float64()$name)       [1] "double"
    // print(arrow::float()$name)         [1] "float"
    // print(arrow::float32()$name)       [1] "float"
    // print(arrow::string()$name)        [1] "utf8"
    // print(arrow::binary()$name)        [1] "binary"
    // print(arrow::large_binary()$name)  [1] "large_binary"
    // print(arrow::bool()$name)          [1] "bool"
    // print(arrow::boolean()$name)       [1] "bool"
    // print(arrow::date64()$name)        [1] "date64"
    // print(arrow::date32()$name)        [1] "date32"
    // print(arrow::time32()$name)        [1] "time32"
    // print(arrow::time64()$name)        [1] "time64"

    check_xptr_tag<tdbs::QueryCondition>(query_cond);
    tiledb_query_condition_op_t op = _op_name_to_tdb_op(cond_op_string);

    if (arrow_type_name == "int64" || arrow_type_name == "uint64") {
        int64_t v = Rcpp::fromInteger64(Rcpp::as<double>(condition_value));
        uint64_t cond_val_size = sizeof(int64_t);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);

    } else if (arrow_type_name == "int32" || arrow_type_name == "uint32") {
        int v = Rcpp::as<int>(condition_value);
        uint64_t cond_val_size = sizeof(int);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);

    } else if (arrow_type_name == "int16" || arrow_type_name == "uint16") {
        int v = Rcpp::as<int>(condition_value);
        uint64_t cond_val_size = sizeof(int16_t);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);

    } else if (arrow_type_name == "int8" || arrow_type_name == "uint8") {
        int v = Rcpp::as<int>(condition_value);
        uint64_t cond_val_size = sizeof(int8_t);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);

    } else if (arrow_type_name == "double") {
        double v = Rcpp::as<double>(condition_value);
        uint64_t cond_val_size = sizeof(double);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);

    } else if (arrow_type_name == "float") {
        float v = static_cast<float>(Rcpp::as<double>(condition_value));
        uint64_t cond_val_size = sizeof(float);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);

    } else if (
        arrow_type_name == "string" || arrow_type_name == "ascii" || arrow_type_name == "utf8" ||
        arrow_type_name == "large_utf8") {
        std::string v = Rcpp::as<std::string>(condition_value);
        query_cond->init(attr_name, v, op);

    } else if (arrow_type_name == "bool") {
        bool v = Rcpp::as<bool>(condition_value);
        uint64_t cond_val_size = sizeof(bool);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);

    } else if (arrow_type_name == "timestamp_s") {
        int64_t v = static_cast<int64_t>(Rcpp::as<double>(condition_value));
        std::stringstream ss;
        ss << "ts3 " << v;
        tdbs::common::logging::LOG_DEBUG(ss.str());
        uint64_t cond_val_size = sizeof(int64_t);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);

    } else if (arrow_type_name == "timestamp_ms") {
        int64_t v = static_cast<int64_t>(Rcpp::as<double>(condition_value) * 1000);
        uint64_t cond_val_size = sizeof(int64_t);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);

    } else if (arrow_type_name == "timestamp_us") {
        int64_t v = static_cast<int64_t>(Rcpp::as<double>(condition_value) * 1e6);
        uint64_t cond_val_size = sizeof(int64_t);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);

    } else if (arrow_type_name == "timestamp_ns") {
        // XXX nanotime ...
        int64_t v = static_cast<int64_t>(Rcpp::as<double>(condition_value) * 1e9);
        uint64_t cond_val_size = sizeof(int64_t);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);

    } else if (arrow_type_name == "date32") {
        // Arrow date32 TileDB DATETIME_DAY
        int64_t v = static_cast<int64_t>(Rcpp::as<double>(condition_value));
        uint64_t cond_val_size = sizeof(int64_t);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);

    } else {
        Rcpp::stop("tiledbsoma query condition: currently unsupported type \"%s\"", arrow_type_name);
    }
}

// [[Rcpp::export]]
Rcpp::XPtr<tdbs::QueryCondition> libtiledbsoma_query_condition_combine(
    Rcpp::XPtr<tdbs::QueryCondition> lhs, Rcpp::XPtr<tdbs::QueryCondition> rhs, const std::string& str) {
    check_xptr_tag<tdbs::QueryCondition>(lhs);
    check_xptr_tag<tdbs::QueryCondition>(lhs);
    tiledb_query_condition_combination_op_t op = _tiledb_query_string_to_condition_combination_op(str);
    tdbs::QueryCondition res = lhs->combine(*rhs.get(), op);
    return make_xptr<tdbs::QueryCondition>(new tdbs::QueryCondition(res));
}

// [[Rcpp::export]]
Rcpp::XPtr<tdbs::QueryCondition> libtiledbsoma_query_condition_in_nin(
    Rcpp::XPtr<somactx_wrap_t> ctxxp, const std::string& attr_name, const std::string& op_name, SEXP values) {
    // Shared pointer to SOMAContext from external pointer wrapper:
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // Shared pointer to TileDB Context from SOMAContext:
    std::shared_ptr<tiledb::Context> ctx = sctx->tiledb_ctx();

    tiledb_query_condition_op_t op = _op_name_to_tdb_op(op_name);

    if (TYPEOF(values) == INTSXP) {
        std::vector<int32_t> iv = Rcpp::as<std::vector<int32_t>>(values);
        auto qc = tdbs::QueryConditionExperimental::create<int32_t>(*ctx.get(), attr_name, iv, op);
        return make_xptr<tdbs::QueryCondition>(new tdbs::QueryCondition(qc));

    } else if (TYPEOF(values) == REALSXP) {
        if (Rcpp::isInteger64(values)) {
            std::vector<int64_t> dv = Rcpp::fromInteger64(Rcpp::NumericVector(values));
            auto qc = tdbs::QueryConditionExperimental::create<int64_t>(*ctx.get(), attr_name, dv, op);
            return make_xptr<tdbs::QueryCondition>(new tdbs::QueryCondition(qc));
        } else {
            std::vector<double> dv = Rcpp::as<std::vector<double>>(values);
            auto qc = tdbs::QueryConditionExperimental::create<double>(*ctx.get(), attr_name, dv, op);
            return make_xptr<tdbs::QueryCondition>(new tdbs::QueryCondition(qc));
        }

    } else if (TYPEOF(values) == STRSXP) {
        std::vector<std::string> sv = Rcpp::as<std::vector<std::string>>(values);
        auto qc = tdbs::QueryConditionExperimental::create(*ctx.get(), attr_name, sv, op);
        return make_xptr<tdbs::QueryCondition>(new tdbs::QueryCondition(qc));

    } else {
        Rcpp::stop("No support (yet) for type '%s'.", Rcpp::type2name(values));
    }
}
