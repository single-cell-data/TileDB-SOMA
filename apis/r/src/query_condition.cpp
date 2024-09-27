#include <Rcpp.h>                   // for R interface to C++
#include <nanoarrow/r.h>            // for C interface to Arrow (via R package)
#include <RcppInt64>                // for fromInteger64
#include <nanoarrow/nanoarrow.hpp>  // for C/C++ interface to Arrow

// we currently get deprecation warnings by default which are noisy
#ifndef TILEDB_NO_API_DEPRECATION_WARNINGS
#define TILEDB_NO_API_DEPRECATION_WARNINGS
#endif

#include "rutilities.h"  // local declarations
#include "xptr-utils.h"  // xptr taggging utilities

// [[Rcpp::export]]
Rcpp::XPtr<tiledb::Query> libtiledbsoma_query_set_condition(
    Rcpp::XPtr<tiledb::Query> query, Rcpp::XPtr<tiledb::QueryCondition> query_cond) {
    check_xptr_tag<tiledb::Query>(query);
    query->set_condition(*query_cond.get());
    return query;
}

/**
 * Query Condition
 */
const char* _tiledb_query_condition_op_to_string(
    tiledb_query_condition_op_t op) {
    switch (op) {
        case TILEDB_LT:
            return "LT";
        case TILEDB_LE:
            return "LE";
        case TILEDB_GT:
            return "GT";
        case TILEDB_GE:
            return "GE";
        case TILEDB_EQ:
            return "EQ";
        case TILEDB_NE:
            return "NE";
        case TILEDB_IN:
            return "IN";
        case TILEDB_NOT_IN:
            return "NOT_IN";
        default:
            Rcpp::stop("Unknown condition op (%d)", op);
    }
}

tiledb_query_condition_op_t _tiledb_query_string_to_condition_op(
    const std::string& opstr) {
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

const char* _tiledb_query_condition_combination_op_to_string(
    tiledb_query_condition_combination_op_t op) {
    switch (op) {
        case TILEDB_AND:
            return "AND";
        case TILEDB_OR:
            return "OR";
        case TILEDB_NOT:
            return "NOT";
        default:
            Rcpp::stop("Unknown condition combination op (%d)", op);
    }
}

tiledb_query_condition_combination_op_t
_tiledb_query_string_to_condition_combination_op(const std::string& opstr) {
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

// [[Rcpp::export]]
Rcpp::XPtr<tiledb::QueryCondition> libtiledbsoma_query_condition(
    Rcpp::XPtr<tiledb::Context> ctx) {
    check_xptr_tag<tiledb::Context>(ctx);
    auto ptr = make_xptr<tiledb::QueryCondition>(
        new tiledb::QueryCondition(*ctx.get()));
    return ptr;
}

// [[Rcpp::export]]
void libtiledbsoma_query_condition_init(
    Rcpp::XPtr<tiledb::QueryCondition> query_cond,
    const std::string& attr_name,
    SEXP condition_value,
    const std::string& cond_val_type,
    const std::string& cond_op_string) {
    check_xptr_tag<tiledb::QueryCondition>(query_cond);
    tiledb_query_condition_op_t op = _tiledb_query_string_to_condition_op(
        cond_op_string);
    if (cond_val_type == "INT32" || cond_val_type == "UINT32") {
        int v = Rcpp::as<int>(condition_value);
        uint64_t cond_val_size = sizeof(int);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);
    } else if (cond_val_type == "FLOAT64") {
        double v = Rcpp::as<double>(condition_value);
        uint64_t cond_val_size = sizeof(double);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);
    } else if (cond_val_type == "INT64" || cond_val_type == "UINT64") {
        int64_t v = Rcpp::fromInteger64(Rcpp::as<double>(condition_value));
        uint64_t cond_val_size = sizeof(int64_t);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);
    } else if (cond_val_type == "INT8" || cond_val_type == "UINT8") {
        int v = Rcpp::as<int>(condition_value);
        uint64_t cond_val_size = sizeof(int8_t);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);
    } else if (cond_val_type == "INT16" || cond_val_type == "UINT16") {
        int v = Rcpp::as<int>(condition_value);
        uint64_t cond_val_size = sizeof(int16_t);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);
    } else if (cond_val_type == "FLOAT32") {
        float v = static_cast<float>(Rcpp::as<double>(condition_value));
        uint64_t cond_val_size = sizeof(float);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);
    } else if (cond_val_type == "ASCII" || cond_val_type == "UTF8") {
        std::string v = Rcpp::as<std::string>(condition_value);
        query_cond->init(attr_name, v, op);
    } else if (cond_val_type == "BOOL") {
        bool v = Rcpp::as<bool>(condition_value);
        uint64_t cond_val_size = sizeof(bool);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);
    } else if (cond_val_type == "DATETIME_MS") {
        int64_t v = static_cast<int64_t>(Rcpp::as<double>(condition_value) * 1000);
        uint64_t cond_val_size = sizeof(int64_t);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);
    } else if (cond_val_type == "DATETIME_DAY") {
        int64_t v = static_cast<int64_t>(Rcpp::as<double>(condition_value));
        uint64_t cond_val_size = sizeof(int64_t);
        query_cond->init(attr_name, (void*)&v, cond_val_size, op);
    } else {
        Rcpp::stop("Currently unsupported type: %s", cond_val_type);
    }
}

// [[Rcpp::export]]
Rcpp::XPtr<tiledb::QueryCondition> libtiledbsoma_query_condition_combine(
    Rcpp::XPtr<tiledb::QueryCondition> lhs,
    Rcpp::XPtr<tiledb::QueryCondition> rhs,
    const std::string& str) {
    check_xptr_tag<tiledb::QueryCondition>(lhs);
    check_xptr_tag<tiledb::QueryCondition>(lhs);
    tiledb_query_condition_combination_op_t
        op = _tiledb_query_string_to_condition_combination_op(str);
    tiledb::QueryCondition res = lhs->combine(*rhs.get(), op);
    auto query_cond = make_xptr<tiledb::QueryCondition>(
        new tiledb::QueryCondition(res));
    return query_cond;
}

// [[Rcpp::export]]
void libtiledbsoma_query_condition_set_use_enumeration(
    Rcpp::XPtr<tiledb::Context> ctx,
    Rcpp::XPtr<tiledb::QueryCondition> cond,
    bool use_enumeration) {
    check_xptr_tag<tiledb::Context>(ctx);
    check_xptr_tag<tiledb::QueryCondition>(cond);
    tiledb::QueryConditionExperimental::set_use_enumeration(
        *ctx.get(), *cond.get(), use_enumeration);
}

// [[Rcpp::export]]
Rcpp::XPtr<tiledb::QueryCondition> libtiledbsoma_query_condition_create(
    Rcpp::XPtr<tiledb::Context> ctx,
    const std::string& name,
    SEXP vec,
    const std::string& cond_op_string) {
    check_xptr_tag<tiledb::Context>(ctx);
    tiledb_query_condition_op_t op = _tiledb_query_string_to_condition_op(
        cond_op_string);
    // consider three cases of 'vec' based on R types:  int, double and
    // int64-as-double
    if (TYPEOF(vec) == INTSXP) {
        std::vector<int32_t> iv = Rcpp::as<std::vector<int32_t>>(vec);
        auto qc = tiledb::QueryConditionExperimental::create<int32_t>(
            *ctx.get(), name, iv, op);
        return make_xptr<tiledb::QueryCondition>(
            new tiledb::QueryCondition(qc));
    } else if (TYPEOF(vec) == REALSXP) {
        if (Rcpp::isInteger64(vec)) {
            std::vector<int64_t> dv = Rcpp::fromInteger64(
                Rcpp::NumericVector(vec));
            auto qc = tiledb::QueryConditionExperimental::create<int64_t>(
                *ctx.get(), name, dv, op);
            return make_xptr<tiledb::QueryCondition>(
                new tiledb::QueryCondition(qc));
        } else {
            std::vector<double> dv = Rcpp::as<std::vector<double>>(vec);
            auto qc = tiledb::QueryConditionExperimental::create<double>(
                *ctx.get(), name, dv, op);
            return make_xptr<tiledb::QueryCondition>(
                new tiledb::QueryCondition(qc));
        }
    } else if (TYPEOF(vec) == STRSXP) {
        std::vector<std::string> sv = Rcpp::as<std::vector<std::string>>(vec);
        auto qc = tiledb::QueryConditionExperimental::create(
            *ctx.get(), name, sv, op);
        return make_xptr<tiledb::QueryCondition>(
            new tiledb::QueryCondition(qc));
    } else {
        Rcpp::stop("No support (yet) for type '%s'.", Rcpp::type2name(vec));
    }
    return make_xptr<tiledb::QueryCondition>(R_NilValue);
}
