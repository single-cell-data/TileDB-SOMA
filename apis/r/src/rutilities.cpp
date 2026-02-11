
// we currently get deprecation warnings by default which are noisy
#ifndef TILEDB_NO_API_DEPRECATION_WARNINGS
#define TILEDB_NO_API_DEPRECATION_WARNINGS
#endif

#include <Rcpp.h>                 // for R interface to C++
#include <nanoarrow/nanoarrow.h>  // for C interface to Arrow
#include <nanoarrow/r.h>          // for C interface to Arrow (via R package)
#include <tiledbsoma/reindexer/reindexer.h>
#include <RcppInt64>  // for fromInteger64
#include <sstream>
#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"  // local declarations
#include "xptr-utils.h"  // xptr taggging utilitie

namespace tdbs = tiledbsoma;

void apply_dim_points(
    tdbs::common::ManagedQuery* mq,
    std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>>& name2dim,
    Rcpp::List lst) {
    std::vector<std::string> colnames = lst.attr("names");
    for (auto& nm : colnames) {
        auto dm = name2dim[nm];
        auto tp = dm->type();
        bool suitable = false;
        if (tp == TILEDB_UINT64) {
            Rcpp::NumericVector payload = lst[nm];
            std::vector<int64_t> iv = Rcpp::fromInteger64(payload, false);
            std::vector<uint64_t> uv(iv.size());
            const std::pair<uint64_t, uint64_t> pr = dm->domain<uint64_t>();
            for (size_t i = 0; i < iv.size(); i++) {
                uv[i] = static_cast<uint64_t>(iv[i]);
                if (uv[i] >= pr.first && uv[i] <= pr.second) {
                    mq->select_point<uint64_t>(nm, uv[i]);  // bonked when use with vector
                    std::stringstream ss;
                    ss << "[apply_dim_points] Applying dim point " << uv[i] << " on " << nm;
                    tdbs::common::logging::LOG_DEBUG(ss.str());
                    suitable = true;
                }
            }
        } else if (tp == TILEDB_INT64) {
            Rcpp::NumericVector payload = lst[nm];
            std::vector<int64_t> iv = Rcpp::fromInteger64(payload, false);
            const std::pair<int64_t, int64_t> pr = dm->domain<int64_t>();
            for (size_t i = 0; i < iv.size(); i++) {
                if (iv[i] >= pr.first && iv[i] <= pr.second) {
                    mq->select_point<int64_t>(nm, iv[i]);
                    std::stringstream ss;
                    ss << "[apply_dim_points] Applying dim point " << iv[i] << " on " << nm;
                    tdbs::common::logging::LOG_DEBUG(ss.str());
                    suitable = true;
                }
            }
        } else if (tp == TILEDB_FLOAT32) {
            Rcpp::NumericVector payload = lst[nm];
            const std::pair<float, float> pr = dm->domain<float>();
            for (R_xlen_t i = 0; i < payload.size(); i++) {
                float v = static_cast<float>(payload[i]);
                if (v >= pr.first && v <= pr.second) {
                    mq->select_point<float>(nm, v);
                    std::stringstream ss;
                    ss << "[apply_dim_points] Applying dim point " << v << " on " << nm;
                    tdbs::common::logging::LOG_DEBUG(ss.str());
                    suitable = true;
                }
            }
        } else if (tp == TILEDB_FLOAT64) {
            Rcpp::NumericVector payload = lst[nm];
            const std::pair<double, double> pr = dm->domain<double>();
            for (R_xlen_t i = 0; i < payload.size(); i++) {
                if (payload[i] >= pr.first && payload[i] <= pr.second) {
                    mq->select_point<double>(nm, payload[i]);
                    std::stringstream ss;
                    ss << "[apply_dim_points] Applying dim point " << payload[i] << " on " << nm;
                    tdbs::common::logging::LOG_DEBUG(ss.str());
                    suitable = true;
                }
            }
        } else if (tp == TILEDB_INT32) {
            Rcpp::IntegerVector payload = lst[nm];
            const std::pair<int32_t, int32_t> pr = dm->domain<int32_t>();
            for (R_xlen_t i = 0; i < payload.size(); i++) {
                if (payload[i] >= pr.first && payload[i] <= pr.second) {
                    mq->select_point<int32_t>(nm, payload[i]);
                    std::stringstream ss;
                    ss << "[apply_dim_points] Applying dim point " << payload[i] << " on " << nm;
                    tdbs::common::logging::LOG_DEBUG(ss.str());
                    suitable = true;
                }
            }
        } else {
            Rcpp::stop("Currently unsupported type: ", tiledb::impl::to_str(tp));
        }
        if (!suitable) {
            Rcpp::stop("Unsuitable dim points on dimension '%s' with domain %s", nm, dm->domain_to_str());
        }
    }
}

void apply_dim_ranges(
    tdbs::common::ManagedQuery* mq,
    std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>>& name2dim,
    Rcpp::List lst) {
    std::vector<std::string> colnames = lst.attr("names");
    for (auto& nm : colnames) {
        auto dm = name2dim[nm];
        auto tp = dm->type();
        bool suitable = false;
        if (tp == TILEDB_UINT64) {
            Rcpp::NumericMatrix mm = lst[nm];
            Rcpp::NumericMatrix::Column lo = mm.column(0);  // works as proxy for int and float types
            Rcpp::NumericMatrix::Column hi = mm.column(1);  // works as proxy for int and float types
            std::vector<std::pair<uint64_t, uint64_t>> vp(mm.nrow());
            const std::pair<uint64_t, uint64_t> pr = dm->domain<uint64_t>();
            for (int i = 0; i < mm.nrow(); i++) {
                uint64_t l = static_cast<uint64_t>(Rcpp::fromInteger64(lo[i]));
                uint64_t h = static_cast<uint64_t>(Rcpp::fromInteger64(hi[i]));
                vp[i] = std::make_pair(std::max(l, pr.first), std::min(h, pr.second));
                std::stringstream ss;
                ss << "[apply_dim_ranges] Applying dim point " << i << " on " << nm << " with " << l << " - " << h;
                tdbs::common::logging::LOG_DEBUG(ss.str());
                suitable = l < pr.second && h > pr.first;  // lower must be less than max, higher
                                                           // more than min
            }
            if (suitable)
                mq->select_ranges<uint64_t>(nm, vp);
        } else if (tp == TILEDB_INT64) {
            Rcpp::NumericMatrix mm = lst[nm];
            std::vector<int64_t> lo = Rcpp::fromInteger64(mm.column(0), false);
            std::vector<int64_t> hi = Rcpp::fromInteger64(mm.column(1), false);
            std::vector<std::pair<int64_t, int64_t>> vp(mm.nrow());
            const std::pair<int64_t, int64_t> pr = dm->domain<int64_t>();
            for (int i = 0; i < mm.nrow(); i++) {
                vp[i] = std::make_pair(std::max(lo[i], pr.first), std::min(hi[i], pr.second));
                std::stringstream ss;
                ss << "[apply_dim_ranges] Applying dim point " << i << " on " << nm << " with " << lo[i] << " - "
                   << hi[i];
                tdbs::common::logging::LOG_DEBUG(ss.str());
                suitable = lo[i] < pr.second && hi[i] > pr.first;  // lower must be less than max,
                                                                   // higher more than min
            }
            if (suitable)
                mq->select_ranges<int64_t>(nm, vp);
        } else if (tp == TILEDB_FLOAT32) {
            Rcpp::NumericMatrix mm = lst[nm];
            Rcpp::NumericMatrix::Column lo = mm.column(0);  // works as proxy for int and float types
            Rcpp::NumericMatrix::Column hi = mm.column(1);  // works as proxy for int and float types
            std::vector<std::pair<float, float>> vp(mm.nrow());
            const std::pair<float, float> pr = dm->domain<float>();
            for (int i = 0; i < mm.nrow(); i++) {
                float l = static_cast<float>(lo[i]);
                float h = static_cast<float>(hi[i]);
                vp[i] = std::make_pair(std::max(l, pr.first), std::min(h, pr.second));
                std::stringstream ss;
                ss << "[apply_dim_ranges] Applying dim point " << i << " on " << nm << " with " << l << " - " << h;
                tdbs::common::logging::LOG_DEBUG(ss.str());
                suitable = l < pr.second && h > pr.first;  // lower must be less than max, higher
                                                           // more than min
            }
            if (suitable)
                mq->select_ranges<float>(nm, vp);
        } else if (tp == TILEDB_FLOAT64) {
            Rcpp::NumericMatrix mm = lst[nm];
            Rcpp::NumericMatrix::Column lo = mm.column(0);  // works as proxy for int and float types
            Rcpp::NumericMatrix::Column hi = mm.column(1);  // works as proxy for int and float types
            std::vector<std::pair<double, double>> vp(mm.nrow());
            const std::pair<double, double> pr = dm->domain<double>();
            for (int i = 0; i < mm.nrow(); i++) {
                vp[i] = std::make_pair(std::max(lo[i], pr.first), std::min(hi[i], pr.second));
                std::stringstream ss;
                ss << "[apply_dim_ranges] Applying dim point " << i << " on " << nm << " with " << lo[i] << " - "
                   << hi[i];
                tdbs::common::logging::LOG_DEBUG(ss.str());
                suitable = lo[i] < pr.second && hi[i] > pr.first;  // lower must be less than max,
                                                                   // higher more than min
            }
            if (suitable)
                mq->select_ranges<double>(nm, vp);
        } else if (tp == TILEDB_INT32) {
            Rcpp::IntegerMatrix mm = lst[nm];
            Rcpp::IntegerMatrix::Column lo = mm.column(0);  // works as proxy for int and float types
            Rcpp::IntegerMatrix::Column hi = mm.column(1);  // works as proxy for int and float types
            std::vector<std::pair<int32_t, int32_t>> vp(mm.nrow());
            const std::pair<int32_t, int32_t> pr = dm->domain<int32_t>();
            for (int i = 0; i < mm.nrow(); i++) {
                vp[i] = std::make_pair(std::max(lo[i], pr.first), std::min(hi[i], pr.second));
                std::stringstream ss;
                ss << "[apply_dim_ranges] Applying dim point " << i << " on " << nm << " with " << lo[i] << " - "
                   << hi[i];
                tdbs::common::logging::LOG_DEBUG(ss.str());
                suitable = lo[i] < pr.second && hi[i] > pr.first;  // lower must be less than max,
                                                                   // higher more than min
            }
            if (suitable)
                mq->select_ranges<int32_t>(nm, vp);
        } else {
            Rcpp::stop("Currently unsupported type: ", tiledb::impl::to_str(tp));
        }
        if (!suitable) {
            Rcpp::stop("Unsuitable dim ranges on dimension '%s' with domain %s", nm, dm->domain_to_str());
        }
    }
}

// initialize arrow schema and array, respectively
Rcpp::XPtr<ArrowSchema> schema_setup_struct(Rcpp::XPtr<ArrowSchema> schxp, int64_t n_children) {
    ArrowSchema* schema = schxp.get();
    auto type = NANOARROW_TYPE_STRUCT;

    ArrowSchemaInit(schema);  // modified from ArrowSchemaInitFromType()
    int result = ArrowSchemaSetType(schema, type);
    if (result != NANOARROW_OK) {
        schema->release(schema);
        Rcpp::stop("Error setting struct schema");
    }

    // now adapted from ArrowSchemaAllocateChildren
    if (schema->children != NULL)
        Rcpp::stop("Error allocation as children not null");

    if (n_children > 0) {
        auto ptr = (struct ArrowSchema**)ArrowMalloc(n_children * sizeof(struct ArrowSchema*));
        Rcpp::XPtr<ArrowSchema*> schema_ptrxp = make_xptr(ptr, false);
        schema->children = schema_ptrxp.get();
        if (schema->children == NULL)
            Rcpp::stop("Failed to allocate ArrowSchema*");

        schema->n_children = n_children;
        memset(schema->children, 0, n_children * sizeof(struct ArrowSchema*));

        for (int64_t i = 0; i < n_children; i++) {
            schema->children[i] = schema_owning_xptr();
            if (schema->children[i] == NULL)
                Rcpp::stop("Error allocation schema child %ld", i);
            schema->children[i]->release = NULL;
        }
    }
    return schxp;
}

std::vector<int64_t> i64_from_rcpp_numeric(const Rcpp::NumericVector& input) {
    auto ndim = input.size();
    std::vector<int64_t> output(ndim);
    for (auto i = 0; i < ndim; i++) {
        output[i] = input[i];
    }

    return output;
}

// formerly stats.cpp

//' TileDB SOMA Statistics
//'
//' These functions expose the TileDB Core functionality for performance
//'  measurements and statistics
//'
//' \itemize{
//'  \item \code{tiledbsoma_stats_enable()}/\code{tiledbsoma_stats_disable()}:
//'   Enable and disable TielDB's internal statistics
//'  \item \code{tiledbsoma_stats_reset()}: Reset all statistics to \code{0}
//'  \item \code{tiledbsoma_stats_dump()}: Dump all statistics as a JSON string
//'  \item \code{tiledbsoma_stats_show()}: Pretty-print the JSON statistics
//' }
//'
//' @return \code{tiledbsoma_stats_show()}: a single-length character vector
//' with the TileDB statistics encoded in JSON format
//'
//' @return All other functions invisibly return \code{NULL}
//'
//' @name tiledbsoma_stats
//'
//' @export
//'
// [[Rcpp::export]]
void tiledbsoma_stats_enable() {
    tiledbsoma::stats::enable();
}

//' @rdname tiledbsoma_stats
//'
//' @export
//'
// [[Rcpp::export]]
void tiledbsoma_stats_disable() {
    tiledbsoma::stats::disable();
}

//' @rdname tiledbsoma_stats
//'
//' @export
//'
// [[Rcpp::export]]
void tiledbsoma_stats_reset() {
    tiledbsoma::stats::reset();
}

//' @rdname tiledbsoma_stats
//'
//' @export
//'
// [[Rcpp::export]]
std::string tiledbsoma_stats_dump() {
    return tiledbsoma::stats::dump();
}

// formerly version.cpp

// libtiledbsoma version
//
// Returns a string with version information for libtiledbsoma and the linked
// TileDB Embedded library. If argument `compact` is set to `TRUE`, a shorter
// version of just the TileDB Embedded library version is returned.
//
// [[Rcpp::export]]
std::string libtiledbsoma_version(const bool compact = false, const bool major_minor_only = false) {
    if (compact) {
        auto v = tiledbsoma::version::embedded_version_triple();
        std::ostringstream txt;
        if (major_minor_only) {
            txt << std::get<0>(v) << "." << std::get<1>(v);
        } else {
            txt << std::get<0>(v) << "." << std::get<1>(v) << "." << std::get<2>(v);
        }
        return txt.str();
    } else {
        return tiledbsoma::version::as_string();
    }
}

// TileDB embedded version
//
// Gets the version of the TileDB Embedded library that is currently in use.
//
// [[Rcpp::export]]
Rcpp::IntegerVector tiledb_embedded_version() {
    std::tuple<int, int, int> triple = tiledbsoma::version::embedded_version_triple();
    return Rcpp::IntegerVector::create(std::get<0>(triple), std::get<1>(triple), std::get<2>(triple));
}

// Also present in tiledb-r but only after 0.23.0 so this can be removed (and
// the call to it updated) once we hit a new tiledb-r release 0.24.0 (or 0.23.1)
//
// [[Rcpp::export]]
size_t tiledb_datatype_max_value(const std::string& datatype) {
    if (datatype == "INT8")
        return std::numeric_limits<int8_t>::max();
    else if (datatype == "UINT8")
        return std::numeric_limits<uint8_t>::max();
    else if (datatype == "INT16")
        return std::numeric_limits<int16_t>::max();
    else if (datatype == "UINT16")
        return std::numeric_limits<uint16_t>::max();
    else if (datatype == "INT32")
        return std::numeric_limits<int32_t>::max();
    else if (datatype == "UINT32")
        return std::numeric_limits<uint32_t>::max();
    else if (datatype == "INT64")
        return std::numeric_limits<int64_t>::max();
    else if (datatype == "UINT64")
        return std::numeric_limits<uint64_t>::max();
    else
        Rcpp::stop("currently unsupported datatype (%s)", datatype);
}

// Make (optional) TimestampRange from (nullable, two-element)
// Rcpp::DatetimeVector
std::optional<tdbs::TimestampRange> makeTimestampRange(Rcpp::Nullable<Rcpp::DatetimeVector> tsvec) {
    // optional timestamp, defaults to 'none' aka std::nullopt
    std::optional<tdbs::TimestampRange> tsrng = std::nullopt;

    if (tsvec.isNotNull()) {
        // an Rcpp 'Nullable' is a decent compromise between adhering to SEXP
        // semantics and having 'optional' behaviour -- but when there is a
        // value we need to be explicit
        Rcpp::DatetimeVector vec(tsvec);  // vector of Rcpp::Datetime ie POSIXct
                                          // w/ (fract.) secs since epoch
        if (vec.size() == 1) {
            tsrng = std::make_pair<uint64_t>(
                0, static_cast<uint64_t>(Rcpp::Datetime(vec[0]).getFractionalTimestamp() * 1000));
        } else if (vec.size() == 2) {
            tsrng = std::make_pair<uint64_t>(
                static_cast<uint64_t>(Rcpp::Datetime(vec[0]).getFractionalTimestamp() * 1000),
                static_cast<uint64_t>(Rcpp::Datetime(vec[1]).getFractionalTimestamp() * 1000));
        } else {
            Rcpp::stop("TimestampRange must be a one or two-element vector");
        }
    }

    return tsrng;
}

SEXP convert_domainish(const tdbs::common::arrow::ArrowTable& arrow_table) {
    ArrowArray* arrow_array = arrow_table.first.get();
    ArrowSchema* arrow_schema = arrow_table.second.get();

    auto schemaxp = nanoarrow_schema_owning_xptr();
    auto sch = nanoarrow_output_schema_from_xptr(schemaxp);
    exitIfError(ArrowSchemaInitFromType(sch, NANOARROW_TYPE_STRUCT), "Bad schema init");
    exitIfError(ArrowSchemaSetName(sch, ""), "Bad schema name");
    exitIfError(ArrowSchemaAllocateChildren(sch, arrow_schema->n_children), "Bad schema children alloc");

    if (arrow_array->n_children != arrow_schema->n_children) {
        Rcpp::stop(
            "schema/data column mismatch %d != %d\n", (int)arrow_array->n_children, (int)arrow_schema->n_children);
    }
    auto ncol = arrow_schema->n_children;

    auto arrayxp = nanoarrow_array_owning_xptr();
    auto arr = nanoarrow_output_array_from_xptr(arrayxp);
    exitIfError(ArrowArrayInitFromType(arr, NANOARROW_TYPE_STRUCT), "Bad array init");
    exitIfError(ArrowArrayAllocateChildren(arr, arrow_array->n_children), "Bad array children alloc");

    for (size_t i = 0; i < static_cast<size_t>(ncol); i++) {
        if (arrow_array->children[i]->n_buffers == 3) {
            // Arrow semantics: variable-length: buffers 0,1,2 are validity,
            // offsets, data
            std::vector<std::string> lohi = tiledbsoma::ArrowAdapter::get_array_string_column(
                arrow_array->children[i], arrow_schema->children[i]);
            std::stringstream ss;
            ss << "[domainish] name {} format {} length {} lo {} hi {}" << std::string(arrow_schema->children[i]->name)
               << " format " << std::string(arrow_schema->children[i]->format) << " length "
               << arrow_array->children[i]->length << "lo " << lohi[0] << " hi " << lohi[1];
            tdbs::common::logging::LOG_DEBUG(ss.str());
        } else {
            // Arrow semantics: non-variable-length: buffers 0,1 are validity &
            // data
            std::stringstream ss;
            ss << "[domainish] name " << std::string(arrow_schema->children[i]->name) << " format "
               << std::string(arrow_schema->children[i]->format) << " length " << arrow_array->children[i]->length;
            tdbs::common::logging::LOG_DEBUG(ss.str());
        }

        ArrowArrayMove(arrow_array->children[i], arr->children[i]);
        ArrowSchemaMove(arrow_schema->children[i], sch->children[i]);
    }

    // Nanoarrow special: stick schema into xptr tag to return single SEXP
    array_xptr_set_schema(arrayxp, schemaxp);  // embed schema in array

    return arrayxp;
}

static std::map<std::string, std::string> _type_name_remap = {
    {"int8", "c"},
    {"int16", "s"},
    {"int32", "i"},
    {"int64", "l"},
    {"uint8", "C"},
    {"uint16", "S"},
    {"uint32", "I"},
    {"uint64", "L"},
    {"utf8", "u"},
    {"large_utf8", "U"},
    {"bool", "b"},
    {"float", "f"},
    {"double", "g"}};

std::string remap_arrow_type_code_r_to_c(std::string input) {
    auto it = _type_name_remap.find(input);
    if (it == _type_name_remap.end()) {
        return input;
    } else {
        return it->second;
    }
}

// Taken from tiledb-r
// https://github.com/TileDB-Inc/TileDB-R/blob/525bdfc0f34aadb74a312a5d8428bd07819a8f83/src/libtiledb.cpp#L246C1-L261C2
const char* _tiledb_layout_to_string(tiledb_layout_t layout) {
    switch (layout) {
        case TILEDB_ROW_MAJOR:
            return "ROW_MAJOR";
        case TILEDB_COL_MAJOR:
            return "COL_MAJOR";
        case TILEDB_GLOBAL_ORDER:
            return "GLOBAL_ORDER";
        case TILEDB_UNORDERED:
            return "UNORDERED";
        case TILEDB_HILBERT:
            return "HILBERT";
        default:
            Rcpp::stop("unknown tiledb_layout_t (%d)", layout);
    }
}

// internal helper function for
// `_get_filter_options()`
// <Numeric> should be `int32_t` or `double`
// TODO: replace with safer API from core when available
template <typename Numeric>
Numeric _get_filter_option(Rcpp::XPtr<tiledb::Filter> filter, tiledb_filter_option_t option) {
    try {
        switch (option) {
            case TILEDB_BIT_WIDTH_MAX_WINDOW:
            case TILEDB_POSITIVE_DELTA_MAX_WINDOW:
                return static_cast<Numeric>(filter->get_option<uint32_t>(option));
            case TILEDB_SCALE_FLOAT_BYTEWIDTH:
                return static_cast<Numeric>(filter->get_option<uint64_t>(option));
            default:
                return filter->get_option<Numeric>(option);
        }
    } catch (const tiledb::TileDBError& e) {
        const std::string msg = e.what();
        if (msg.find("unknown option") != std::string::npos) {
            if (std::is_same<Numeric, double>::value) {
                return R_NaReal;
            }
            return R_NaInt;
        }
        if (msg.find("Filter does not support options") != std::string::npos) {
            if (std::is_same<Numeric, double>::value) {
                return R_NaReal;
            }
            return R_NaInt;
        }
        throw(e);
    }
}

// filter options taken from tiledb-r
// https://github.com/TileDB-Inc/TileDB-R/blob/525bdfc0f34aadb74a312a5d8428bd07819a8f83/src/libtiledb.cpp#L369-L388
Rcpp::List _get_filter_options(Rcpp::XPtr<tiledb::Filter> filter) {
    return Rcpp::List::create(
        Rcpp::Named("filter_type") = tiledb::Filter::to_str(filter->filter_type()),
        Rcpp::Named("compression_level") = _get_filter_option<int32_t>(filter, TILEDB_COMPRESSION_LEVEL),
        Rcpp::Named("bit_width") = _get_filter_option<int32_t>(filter, TILEDB_BIT_WIDTH_MAX_WINDOW),
        Rcpp::Named("positive_delta") = _get_filter_option<int32_t>(filter, TILEDB_POSITIVE_DELTA_MAX_WINDOW),
        Rcpp::Named("float_bytewidth") = _get_filter_option<double>(filter, TILEDB_SCALE_FLOAT_BYTEWIDTH),
        Rcpp::Named("float_factor") = _get_filter_option<double>(filter, TILEDB_SCALE_FLOAT_FACTOR),
        Rcpp::Named("float_offset") = _get_filter_option<double>(filter, TILEDB_SCALE_FLOAT_OFFSET));
}

// adapted from tiledb-r
// https://github.com/TileDB-Inc/TileDB-R/blob/57d3d929a5cd6d6e9bc4bf5bda0036c736ca46dc/src/libtiledb.cpp#L927-L1039
SEXP _get_dim_domain(Rcpp::XPtr<tiledb::Dimension> dim) {
    auto dim_type = dim->type();
    switch (dim_type) {
        case TILEDB_FLOAT32: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_FLOAT32>::type;
            return Rcpp::NumericVector({dim->domain<DataType>().first, dim->domain<DataType>().second});
        }
        case TILEDB_FLOAT64: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_FLOAT64>::type;
            auto d1 = dim->domain<DataType>().first;
            auto d2 = dim->domain<DataType>().second;
            if (d1 == R_NaReal || d2 == R_NaReal) {
                Rcpp::stop("tiledb_dim domain FLOAT64 value not representable as an R double");
            }
            return Rcpp::NumericVector({d1, d2});
        }
        case TILEDB_INT8: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_INT8>::type;
            return Rcpp::IntegerVector({dim->domain<DataType>().first, dim->domain<DataType>().second});
        }
        case TILEDB_UINT8: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_UINT8>::type;
            return Rcpp::IntegerVector({dim->domain<DataType>().first, dim->domain<DataType>().second});
        }
        case TILEDB_INT16: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_INT16>::type;
            return Rcpp::IntegerVector({dim->domain<DataType>().first, dim->domain<DataType>().second});
        }
        case TILEDB_UINT16: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_UINT16>::type;
            return Rcpp::IntegerVector({dim->domain<DataType>().first, dim->domain<DataType>().second});
        }
        case TILEDB_INT32: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_INT32>::type;
            auto d1 = dim->domain<DataType>().first;
            auto d2 = dim->domain<DataType>().second;
            if (d1 == R_NaInt || d2 == R_NaInt) {
                Rcpp::stop("tiledb_dim domain INT32 value not representable as an R integer");
            }
            return Rcpp::IntegerVector({d1, d2});
        }
        case TILEDB_UINT32: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_UINT32>::type;
            auto d1 = dim->domain<DataType>().first;
            auto d2 = dim->domain<DataType>().second;
            auto uint_max = std::numeric_limits<uint32_t>::max();
            if (d1 > uint_max || d2 > uint_max) {
                Rcpp::stop("tiledb_dim domain UINT32 value not representable as an R integer64 type");
            }
            return Rcpp::NumericVector({static_cast<double>(d1), static_cast<double>(d2)});
        }
        case TILEDB_INT64: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_INT64>::type;
            auto d1 = dim->domain<DataType>().first;
            auto d2 = dim->domain<DataType>().second;
            return Rcpp::NumericVector({static_cast<double>(d1), static_cast<double>(d2)});
        }
        case TILEDB_UINT64: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_UINT64>::type;
            auto d1 = dim->domain<DataType>().first;
            auto d2 = dim->domain<DataType>().second;
            return Rcpp::NumericVector({static_cast<double>(d1), static_cast<double>(d2)});
        }
        case TILEDB_DATETIME_YEAR:
        case TILEDB_DATETIME_MONTH:
        case TILEDB_DATETIME_WEEK:
        case TILEDB_DATETIME_DAY:
        case TILEDB_DATETIME_HR:
        case TILEDB_DATETIME_MIN:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_DATETIME_PS:
        case TILEDB_DATETIME_FS:
        case TILEDB_DATETIME_AS: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_INT64>::type;
            auto d1 = dim->domain<DataType>().first;
            auto d2 = dim->domain<DataType>().second;
            auto int32_max = std::numeric_limits<int32_t>::max();
            if (d1 <= R_NaInt || d1 > int32_max || d2 <= R_NaInt || d2 > int32_max) {
                return Rcpp::NumericVector({static_cast<double>(d1), static_cast<double>(d2)});
            }
            return Rcpp::IntegerVector({static_cast<int32_t>(d1), static_cast<int32_t>(d2)});
        }
        default:
            Rcpp::stop("invalid tiledb_dim domain type (%s)", tiledb::impl::to_str(dim_type));
    }
}

// adapted from tiledb-r
// https://github.com/TileDB-Inc/TileDB-R/blob/57d3d929a5cd6d6e9bc4bf5bda0036c736ca46dc/src/libtiledb.cpp#L1042-L1137
SEXP _get_dim_tile(Rcpp::XPtr<tiledb::Dimension> dim) {
    auto dim_type = dim->type();
    switch (dim_type) {
        case TILEDB_FLOAT32: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_FLOAT32>::type;
            return Rcpp::wrap(static_cast<double>(dim->tile_extent<DataType>()));
        }
        case TILEDB_FLOAT64: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_FLOAT64>::type;
            auto t = dim->tile_extent<DataType>();
            if (t == R_NaReal) {
                Rcpp::stop("tiledb_dim tile FLOAT64 value not representable as an R double");
            }
            return Rcpp::wrap(static_cast<double>(t));
        }
        case TILEDB_INT8: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_INT8>::type;
            return Rcpp::wrap(static_cast<int32_t>(dim->tile_extent<DataType>()));
        }
        case TILEDB_UINT8: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_UINT8>::type;
            return Rcpp::wrap(static_cast<int32_t>(dim->tile_extent<DataType>()));
        }
        case TILEDB_INT16: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_INT16>::type;
            return Rcpp::wrap(static_cast<int32_t>(dim->tile_extent<DataType>()));
        }
        case TILEDB_UINT16: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_UINT16>::type;
            return Rcpp::wrap(static_cast<int32_t>(dim->tile_extent<DataType>()));
        }
        case TILEDB_INT32: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_INT32>::type;
            auto t = dim->tile_extent<DataType>();
            if (t == R_NaInt) {
                Rcpp::stop("tiledb_dim tile INT32 value not representable as an R integer");
            }
            return Rcpp::wrap(static_cast<int32_t>(t));
        }
        case TILEDB_UINT32: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_UINT32>::type;
            auto t = dim->tile_extent<DataType>();
            if (t > std::numeric_limits<int32_t>::max()) {
                Rcpp::warning("tiledb_dim tile UINT32 value not representable as an R integer, returning double");
                return Rcpp::wrap(static_cast<double>(t));
            }
            return Rcpp::wrap(static_cast<int32_t>(t));
        }
        case TILEDB_DATETIME_YEAR:
        case TILEDB_DATETIME_MONTH:
        case TILEDB_DATETIME_WEEK:
        case TILEDB_DATETIME_DAY:
        case TILEDB_DATETIME_HR:
        case TILEDB_DATETIME_MIN:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_DATETIME_PS:
        case TILEDB_DATETIME_FS:
        case TILEDB_DATETIME_AS:
        case TILEDB_INT64: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_INT64>::type;
            return Rcpp::wrap(static_cast<double>(dim->tile_extent<DataType>()));
        }
        case TILEDB_UINT64: {
            using DataType = tiledb::impl::tiledb_to_type<TILEDB_UINT64>::type;
            return Rcpp::wrap(static_cast<double>(dim->tile_extent<DataType>()));
        }
        default:
            Rcpp::stop("invalid tiledb_dim domain type (%s)", tiledb::impl::to_str(dim_type));
    }
}

Rcpp::List metadata_as_rlist(std::map<std::string, tiledbsoma::MetadataValue>& mvmap) {
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
            auto str = std::string((char*)ptr, len);
            lst.push_back(str);
        } else if (dtype == TILEDB_INT64) {
            std::vector<int64_t> v(len);
            std::memcpy(&(v[0]), ptr, len * sizeof(int64_t));
            lst.push_back(Rcpp::toInteger64(v));
        } else if (dtype == TILEDB_INT32) {
            Rcpp::IntegerVector v(len);
            std::memcpy(v.begin(), ptr, len * sizeof(int32_t));
        } else {
            auto txt = tiledb::impl::type_to_str(dtype);
            Rcpp::stop("Currently unsupported type '%s'", txt.c_str());
        }
    }
    lst.attr("names") = Rcpp::CharacterVector(namvec.begin(), namvec.end());
    return lst;
};
