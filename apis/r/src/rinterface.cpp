#include <Rcpp.h>                   // for R interface to C++
#include <nanoarrow/r.h>            // for C interface to Arrow (via R package)
#include <RcppInt64>                // for fromInteger64
#include <nanoarrow/nanoarrow.hpp>  // for C/C++ interface to Arrow
#include <sstream>
#include <type_traits>

// we currently get deprecation warnings by default which are noisy
#ifndef TILEDB_NO_API_DEPRECATION_WARNINGS
#define TILEDB_NO_API_DEPRECATION_WARNINGS
#endif

// We get these via nanoarrow and must cannot include carrow.h again
#define ARROW_SCHEMA_AND_ARRAY_DEFINED 1
#include <tiledbsoma/reindexer/reindexer.h>
#include <tiledbsoma/tiledbsoma>
#if TILEDB_VERSION_MAJOR == 2 && TILEDB_VERSION_MINOR >= 4
#include <tiledb/tiledb_experimental>
#endif

#include "rutilities.h"  // local declarations
#include "xptr-utils.h"  // xptr taggging utilities

// (Adapted) helper functions from nanoarrow
//
// Create an external pointer with the proper class and that will release any
// non-null, non-released pointer when garbage collected. We use a tagged XPtr,
// but do not set an XPtr finalizer
Rcpp::XPtr<ArrowSchema> schema_owning_xptr(void) {
    struct ArrowSchema* schema = (struct ArrowSchema*)ArrowMalloc(sizeof(struct ArrowSchema));
    if (schema == NULL)
        Rcpp::stop("Failed to allocate ArrowSchema");
    schema->release = NULL;
    Rcpp::XPtr<ArrowSchema> schema_xptr = make_xptr(schema, false);
    return schema_xptr;
}

// Create an external pointer with the proper class and that will release any
// non-null, non-released pointer when garbage collected. We use a tagged XPtr,
// but do not set an XPtr finalizer
Rcpp::XPtr<ArrowArray> array_owning_xptr(void) {
    struct ArrowArray* array = (struct ArrowArray*)ArrowMalloc(sizeof(struct ArrowArray));
    if (array == NULL)
        Rcpp::stop("Failed to allocate ArrowArray");
    array->release = NULL;
    Rcpp::XPtr<ArrowArray> array_xptr = make_xptr(array, false);
    return array_xptr;
}

namespace tdbs = tiledbsoma;

SEXP soma_array_read_impl(
    tiledbsoma::SOMAArray* soma_array,
    Rcpp::Nullable<Rcpp::CharacterVector> colnames = R_NilValue,
    Rcpp::Nullable<Rcpp::XPtr<tiledb::QueryCondition>> qc = R_NilValue,
    Rcpp::Nullable<Rcpp::List> dim_points = R_NilValue,
    Rcpp::Nullable<Rcpp::List> dim_ranges = R_NilValue,
    std::string batch_size = "auto",
    std::string result_order = "auto",
    const std::string& loglevel = "auto") {
    if (loglevel != "auto") {
        tdbs::common::logging::LOG_SET_LEVEL(loglevel);
    }

    std::stringstream ss;
    ss << "[soma_array_reader] Reading from " << soma_array->uri();
    tdbs::common::logging::LOG_DEBUG(ss.str());

    auto mq = soma_array->create_managed_query();

    auto tdb_result_order = get_tdb_result_order(result_order);
    mq.set_layout(tdb_result_order);

    std::vector<std::string> column_names = {};
    if (!colnames.isNull()) {  // If we have column names, select them
        column_names = Rcpp::as<std::vector<std::string>>(colnames);
        std::stringstream ss;
        ss << "[soma_array_reader] Selecting " << column_names.size() << " columns";
        tdbs::common::logging::LOG_DEBUG(ss.str());
        if (!column_names.empty()) {
            mq.select_columns(column_names);
        }
    }

    std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>> name2dim;
    std::shared_ptr<tiledb::ArraySchema> schema = soma_array->tiledb_schema();
    tiledb::Domain domain = schema->domain();
    std::vector<tiledb::Dimension> dims = domain.dimensions();
    for (auto& dim : dims) {
        std::stringstream ss;
        ss << "[soma_array_reader] Dimension " << dim.name() << " type " << tiledb::impl::to_str(dim.type())
           << " extent " << dim.tile_extent_to_str();
        tdbs::common::logging::LOG_DEBUG(ss.str());
        name2dim.emplace(std::make_pair(dim.name(), std::make_shared<tiledb::Dimension>(dim)));
    }

    // If we have a query condition, apply it
    if (!qc.isNull()) {
        tdbs::common::logging::LOG_DEBUG("[soma_array_reader_impl] Applying query condition");
        Rcpp::XPtr<tiledb::QueryCondition> qcxp(qc);
        mq.set_condition(*qcxp);
    }

    // If we have dimension points, apply them
    // The interface is named list, where each (named) list elements is one
    // (named) dimesion The List element is a simple vector of points and each
    // point is applied to the named dimension
    if (!dim_points.isNull()) {
        Rcpp::List lst(dim_points);
        apply_dim_points(&mq, name2dim, lst);
    }

    // If we have a dimension points, apply them
    if (!dim_ranges.isNull()) {
        Rcpp::List lst(dim_ranges);
        apply_dim_ranges(&mq, name2dim, lst);
    }

    // Getting next batch:  std::optional<std::shared_ptr<ArrayBuffers>>
    auto sr_data = mq.read_next();
    if (!mq.results_complete()) {
        Rcpp::stop(
            "Read of '%s' is incomplete.\nConsider increasing the memory allocation via the configuration\noption "
            "'soma.init_buffer_bytes', or using iterated partial reads.",
            soma_array->uri());
    }
    {
        std::stringstream ss;
        ss << "[soma_array_reader] Read complete with " << sr_data->get()->num_rows() << " rows and "
           << sr_data->get()->names().size() << " cols";
        tdbs::common::logging::LOG_DEBUG(ss.str());
    }
    const std::vector<std::string> names = sr_data->get()->names();
    auto ncol = names.size();
    // Schema first
    auto schemaxp = nanoarrow_schema_owning_xptr();
    auto sch = nanoarrow_output_schema_from_xptr(schemaxp);
    exitIfError(ArrowSchemaInitFromType(sch, NANOARROW_TYPE_STRUCT), "Bad schema init");
    exitIfError(ArrowSchemaSetName(sch, ""), "Bad schema name");
    exitIfError(ArrowSchemaAllocateChildren(sch, ncol), "Bad schema children alloc");

    // Array second
    auto arrayxp = nanoarrow_array_owning_xptr();
    auto arr = nanoarrow_output_array_from_xptr(arrayxp);
    exitIfError(ArrowArrayInitFromType(arr, NANOARROW_TYPE_STRUCT), "Bad array init");
    exitIfError(ArrowArrayAllocateChildren(arr, ncol), "Bad array children alloc");

    auto pp_vector = tdbs::ArrowAdapter::buffer_to_arrow(sr_data.value(), true);

    arr->length = 0;  // initial value
    for (size_t i = 0; i < ncol; i++) {
        {
            std::stringstream ss;
            ss << "[soma_array_reader] Accessing '" << names[i] << "' at pos " << i;
            tdbs::common::logging::LOG_DEBUG(ss.str());
        }

        auto& pp = pp_vector[i];  // this is pair of array and schema pointer

        // pp.first.get(), sizeof(ArrowArray));
        ArrowArrayMove(pp.first.get(), arr->children[i]);
        ArrowSchemaMove(pp.second.get(), sch->children[i]);
        {
            std::stringstream ss;
            ss << "[soma_array_reader] Incoming name " << std::string(pp.second->name) << " length "
               << pp.first->length;
            tdbs::common::logging::LOG_DEBUG(ss.str());
        }
        if (pp.first->length > arr->length) {
            std::stringstream ss;
            ss << "[soma_array_reader] Setting array length to " << pp.first->length;
            tdbs::common::logging::LOG_DEBUG(ss.str());
            arr->length = pp.first->length;
        }
    }

    // Nanoarrow special: stick schema into xptr tag to return single SEXP
    array_xptr_set_schema(arrayxp, schemaxp);  // embed schema in array
    return arrayxp;
}

// [[Rcpp::export]]
SEXP soma_array_read(
    Rcpp::XPtr<tiledbsoma::SOMAArray> soma_array,
    Rcpp::Nullable<Rcpp::CharacterVector> colnames = R_NilValue,
    Rcpp::Nullable<Rcpp::XPtr<tiledb::QueryCondition>> qc = R_NilValue,
    Rcpp::Nullable<Rcpp::List> dim_points = R_NilValue,
    Rcpp::Nullable<Rcpp::List> dim_ranges = R_NilValue,
    std::string batch_size = "auto",
    std::string result_order = "auto",
    const std::string& loglevel = "auto",
    Rcpp::Nullable<Rcpp::DatetimeVector> timestamprange = R_NilValue) {
    return soma_array_read_impl(
        soma_array.get(), colnames, qc, dim_points, dim_ranges, batch_size, result_order, loglevel);
}

// [[Rcpp::export]]
SEXP soma_array_reader_impl(
    const std::string& uri,
    Rcpp::XPtr<somactx_wrap_t> ctxxp,
    Rcpp::Nullable<Rcpp::CharacterVector> colnames = R_NilValue,
    Rcpp::Nullable<Rcpp::XPtr<tiledb::QueryCondition>> qc = R_NilValue,
    Rcpp::Nullable<Rcpp::List> dim_points = R_NilValue,
    Rcpp::Nullable<Rcpp::List> dim_ranges = R_NilValue,
    std::string batch_size = "auto",
    std::string result_order = "auto",
    const std::string& loglevel = "auto",
    Rcpp::Nullable<Rcpp::DatetimeVector> timestamprange = R_NilValue) {
    if (loglevel != "auto") {
        tdbs::common::logging::LOG_SET_LEVEL(loglevel);
    }

    std::optional<tdbs::TimestampRange> tsrng = makeTimestampRange(timestamprange);
    if (timestamprange.isNotNull()) {
        Rcpp::DatetimeVector vec(timestamprange);
        std::stringstream ss;
        ss << "[soma_array_reader] timestamp range (" << vec[0] << ", " << vec[1] << ")";
        tdbs::common::logging::LOG_DEBUG(ss.str());
    }

    auto soma_array = tiledbsoma::SOMAArray::open(OpenMode::soma_read, uri, ctxxp->ctxptr, tsrng);

    auto retval = soma_array_read_impl(
        soma_array.get(), colnames, qc, dim_points, dim_ranges, batch_size, result_order, loglevel);
    soma_array->close();
    return retval;
}

//' Set TileDB-SOMA Logging Level
//'
//' Set the logging level for the \R package and underlying C++ library
//'
//' @param level A character value with logging level. May be \dQuote{trace}, \dQuote{debug}, \dQuote{info},
//' or \dQuote{warn}
//'
//' @return Invisibly returns \code{NULL}
//'
//' @export
//'
// [[Rcpp::export]]
void set_log_level(const std::string& level) {
    tdbs::common::logging::LOG_SET_LEVEL(level);
}

//' Set a trace message
//'
//' @param msg The message to set with level 'trace'
//'
//' @return Invisibly returns \code{NULL}
//'
//' @noRd
//'
// [[Rcpp::export]]
void soma_trace(const std::string& msg) {
    tdbs::common::logging::LOG_TRACE(msg);
}

//' Set a debug message
//'
//' @param msg The message to set with level 'debug'
//'
//' @return Invisibly returns \code{NULL}
//'
//' @noRd
//'
// [[Rcpp::export]]
void soma_debug(const std::string& msg) {
    tdbs::common::logging::LOG_DEBUG(msg);
}

//' Set a info message
//'
//' @param msg The message to set with level 'info'
//'
//' @return Invisibly returns \code{NULL}
//'
//' @noRd
//'
// [[Rcpp::export]]
void soma_info(const std::string& msg) {
    tdbs::common::logging::LOG_INFO(msg);
}

//' Set a warn message
//'
//' @param msg The message to set with level 'warn'
//'
//' @return Invisibly returns \code{NULL}
//'
//' @noRd
//'
// [[Rcpp::export]]
void soma_warn(const std::string& msg) {
    tdbs::common::logging::LOG_WARN(msg);
}

// [[Rcpp::export]]
double nnz(Rcpp::XPtr<tiledbsoma::SOMAArray> array) {
    if (!array) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    return static_cast<double>(array->nnz());
}

// [[Rcpp::export]]
bool check_arrow_schema_tag(Rcpp::XPtr<ArrowSchema> xp) {
    check_xptr_tag<ArrowSchema>(xp);  // throws if mismatched
    return true;
}

// [[Rcpp::export]]
bool check_arrow_array_tag(Rcpp::XPtr<ArrowArray> xp) {
    check_xptr_tag<ArrowArray>(xp);  // throws if mismatched
    return true;
}

// [[Rcpp::export]]
Rcpp::NumericVector shape(Rcpp::XPtr<tiledbsoma::SOMAArray> array) {
    if (!array) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    return Rcpp::toInteger64(array->shape());
}

// [[Rcpp::export]]
Rcpp::NumericVector maxshape(Rcpp::XPtr<tiledbsoma::SOMAArray> array) {
    if (!array) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    return Rcpp::toInteger64(array->maxshape());
}

// [[Rcpp::export]]
SEXP non_empty_domain(Rcpp::XPtr<tiledbsoma::SOMAArray> array) {
    if (!array) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    auto arrow_table = array->get_non_empty_domain();
    return convert_domainish(arrow_table);
}

// [[Rcpp::export]]
SEXP domain(Rcpp::XPtr<tiledbsoma::SOMADataFrame> dataframe) {
    if (!dataframe) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    return convert_domainish(dataframe->get_soma_domain());
}

// [[Rcpp::export]]
SEXP maxdomain(Rcpp::XPtr<tiledbsoma::SOMADataFrame> dataframe) {
    if (!dataframe) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    return convert_domainish(dataframe->get_soma_maxdomain());
}

/** Only used for testing. */
// [[Rcpp::export]]
Rcpp::NumericVector maybe_soma_joinid_shape(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMADataFrame::open(uri, OpenMode::soma_read, ctxxp->ctxptr);
    auto retval = sr->maybe_soma_joinid_shape();
    sr->close();
    if (retval.has_value()) {
        return Rcpp::toInteger64(retval.value());
    } else {
        return Rcpp::NumericVector::create(NA_REAL);  // one element vector, and is.na() is true
    }
}

/** Only used for testing. */
// [[Rcpp::export]]
Rcpp::NumericVector maybe_soma_joinid_maxshape(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMADataFrame::open(uri, OpenMode::soma_read, ctxxp->ctxptr);
    auto retval = sr->maybe_soma_joinid_maxshape();
    sr->close();
    if (retval.has_value()) {
        return Rcpp::toInteger64(retval.value());
    } else {
        return Rcpp::NumericVector::create(NA_REAL);  // one element vector, and is.na() is true
    }
}

// [[Rcpp::export]]
Rcpp::LogicalVector has_current_domain(Rcpp::XPtr<tiledbsoma::SOMAArray> array) {
    if (!array) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    return Rcpp::LogicalVector(array->has_current_domain());
}

// [[Rcpp::export]]
Rcpp::NumericVector ndim(Rcpp::XPtr<tiledbsoma::SOMAArray> array) {
    if (!array) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    return Rcpp::NumericVector::create(array->ndim());
}

// [[Rcpp::export]]
Rcpp::CharacterVector c_dimnames(Rcpp::XPtr<tiledbsoma::SOMAArray> array) {
    if (!array) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    auto lib_retval = array->dimension_names();

    size_t n = lib_retval.size();
    Rcpp::CharacterVector retval(n);
    for (size_t i = 0; i < n; i++) {
        retval[i] = lib_retval[i];
    }
    return retval;
}

// [[Rcpp::export]]
Rcpp::CharacterVector c_attrnames(Rcpp::XPtr<tiledbsoma::SOMAArray> array) {
    if (!array) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    auto lib_retval = array->attribute_names();

    size_t n = lib_retval.size();
    Rcpp::CharacterVector retval(n);
    for (size_t i = 0; i < n; i++) {
        retval[i] = lib_retval[i];
    }
    return retval;
}

// [[Rcpp::export]]
SEXP c_schema(Rcpp::XPtr<tiledbsoma::SOMAArray> array) {
    if (!array) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    auto lib_retval = array->arrow_schema(true);

    auto schemaxp = nanoarrow_schema_owning_xptr();
    auto sch = nanoarrow_output_schema_from_xptr(schemaxp);
    exitIfError(ArrowSchemaInitFromType(sch, NANOARROW_TYPE_STRUCT), "Bad schema init");
    exitIfError(ArrowSchemaSetName(sch, ""), "Bad schema name");
    exitIfError(ArrowSchemaAllocateChildren(sch, lib_retval->n_children), "Bad schema children alloc");

    for (size_t i = 0; i < static_cast<size_t>(lib_retval->n_children); i++) {
        std::stringstream ss;
        ss << "[c_schema] Accessing name '" << std::string(lib_retval->children[i]->name) << "' format '"
           << std::string(lib_retval->children[i]->format) << "' at position " << i;
        tdbs::common::logging::LOG_DEBUG(ss.str());
        ArrowSchemaMove(lib_retval->children[i], sch->children[i]);
    }

    return schemaxp;
}

// [[Rcpp::export]]
bool c_is_sparse(Rcpp::XPtr<tiledbsoma::SOMAArray> array) {
    if (!array) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    return array->tiledb_schema()->array_type() == TILEDB_SPARSE;
}

// [[Rcpp::export]]
bool c_allows_dups(Rcpp::XPtr<tiledbsoma::SOMAArray> array) {
    if (!array) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    return array->tiledb_schema()->allows_dups();
}

// [[Rcpp::export]]
double c_capacity(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::soma_read, uri, ctxxp->ctxptr);
    std::shared_ptr<tiledb::ArraySchema> sch = sr->tiledb_schema();
    sr->close();

    uint64_t cap = sch->capacity();
    return static_cast<double>(cap);
}

// [[Rcpp::export]]
std::string c_tile_order(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::soma_read, uri, ctxxp->ctxptr);
    std::shared_ptr<tiledb::ArraySchema> sch = sr->tiledb_schema();
    sr->close();

    auto order = sch->tile_order();
    return _tiledb_layout_to_string(order);
}

// [[Rcpp::export]]
std::string c_cell_order(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::soma_read, uri, ctxxp->ctxptr);
    std::shared_ptr<tiledb::ArraySchema> sch = sr->tiledb_schema();
    sr->close();

    auto order = sch->cell_order();
    return _tiledb_layout_to_string(order);
}

// [[Rcpp::export]]
Rcpp::List c_schema_filters(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::soma_read, uri, ctxxp->ctxptr);
    std::shared_ptr<tiledb::ArraySchema> sch = sr->tiledb_schema();
    sr->close();

    Rcpp::List filter_list, coord_filters, offset_filters, validity_filters;

    auto coords_filter_list = make_xptr<tiledb::FilterList>(new tiledb::FilterList(sch->coords_filter_list()));
    int ncoords_filters = static_cast<int32_t>(coords_filter_list->nfilters());
    for (int i = 0; i < ncoords_filters; i++) {
        auto filter = make_xptr<tiledb::Filter>(new tiledb::Filter(coords_filter_list->filter(i)));
        auto filter_type = tiledb::Filter::to_str(filter->filter_type());
        coord_filters[filter_type] = _get_filter_options(filter);
    }
    filter_list["coords"] = coord_filters;

    auto offset_filter_list = make_xptr<tiledb::FilterList>(new tiledb::FilterList(sch->offsets_filter_list()));
    int noffset_filters = static_cast<int32_t>(offset_filter_list->nfilters());
    for (int i = 0; i < noffset_filters; i++) {
        auto filter = make_xptr<tiledb::Filter>(new tiledb::Filter(offset_filter_list->filter(i)));
        auto filter_type = tiledb::Filter::to_str(filter->filter_type());
        offset_filters[filter_type] = _get_filter_options(filter);
    }
    filter_list["offsets"] = offset_filters;

    auto validity_filter_list = make_xptr<tiledb::FilterList>(new tiledb::FilterList(sch->validity_filter_list()));
    int nvalidity_filters = static_cast<int32_t>(validity_filter_list->nfilters());
    for (int i = 0; i < nvalidity_filters; i++) {
        auto filter = make_xptr<tiledb::Filter>(new tiledb::Filter(validity_filter_list->filter(i)));
        auto filter_type = tiledb::Filter::to_str(filter->filter_type());
        validity_filters[filter_type] = _get_filter_options(filter);
    }
    filter_list["validity"] = validity_filters;

    return filter_list;
}

// Taken from tiledb-r
// https://github.com/TileDB-Inc/TileDB-R/blob/525bdfc0f34aadb74a312a5d8428bd07819a8f83/src/libtiledb.cpp#L31C1-L110C2
const char* _tiledb_datatype_to_string(tiledb_datatype_t dtype) {
    switch (dtype) {
        case TILEDB_INT8:
            return "INT8";
        case TILEDB_UINT8:
            return "UINT8";
        case TILEDB_INT16:
            return "INT16";
        case TILEDB_UINT16:
            return "UINT16";
        case TILEDB_INT32:
            return "INT32";
        case TILEDB_UINT32:
            return "UINT32";
        case TILEDB_INT64:
            return "INT64";
        case TILEDB_UINT64:
            return "UINT64";
        case TILEDB_FLOAT32:
            return "FLOAT32";
        case TILEDB_FLOAT64:
            return "FLOAT64";
        case TILEDB_CHAR:
            return "CHAR";
        case TILEDB_STRING_ASCII:
            return "ASCII";
        case TILEDB_STRING_UTF8:
            return "UTF8";
        case TILEDB_STRING_UTF16:
            return "UTF16";
        case TILEDB_STRING_UTF32:
            return "UTF32";
        case TILEDB_STRING_UCS2:
            return "UCS2";
        case TILEDB_STRING_UCS4:
            return "UCS4";
        case TILEDB_ANY:
            return "ANY";
        case TILEDB_DATETIME_YEAR:
            return "DATETIME_YEAR";
        case TILEDB_DATETIME_MONTH:
            return "DATETIME_MONTH";
        case TILEDB_DATETIME_WEEK:
            return "DATETIME_WEEK";
        case TILEDB_DATETIME_DAY:
            return "DATETIME_DAY";
        case TILEDB_DATETIME_HR:
            return "DATETIME_HR";
        case TILEDB_DATETIME_MIN:
            return "DATETIME_MIN";
        case TILEDB_DATETIME_SEC:
            return "DATETIME_SEC";
        case TILEDB_DATETIME_MS:
            return "DATETIME_MS";
        case TILEDB_DATETIME_US:
            return "DATETIME_US";
        case TILEDB_DATETIME_NS:
            return "DATETIME_NS";
        case TILEDB_DATETIME_PS:
            return "DATETIME_PS";
        case TILEDB_DATETIME_FS:
            return "DATETIME_FS";
        case TILEDB_DATETIME_AS:
            return "DATETIME_AS";
        case TILEDB_BLOB:
            return "BLOB";
        case TILEDB_BOOL:
            return "BOOL";
        case TILEDB_GEOM_WKB:
            return "GEOM_WKB";
        case TILEDB_GEOM_WKT:
            return "GEOM_WKT";
        default:
            Rcpp::stop("unknown tiledb_datatype_t (%d)", dtype);
    }
}

// identify ncells
// taken from tiledb-r
// https://github.com/TileDB-Inc/TileDB-R/blob/525bdfc0f34aadb74a312a5d8428bd07819a8f83/src/libtiledb.cpp#L1590-L1599
template <typename AttrOrDim>
int _get_ncells(AttrOrDim x) {
    int ncells;
    if (x->cell_val_num() == TILEDB_VAR_NUM) {
        ncells = R_NaInt;
    } else if (x->cell_val_num() > std::numeric_limits<int32_t>::max()) {
        Rcpp::stop("tiledb_attr ncells value not representable as an R integer");
    } else {
        ncells = static_cast<int32_t>(x->cell_val_num());
    }
    return ncells;
}

// [[Rcpp::export]]
Rcpp::List c_attributes(Rcpp::XPtr<tiledbsoma::SOMAArray> array) {
    if (!array) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    auto sch = array->tiledb_schema();

    Rcpp::List result;
    int nattr = sch->attribute_num();
    for (auto i = 0; i < nattr; i++) {
        auto attr = make_xptr<tiledb::Attribute>(new tiledb::Attribute(sch->attribute(i)));
        auto name = attr->name();

        // identify the filters
        // filter options taken from tiledb-r
        // https://github.com/TileDB-Inc/TileDB-R/blob/525bdfc0f34aadb74a312a5d8428bd07819a8f83/src/libtiledb.cpp#L369-L388
        Rcpp::List filters;
        auto filter_list = make_xptr<tiledb::FilterList>(new tiledb::FilterList(attr->filter_list()));
        int nfilters = static_cast<int32_t>(filter_list->nfilters());
        for (auto j = 0; j < nfilters; j++) {
            auto filter = make_xptr<tiledb::Filter>(new tiledb::Filter(filter_list->filter(j)));
            auto filter_type = tiledb::Filter::to_str(filter->filter_type());
            filters[filter_type] = _get_filter_options(filter);
        }

        // assemble the attribute information list
        result[name] = Rcpp::List::create(
            Rcpp::Named("name") = name,
            Rcpp::Named("type") = _tiledb_datatype_to_string(attr->type()),
            Rcpp::Named("ncells") = _get_ncells<Rcpp::XPtr<tiledb::Attribute>>(attr),
            Rcpp::Named("nullable") = attr->nullable(),
            Rcpp::Named("filter_list") = filters);
    }

    return result;
}

// adapted from tiledb-r
// https://github.com/TileDB-Inc/TileDB-R/blob/525bdfc0f34aadb74a312a5d8428bd07819a8f83/src/libtiledb.cpp#L3027-L3042
// [[Rcpp::export]]
Rcpp::LogicalVector c_attributes_enumerated(Rcpp::XPtr<tiledbsoma::SOMADataFrame> dataframe) {
    if (!dataframe) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    auto sch = dataframe->tiledb_schema();

    int nattrs = sch->attribute_num();
    Rcpp::LogicalVector has_enum = Rcpp::LogicalVector(nattrs);
    Rcpp::CharacterVector names = Rcpp::CharacterVector(nattrs);
    for (int i = 0; i < nattrs; i++) {
        auto attr = sch->attribute(i);
        auto enmr = tiledb::AttributeExperimental::get_enumeration_name(*(dataframe->ctx()->tiledb_ctx()), attr);
        has_enum(i) = enmr != std::nullopt;
        names(i) = attr.name();
    }

    has_enum.attr("names") = names;
    return has_enum;
}

// [[Rcpp::export]]
Rcpp::CharacterVector c_attribute_enumeration_levels(
    Rcpp::XPtr<tiledbsoma::SOMADataFrame> dataframe, const std::string& name) {
    if (!dataframe) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    std::pair<ArrowArray*, ArrowSchema*> enum_values = dataframe->get_enumeration_values_for_column(name);

    if (enum_values.first->length > std::numeric_limits<int32_t>::max()) {
        Rcpp::stop("too many enumeration levels for R");
    }

    nanoarrow::UniqueArrayView enum_view;
    ArrowArrayViewInitFromType(enum_view.get(), NANOARROW_TYPE_LARGE_STRING);
    NANOARROW_RETURN_NOT_OK(ArrowArrayViewSetArray(enum_view.get(), enum_values.first, nullptr));

    int nlevels = static_cast<int32_t>(enum_values.first->length);
    Rcpp::CharacterVector enumerations = Rcpp::CharacterVector(nlevels);
    for (int i = 0; i < nlevels; i++) {
        if (ArrowArrayViewIsNull(enum_view.get(), i)) {
            enumerations(i) = R_NaString;
        } else {
            ArrowStringView item = ArrowArrayViewGetStringUnsafe(enum_view.get(), i);
            enumerations(i) = std::string(item.data, item.size_bytes);
        }
    }

    return enumerations;
}

// [[Rcpp::export]]
Rcpp::List c_domain(Rcpp::XPtr<tiledbsoma::SOMAArray> array) {
    if (!array) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    auto sch = array->tiledb_schema();

    Rcpp::List result;
    auto domain = make_xptr<tiledb::Domain>(new tiledb::Domain(sch->domain()));
    // adapted from tiledb-r
    // https://github.com/TileDB-Inc/TileDB-R/blob/525bdfc0f34aadb74a312a5d8428bd07819a8f83/src/libtiledb.cpp#L1262C3-L1264C4
    uint32_t rank = domain->ndim();
    if (rank > std::numeric_limits<int32_t>::max()) {
        Rcpp::stop("tiledb::Domain rank is not representable by an R integer");
    }
    int ndim = static_cast<int32_t>(rank);
    for (int i = 0; i < ndim; i++) {
        auto dim = make_xptr<tiledb::Dimension>(new tiledb::Dimension(domain->dimension(i)));
        auto name = dim->name();

        // identify the filters
        Rcpp::List filters;
        auto filter_list = make_xptr<tiledb::FilterList>(new tiledb::FilterList(dim->filter_list()));
        int nfilters = static_cast<int32_t>(filter_list->nfilters());
        for (auto j = 0; j < nfilters; j++) {
            auto filter = make_xptr<tiledb::Filter>(new tiledb::Filter(filter_list->filter(j)));
            auto filter_type = tiledb::Filter::to_str(filter->filter_type());
            filters[filter_type] = _get_filter_options(filter);
        }

        // assemble the dimension information list
        result[name] = Rcpp::List::create(
            Rcpp::Named("name") = name,
            Rcpp::Named("type") = _tiledb_datatype_to_string(dim->type()),
            Rcpp::Named("ncells") = _get_ncells<Rcpp::XPtr<tiledb::Dimension>>(dim),
            Rcpp::Named("domain") = _get_dim_domain(dim),
            Rcpp::Named("tile") = _get_dim_tile(dim),
            Rcpp::Named("filters") = filters);
    }

    return result;
}

// [[Rcpp::export]]
std::string resize(
    Rcpp::XPtr<tiledbsoma::SOMAArray> ndarray,
    Rcpp::NumericVector new_shape,
    std::string function_name_for_messages,
    bool check_only) {
    if (!ndarray) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    // This function is solely for SparseNDArray and DenseNDArray for which the
    // dims are required by the SOMA spec to be of type int64. Domain-resize for
    // variant-indexed dataframes is via upgrade_domain and change_domain.
    std::vector<int64_t> new_shape_i64 = i64_from_rcpp_numeric(new_shape);

    std::string retval = "";
    if (check_only) {
        auto status_and_reason = ndarray->can_resize(new_shape_i64, function_name_for_messages);
        retval = status_and_reason.second;
    } else {
        ndarray->resize(new_shape_i64, function_name_for_messages);
    }
    return retval;
}

// [[Rcpp::export]]
void resize_soma_joinid_shape(
    Rcpp::XPtr<tiledbsoma::SOMADataFrame> dataframe,
    Rcpp::NumericVector new_shape,
    std::string function_name_for_messages) {
    if (!dataframe) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    std::vector<int64_t> new_shape_i64 = i64_from_rcpp_numeric(new_shape);
    dataframe->resize_soma_joinid_shape(new_shape_i64[0], function_name_for_messages);
}

// [[Rcpp::export]]
std::string tiledbsoma_upgrade_shape(
    Rcpp::XPtr<tiledbsoma::SOMAArray> ndarray,
    Rcpp::NumericVector new_shape,
    std::string function_name_for_messages,
    bool check_only) {
    if (!ndarray) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    std::vector<int64_t> new_shape_i64 = i64_from_rcpp_numeric(new_shape);

    std::string retval = "";
    if (check_only) {
        auto status_and_reason = ndarray->can_upgrade_shape(new_shape_i64, function_name_for_messages);
        retval = status_and_reason.second;
    } else {
        ndarray->upgrade_shape(new_shape_i64, function_name_for_messages);
    }

    return retval;
}

// [[Rcpp::export]]
std::string upgrade_or_change_domain(
    Rcpp::XPtr<tiledbsoma::SOMADataFrame> dataframe,
    bool is_change_domain,
    naxpArray nadimap,
    naxpSchema nadimsp,
    std::string function_name_for_messages,
    bool check_only) {
    if (!dataframe) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    // This is pointer manipulation from R -> Rcpp SEXP -> libtiledbsoma:
    nanoarrow::UniqueArray apdim{nanoarrow_array_from_xptr(nadimap)};
    nanoarrow::UniqueSchema spdim{nanoarrow_schema_from_xptr(nadimsp)};

    auto dimarr = tdbs::common::arrow::make_managed_unique<ArrowArray>();
    auto dimsch = tdbs::common::arrow::make_managed_unique<ArrowSchema>();

    apdim.move(dimarr.get());
    spdim.move(dimsch.get());

    tdbs::common::arrow::ArrowTable arrow_table(std::move(dimarr), std::move(dimsch));

    // Now call libtiledbsoma
    std::string reason_string = "";
    if (is_change_domain) {
        if (check_only) {
            auto status_and_reason = dataframe->can_change_domain(arrow_table, function_name_for_messages);
            reason_string = status_and_reason.second;
        } else {
            dataframe->change_domain(arrow_table, function_name_for_messages);
        }
    } else {
        if (check_only) {
            auto status_and_reason = dataframe->can_upgrade_domain(arrow_table, function_name_for_messages);
            reason_string = status_and_reason.second;
        } else {
            dataframe->upgrade_domain(arrow_table, function_name_for_messages);
        }
    }
    return reason_string;
}

// [[Rcpp::export]]
void c_update_dataframe_schema(
    Rcpp::XPtr<tiledbsoma::SOMADataFrame> dataframe,
    Rcpp::CharacterVector column_names_to_drop,
    Rcpp::List add_cols_types,
    Rcpp::List add_cols_enum_value_types,
    Rcpp::List add_cols_enum_ordered) {
    if (!dataframe) {
        Rcpp::exception("Internal error: SOMAObject handle is not initialized.");
    }
    // Drop columns is just a list of column names: it goes right through
    // from R to C++.
    std::vector<std::string> drop_attrs = Rcpp::as<std::vector<std::string>>(column_names_to_drop);

    // For add columns: coming from R we have a named list from attr name to:
    // * for non-enum attrs: the datatype of the attr
    // * for enum attrs: the index type (e.g. int8) of the enumeration attr
    std::map<std::string, std::string> add_attrs;
    int n_add = add_cols_types.length();

    if (n_add > 0) {
        // Calling .names on empty list results in:
        // Not compatible with STRSXP: [type=NULL].Abort trap: 6
        Rcpp::CharacterVector add_col_names = add_cols_types.names();
        for (int i = 0; i < n_add; i++) {
            std::string type_name = Rcpp::as<std::string>(add_cols_types[i]);

            // Map type names like "int8" in the R Arrow API to type names like
            // "c" in the C NanoArrow API. I looked and didn't find an R
            // accessor for this; no big deal. Here we remap what we know about,
            // and if there's still an unrecognized type name (which is
            // developer error, not user error), we will let libtiledbsoma
            // throw.
            type_name = remap_arrow_type_code_r_to_c(type_name);

            add_attrs.emplace(add_col_names[i], type_name);
        }
    }

    // For enum columns, two more things: value type (e.g. string) and
    // is-ordered-enum boolean. These come into us as separate lists but
    // we will reshape them into a map from enum-attr name to tuple of
    // (value_type, ordered).

    // First do integrity checks.
    if (add_cols_enum_value_types.length() != add_cols_enum_ordered.length()) {
        // This isn't user error
        throw Rcpp::exception("c_update_dataframe_schema: internal coding error");
    }

    std::map<std::string, std::pair<std::string, bool>> add_enmrs;
    int n_add_enum = add_cols_enum_value_types.length();
    if (n_add_enum > 0) {
        // Calling .names on empty list results in:
        // Not compatible with STRSXP: [type=NULL].Abort trap: 6
        Rcpp::CharacterVector add_enum_col_names = add_cols_enum_value_types.names();
        Rcpp::CharacterVector other_names = add_cols_enum_ordered.names();
        for (int i = 0; i < n_add_enum; i++) {
            if (add_enum_col_names[i] != other_names[i]) {
                // This also isn't user error
                throw Rcpp::exception("c_update_dataframe_schema: internal coding error");
            }
        }

        for (int i = 0; i < n_add_enum; i++) {
            std::string key = Rcpp::as<std::string>(add_enum_col_names[i]);
            std::string type_name = Rcpp::as<std::string>(add_cols_enum_value_types[i]);
            type_name = remap_arrow_type_code_r_to_c(type_name);
            bool ordered = Rcpp::as<bool>(add_cols_enum_ordered[i]);
            add_enmrs.emplace(key, std::make_pair(type_name, ordered));
        }
    }

    dataframe->update_dataframe_schema(drop_attrs, add_attrs, add_enmrs);
}
