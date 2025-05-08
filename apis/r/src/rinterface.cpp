#include <Rcpp.h>                   // for R interface to C++
#include <nanoarrow/r.h>            // for C interface to Arrow (via R package)
#include <RcppInt64>                // for fromInteger64
#include <nanoarrow/nanoarrow.hpp>  // for C/C++ interface to Arrow
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
    struct ArrowSchema* schema = (struct ArrowSchema*)ArrowMalloc(
        sizeof(struct ArrowSchema));
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
    struct ArrowArray* array = (struct ArrowArray*)ArrowMalloc(
        sizeof(struct ArrowArray));
    if (array == NULL)
        Rcpp::stop("Failed to allocate ArrowArray");
    array->release = NULL;
    Rcpp::XPtr<ArrowArray> array_xptr = make_xptr(array, false);
    return array_xptr;
}

namespace tdbs = tiledbsoma;

//' @noRd
// [[Rcpp::export(soma_array_reader_impl)]]
SEXP soma_array_reader(
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
        spdl::set_level(loglevel);
        tdbs::LOG_SET_LEVEL(loglevel);
    }

    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> somactx = ctxxp->ctxptr;

    tdbs::LOG_INFO(fmt::format("[soma_array_reader] Reading from {}", uri));

    std::vector<std::string> column_names = {};
    if (!colnames.isNull()) {  // If we have column names, select them
        column_names = Rcpp::as<std::vector<std::string>>(colnames);
        tdbs::LOG_DEBUG(fmt::format("[soma_array_reader] Selecting {} columns", column_names.size()));
    }

    auto tdb_result_order = get_tdb_result_order(result_order);

    // optional timestamp range
    std::optional<tdbs::TimestampRange> tsrng = makeTimestampRange(
        timestamprange);
    if (timestamprange.isNotNull()) {
        Rcpp::DatetimeVector vec(timestamprange);
        tdbs::LOG_DEBUG(fmt::format("[soma_array_reader] timestamprange ({},{})", vec[0], vec[1]));
    }

    // Read selected columns from the uri (return is unique_ptr<SOMAArray>)
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, somactx, tsrng);

    auto mq = tdbs::ManagedQuery(*sr, somactx->tiledb_ctx(), "unnamed");
    mq.set_layout(tdb_result_order);
    if(!column_names.empty()){
        mq.select_columns(column_names);
    }

    std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>>
        name2dim;
    std::shared_ptr<tiledb::ArraySchema> schema = sr->tiledb_schema();
    tiledb::Domain domain = schema->domain();
    std::vector<tiledb::Dimension> dims = domain.dimensions();
    for (auto& dim : dims) {
        tdbs::LOG_INFO(fmt::format("[soma_array_reader] Dimension {} type {} domain {} extent {}",
            dim.name(),
            tiledb::impl::to_str(dim.type()),
            dim.domain_to_str(),
            dim.tile_extent_to_str()));
        name2dim.emplace(std::make_pair(
            dim.name(), std::make_shared<tiledb::Dimension>(dim)));
    }

    // If we have a query condition, apply it
    if (!qc.isNull()) {
        tdbs::LOG_INFO("[soma_array_reader_impl] Applying query condition");
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
            "Read of '%s' is incomplete.\nConsider increasing the memory "
            "allocation via the configuration\noption "
            "'soma.init_buffer_bytes', "
            "or using iterated partial reads.",
            uri);
    }
    tdbs::LOG_INFO(fmt::format("[soma_array_reader] Read complete with {} rows and {} cols",
        sr_data->get()->num_rows(),
        sr_data->get()->names().size()));

    const std::vector<std::string> names = sr_data->get()->names();
    auto ncol = names.size();
    // Schema first
    auto schemaxp = nanoarrow_schema_owning_xptr();
    auto sch = nanoarrow_output_schema_from_xptr(schemaxp);
    exitIfError(
        ArrowSchemaInitFromType(sch, NANOARROW_TYPE_STRUCT), "Bad schema init");
    exitIfError(ArrowSchemaSetName(sch, ""), "Bad schema name");
    exitIfError(
        ArrowSchemaAllocateChildren(sch, ncol), "Bad schema children alloc");

    // Array second
    auto arrayxp = nanoarrow_array_owning_xptr();
    auto arr = nanoarrow_output_array_from_xptr(arrayxp);
    exitIfError(
        ArrowArrayInitFromType(arr, NANOARROW_TYPE_STRUCT), "Bad array init");
    exitIfError(
        ArrowArrayAllocateChildren(arr, ncol), "Bad array children alloc");

    arr->length = 0;  // initial value
    for (size_t i = 0; i < ncol; i++) {
        tdbs::LOG_INFO(fmt::format("[soma_array_reader] Accessing '{}' at pos {}", names[i], i));

        // now buf is a shared_ptr to ColumnBuffer
        auto buf = sr_data->get()->at(names[i]);

        // this is pair of array and schema pointer
        auto pp = tdbs::ArrowAdapter::to_arrow(buf);

        // memcpy((void*) sch->children[i], pp.second.get(),
        // sizeof(ArrowSchema)); memcpy((void*) arr->children[i],
        // pp.first.get(), sizeof(ArrowArray));
        ArrowArrayMove(pp.first.get(), arr->children[i]);
        ArrowSchemaMove(pp.second.get(), sch->children[i]);
        tdbs::LOG_INFO(fmt::format("[soma_array_reader] Incoming name {} length {}",
            std::string(pp.second->name),
            pp.first->length));

        if (pp.first->length > arr->length) {
            tdbs::LOG_DEBUG(fmt::format("[soma_array_reader] Setting array length to {}",
                pp.first->length));
            arr->length = pp.first->length;
        }
    }

    // Nanoarrow special: stick schema into xptr tag to return single SEXP
    array_xptr_set_schema(arrayxp, schemaxp);  // embed schema in array
    sr->close();
    return arrayxp;
}

//' Set the logging level for the R package and underlying C++ library
//'
//' @param level A character value with logging level understood by
//' \sQuote{spdlog} such as \dQuote{trace}, \dQuote{debug}, \dQuote{info}, or
//' \dQuote{warn}.
//' @return Nothing is returned as the function is invoked for
//' the side-effect.
//' @export
// [[Rcpp::export]]
void set_log_level(const std::string& level) {
    spdl::setup("R", level);
    tdbs::LOG_SET_LEVEL(level);
}

//' @noRd
// [[Rcpp::export]]
Rcpp::CharacterVector get_column_types(
    const std::string& uri, const std::vector<std::string>& colnames) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri);
    auto mq = tdbs::ManagedQuery(*sr, sr->ctx()->tiledb_ctx());
    auto sr_data = mq.read_next();
    size_t n = colnames.size();
    Rcpp::CharacterVector vs(n);
    for (size_t i = 0; i < n; i++) {
        auto datatype = sr_data->get()->at(colnames[i])->type();
        vs[i] = std::string(tiledb::impl::to_str(datatype));
    }
    vs.attr("names") = colnames;
    sr->close();
    return vs;
}

// [[Rcpp::export]]
double nnz(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    auto retval = static_cast<double>(sr->nnz());
    sr->close();
    return retval;
}

//' @noRd
// [[Rcpp::export]]
bool check_arrow_schema_tag(Rcpp::XPtr<ArrowSchema> xp) {
    check_xptr_tag<ArrowSchema>(xp);  // throws if mismatched
    return true;
}

//' @noRd
// [[Rcpp::export]]
bool check_arrow_array_tag(Rcpp::XPtr<ArrowArray> xp) {
    check_xptr_tag<ArrowArray>(xp);  // throws if mismatched
    return true;
}

// [[Rcpp::export]]
Rcpp::NumericVector shape(
    const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    auto retval = Rcpp::toInteger64(sr->shape());
    sr->close();
    return retval;
}

// [[Rcpp::export]]
Rcpp::NumericVector maxshape(
    const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    auto retval = Rcpp::toInteger64(sr->maxshape());
    sr->close();
    return retval;
}

// [[Rcpp::export]]
SEXP non_empty_domain(
    const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sdf = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    tdbs::ArrowTable arrow_table = sdf->get_non_empty_domain();
    SEXP retval = convert_domainish(arrow_table);
    sdf->close();
    return retval;
}

// [[Rcpp::export]]
SEXP domain(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sdf = tdbs::SOMADataFrame::open(uri, OpenMode::read, ctxxp->ctxptr);
    tdbs::ArrowTable arrow_table = sdf->get_soma_domain();
    SEXP retval = convert_domainish(arrow_table);
    sdf->close();
    return retval;
}

// [[Rcpp::export]]
SEXP maxdomain(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sdf = tdbs::SOMADataFrame::open(uri, OpenMode::read, ctxxp->ctxptr);
    tdbs::ArrowTable arrow_table = sdf->get_soma_maxdomain();
    SEXP retval = convert_domainish(arrow_table);
    sdf->close();
    return retval;
}

// [[Rcpp::export]]
Rcpp::NumericVector maybe_soma_joinid_shape(
    const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // Pro-tip:
    // * Open with mode and uri gives a SOMAArray.
    // * Open with uri and mode gives a SOMADataFrame.
    // This was done intentionally to resolve an ambiguous-overload compiler
    // error. ^ Unsure. This is C++, and it is typed so member functions return
    // objects of their class.
    auto sr = tdbs::SOMADataFrame::open(uri, OpenMode::read, ctxxp->ctxptr);
    auto retval = sr->maybe_soma_joinid_shape();
    sr->close();
    if (retval.has_value()) {
        return Rcpp::toInteger64(retval.value());
    } else {
        return Rcpp::NumericVector::create(
            NA_REAL);  // one element vector, and is.na() is true
    }
}

// [[Rcpp::export]]
Rcpp::NumericVector maybe_soma_joinid_maxshape(
    const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMADataFrame::open(uri, OpenMode::read, ctxxp->ctxptr);
    auto retval = sr->maybe_soma_joinid_maxshape();
    sr->close();
    if (retval.has_value()) {
        return Rcpp::toInteger64(retval.value());
    } else {
        return Rcpp::NumericVector::create(
            NA_REAL);  // one element vector, and is.na() is true
    }
}

// [[Rcpp::export]]
Rcpp::LogicalVector has_current_domain(
    const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    auto retval = Rcpp::LogicalVector(sr->has_current_domain());
    sr->close();
    return retval;
}

// [[Rcpp::export]]
Rcpp::NumericVector ndim(
    const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    auto lib_retval = sr->ndim();
    sr->close();

    return Rcpp::NumericVector::create(lib_retval);
}

// [[Rcpp::export]]
Rcpp::CharacterVector c_dimnames(
    const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    auto lib_retval = sr->dimension_names();
    sr->close();

    size_t n = lib_retval.size();
    Rcpp::CharacterVector retval(n);
    for (size_t i = 0; i < n; i++) {
        retval[i] = lib_retval[i];
    }
    return retval;
}

// [[Rcpp::export]]
Rcpp::CharacterVector c_attrnames(
    const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    auto lib_retval = sr->attribute_names();
    sr->close();

    size_t n = lib_retval.size();
    Rcpp::CharacterVector retval(n);
    for (size_t i = 0; i < n; i++) {
        retval[i] = lib_retval[i];
    }
    return retval;
}

// [[Rcpp::export]]
SEXP c_schema(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    std::unique_ptr<ArrowSchema> lib_retval = sr->arrow_schema();
    sr->close();

    auto schemaxp = nanoarrow_schema_owning_xptr();
    auto sch = nanoarrow_output_schema_from_xptr(schemaxp);
    exitIfError(
        ArrowSchemaInitFromType(sch, NANOARROW_TYPE_STRUCT), "Bad schema init");
    exitIfError(ArrowSchemaSetName(sch, ""), "Bad schema name");
    exitIfError(
        ArrowSchemaAllocateChildren(sch, lib_retval->n_children),
        "Bad schema children alloc");

    for (size_t i = 0; i < static_cast<size_t>(lib_retval->n_children); i++) {
        tdbs::LOG_INFO(fmt::format("[c_schema] Accessing name '{}' format '{}' at position {}",
            std::string(lib_retval->children[i]->name),
            std::string(lib_retval->children[i]->format),
            i));

        ArrowSchemaMove(lib_retval->children[i], sch->children[i]);
    }

    return schemaxp;
}

// [[Rcpp::export]]
bool c_is_sparse(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    std::shared_ptr<tiledb::ArraySchema> sch = sr->tiledb_schema();
    sr->close();

    return sch->array_type() == TILEDB_SPARSE;
}

// [[Rcpp::export]]
bool c_allows_dups(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    std::shared_ptr<tiledb::ArraySchema> sch = sr->tiledb_schema();
    sr->close();

    return sch->allows_dups();
}

// [[Rcpp::export]]
double c_capacity(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    std::shared_ptr<tiledb::ArraySchema> sch = sr->tiledb_schema();
    sr->close();

    uint64_t cap = sch->capacity();
    return static_cast<double>(cap);
}

// [[Rcpp::export]]
std::string c_tile_order(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    std::shared_ptr<tiledb::ArraySchema> sch = sr->tiledb_schema();
    sr->close();

    auto order = sch->tile_order();
    return _tiledb_layout_to_string(order);
}

// [[Rcpp::export]]
std::string c_cell_order(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    std::shared_ptr<tiledb::ArraySchema> sch = sr->tiledb_schema();
    sr->close();

    auto order = sch->cell_order();
    return _tiledb_layout_to_string(order);
}

// [[Rcpp::export]]
Rcpp::List c_schema_filters(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
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
const char *_tiledb_datatype_to_string(tiledb_datatype_t dtype) {
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
Rcpp::List c_attributes(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    std::shared_ptr<tiledb::ArraySchema> sch = sr->tiledb_schema();
    sr->close();

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
            Rcpp::Named("filter_list") = filters
        );
    }

    return result;
}

// adapted from tiledb-r
// https://github.com/TileDB-Inc/TileDB-R/blob/525bdfc0f34aadb74a312a5d8428bd07819a8f83/src/libtiledb.cpp#L3027-L3042
// [[Rcpp::export]]
Rcpp::LogicalVector c_attributes_enumerated(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    std::shared_ptr<tiledb::ArraySchema> sch = sr->tiledb_schema();
    sr->close();

    int nattrs = sch->attribute_num();
    Rcpp::LogicalVector has_enum = Rcpp::LogicalVector(nattrs);
    Rcpp::CharacterVector names = Rcpp::CharacterVector(nattrs);
    for (int i = 0; i < nattrs; i++) {
        auto attr = make_xptr<tiledb::Attribute>(new tiledb::Attribute(sch->attribute(i)));
        auto enmr = tiledb::AttributeExperimental::get_enumeration_name(
            *(ctxxp->ctxptr->tiledb_ctx()),
            *attr.get()
        );
        has_enum(i) = enmr != std::nullopt;
        names(i) = attr->name();
    }

    has_enum.attr("names") = names;
    return has_enum;
}

// [[Rcpp::export]]
Rcpp::CharacterVector c_attribute_enumeration_levels(
        const std::string& uri,
        Rcpp::XPtr<somactx_wrap_t> ctxxp,
        const std::string& name
) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    std::pair<ArrowArray*, ArrowSchema*> enum_values = sr->get_enumeration_values_for_column(name);
    sr->close();

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
Rcpp::List c_domain(const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, ctxxp->ctxptr);
    std::shared_ptr<tiledb::ArraySchema> sch = sr->tiledb_schema();
    sr->close();

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
            Rcpp::Named("filters") = filters
        );
    }

    return result;
}

// [[Rcpp::export]]
std::string resize(
    const std::string& uri,
    Rcpp::NumericVector new_shape,
    std::string function_name_for_messages,
    bool check_only,
    Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // This function is solely for SparseNDArray and DenseNDArray for which the
    // dims are required by the SOMA spec to be of type int64. Domain-resize for
    // variant-indexed dataframes is via upgrade_domain and change_domain.
    auto sr = tdbs::SOMAArray::open(OpenMode::write, uri, ctxxp->ctxptr);
    std::vector<int64_t> new_shape_i64 = i64_from_rcpp_numeric(new_shape);

    std::string retval = "";
    if (check_only) {
        auto status_and_reason = sr->can_resize(
            new_shape_i64, function_name_for_messages);
        retval = status_and_reason.second;
    } else {
        sr->resize(new_shape_i64, function_name_for_messages);
    }

    sr->close();
    return retval;
}

// [[Rcpp::export]]
void resize_soma_joinid_shape(
    const std::string& uri,
    Rcpp::NumericVector new_shape,
    std::string function_name_for_messages,
    Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // This function is solely for SOMADataFrame.
    auto sr = tdbs::SOMADataFrame::open(uri, OpenMode::write, ctxxp->ctxptr);
    std::vector<int64_t> new_shape_i64 = i64_from_rcpp_numeric(new_shape);
    sr->resize_soma_joinid_shape(new_shape_i64[0], function_name_for_messages);
    sr->close();
}

// [[Rcpp::export]]
std::string tiledbsoma_upgrade_shape(
    const std::string& uri,
    Rcpp::NumericVector new_shape,
    std::string function_name_for_messages,
    bool check_only,
    Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // This function is solely for SparseNDArray and DenseNDArray for which the
    // dims are required by the SOMA spec to be of type int64. Domain-resize for
    // variant-indexed dataframes is via upgrade_domain and change_domain.
    auto sr = tdbs::SOMAArray::open(OpenMode::write, uri, ctxxp->ctxptr);
    std::vector<int64_t> new_shape_i64 = i64_from_rcpp_numeric(new_shape);

    std::string retval = "";
    if (check_only) {
        auto status_and_reason = sr->can_upgrade_shape(
            new_shape_i64, function_name_for_messages);
        retval = status_and_reason.second;
    } else {
        sr->upgrade_shape(new_shape_i64, function_name_for_messages);
    }

    sr->close();
    return retval;
}

// [[Rcpp::export]]
std::string upgrade_or_change_domain(
    const std::string& uri,
    bool is_change_domain,
    naxpArray nadimap,
    naxpSchema nadimsp,
    std::string function_name_for_messages,
    bool check_only,
    Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // This is pointer manipulation from R -> Rcpp SEXP -> libtiledbsoma:
    nanoarrow::UniqueArray apdim{nanoarrow_array_from_xptr(nadimap)};
    nanoarrow::UniqueSchema spdim{nanoarrow_schema_from_xptr(nadimsp)};

    auto dimarr = std::make_unique<ArrowArray>();
    auto dimsch = std::make_unique<ArrowSchema>();

    apdim.move(dimarr.get());
    spdim.move(dimsch.get());

    tdbs::ArrowTable arrow_table(std::move(dimarr), std::move(dimsch));

    // Now call libtiledbsoma
    std::string reason_string = "";
    auto sr = tdbs::SOMADataFrame::open(uri, OpenMode::write, ctxxp->ctxptr);
    if (is_change_domain) {
        if (check_only) {
            auto status_and_reason = sr->can_change_domain(
                arrow_table, function_name_for_messages);
            reason_string = status_and_reason.second;
        } else {
            sr->change_domain(arrow_table, function_name_for_messages);
        }
    } else {
        if (check_only) {
            auto status_and_reason = sr->can_upgrade_domain(
                arrow_table, function_name_for_messages);
            reason_string = status_and_reason.second;
        } else {
            sr->upgrade_domain(arrow_table, function_name_for_messages);
        }
    }
    sr->close();
    return reason_string;
}

// [[Rcpp::export]]
void c_update_dataframe_schema(
    const std::string& uri,
    Rcpp::XPtr<somactx_wrap_t> ctxxp,
    Rcpp::CharacterVector column_names_to_drop,
    Rcpp::List add_cols_types,
    Rcpp::List add_cols_enum_value_types,
    Rcpp::List add_cols_enum_ordered) {
    // Drop columns is just a list of column names: it goes right through
    // from R to C++.
    std::vector<std::string> drop_attrs = Rcpp::as<std::vector<std::string>>(
        column_names_to_drop);

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
        throw Rcpp::exception(
            "c_update_dataframe_schema: internal coding error");
    }

    std::map<std::string, std::pair<std::string, bool>> add_enmrs;
    int n_add_enum = add_cols_enum_value_types.length();
    if (n_add_enum > 0) {
        // Calling .names on empty list results in:
        // Not compatible with STRSXP: [type=NULL].Abort trap: 6
        Rcpp::CharacterVector add_enum_col_names = add_cols_enum_value_types
                                                       .names();
        Rcpp::CharacterVector other_names = add_cols_enum_ordered.names();
        for (int i = 0; i < n_add_enum; i++) {
            if (add_enum_col_names[i] != other_names[i]) {
                // This also isn't user error
                throw Rcpp::exception(
                    "c_update_dataframe_schema: internal coding error");
            }
        }

        for (int i = 0; i < n_add_enum; i++) {
            std::string key = Rcpp::as<std::string>(add_enum_col_names[i]);
            std::string type_name = Rcpp::as<std::string>(
                add_cols_enum_value_types[i]);
            type_name = remap_arrow_type_code_r_to_c(type_name);
            bool ordered = Rcpp::as<bool>(add_cols_enum_ordered[i]);
            add_enmrs.emplace(key, std::make_pair(type_name, ordered));
        }
    }

    auto sdf = tdbs::SOMADataFrame::open(uri, OpenMode::write, ctxxp->ctxptr);
    sdf->update_dataframe_schema(drop_attrs, add_attrs, add_enmrs);
    sdf->close();
}
