#include <Rcpp.h>                   // for R interface to C++
#include <nanoarrow/r.h>            // for C interface to Arrow (via R package)
#include <RcppInt64>                // for fromInteger64
#include <nanoarrow/nanoarrow.hpp>  // for C/C++ interface to Arrow

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

    spdl::info("[soma_array_reader] Reading from {}", uri);

    std::vector<std::string> column_names = {};
    if (!colnames.isNull()) {  // If we have column names, select them
        column_names = Rcpp::as<std::vector<std::string>>(colnames);
        spdl::debug(
            "[soma_array_reader] Selecting {} columns", column_names.size());
    }

    auto tdb_result_order = get_tdb_result_order(result_order);

    // optional timestamp range
    std::optional<tdbs::TimestampRange> tsrng = makeTimestampRange(
        timestamprange);
    if (timestamprange.isNotNull()) {
        Rcpp::DatetimeVector vec(timestamprange);
        spdl::debug(
            "[soma_array_reader] timestamprange ({},{})", vec[0], vec[1]);
    }

    // Read selected columns from the uri (return is unique_ptr<SOMAArray>)
    auto sr = tdbs::SOMAArray::open(
        OpenMode::read,
        uri,
        somactx,
        "unnamed",  // name parameter could be added
        column_names,
        batch_size,
        tdb_result_order,
        tsrng);

    std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>>
        name2dim;
    std::shared_ptr<tiledb::ArraySchema> schema = sr->tiledb_schema();
    tiledb::Domain domain = schema->domain();
    std::vector<tiledb::Dimension> dims = domain.dimensions();
    for (auto& dim : dims) {
        spdl::info(
            "[soma_array_reader] Dimension {} type {} domain {} extent {}",
            dim.name(),
            tiledb::impl::to_str(dim.type()),
            dim.domain_to_str(),
            dim.tile_extent_to_str());
        name2dim.emplace(std::make_pair(
            dim.name(), std::make_shared<tiledb::Dimension>(dim)));
    }

    // If we have a query condition, apply it
    if (!qc.isNull()) {
        spdl::info("[soma_array_reader_impl] Applying query condition");
        Rcpp::XPtr<tiledb::QueryCondition> qcxp(qc);
        sr->set_condition(*qcxp);
    }

    // If we have dimension points, apply them
    // The interface is named list, where each (named) list elements is one
    // (named) dimesion The List element is a simple vector of points and each
    // point is applied to the named dimension
    if (!dim_points.isNull()) {
        Rcpp::List lst(dim_points);
        apply_dim_points(sr.get(), name2dim, lst);
    }

    // If we have a dimension points, apply them
    if (!dim_ranges.isNull()) {
        Rcpp::List lst(dim_ranges);
        apply_dim_ranges(sr.get(), name2dim, lst);
    }

    // Getting next batch:  std::optional<std::shared_ptr<ArrayBuffers>>
    auto sr_data = sr->read_next();
    if (!sr->results_complete()) {
        Rcpp::stop(
            "Read of '%s' is incomplete.\nConsider increasing the memory "
            "allocation via the configuration\noption "
            "'soma.init_buffer_bytes', "
            "or using iterated partial reads.",
            uri);
    }
    spdl::info(
        "[soma_array_reader] Read complete with {} rows and {} cols",
        sr_data->get()->num_rows(),
        sr_data->get()->names().size());

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
        spdl::info("[soma_array_reader] Accessing '{}' at pos {}", names[i], i);

        // now buf is a shared_ptr to ColumnBuffer
        auto buf = sr_data->get()->at(names[i]);

        // this is pair of array and schema pointer
        auto pp = tdbs::ArrowAdapter::to_arrow(buf);

        // memcpy((void*) sch->children[i], pp.second.get(),
        // sizeof(ArrowSchema)); memcpy((void*) arr->children[i],
        // pp.first.get(), sizeof(ArrowArray));
        ArrowArrayMove(pp.first.get(), arr->children[i]);
        ArrowSchemaMove(pp.second.get(), sch->children[i]);

        spdl::info(
            "[soma_array_reader] Incoming name {} length {}",
            std::string(pp.second->name),
            pp.first->length);

        if (pp.first->length > arr->length) {
            spdl::debug(
                "[soma_array_reader] Setting array length to {}",
                pp.first->length);
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
    auto sr_data = sr->read_next();
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
SEXP non_empty_domain_new(
    const std::string& uri, Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    auto sdf = tdbs::SOMADataFrame::open(uri, OpenMode::read, ctxxp->ctxptr);
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

    for (size_t i = 0; i < lib_retval->n_children; i++) {
        spdl::info(
            "[c_schema] Accessing name '{}' format '{}' at position {}",
            std::string(lib_retval->children[i]->name),
            std::string(lib_retval->children[i]->format),
            i);

        ArrowSchemaMove(lib_retval->children[i], sch->children[i]);
    }

    return schemaxp;
}

// [[Rcpp::export]]
void resize(
    const std::string& uri,
    Rcpp::NumericVector new_shape,
    Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // This function is solely for SparseNDArray and DenseNDArray for which the
    // dims are required by the SOMA spec to be of type int64. Domain-resize for
    // variant-indexed dataframes will be separate work as tracked on
    // https://github.com/single-cell-data/TileDB-SOMA/issues/2407.
    auto sr = tdbs::SOMAArray::open(OpenMode::write, uri, ctxxp->ctxptr);
    std::vector<int64_t> new_shape_i64 = i64_from_rcpp_numeric(new_shape);
    sr->resize(new_shape_i64);
    sr->close();
}

// [[Rcpp::export]]
void resize_soma_joinid(
    const std::string& uri,
    Rcpp::NumericVector new_shape,
    Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // This function is solely for SOMADataFrame.
    auto sr = tdbs::SOMADataFrame::open(uri, OpenMode::write, ctxxp->ctxptr);
    std::vector<int64_t> new_shape_i64 = i64_from_rcpp_numeric(new_shape);
    sr->resize_soma_joinid(new_shape_i64[0]);
    sr->close();
}

// [[Rcpp::export]]
void tiledbsoma_upgrade_shape(
    const std::string& uri,
    Rcpp::NumericVector new_shape,
    Rcpp::XPtr<somactx_wrap_t> ctxxp) {
    // This function is solely for SparseNDArray and DenseNDArray for which the
    // dims are required by the SOMA spec to be of type int64. Domain-resize for
    // variant-indexed dataframes will be separate work as tracked on
    // https://github.com/single-cell-data/TileDB-SOMA/issues/2407.
    auto sr = tdbs::SOMAArray::open(OpenMode::write, uri, ctxxp->ctxptr);
    std::vector<int64_t> new_shape_i64 = i64_from_rcpp_numeric(new_shape);
    sr->upgrade_shape(new_shape_i64);
    sr->close();
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
    // we will reshape them into a map from enum-attr name to pair of
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

            add_enmrs.emplace(key, std::pair(type_name, ordered));
        }
    }

    tdbs::SOMADataFrame::update_dataframe_schema(
        uri, ctxxp->ctxptr, drop_attrs, add_attrs, add_enmrs);
}
