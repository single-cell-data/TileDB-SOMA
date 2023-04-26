#include <Rcpp.h>               // for R interface to C++
#include <nanoarrow.h>          // for C interface to Arrow

// We get these via nanoarrow and must cannot include carrow.h again
#define ARROW_SCHEMA_AND_ARRAY_DEFINED 1
#include <tiledbsoma/tiledbsoma>
#if TILEDB_VERSION_MAJOR == 2 && TILEDB_VERSION_MINOR >= 4
#include <tiledb/tiledb_experimental>
#endif

#include "rutilities.h"         // local declarations
#include "xptr-utils.h"         // xptr taggging utilities

// (Adapted) helper functions from nanoarrow
//
// Create an external pointer with the proper class and that will release any
// non-null, non-released pointer when garbage collected. We use a tagged XPtr,
// but do not set an XPtr finalizer
Rcpp::XPtr<ArrowSchema> schema_owning_xptr(void) {
  struct ArrowSchema* schema = (struct ArrowSchema*)ArrowMalloc(sizeof(struct ArrowSchema));
  if (schema == NULL) Rcpp::stop("Failed to allocate ArrowSchema");
  schema->release = NULL;
  Rcpp::XPtr<ArrowSchema> schema_xptr = make_xptr(schema, false);
  return schema_xptr;
}
// Create an external pointer with the proper class and that will release any
// non-null, non-released pointer when garbage collected. We use a tagged XPtr,
// but do not set an XPtr finalizer
Rcpp::XPtr<ArrowArray> array_owning_xptr(void) {
  struct ArrowArray* array = (struct ArrowArray*)ArrowMalloc(sizeof(struct ArrowArray));
  if (array == NULL) Rcpp::stop("Failed to allocate ArrowArray");
  array->release = NULL;
  Rcpp::XPtr<ArrowArray> array_xptr = make_xptr(array, false);
  return array_xptr;
}

namespace tdbs = tiledbsoma;

//' Read SOMA Data From a Given URI
//'
//' This functions access a given SOMA URI and returns a complete data.frame. It does
//' not iterate; if your data is large than the initial read size consider the \code{sr_*}
//' functions.
//'
//' @param uri Character value with URI path to a SOMA data set
//' @param colnames Optional vector of character value with the name of the columns to retrieve
//' @param qc Optional external Pointer object to TileDB Query Condition, defaults to \sQuote{NULL} i.e.
//' no query condition
//' @param dim_points Optional named list with vector of data points to select on the given
//' dimension(s). Each dimension can be one entry in the list.
//' @param dim_ranges Optional named list with two-column matrix where each row select a range
//' for the given dimension. Each dimension can be one entry in the list.
//' @param batch_size Character value with the desired batch size, defaults to \sQuote{auto}
//' @param result_order Character value with the desired result order, defaults to \sQuote{auto}
//' @param loglevel Character value with the desired logging level, defaults to \sQuote{auto}
//' which lets prior setting prevail, any other value is set as new logging level.
//' @param config Optional character vector containing TileDB config.
//' @return A List object with two pointers to Arrow array data and schema is returned
//' @examples
//' \dontrun{
//' uri <- "test/soco/pbmc3k_processed/obs"
//' z <- soma_array_reader(uri)
//' tb <- as_arrow_table(z)
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List soma_array_reader(const std::string& uri,
                             Rcpp::Nullable<Rcpp::CharacterVector> colnames = R_NilValue,
                             Rcpp::Nullable<Rcpp::XPtr<tiledb::QueryCondition>> qc = R_NilValue,
                             Rcpp::Nullable<Rcpp::List> dim_points = R_NilValue,
                             Rcpp::Nullable<Rcpp::List> dim_ranges = R_NilValue,
                             std::string batch_size = "auto",
                             std::string result_order = "auto",
                             const std::string& loglevel = "auto",
                             Rcpp::Nullable<Rcpp::CharacterVector> config = R_NilValue) {

    if (loglevel != "auto") {
        spdl::set_level(loglevel);
        tdbs::LOG_SET_LEVEL(loglevel);
    }

    spdl::info("[soma_array_reader] Reading from {}", uri);

    std::map<std::string, std::string> platform_config = config_vector_to_map(config);

    std::vector<std::string> column_names = {};
    if (!colnames.isNull()) {    // If we have column names, select them
        column_names = Rcpp::as<std::vector<std::string>>(colnames);
        spdl::debug("[soma_array_reader] Selecting {} columns", column_names.size());
    }

    // Read selected columns from the uri (return is unique_ptr<SOMAArray>)
    auto sr = tdbs::SOMAArray::open(uri,
                                          "unnamed",         // name parameter could be added
                                          platform_config,   // to add, done in iterated reader
                                          column_names,
                                          batch_size,
                                          result_order);

    std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>> name2dim;
    std::shared_ptr<tiledb::ArraySchema> schema = sr->schema();
    tiledb::Domain domain = schema->domain();
    std::vector<tiledb::Dimension> dims = domain.dimensions();
    for (auto& dim: dims) {
        spdl::info("[soma_array_reader] Dimension {} type {} domain {} extent {}",
                   dim.name(), tiledb::impl::to_str(dim.type()),
                   dim.domain_to_str(), dim.tile_extent_to_str());
        name2dim.emplace(std::make_pair(dim.name(), std::make_shared<tiledb::Dimension>(dim)));
    }

    // If we have a query condition, apply it
    if (!qc.isNull()) {
        spdl::info("[soma_array_reader] Applying query condition");
        Rcpp::XPtr<tiledb::QueryCondition> qcxp(qc);
        sr->set_condition(*qcxp);
    }

    // If we have dimension points, apply them
    // The interface is named list, where each (named) list elements is one (named) dimesion
    // The List element is a simple vector of points and each point is applied to the named dimension
    if (!dim_points.isNull()) {
        Rcpp::List lst(dim_points);
        apply_dim_points(sr.get(), name2dim, lst);
    }

    // If we have a dimension points, apply them
    if (!dim_ranges.isNull()) {
        Rcpp::List lst(dim_ranges);
        apply_dim_ranges(sr.get(), name2dim, lst);
    }

    sr->submit();

    // Getting next batch:  std::optional<std::shared_ptr<ArrayBuffers>>
    auto sr_data = sr->read_next();
    if (!sr->results_complete()) {
        Rcpp::stop("Read of '%s' is incomplete.\nConsider increasing the memory "
                   "allocation via the configuration\noption 'soma.init_buffer_bytes', "
                   "or using iterated partial reads.", uri);
    }
    spdl::info("[soma_array_reader] Read complete with {} rows and {} cols",
               sr_data->get()->num_rows(), sr_data->get()->names().size());

    const std::vector<std::string> names = sr_data->get()->names();
    auto ncol = names.size();
    Rcpp::XPtr<ArrowSchema> schemaxp = schema_owning_xptr();
    Rcpp::XPtr<ArrowArray> arrayxp = array_owning_xptr();
    ArrowSchemaInitFromType((ArrowSchema*)R_ExternalPtrAddr(schemaxp), NANOARROW_TYPE_STRUCT);
    ArrowSchemaAllocateChildren((ArrowSchema*)R_ExternalPtrAddr(schemaxp), ncol);
    ArrowArrayInitFromType((ArrowArray*)R_ExternalPtrAddr(arrayxp), NANOARROW_TYPE_STRUCT);
    ArrowArrayAllocateChildren((ArrowArray*)R_ExternalPtrAddr(arrayxp), ncol);

    arrayxp->length = 0;

    for (size_t i=0; i<ncol; i++) {
        // this allocates, and properly wraps as external pointers controlling lifetime
        Rcpp::XPtr<ArrowSchema> chldschemaxp = schema_owning_xptr();
        Rcpp::XPtr<ArrowArray> chldarrayxp = array_owning_xptr();

        spdl::info("[soma_array_reader] Accessing {} at {}", names[i], i);

        // now buf is a shared_ptr to ColumnBuffer
        auto buf = sr_data->get()->at(names[i]);

        // this is pair of array and schema pointer
        auto pp = tdbs::ArrowAdapter::to_arrow(buf);

        memcpy((void*) chldschemaxp, pp.second.get(), sizeof(ArrowSchema));
        memcpy((void*) chldarrayxp, pp.first.get(), sizeof(ArrowArray));

        spdl::info("[soma_array_reader] Incoming name {} length {}", std::string(pp.second->name), pp.first->length);

        schemaxp->children[i] = chldschemaxp;
        arrayxp->children[i] = chldarrayxp;

        if (pp.first->length > arrayxp->length) {
            spdl::debug("[soma_array_reader] Setting array length to {}", pp.first->length);
            arrayxp->length = pp.first->length;
        }
    }

    Rcpp::List as = Rcpp::List::create(Rcpp::Named("array_data") = arrayxp,
                                       Rcpp::Named("schema") = schemaxp);
    return as;
}

//' @noRd
// [[Rcpp::export]]
void set_log_level(const std::string& level) {
    spdl::set_level(level);
    tdbs::LOG_SET_LEVEL(level);
}

//' @noRd
// [[Rcpp::export]]
Rcpp::CharacterVector get_column_types(const std::string& uri,
                                       const std::vector<std::string>& colnames) {

    auto sr = tdbs::SOMAArray::open(uri);
    sr->submit();
    auto sr_data = sr->read_next();
    size_t n = colnames.size();
    Rcpp::CharacterVector vs(n);
    for (size_t i=0; i<n; i++) {
        auto datatype = sr_data->get()->at(colnames[i])->type();
        vs[i] = std::string(tiledb::impl::to_str(datatype));
    }
    vs.attr("names") = colnames;
    return vs;
}

//' @rdname soma_array_reader
//' @export
// [[Rcpp::export]]
double nnz(const std::string& uri, Rcpp::Nullable<Rcpp::CharacterVector> config = R_NilValue) {
    auto sr = tdbs::SOMAArray::open(uri, "unnamed", config_vector_to_map(config));
    return static_cast<double>(sr->nnz());
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

//' @rdname soma_array_reader
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector shape(const std::string& uri,
                          Rcpp::Nullable<Rcpp::CharacterVector> config = R_NilValue) {
    auto sr = tdbs::SOMAArray::open(uri, "unnamed", config_vector_to_map(Rcpp::wrap(config)));
    return makeInteger64(sr->shape());
}
