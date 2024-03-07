#include <Rcpp.h>               // for R interface to C++
#include <nanoarrow/r.h>        // for C interface to Arrow
#include <nanoarrow.hpp>        // for C/C++ interface to Arrow
#include <RcppInt64>            // for fromInteger64

// we currently get deprecation warnings by default which are noisy
#ifndef TILEDB_NO_API_DEPRECATION_WARNINGS
#define TILEDB_NO_API_DEPRECATION_WARNINGS
#endif

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

Rcpp::XPtr<ArrowSchema> schema_setup_struct(Rcpp::XPtr<ArrowSchema> schxp, int64_t n_children);
Rcpp::XPtr<ArrowArray> array_setup_struct(Rcpp::XPtr<ArrowArray> arrxp, int64_t n_children);


//' @noRd
// [[Rcpp::export(soma_array_reader_impl)]]
nanoarrowXPtr soma_array_reader(const std::string& uri,
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
    // to create a Context object:
    //    std::make_shared<Context>(Config(platform_config)),

    std::vector<std::string> column_names = {};
    if (!colnames.isNull()) {    // If we have column names, select them
        column_names = Rcpp::as<std::vector<std::string>>(colnames);
        spdl::debug("[soma_array_reader] Selecting {} columns", column_names.size());
    }

    auto tdb_result_order = get_tdb_result_order(result_order);

    // Read selected columns from the uri (return is unique_ptr<SOMAArray>)
    auto sr = tdbs::SOMAArray::open(OpenMode::read,
                                    uri,
                                    "unnamed",         // name parameter could be added
                                    platform_config,
                                    column_names,
                                    batch_size,
                                    tdb_result_order);

    std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>> name2dim;
    std::shared_ptr<tiledb::ArraySchema> schema = sr->tiledb_schema();
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

    arr->length = 0;             // initial value

    for (size_t i=0; i<ncol; i++) {
        spdl::info("[soma_array_reader] Accessing {} at {}", names[i], i);

        // now buf is a shared_ptr to ColumnBuffer
        auto buf = sr_data->get()->at(names[i]);

        // this is pair of array and schema pointer
        auto pp = tdbs::ArrowAdapter::to_arrow(buf);

        //memcpy((void*) sch->children[i], pp.second.get(), sizeof(ArrowSchema));
        //memcpy((void*) arr->children[i], pp.first.get(), sizeof(ArrowArray));
        ArrowArrayMove(pp.first.get(),  arr->children[i]);
        ArrowSchemaMove(pp.second.get(), sch->children[i]);

        spdl::info("[soma_array_reader] Incoming name {} length {}",
                   std::string(pp.second->name), pp.first->length);

        if (pp.first->length > arr->length) {
            spdl::debug("[soma_array_reader] Setting array length to {}", pp.first->length);
            arr->length = pp.first->length;
        }
    }

   // Nanoarrow special: stick schema into xptr tag to return single SEXP
   array_xptr_set_schema(arrayxp, schemaxp); 			// embed schema in array
   return arrayxp;
}

//' Set the logging level for the R package and underlying C++ library
//'
//' @param level A character value with logging level understood by \sQuote{spdlog}
//' such as \dQuote{trace}, \dQuote{debug}, \dQuote{info}, or \dQuote{warn}.
//' @return Nothing is returned as the function is invoked for the side-effect.
//' @export
// [[Rcpp::export]]
void set_log_level(const std::string& level) {
    spdl::set_level(level);
    tdbs::LOG_SET_LEVEL(level);
}

//' @noRd
// [[Rcpp::export]]
Rcpp::CharacterVector get_column_types(const std::string& uri,
                                       const std::vector<std::string>& colnames) {

    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri);
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

// [[Rcpp::export]]
double nnz(const std::string& uri, Rcpp::Nullable<Rcpp::CharacterVector> config = R_NilValue) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, "unnamed", config_vector_to_map(config));
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

// [[Rcpp::export]]
Rcpp::NumericVector shape(const std::string& uri,
                          Rcpp::Nullable<Rcpp::CharacterVector> config = R_NilValue) {
    auto sr = tdbs::SOMAArray::open(OpenMode::read, uri, "unnamed", config_vector_to_map(Rcpp::wrap(config)));
    return Rcpp::toInteger64(sr->shape());
}
