// Work in progress

#include <Rcpp.h>

#include <tiledb/tiledb>
#if TILEDB_VERSION_MAJOR == 2 && TILEDB_VERSION_MINOR >= 4
#include <tiledb/tiledb_experimental>
#endif

#include <tiledbsoma/tiledbsoma>
#include "pyapi/arrow_adapter.h"
#include <archAPI.h>
#include "rutilities.h"

namespace tdbs = tiledbsoma;

// enum for TileDB XPtr Object type using int32_t payload (for R)
enum tiledb_xptr_object : int32_t {};
const tiledb_xptr_object tiledb_xptr_default                     { 0 };
const tiledb_xptr_object tiledb_xptr_object_array                { 10 };
const tiledb_xptr_object tiledb_xptr_object_arrayschema          { 20 };
const tiledb_xptr_object tiledb_xptr_object_arrayschemaevolution { 30 };
const tiledb_xptr_object tiledb_xptr_object_attribute            { 40 };
const tiledb_xptr_object tiledb_xptr_object_config               { 50 };
const tiledb_xptr_object tiledb_xptr_object_context              { 60 };
const tiledb_xptr_object tiledb_xptr_object_dimension            { 70 };
const tiledb_xptr_object tiledb_xptr_object_domain               { 80 };
const tiledb_xptr_object tiledb_xptr_object_filter               { 90 };
const tiledb_xptr_object tiledb_xptr_object_filterlist           { 100 };
const tiledb_xptr_object tiledb_xptr_object_fragmentinfo         { 110 };
const tiledb_xptr_object tiledb_xptr_object_group                { 120 };
const tiledb_xptr_object tiledb_xptr_object_query                { 130 };
const tiledb_xptr_object tiledb_xptr_object_querycondition       { 140 };
const tiledb_xptr_object tiledb_xptr_object_vfs                  { 150 };
const tiledb_xptr_object tiledb_xptr_vfs_fh_t              		 { 160 };
const tiledb_xptr_object tiledb_xptr_vlc_buf_t                   { 170 };
const tiledb_xptr_object tiledb_xptr_vlv_buf_t                   { 180 };
const tiledb_xptr_object tiledb_xptr_query_buf_t                 { 190 };

// the definitions above are internal to tiledb-r but we need a new value here if we want tag the external pointer
const tiledb_xptr_object tiledb_soma_reader_t                    { 500 };

// templated checkers for external pointer tags
template <typename T> const int32_t XPtrTagType                            = tiledb_xptr_default; // clang++ wants a value
template <> inline const int32_t XPtrTagType<tiledb::Array>                = tiledb_xptr_object_array;
template <> inline const int32_t XPtrTagType<tiledb::ArraySchema>          = tiledb_xptr_object_arrayschema;
template <> inline const int32_t XPtrTagType<tiledb::ArraySchemaEvolution> = tiledb_xptr_object_arrayschemaevolution;
template <> inline const int32_t XPtrTagType<tiledb::Attribute>            = tiledb_xptr_object_attribute;
template <> inline const int32_t XPtrTagType<tiledb::Config>               = tiledb_xptr_object_config;
template <> inline const int32_t XPtrTagType<tiledb::Context>              = tiledb_xptr_object_context;
template <> inline const int32_t XPtrTagType<tiledb::Dimension>            = tiledb_xptr_object_dimension;
template <> inline const int32_t XPtrTagType<tiledb::Domain>               = tiledb_xptr_object_domain;
template <> inline const int32_t XPtrTagType<tiledb::Filter>               = tiledb_xptr_object_filter;
template <> inline const int32_t XPtrTagType<tiledb::FilterList>           = tiledb_xptr_object_filterlist;
template <> inline const int32_t XPtrTagType<tiledb::FragmentInfo>         = tiledb_xptr_object_fragmentinfo;
template <> inline const int32_t XPtrTagType<tiledb::Group>                = tiledb_xptr_object_group;
template <> inline const int32_t XPtrTagType<tiledb::Query>                = tiledb_xptr_object_query;
template <> inline const int32_t XPtrTagType<tiledb::QueryCondition>       = tiledb_xptr_object_query;
template <> inline const int32_t XPtrTagType<tiledb::VFS>                  = tiledb_xptr_object_vfs;
// this need the C API for which we do not include a header
// template <> inline const int32_t XPtrTagType<vfs_fh_t>                     = tiledb_xptr_vfs_fh_t;
// template <> inline const int32_t XPtrTagType<vlc_buf_t>                    = tiledb_xptr_vlc_buf_t;
// template <> inline const int32_t XPtrTagType<vlv_buf_t>                    = tiledb_xptr_vlv_buf_t;
// template <> inline const int32_t XPtrTagType<query_buf_t>                  = tiledb_xptr_query_buf_t;

template <> inline const int32_t XPtrTagType<tdbs::SOMAReader>             = tiledb_xptr_query_buf_t;

template <typename T> Rcpp::XPtr<T> make_xptr(T* p) {
    return Rcpp::XPtr<T>(p, true, Rcpp::wrap(XPtrTagType<T>), R_NilValue);
}

template <typename T> Rcpp::XPtr<T> make_xptr(SEXP p) {
    return Rcpp::XPtr<T>(p); 	// the default XPtr ctor with deleter on and tag and prot nil
}

template<typename T> void check_xptr_tag(Rcpp::XPtr<T> ptr) {
    if (R_ExternalPtrTag(ptr) == R_NilValue) {
        Rcpp::stop("External pointer without tag, expected tag %d\n", XPtrTagType<T>);
    }
    if (R_ExternalPtrTag(ptr) != R_NilValue) {
        int32_t tag = Rcpp::as<int32_t>(R_ExternalPtrTag(ptr));
        if (XPtrTagType<T> != tag) {
            Rcpp::stop("Wrong tag type: expected %d but received %d\n", XPtrTagType<T>, tag);
        }
    }
}

//' Iterator-Style Access to SOMA Array via SOMAReader
//'
//' The `sr_*` functions provide low-level access to an instance of the SOMAReader
//' class so that iterative access over parts of a (large) array is possible.
//' \describe{
//'   \item{\code{sr_setup}}{instantiates and by default also submits a query}
//'   \item{\code{sr_complete}}{checks is more data is available}
//'   \code{\code{sr_next}} returns the next chunk.
//' }
//'
//' @param ctx An external pointer to a TileDB Context object
//' @param uri Character value with URI path to a SOMA data set
//' @param loglevel Character value with the desired logging level, defaults to \sQuote{warn}
//' @param sr An external pointer to a TileDB SOMAReader object
//'
//' @return \code{sr_setup} returns an external pointer to a SOMAReader. \code{sr_complete}
//' returns a boolean, and \code{sr_next} returns an Arrow array helper object.
//'
//' @examples
//' \dontrun{
//' ctx <- tiledb_ctx()
//' sr <- sr_setup(ctx@ptr, uri, "warn")
//' rl <- data.frame()
//' while (!tiledbsoma:::sr_complete(sr)) {
//'     dat <- tiledbsoma:::sr_next(sr)
//'     dat |>
//'         arch::from_arch_array(arrow::RecordBatch) |>
//'         arrow::as_arrow_table() |>
//'         collect() |>
//'         as.data.frame() |>
//'         data.table() -> D
//'     rl <- rbind(rl, D)
//' }
//' summary(rl)
//' }
//' @export
// [[Rcpp::export]]
Rcpp::XPtr<tdbs::SOMAReader> sr_setup(Rcpp::XPtr<tiledb::Context> ctx,
                                      const std::string& uri,
                                      const std::string& loglevel = "warn") {
    check_xptr_tag<tiledb::Context>(ctx);
    tdbs::LOG_SET_LEVEL(loglevel);

    tdbs::LOG_INFO(fmt::format("[sr_setup] Setting up {}", uri));

    //std::map<std::string, std::string> platform_config;
    std::string_view name = "unnamed";
    std::vector<std::string> column_names = {};
    std::string_view batch_size = "auto";
    std::string_view result_order = "auto";

    tiledb::Config cfg{ctx.get()->config()}; // get config in order to make shared_ptr
    std::shared_ptr<tiledb::Context> ctxptr = std::make_shared<tiledb::Context>(cfg);

    auto ptr = new tdbs::SOMAReader(uri, name, ctxptr, column_names, batch_size, result_order);
    ptr->submit();
    Rcpp::XPtr<tdbs::SOMAReader> xptr = make_xptr<tdbs::SOMAReader>(ptr);
    return xptr;
}
//' @rdname sr_setup
//' @export
// [[Rcpp::export]]
bool sr_complete(Rcpp::XPtr<tdbs::SOMAReader> sr) {
   check_xptr_tag<tdbs::SOMAReader>(sr);
   tdbs::LOG_INFO(fmt::format("[sr_complete] Complete test is {}", sr->is_complete()));
   return sr->is_complete();
}

//' @rdname sr_setup
//' @export
// [[Rcpp::export]]
Rcpp::List sr_next(Rcpp::XPtr<tdbs::SOMAReader> sr) {
   check_xptr_tag<tdbs::SOMAReader>(sr);

   auto sr_data = sr->read_next();
   if (!sr->results_complete()) {
       tdbs::LOG_TRACE(fmt::format("[sr_next] Read is incomplete"));
   }
   tdbs::LOG_INFO(fmt::format("[sr_next] Read {} rows and {} cols",
                              sr_data->get()->num_rows(), sr_data->get()->names().size()));

   const std::vector<std::string> names = sr_data->get()->names();
   auto ncol = names.size();
   Rcpp::List schlst(ncol), arrlst(ncol);

   for (size_t i=0; i<ncol; i++) {
       // this allocates, and properly wraps as external pointers controlling lifetime
       SEXP schemaxp = arch_c_allocate_schema();
       SEXP arrayxp = arch_c_allocate_array_data();

       tdbs::LOG_TRACE(fmt::format("[sr_next] Accessing {} at {}", names[i], i));

       // now buf is a shared_ptr to ColumnBuffer
       auto buf = sr_data->get()->at(names[i]);

       // this is pair of array and schema pointer
       auto pp = tdbs::ArrowAdapter::to_arrow(buf);

       memcpy((void*) R_ExternalPtrAddr(schemaxp), pp.second.get(), sizeof(ArrowSchema));
       memcpy((void*) R_ExternalPtrAddr(arrayxp), pp.first.get(), sizeof(ArrowArray));

       schlst[i] = schemaxp;
       arrlst[i] = arrayxp;
   }

   struct ArrowArray* array_data_tmp = (struct ArrowArray*) R_ExternalPtrAddr(arrlst[0]);
   int rows = static_cast<int>(array_data_tmp->length);
   SEXP sxp = arch_c_schema_xptr_new(Rcpp::wrap("+s"), 	// format
                                     Rcpp::wrap(""),   	// name
                                     Rcpp::List(),       	// metadata
                                     Rcpp::wrap(2),      	// flags, 2 == unordered, nullable, no sorted map keys
                                     schlst, 	        	// children
                                     R_NilValue);        	// dictionary
   SEXP axp = arch_c_array_from_sexp(Rcpp::List::create(Rcpp::Named("")=R_NilValue), // buffers
                                     Rcpp::wrap(rows), 	// length
                                     Rcpp::wrap(-1), 	    // null count, -1 means not determined
                                     Rcpp::wrap(0),    	// offset (in bytes)
                                     arrlst,               // children
                                     R_NilValue);          // dictionary
   Rcpp::List as = Rcpp::List::create(Rcpp::Named("schema") = sxp,
                                      Rcpp::Named("array_data") = axp);
   as.attr("class") = "arch_array";
   return as;
}
