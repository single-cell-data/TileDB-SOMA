#include <Rcpp/Lighter>                         // for R interface to C++
#include <nanoarrow/r.h>                        // for C/C++ interface to Arrow (via header exported from the R package)
#include <nanoarrow/nanoarrow.hpp>              // for C/C++ interface to Arrow (vendored)
#include <RcppInt64>                            // for fromInteger64

#include <tiledbsoma/tiledbsoma>
#include <tiledbsoma/reindexer/reindexer.h>

#include "rutilities.h"         				// local declarations
#include "xptr-utils.h"         				// xptr taggging utilities

namespace tdbs = tiledbsoma;

void _show_content(const nanoarrow::UniqueArray& ap, const nanoarrow::UniqueSchema& sp) {
    int n = sp.get()->n_children;
    ArrowError ec;
    for (auto i=0; i<n; i++) {
        Rcpp::Rcout << " " << sp.get()->children[i]->name << " : ";

        nanoarrow::UniqueArrayView col;
        ArrowArrayViewInitFromSchema(col.get(), sp.get()->children[i], &ec);
        ArrowArrayViewSetArray(col.get(), ap.get()->children[i], &ec);

        int m = col.get()->length;
        for (auto j=0; j<m; j++)
            Rcpp::Rcout << ArrowArrayViewGetDoubleUnsafe(col.get(), j) << " ";
        Rcpp::Rcout << std::endl;
    }
}

//  ctx_wrap_t* ctxwrap_p = new ContextWrapper(ctxptr);
//  Rcpp::XPtr<ctx_wrap_t> ctx_wrap_xptr = make_xptr<ctx_wrap_t>(ctxwrap_p, false);

// [[Rcpp::export]]
void createSchemaFromArrow(const std::string& uri,
                           naxpSchema nasp, naxpArray nadimap, naxpSchema nadimsp,
                           bool sparse, std::string datatype,
                           Rcpp::List pclst,
                           Rcpp::Nullable<Rcpp::XPtr<ctx_wrap_t>> ctxptr = R_NilValue) {
    //struct ArrowArray* ap = (struct ArrowArray*) R_ExternalPtrAddr(naap);
    //struct ArrowSchema* sp = (struct ArrowSchema*) R_ExternalPtrAddr(nasp);
    //
    // or
    // auto ap = nanoarrow_array_from_xptr(naap);
    // auto sp = nanoarrow_schema_from_xptr(nasp);
    //
    // or:
    //nanoarrow::UniqueArray ap{nanoarrow_array_from_xptr(naap)};
    nanoarrow::UniqueSchema sp{nanoarrow_schema_from_xptr(nasp)};
    //_show_content(ap, sp);
    nanoarrow::UniqueArray apdim{nanoarrow_array_from_xptr(nadimap)};
    nanoarrow::UniqueSchema spdim{nanoarrow_schema_from_xptr(nadimsp)};
    //_show_content(apdim, spdim);

    auto schema = std::make_unique<ArrowSchema>();
    sp.move(schema.get());

    auto dimsch = std::make_unique<ArrowSchema>();
    spdim.move(dimsch.get());
    auto dimarr = std::make_unique<ArrowArray>();
    apdim.move(dimarr.get());

    tdbs::PlatformConfig pltcfg;
    pltcfg.dataframe_dim_zstd_level       = Rcpp::as<int>(pclst["dataframe_dim_zstd_level"]);
    pltcfg.sparse_nd_array_dim_zstd_level = Rcpp::as<int>(pclst["sparse_nd_array_dim_zstd_level"]);
    pltcfg.dense_nd_array_dim_zstd_level  = Rcpp::as<int>(pclst["dense_nd_array_dim_zstd_level"]);
    pltcfg.write_X_chunked                = Rcpp::as<bool>(pclst["write_X_chunked"]);
    pltcfg.goal_chunk_nnz                 = Rcpp::as<double>(pclst["goal_chunk_nnz"]);
    pltcfg.capacity                       = Rcpp::as<double>(pclst["capacity"]);
    pltcfg.offsets_filters                = Rcpp::as<std::string>(pclst["offsets_filters"]);
    pltcfg.validity_filters               = Rcpp::as<std::string>(pclst["validity_filters"]);
    pltcfg.allows_duplicates              = Rcpp::as<bool>(pclst["allows_duplicates"]);
    pltcfg.cell_order                     = Rcpp::as<std::string>(pclst["cell_order"]);
    pltcfg.tile_order                     = Rcpp::as<std::string>(pclst["tile_order"]);
    pltcfg.attrs                          = Rcpp::as<std::string>(pclst["attrs"]);
    pltcfg.dims                           = Rcpp::as<std::string>(pclst["dims"]);

    std::shared_ptr<tiledb::Context> ctxsp;
    if (ctxptr.isNotNull()) {   				// if optional context wrapper pointer was present
        Rcpp::XPtr<ctx_wrap_t> ctxxp(ctxptr); 	// Rcpp::Nullable<> needs instantiation
        ctxsp = (*ctxxp.get()).ctxptr; 			// access shared pointer to contexct from struct
    } else {
        //std::map<std::string, std::string> pltfrmcfgmap = pltcfg.to_list();
        tiledb::Config cfg; //pltcfg); 			// create plain config
        										// create shared pointer to context given config
        ctxsp = std::make_shared<tiledb::Context>(cfg);
        //ctx_wrap_t* ctxwrap_p = new ContextWrapper(ctxsp); 	// create wrapper struct
        //ctxptr = make_xptr<ctx_wrap_t>(ctxwrap_p, false);   // and create and assign extptr
    }
    // create the ArraySchema
    auto as = tdbs::ArrowAdapter::tiledb_schema_from_arrow_schema(ctxsp, std::move(schema),
                                                                  std::pair(std::move(dimarr),
                                                                            std::move(dimsch)),
                                                                  datatype, sparse,
                                                                  pltcfg);

    // We can inspect the (TileDB) ArraySchema via a simple helper:  as.dump();
    tiledb::VFS vfs(*ctxsp);
    //spdl::warn("[createSchemaFromArrow] About to create {}, dir exists ? {}", uri, vfs.is_dir(uri));

    if (vfs.is_dir(uri)) {
        Rcpp::stop(tfm::format("Error: Array '%s' already exists", uri));
    } else {
        // Create the schema at the given URI
        tiledb::Array::create(uri, as);
    }
}


// [[Rcpp::export]]
void writeArrayFromArrow(const std::string& uri, naxpArray naap, naxpSchema nasp,
                         const std::string arraytype = "",
                         Rcpp::Nullable<Rcpp::CharacterVector> config = R_NilValue) {

    //struct ArrowArray* ap = (struct ArrowArray*) R_ExternalPtrAddr(naap);
    //struct ArrowSchema* sp = (struct ArrowSchema*) R_ExternalPtrAddr(nasp);
    //
    // or
    // auto ap = nanoarrow_array_from_xptr(naap);
    // auto sp = nanoarrow_schema_from_xptr(nasp);
    //
    // or:
    nanoarrow::UniqueArray ap{nanoarrow_array_from_xptr(naap)};
    nanoarrow::UniqueSchema sp{nanoarrow_schema_from_xptr(nasp)};

    // now move nanoarrow unique arrays (created from objects handed from R) into
    // proper unique pointers to arrow schema and array
    auto schema = std::make_unique<ArrowSchema>();
    sp.move(schema.get());
    auto array = std::make_unique<ArrowArray>();
    ap.move(array.get());

    // if we hae a coonfig, use it
    std::shared_ptr<tdbs::SOMAContext> somactx;
    if (config.isNotNull()) {
        std::map<std::string, std::string> smap;
        auto config_vec = config.as();
        auto config_names = Rcpp::as<Rcpp::CharacterVector>(config_vec.names());
        for (auto &name : config_names) {
            std::string param = Rcpp::as<std::string>(name);
            std::string value = Rcpp::as<std::string>(config_vec[param]);
            smap[param] = value;
        }
        somactx = std::make_shared<tdbs::SOMAContext>(smap);
    } else {
        somactx = std::make_shared<tdbs::SOMAContext>();
    }

    std::shared_ptr<tdbs::SOMAArray> arrup;
    if (arraytype == "SOMADataFrame") {
        arrup = tdbs::SOMADataFrame::open(OpenMode::write, uri, somactx);
    } else if (arraytype == "SOMADenseNDArray") {
        arrup = tdbs::SOMADenseNDArray::open(OpenMode::write, uri, somactx,
                                             "unnamed", {}, "auto", ResultOrder::colmajor);
    }

    arrup.get()->set_array_data(std::move(schema), std::move(array));
    arrup.get()->write();
    arrup.get()->close();

}
