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

// [[Rcpp::export]]
Rcpp::XPtr<somactx_wrap_t> createSOMAContext(Rcpp::Nullable<Rcpp::CharacterVector> config = R_NilValue) {

    // if we hae a config, use it
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

    auto ptr = new somactx_wrap_t(somactx);
    auto xp = make_xptr<somactx_wrap_t>(ptr);
    return xp;
}


//  ctx_wrap_t* ctxwrap_p = new ContextWrapper(ctxptr);
//  Rcpp::XPtr<ctx_wrap_t> ctx_wrap_xptr = make_xptr<ctx_wrap_t>(ctxwrap_p, false);

// [[Rcpp::export]]
void createSchemaFromArrow(const std::string& uri, naxpSchema nasp, naxpArray nadimap, naxpSchema nadimsp,
                           bool sparse, std::string datatype, Rcpp::List pclst, Rcpp::XPtr<somactx_wrap_t> ctxxp) {


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

    // shared pointer to SOMAContext from external pointer wrapper
    std::shared_ptr<tdbs::SOMAContext> sctx = ctxxp->ctxptr;
    // shared pointer to TileDB Context from SOMAContext
    std::shared_ptr<tiledb::Context> ctx = sctx->tiledb_ctx();

    bool exists = false;
    if (datatype == "SOMADataFrame") {
        exists = tdbs::SOMADataFrame::exists(uri, sctx);
    } else if (datatype == "SOMASparseNDArray") {
        exists = tdbs::SOMASparseNDArray::exists(uri, sctx);
    } else if (datatype == "SOMADenseNDArray") {
        exists = tdbs::SOMADenseNDArray::exists(uri, sctx);
    } else {
        Rcpp::stop(tfm::format("Error: Invalid SOMA type_argument '%s'", datatype));
    }

    if (exists) {
        Rcpp::stop(tfm::format("Error: Array '%s' already exists", uri));
    }

    if (datatype == "SOMADataFrame") {
        tdbs::SOMADataFrame::create(uri, std::move(schema),
                                    std::pair(std::move(dimarr), std::move(dimsch)),
                                    sctx, pltcfg);
    } else if (datatype == "SOMASparseNDArray") {
        // for arrays n_children will be three as we have two dims and a data col
        std::string datacoltype = sp->children[sp->n_children-1]->format;
        tdbs::SOMASparseNDArray::create(uri, datacoltype,
                                        std::pair(std::move(dimarr), std::move(dimsch)),
                                        sctx, pltcfg);
    } else if (datatype == "SOMADenseNDArray") {
        // for arrays n_children will be three as we have two dims and a data col
        std::string datacoltype = sp->children[sp->n_children-1]->format;
        tdbs::SOMADenseNDArray::create(uri, datacoltype,
                                       std::pair(std::move(dimarr), std::move(dimsch)),
                                       sctx, pltcfg);
    } else {
        Rcpp::stop(tfm::format("Error: Invalid SOMA type_argument '%s'", datatype));
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
    } else if (arraytype == "SOMASparseNDArray") {
        arrup = tdbs::SOMASparseNDArray::open(OpenMode::write, uri, somactx);
    } else {  // not reached
        Rcpp::stop(tfm::format("Unexpected array type '%s'", arraytype));
    }

    arrup.get()->set_array_data(std::move(schema), std::move(array));
    arrup.get()->write();
    arrup.get()->close();
}
