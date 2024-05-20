#include <Rcpp/Lighter>                         // for R interface to C++
#include <nanoarrow/r.h>                        // for C/C++ interface to Arrow (via header export from the R package)
#include <nanoarrow/nanoarrow.hpp>              // for C/C++ interface to Arrow (vendored)
#include <RcppInt64>                            // for fromInteger64

#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"         // local declarations
#include "xptr-utils.h"         // xptr taggging utilities

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
void createSchemaFromArrow(const std::string& uri,
                           SEXP nasp,
                           SEXP nadimap, SEXP nadimsp,
                           Rcpp::List pclst) {
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

    auto ctx = tiledb::Context();
    auto vfs = tiledb::VFS(ctx);

    auto ctxsp = std::make_shared<tiledb::Context>(ctx);
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

    // create the ArraySchema
    auto as = tdbs::ArrowAdapter::tiledb_schema_from_arrow_schema(ctxsp, std::move(schema),
                                                                  std::pair(std::move(dimarr),
                                                                            std::move(dimsch)),
                                                                  "SOMADataFrame",
                                                                  true, // sparse
                                                                  pltcfg);
    tiledb::Array::create(uri, as);
}

// [[Rcpp::export]]
void writeArrayFromArrow(const std::string& uri, SEXP naap, SEXP nasp) {
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

    auto somactx = std::make_shared<tdbs::SOMAContext>();
    auto arrup = tdbs::SOMADataFrame::open(OpenMode::write, uri, somactx);

    auto schema = std::make_unique<ArrowSchema>();
    sp.move(schema.get());
    auto array = std::make_unique<ArrowArray>();
    ap.move(array.get());

    arrup.get()->set_array_data(std::move(schema), std::move(array));
    arrup.get()->write();
    arrup.get()->close();
}
