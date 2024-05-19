#include <Rcpp/Lighter>                         // for R interface to C++
#include <nanoarrow/r.h>                        // for C/C++ interface to Arrow (via header export from the R package)
#include <nanoarrow/nanoarrow.hpp>              // for C/C++ interface to Arrow (vendored)
#include <RcppInt64>                            // for fromInteger64

#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"         // local declarations
#include "xptr-utils.h"         // xptr taggging utilities

namespace tdbs = tiledbsoma;

// we should not have to add this ...
inline ArrowType format2type(std::string fmt) {
    if      (fmt == "n") return NANOARROW_TYPE_NA;
    else if (fmt == "b") return NANOARROW_TYPE_BOOL;
    else if (fmt == "C") return NANOARROW_TYPE_UINT8;
    else if (fmt == "c") return NANOARROW_TYPE_INT8;
    else if (fmt == "S") return NANOARROW_TYPE_UINT16;
    else if (fmt == "s") return NANOARROW_TYPE_INT16;
    else if (fmt == "I") return NANOARROW_TYPE_UINT32;
    else if (fmt == "i") return NANOARROW_TYPE_INT32;
    else if (fmt == "L") return NANOARROW_TYPE_UINT64;
    else if (fmt == "l") return NANOARROW_TYPE_INT64;
    else if (fmt == "f") return NANOARROW_TYPE_FLOAT;
    else if (fmt == "g") return NANOARROW_TYPE_DOUBLE;
    else if (fmt == "u") return NANOARROW_TYPE_STRING;
    else if (fmt == "U") return NANOARROW_TYPE_LARGE_STRING;
    else if (fmt == "z") return NANOARROW_TYPE_BINARY;
    else if (fmt == "Z") return NANOARROW_TYPE_LARGE_BINARY;
    else Rcpp::stop("Unsupported format string: ", fmt);
}

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

// void _fill_in_helper(ArrowArray& dom, ArrowArray&ext, int pos, ArrowType type, int domlo, int domhi, int extend) {
//     ArrowArrayInitFromType(dom.children[pos], type);
//     ArrowArrayStartAppending(dom.children[pos]);
//     ArrowArrayAppendInt(dom.children[pos],domlo);
//     ArrowArrayAppendInt(dom.children[pos],domhi);
//     ArrowArrayFinishBuildingDefault(dom.children[pos], nullptr);

//     ArrowArrayInitFromType(ext.children[pos], type);
//     ArrowArrayStartAppending(ext.children[pos]);
//     ArrowArrayAppendInt(ext.children[pos],extend);
//     ArrowArrayFinishBuildingDefault(ext.children[pos], nullptr);
// }


const char* _tiledb_filter_to_string(tiledb_filter_type_t filter) {
  switch(filter) {
    case TILEDB_FILTER_NONE:
     return "NONE";
    case TILEDB_FILTER_GZIP:
      return "GZIP";
    case TILEDB_FILTER_ZSTD:
      return "ZSTD";
    case TILEDB_FILTER_LZ4:
      return "LZ4";
    case TILEDB_FILTER_RLE:
      return "RLE";
    case TILEDB_FILTER_BZIP2:
      return "BZIP2";
    case TILEDB_FILTER_DOUBLE_DELTA:
      return "DOUBLE_DELTA";
    case TILEDB_FILTER_BIT_WIDTH_REDUCTION:
      return "BIT_WIDTH_REDUCTION";
    case TILEDB_FILTER_BITSHUFFLE:
      return "BITSHUFFLE";
    case TILEDB_FILTER_BYTESHUFFLE:
      return "BYTESHUFFLE";
    case TILEDB_FILTER_POSITIVE_DELTA:
      return "POSITIVE_DELTA";
    case TILEDB_FILTER_CHECKSUM_MD5:
      return "CHECKSUM_MD5";
    case TILEDB_FILTER_CHECKSUM_SHA256:
      return "CHECKSUM_SHA256";
#if TILEDB_VERSION >= TileDB_Version(2,9,0)
    case TILEDB_FILTER_DICTIONARY:
      return "DICTIONARY_ENCODING";
#endif
#if TILEDB_VERSION >= TileDB_Version(2,11,0)
    case TILEDB_FILTER_SCALE_FLOAT:
      return "SCALE_FLOAT";
#endif
#if TILEDB_VERSION >= TileDB_Version(2,12,0)
  case TILEDB_FILTER_XOR:
    return "FILTER_XOR";
#endif
  default: {
      Rcpp::stop("unknown tiledb_filter_t (%d)", filter);
    }
  }
}

std::vector<std::string> _list2vector(Rcpp::List lst) {
    std::vector<std::string> v;
    for (auto i=0; i<lst.size(); i++) {
        Rcpp::S4 s{lst[i]};
        Rcpp::XPtr<tiledb::Filter> fxp = s.slot("ptr");
        std::string str{ _tiledb_filter_to_string( fxp->filter_type() ) };
        v.push_back(str);
    }
    return v;
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

    // compAArrowError ec;
    // std::vector<std::string> dimnam = {"d1", "d2"};
    // ArrowArray dom, ext;
    // int nc = 2;
    // ArrowArrayInitFromType(&dom, NANOARROW_TYPE_STRUCT);
    // ArrowArrayAllocateChildren(&dom, nc);
    // ArrowArrayInitFromType(&ext, NANOARROW_TYPE_STRUCT);
    // ArrowArrayAllocateChildren(&ext, nc);
    // _fill_in_helper(dom, ext, 0, NANOARROW_TYPE_INT16, 0, 100, 100); // could get types and values from payload
    // _fill_in_helper(dom, ext, 1, NANOARROW_TYPE_INT16, 0, 50, 50);
    // ArrowArrayFinishBuildingDefault(&dom, &ec);
    // ArrowArrayFinishBuildingDefault(&ext, &ec);

    //auto domsp = std::make_shared<ArrowArray>(dom);
    //auto extsp = std::make_shared<ArrowArray>(ext);
    //tiledbsoma::ColumnIndexInfo ici = { dimnam, domsp, extsp };

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
    if (vfs.is_dir(uri)) vfs.remove_dir(uri);
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

// [[Rcpp::export]]
SEXP getArrowSchema(const std::string& uri) {
    auto somactx = std::make_shared<tdbs::SOMAContext>();

    auto arrup = tdbs::SOMAArray::open(OpenMode::read, uri, somactx);
    auto schup = arrup.get()->arrow_schema();

    auto schxp = nanoarrow_schema_owning_xptr();
    ArrowSchemaMove(schup.get(), (ArrowSchema*) xptr_addr(schxp));
    return schxp;
}
