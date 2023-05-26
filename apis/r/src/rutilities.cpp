
// we currently get deprecation warnings by default which are noisy
#ifndef TILEDB_NO_API_DEPRECATION_WARNINGS
#define TILEDB_NO_API_DEPRECATION_WARNINGS
#endif

#include <Rcpp.h>               // for R interface to C++
#include <nanoarrow.h>          // for C interface to Arrow

// We get these via nanoarrow and must cannot include carrow.h again
#define ARROW_SCHEMA_AND_ARRAY_DEFINED 1
#include <tiledbsoma/tiledbsoma>

#include "rutilities.h"         // local declarations
#include "xptr-utils.h"         // xptr taggging utilitie

namespace tdbs = tiledbsoma;

void apply_dim_points(tdbs::SOMAArray *sr,
                      std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>>& name2dim,
                      Rcpp::List lst) {
    std::vector<std::string> colnames = lst.attr("names");
    for (auto& nm: colnames) {
        auto dm = name2dim[nm];
        auto tp = dm->type();
        if (tp == TILEDB_UINT64) {
            Rcpp::NumericVector payload = lst[nm];
            std::vector<int64_t> iv = getInt64Vector(payload);
            std::vector<uint64_t> uv(iv.size());
            const std::pair<uint64_t,uint64_t> pr = dm->domain<uint64_t>();
            for (size_t i=0; i<iv.size(); i++) {
                uv[i] = static_cast<uint64_t>(iv[i]);
                if (uv[i] >= pr.first && uv[i] <= pr.second) {
                    sr->set_dim_point<uint64_t>(nm, uv[i]);  // bonked when use with vector
                    spdl::info("[apply_dim_points] Applying dim point {} on {}", uv[i], nm);
                }
            }
        } else if (tp == TILEDB_INT64) {
            Rcpp::NumericVector payload = lst[nm];
            std::vector<int64_t> iv = getInt64Vector(payload);
            const std::pair<int64_t,int64_t> pr = dm->domain<int64_t>();
            for (size_t i=0; i<iv.size(); i++) {
                if (iv[i] >= pr.first && iv[i] <= pr.second) {
                    sr->set_dim_point<int64_t>(nm, iv[i]);
                    spdl::info("[apply_dim_points] Applying dim point {} on {}", iv[i], nm);
                }
            }
        } else if (tp == TILEDB_FLOAT32) {
            Rcpp::NumericVector payload = lst[nm];
            const std::pair<float,float> pr = dm->domain<float>();
            for (R_xlen_t i=0; i<payload.size(); i++) {
                float v = static_cast<float>(payload[i]);
                if (v >= pr.first && v <= pr.second) {
                    sr->set_dim_point<float>(nm, v);
                    spdl::info("[apply_dim_points] Applying dim point {} on {}", v, nm);
                }
            }
        } else if (tp == TILEDB_FLOAT64) {
            Rcpp::NumericVector payload = lst[nm];
            const std::pair<double,double> pr = dm->domain<double>();
            for (R_xlen_t i=0; i<payload.size(); i++) {
                if (payload[i] >= pr.first && payload[i] <= pr.second) {
                    sr->set_dim_point<double>(nm,payload[i]);
                    spdl::info("[apply_dim_points] Applying dim point {} on {}", payload[i], nm);
                }
            }
        } else if (tp == TILEDB_INT32) {
            Rcpp::IntegerVector payload = lst[nm];
            const std::pair<int32_t,int32_t> pr = dm->domain<int32_t>();
            for (R_xlen_t i=0; i<payload.size(); i++) {
                if (payload[i] >= pr.first && payload[i] <= pr.second) {
                    sr->set_dim_point<int32_t>(nm,payload[i]);
                    spdl::info("[apply_dim_points] Applying dim point {} on {}", payload[i], nm);
                }
            }
        } else {
            Rcpp::stop("Currently unsupported type: ", tiledb::impl::to_str(tp));
        }
    }
}

void apply_dim_ranges(tdbs::SOMAArray* sr,
                      std::unordered_map<std::string, std::shared_ptr<tiledb::Dimension>>& name2dim,
                      Rcpp::List lst) {
    std::vector<std::string> colnames = lst.attr("names");
    for (auto& nm: colnames) {
        auto dm = name2dim[nm];
        auto tp = dm->type();
        if (tp == TILEDB_UINT64) {
            Rcpp::NumericMatrix mm = lst[nm];
            Rcpp::NumericMatrix::Column lo = mm.column(0); // works as proxy for int and float types
            Rcpp::NumericMatrix::Column hi = mm.column(1); // works as proxy for int and float types
            std::vector<std::pair<uint64_t, uint64_t>> vp(mm.nrow());
            const std::pair<uint64_t,uint64_t> pr = dm->domain<uint64_t>();
            for (int i=0; i<mm.nrow(); i++) {
                uint64_t l = static_cast<uint64_t>(makeScalarInteger64(lo[i]));
                uint64_t h = static_cast<uint64_t>(makeScalarInteger64(hi[i]));
                vp[i] = std::make_pair(std::max(l,pr.first), std::min(h, pr.second));
                spdl::info("[apply_dim_ranges] Applying dim point {} on {} with {} - {}", i, nm, l, h) ;
            }
            sr->set_dim_ranges<uint64_t>(nm, vp);
        } else if (tp == TILEDB_INT64) {
            Rcpp::NumericMatrix mm = lst[nm];
            std::vector<int64_t> lo = getInt64Vector(mm.column(0));
            std::vector<int64_t> hi = getInt64Vector(mm.column(1));
            std::vector<std::pair<int64_t, int64_t>> vp(mm.nrow());
            const std::pair<int64_t,int64_t> pr = dm->domain<int64_t>();
            for (int i=0; i<mm.nrow(); i++) {
                vp[i] = std::make_pair(std::max(lo[i],pr.first), std::min(hi[i], pr.second));
                spdl::info("[apply_dim_ranges] Applying dim point {} on {} with {} - {}", i, nm, lo[i], hi[i]) ;
            }
            sr->set_dim_ranges<int64_t>(nm, vp);
        } else if (tp == TILEDB_FLOAT32) {
            Rcpp::NumericMatrix mm = lst[nm];
            Rcpp::NumericMatrix::Column lo = mm.column(0); // works as proxy for int and float types
            Rcpp::NumericMatrix::Column hi = mm.column(1); // works as proxy for int and float types
            std::vector<std::pair<float, float>> vp(mm.nrow());
            const std::pair<float,float> pr = dm->domain<float>();
            for (int i=0; i<mm.nrow(); i++) {
                float l = static_cast<float>(lo[i]);
                float h = static_cast<float>(hi[i]);
                vp[i] = std::make_pair(std::max(l,pr.first), std::min(h, pr.second));
                spdl::info("[apply_dim_ranges] Applying dim point {} on {} with {} - {}", i, nm, l, h) ;
            }
            sr->set_dim_ranges<float>(nm, vp);
        } else if (tp == TILEDB_FLOAT64) {
            Rcpp::NumericMatrix mm = lst[nm];
            Rcpp::NumericMatrix::Column lo = mm.column(0); // works as proxy for int and float types
            Rcpp::NumericMatrix::Column hi = mm.column(1); // works as proxy for int and float types
            std::vector<std::pair<double, double>> vp(mm.nrow());
            const std::pair<double,double> pr = dm->domain<double>();
            for (int i=0; i<mm.nrow(); i++) {
                vp[i] = std::make_pair(std::max(lo[i],pr.first), std::min(hi[i], pr.second));
                spdl::info("[apply_dim_ranges] Applying dim point {} on {} with {} - {}", i, nm, lo[i], hi[i]) ;
            }
            sr->set_dim_ranges<double>(nm, vp);
        } else if (tp == TILEDB_INT32) {
            Rcpp::IntegerMatrix mm = lst[nm];
            Rcpp::IntegerMatrix::Column lo = mm.column(0); // works as proxy for int and float types
            Rcpp::IntegerMatrix::Column hi = mm.column(1); // works as proxy for int and float types
            std::vector<std::pair<int32_t, int32_t>> vp(mm.nrow());
            const std::pair<int32_t,int32_t> pr = dm->domain<int32_t>();
            for (int i=0; i<mm.nrow(); i++) {
                vp[i] = std::make_pair(std::max(lo[i],pr.first), std::min(hi[i], pr.second));
                spdl::info("[apply_dim_ranges] Applying dim point {} on {} with {} - {}", i, nm[i], lo[i], hi[i]) ;
            }
            sr->set_dim_ranges<int32_t>(nm, vp);
        } else {
            Rcpp::stop("Currently unsupported type: ", tiledb::impl::to_str(tp));
        }
    }
}


// initialize arrow schema and array, respectively
Rcpp::XPtr<ArrowSchema> schema_setup_struct(Rcpp::XPtr<ArrowSchema> schxp, int64_t n_children) {
    ArrowSchema* schema = schxp.get();
    auto type = NANOARROW_TYPE_STRUCT;

    ArrowSchemaInit(schema);    					// modified from ArrowSchemaInitFromType()
    int result = ArrowSchemaSetType(schema, type);
    if (result != NANOARROW_OK) {
        schema->release(schema);
        Rcpp::stop("Error setting struct schema");
    }

    // now adapted from ArrowSchemaAllocateChildren
    if (schema->children != NULL) Rcpp::stop("Error allocation as children not null");

    if (n_children > 0) {
        auto ptr = (struct ArrowSchema**) ArrowMalloc(n_children * sizeof(struct ArrowSchema*));
        Rcpp::XPtr<ArrowSchema*> schema_ptrxp = make_xptr(ptr, false);
        schema->children = schema_ptrxp.get();
        if (schema->children == NULL) Rcpp::stop("Failed to allocate ArrowSchema*");

        schema->n_children = n_children;
        memset(schema->children, 0, n_children * sizeof(struct ArrowSchema*));

        for (int64_t i = 0; i < n_children; i++) {
            schema->children[i] = schema_owning_xptr();
            if (schema->children[i] == NULL) Rcpp::stop("Error allocation schema child %ld", i);
            schema->children[i]->release = NULL;
        }
    }
    return schxp;
}

extern "C" {
    void ArrowArrayRelease(struct ArrowArray *array); 		// made non-static in nanoarrow.c
    ArrowErrorCode ArrowArraySetStorageType(struct ArrowArray* array,	// ditto
                                            enum ArrowType storage_type);
}

Rcpp::XPtr<ArrowArray> array_setup_struct(Rcpp::XPtr<ArrowArray> arrxp, int64_t n_children) {
    ArrowArray* array = arrxp.get();
    auto storage_type = NANOARROW_TYPE_STRUCT;

    array->length = 0;
    array->null_count = 0;
    array->offset = 0;
    array->n_buffers = 0;
    array->n_children = 0;
    array->buffers = NULL;
    array->children = NULL;
    array->dictionary = NULL;
    array->release = &ArrowArrayRelease;
    array->private_data = NULL;

    auto private_data = (struct ArrowArrayPrivateData*) ArrowMalloc(sizeof(struct ArrowArrayPrivateData));
    if (private_data == NULL) {
        array->release = NULL;
        Rcpp::stop("Error allocating array private data");
    }
    ArrowBitmapInit(&private_data->bitmap);
    ArrowBufferInit(&private_data->buffers[0]);
    ArrowBufferInit(&private_data->buffers[1]);
    private_data->buffer_data[0] = NULL;
    private_data->buffer_data[1] = NULL;
    private_data->buffer_data[2] = NULL;
    array->private_data = private_data;
    array->buffers = (const void**)(&private_data->buffer_data);
    int result = ArrowArraySetStorageType(array, storage_type);
    if (result != NANOARROW_OK) {
        array->release(array);
        Rcpp::stop("Error setting array storage type");
    }

    ArrowLayoutInit(&private_data->layout, storage_type);
    // We can only know this not to be true when initializing based on a schema so assume this to be true.
    private_data->union_type_id_is_child_index = 1;


    // remainder from ArrowArrayAllocateChildren()
    if (array->children != NULL) Rcpp::stop("Error allocating array children as pointer not null");

    if (n_children == 0) {
        return arrxp;
    }

    auto ptr = (struct ArrowArray**) ArrowMalloc(n_children * sizeof(struct ArrowArray*));
    Rcpp::XPtr<ArrowArray*> array_ptrxp = make_xptr(ptr, false);
    array->children = array_ptrxp.get();
    if (array->children == NULL) Rcpp::stop("Failed to allocated ArrayArray*");

    memset(array->children, 0, n_children * sizeof(struct ArrowArray*));

    for (int64_t i = 0; i < n_children; i++) {
        array->children[i] = array_owning_xptr();
        if (array->children[i] == NULL) Rcpp::stop("Error allocation array child %ld", i);
        array->children[i]->release = NULL;
    }
    array->n_children = n_children;
    return arrxp;
}
