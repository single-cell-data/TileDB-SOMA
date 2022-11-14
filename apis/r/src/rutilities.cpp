
#include <Rcpp.h>
#include <tiledbsoma/tiledbsoma>
#include <archAPI.h>
#include "rutilities.h"

namespace tdbs = tiledbsoma;

void apply_dim_points(tdbs::SOMAReader *sr,
                      std::unordered_map<std::string, tiledb_datatype_t>& name2type,
                      Rcpp::List lst) {
    std::vector<std::string> colnames = lst.attr("names");
    for (auto& nm: colnames) {
        auto tp = name2type[nm];
        if (tp == TILEDB_UINT64) {
            Rcpp::NumericVector payload = lst[nm];
            std::vector<int64_t> iv = getInt64Vector(payload);
            std::vector<uint64_t> uv(iv.size());
            for (size_t i=0; i<iv.size(); i++) {
                uv[i] = static_cast<uint64_t>(iv[i]);
                sr->set_dim_point<uint64_t>(nm, uv[i]);  // bonked when use with vector
                spdl::info(fmt::format("[export_arrow_array] Applying dim point {} on {}", uv[i], nm));
            }
        } else if (tp == TILEDB_INT64) {
            Rcpp::NumericVector payload = lst[nm];
            std::vector<int64_t> iv = getInt64Vector(payload);
            for (size_t i=0; i<iv.size(); i++) {
                sr->set_dim_point<int64_t>(nm, iv[i]);
                spdl::info(fmt::format("[export_arrow_array] Applying dim point {} on {}", iv[i], nm));
            }
        } else if (tp == TILEDB_FLOAT32) {
            Rcpp::NumericVector payload = lst[nm];
            for (R_xlen_t i=0; i<payload.size(); i++) {
                float v = static_cast<uint64_t>(payload[i]);
                sr->set_dim_point<float>(nm, v);
                spdl::info(fmt::format("[export_arrow_array] Applying dim point {} on {}", v, nm));
            }
        } else if (tp == TILEDB_FLOAT64) {
            Rcpp::NumericVector payload = lst[nm];
            for (R_xlen_t i=0; i<payload.size(); i++) {
                sr->set_dim_point<double>(nm,payload[i]);
                spdl::info(fmt::format("[export_arrow_array] Applying dim point {} on {}", payload[i], nm));
            }
        } else if (tp == TILEDB_INT32) {
            Rcpp::IntegerVector payload = lst[nm];
            for (R_xlen_t i=0; i<payload.size(); i++) {
                sr->set_dim_point<int32_t>(nm,payload[i]);
                spdl::info(fmt::format("[export_arrow_array] Applying dim point {} on {}", payload[i], nm));
            }
        } else {
            Rcpp::stop("Currently unsupported type: ", tiledb::impl::to_str(tp));
        }
    }
}

void apply_dim_ranges(tdbs::SOMAReader* sr,
                      std::unordered_map<std::string, tiledb_datatype_t>& name2type,
                      Rcpp::List lst) {
    std::vector<std::string> colnames = lst.attr("names");
    for (auto& nm: colnames) {
        auto tp = name2type[nm];
        if (tp == TILEDB_UINT64) {
            Rcpp::NumericMatrix mm = lst[nm];
            Rcpp::NumericMatrix::Column lo = mm.column(0); // works as proxy for int and float types
            Rcpp::NumericMatrix::Column hi = mm.column(1); // works as proxy for int and float types
            for (int i=0; i<mm.nrow(); i++) {
                uint64_t l = static_cast<uint64_t>(makeScalarInteger64(lo[i]));
                uint64_t h = static_cast<uint64_t>(makeScalarInteger64(hi[i]));
                std::vector<std::pair<uint64_t, uint64_t>> vp{std::make_pair(l,h)};
                sr->set_dim_ranges<uint64_t>(nm, vp);
                spdl::info(fmt::format("[export_arrow_array] Applying dim point {} on {} with {} - {}", i, nm, l, h));
            }
        } else if (tp == TILEDB_INT64) {
            Rcpp::NumericMatrix mm = lst[nm];
            std::vector<int64_t> lo = getInt64Vector(mm.column(0));
            std::vector<int64_t> hi = getInt64Vector(mm.column(1));
            for (int i=0; i<mm.nrow(); i++) {
                std::vector<std::pair<int64_t, int64_t>> vp{std::make_pair(lo[i], hi[i])};
                sr->set_dim_ranges<int64_t>(nm, vp);
                spdl::info(fmt::format("[export_arrow_array] Applying dim point {} on {} with {} - {}", i, nm, lo[i], hi[i]));
            }
        } else if (tp == TILEDB_FLOAT32) {
            Rcpp::NumericMatrix mm = lst[nm];
            Rcpp::NumericMatrix::Column lo = mm.column(0); // works as proxy for int and float types
            Rcpp::NumericMatrix::Column hi = mm.column(1); // works as proxy for int and float types
            for (int i=0; i<mm.nrow(); i++) {
                float l = static_cast<float_t>(lo[i]);
                float h = static_cast<float_t>(hi[i]);
                std::vector<std::pair<float, float>> vp{std::make_pair(l,h)};
                sr->set_dim_ranges<float>(nm, vp);
                spdl::info(fmt::format("[export_arrow_array] Applying dim point {} on {} with {} - {}", i, nm, l, h));
            }
        } else if (tp == TILEDB_FLOAT64) {
            Rcpp::NumericMatrix mm = lst[nm];
            Rcpp::NumericMatrix::Column lo = mm.column(0); // works as proxy for int and float types
            Rcpp::NumericMatrix::Column hi = mm.column(1); // works as proxy for int and float types
            for (int i=0; i<mm.nrow(); i++) {
                std::vector<std::pair<double, double>> vp{std::make_pair(lo[i],hi[i])};
                sr->set_dim_ranges<double>(nm, vp);
                spdl::info(fmt::format("[export_arrow_array] Applying dim point {} on {} with {} - {}", i, nm, lo[i], hi[i]));
            }
        } else if (tp == TILEDB_INT32) {
            Rcpp::IntegerMatrix mm = lst[nm];
            Rcpp::IntegerMatrix::Column lo = mm.column(0); // works as proxy for int and float types
            Rcpp::IntegerMatrix::Column hi = mm.column(1); // works as proxy for int and float types
            for (int i=0; i<mm.nrow(); i++) {
                int32_t l = static_cast<int32_t>(lo[i]);
                int32_t h = static_cast<int32_t>(hi[i]);
                std::vector<std::pair<int32_t, int32_t>> vp{std::make_pair(l,h)};
                sr->set_dim_ranges<int32_t>(nm, vp);
                spdl::info(fmt::format("[export_arrow_array] Applying dim point {} on {} with {} - {}", i, nm[i], l, h));
            }
        } else {
            Rcpp::stop("Currently unsupported type: ", tiledb::impl::to_str(tp));
        }
    }
}
