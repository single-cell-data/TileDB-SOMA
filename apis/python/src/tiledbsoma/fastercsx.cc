/**
 * @file   fastercsx.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2024 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 * Python bindings for CSX conversion primitives.
 */

#include <variant>

#include <tiledbsoma/utils/fastercsx.h>
#include "common.h"

namespace libtiledbsomacpp {

namespace py = pybind11;
using npy_api = pybind11::detail::npy_api;

using namespace py::literals;
using namespace tiledbsoma;

/**
 * @brief Convert py::tuple<py::array> to std::vector<py::array>.
 *
 * Will throw if the tuple does not contain a Numpy ndarray of the correct type.
 *
 * @tparam T
 * @param tup
 * @return std::vector<T>
 */
template <typename T>
std::vector<T> to_vector(const py::tuple& tup) {
    auto vec = std::vector<T>();
    vec.reserve(tup.size());
    for (py::size_t i = 0; i < tup.size(); ++i) {
        auto v = tup[i];
        if (!py::isinstance<T>(v))
            throw py::type_error(
                std::string("Unable to cast tuple element to '") +
                py::type_id<T>() + "'");
        vec.push_back(v.cast<T>());
    }
    return vec;
}

template <typename T>
tcb::span<T const> make_span(py::array arr) {
    assert(py::isinstance<py::array_t<T>>(arr));
    auto arr_typed = py::cast<py::array_t<T>>(arr);
    return tcb::span<const T>(arr.unchecked<T, 1>().data(0), arr_typed.size());
}

template <typename T>
tcb::span<T> make_mutable_span(py::array arr) {
    assert(py::isinstance<py::array_t<T>>(arr));
    auto arr_typed = py::cast<py::array_t<T>>(arr);
    return tcb::span<T>(
        arr.mutable_unchecked<T, 1>().mutable_data(0), arr_typed.size());
}

/**
 * @brief C++ std::type_identity, which is not available until C++20
 *
 * @tparam T
 */
template <class T>
struct type_identity {
    using type = T;
};

/**
 * @brief Value type - types supported in the sparse matrix value
 */
// Dispatched by width, not actual type (e.g., all 8-bit use uint8_t, etc).
// This reduces template instantiation combinatorics (binary size, compile time,
// etc).
#define USE_VALUE_TYPE_WIDTH 0

#if USE_VALUE_TYPE_WIDTH
// TODO: this approach doesn't work, as numpy complains when it attempts to cast
// types (e.g., float32->uint32). What we need to do is change the templating to
// take a "width" int value, rather than a type, and ripple that through the
// entire fastercsx.h codebase.
using ValueType = std::variant<
    type_identity<uint8_t>,
    type_identity<uint16_t>,
    type_identity<uint32_t>,
    type_identity<uint64_t>>;
#else
using ValueType = std::variant<
    type_identity<int8_t>,
    type_identity<int16_t>,
    type_identity<int32_t>,
    type_identity<int64_t>,
    type_identity<uint8_t>,
    type_identity<uint16_t>,
    type_identity<uint32_t>,
    type_identity<uint64_t>,
    type_identity<float>,
    type_identity<double>>;

#endif

/**
 * @brief Index type - types supported as CSx indices.
 *
 * Small (width) types are primarily for memory savings, so this is not an
 * exhaustive list. int32/int64 are important for SciPy compatibility.
 */
using CsxIndexType = std::variant<
    type_identity<int32_t>,
    type_identity<int64_t>,
    type_identity<uint16_t>,
    type_identity<uint32_t>>;

using CooIndexType =
    std::variant<type_identity<int32_t>, type_identity<int64_t>>;

#if USE_VALUE_TYPE_WIDTH
static const std::unordered_map<int, ValueType> value_type_dispatch = {
    {npy_api::NPY_INT8_, type_identity<uint8_t>{}},
    {npy_api::NPY_INT16_, type_identity<uint16_t>{}},
    {npy_api::NPY_INT32_, type_identity<uint32_t>{}},
    {npy_api::NPY_INT64_, type_identity<uint64_t>{}},
    {npy_api::NPY_UINT8_, type_identity<uint8_t>{}},
    {npy_api::NPY_UINT16_, type_identity<uint16_t>{}},
    {npy_api::NPY_UINT32_, type_identity<uint32_t>{}},
    {npy_api::NPY_UINT64_, type_identity<uint64_t>{}},
    {npy_api::NPY_FLOAT_, type_identity<uint32_t>{}},
    {npy_api::NPY_DOUBLE_, type_identity<uint64_t>{}},
};
#else
static const std::unordered_map<int, ValueType> value_type_dispatch = {
    {npy_api::NPY_INT8_, type_identity<int8_t>{}},
    {npy_api::NPY_INT16_, type_identity<int16_t>{}},
    {npy_api::NPY_INT32_, type_identity<int32_t>{}},
    {npy_api::NPY_INT64_, type_identity<int64_t>{}},
    {npy_api::NPY_UINT8_, type_identity<uint8_t>{}},
    {npy_api::NPY_UINT16_, type_identity<uint16_t>{}},
    {npy_api::NPY_UINT32_, type_identity<uint32_t>{}},
    {npy_api::NPY_UINT64_, type_identity<uint64_t>{}},
    {npy_api::NPY_FLOAT_, type_identity<float>{}},
    {npy_api::NPY_DOUBLE_, type_identity<double>{}},
};
#endif

static const std::unordered_map<int, CsxIndexType> csx_index_type_dispatch = {
    {npy_api::NPY_INT32_, type_identity<int32_t>{}},
    {npy_api::NPY_INT64_, type_identity<int64_t>{}},
    {npy_api::NPY_UINT16_, type_identity<uint16_t>{}},
    {npy_api::NPY_UINT32_, type_identity<uint32_t>{}}};

static const std::unordered_map<int, CooIndexType> coo_index_type_dispatch = {
    {npy_api::NPY_INT32_, type_identity<int32_t>{}},
    {npy_api::NPY_INT64_, type_identity<int64_t>{}}};

void compress_coo_validate_args(
    const int64_t n_row,
    const int64_t n_col,
    std::vector<py::array> Ai,
    std::vector<py::array> Aj,
    std::vector<py::array> Ad,
    py::array Bp,
    py::array Bj,
    py::array Bd) {
    /*
    Checks performed:
    1. all arrays ndim == 1
    2. Bd.dtype==Ad.dtype
    3. sum(Ai.size())==sum(Aj.size())==sum(Ad.size())==Bj.size()==Bd.size()
    (this is nnz).
    4. num chunks/items in Ai/Aj/Ad is same size and type
    5. ensure B* are writable
    6. Ensure each element in A* tuples are same type
    */
    if (n_row < 0 || n_col < 0)
        throw std::runtime_error("n_row and n_col must be >= 0");

    auto n_chunks = Ai.size();
    for (auto& vec : {Ai, Aj, Ad}) {
        if (vec.size() != n_chunks)
            throw std::runtime_error(
                "All COO array tuples must contain same number of chunks.");
        for (auto& arr : vec) {
            if (arr.ndim() != 1)
                throw std::runtime_error(
                    "All arrays must be of dimension rank 1.");
            if (arr.dtype().num() != vec[0].dtype().num())
                throw std::runtime_error(
                    "All chunks of COO arrays must be of same type.");
        }
    }
    if (Bp.ndim() != 1 || Bj.ndim() != 1 || Bd.ndim() != 1)
        throw std::runtime_error("All arrays must be of dimension rank 1.");

    for (auto& arr : Ad)
        if (arr.dtype().num() != Bd.dtype().num())
            throw std::runtime_error("All data arrays must be of same type.");

    uint64_t nnz = Bd.size();
    for (auto& vec : {Ai, Aj, Ad}) {
        std::size_t s = std::transform_reduce(
            vec.cbegin(), vec.cend(), 0ul, std::plus<>{}, [](py::array a) {
                return a.size();
            });
        if (s != nnz)
            throw std::runtime_error(
                "All COO arrays must have same size (nnz).");
    }

    if (static_cast<uint64_t>(Bj.size()) != nnz)
        throw std::runtime_error(
            "All output data arrays must have same size (nnz).");
    if (Bp.size() != (n_row + 1))
        throw std::runtime_error("Pointer array size does not match n_rows.");

    if (!Bp.writeable() || !Bj.writeable() || !Bd.writeable())
        throw std::runtime_error("Output arrays must be writable.");

    if (Ai.size() > 0) {
        if (!Ai[0].dtype().is(Aj[0].dtype()))
            throw std::runtime_error("COO index arrays must have same dtype.");
    }
}

template <typename T>
T lookup_dtype(
    const std::unordered_map<int, T>& index,
    const py::dtype& dtype,
    const std::string& array_name) {
    try {
        return index.at(dtype.num());
    } catch (const std::out_of_range& oor) {
        // will bubble up as a ValueError
        throw std::invalid_argument(
            "Unsupported type: " + array_name + " has an unsupported dtype");
    }
}

void compress_coo(
    std::shared_ptr<tiledbsoma::SOMAContext> ctx,
    const std::pair<int64_t, int64_t>& shape,
    py::tuple Ai_,
    py::tuple Aj_,
    py::tuple Ad_,
    py::array Bp,
    py::array Bj,
    py::array Bd) {
    std::vector<py::array> Ai, Aj, Ad;
    try  // convert py::tuple[py::array] to vector[py::array]
    {
        Ai = to_vector<py::array>(Ai_);
        Aj = to_vector<py::array>(Aj_);
        Ad = to_vector<py::array>(Ad_);
    } catch (const py::cast_error& e) {
        throw std::runtime_error(
            std::string("All COO data must be tuple of ndarray (") + e.what() +
            ")");
    }

    const auto [n_row, n_col] = shape;
    const auto nnz = Bd.size();

    compress_coo_validate_args(n_row, n_col, Ai, Aj, Ad, Bp, Bj, Bd);

    CooIndexType coo_index_type = lookup_dtype(
        coo_index_type_dispatch, Ai[0].dtype(), "COO index (row, col) arrays");
    CsxIndexType csx_major_index_type = lookup_dtype(
        csx_index_type_dispatch, Bp.dtype(), "CSx indptr array");
    CsxIndexType csx_minor_index_type = lookup_dtype(
        csx_index_type_dispatch, Bj.dtype(), "CSx indices array");
    ValueType value_type = lookup_dtype(
        value_type_dispatch, Bd.dtype(), "CSx data array");

    std::visit(
        [&](auto value_type,
            auto csx_major_index_type,
            auto csx_minor_index_type,
            auto coo_index_type) {
            using COO_INDEX = typename decltype(coo_index_type)::type;
            using CSX_MAJOR_INDEX =
                typename decltype(csx_major_index_type)::type;
            using CSX_MINOR_INDEX =
                typename decltype(csx_minor_index_type)::type;
            using VALUE = typename decltype(value_type)::type;

            std::vector<tcb::span<COO_INDEX const>> Ai_views, Aj_views;
            std::vector<tcb::span<VALUE const>> Ad_views;
            for (size_t i = 0; i < Ai.size(); ++i) {
                Ai_views.push_back(make_span<COO_INDEX>(Ai[i]));
                Aj_views.push_back(make_span<COO_INDEX>(Aj[i]));
                Ad_views.push_back(make_span<VALUE>(Ad[i]));
            }
            auto Bp_view = make_mutable_span<CSX_MAJOR_INDEX>(Bp);
            auto Bj_view = make_mutable_span<CSX_MINOR_INDEX>(Bj);
            auto Bd_view = make_mutable_span<VALUE>(Bd);

            py::gil_scoped_release release;
            return fastercsx::compress_coo(
                ctx->thread_pool().get(),
                shape,
                nnz,
                Ai_views,
                Aj_views,
                Ad_views,
                Bp_view,
                Bj_view,
                Bd_view);
        },
        value_type,
        csx_major_index_type,
        csx_minor_index_type,
        coo_index_type);
}

void sort_indices(
    std::shared_ptr<tiledbsoma::SOMAContext> ctx,
    py::array Bp,
    py::array Bj,
    py::array Bd) {
    // Error checks first
    //
    if (Bp.ndim() != 1 || Bj.ndim() != 1 || Bd.ndim() != 1)
        throw std::runtime_error("All arrays must be 1D");

    if (!Bp.writeable() || !Bj.writeable() || !Bd.writeable())
        throw std::runtime_error("Output arrays must be writable.");

    // Get dispatch types (TODO: need to throw meaningful errors if missing)
    CsxIndexType csx_major_index_type = lookup_dtype(
        csx_index_type_dispatch, Bp.dtype(), "CSx indptr array");
    CsxIndexType csx_minor_index_type = lookup_dtype(
        csx_index_type_dispatch, Bj.dtype(), "CSx indices array");
    ValueType value_type = lookup_dtype(
        value_type_dispatch, Bd.dtype(), "CSx data array");

    std::visit(
        [&](auto value_type,
            auto csx_major_index_type,
            auto csx_minor_index_type) {
            using CSX_MAJOR_INDEX =
                typename decltype(csx_major_index_type)::type;
            using CSX_MINOR_INDEX =
                typename decltype(csx_minor_index_type)::type;
            using VALUE = typename decltype(value_type)::type;

            auto n_row = Bp.size() - 1;
            int64_t nnz = py::cast<py::array_t<CSX_MAJOR_INDEX>>(Bp).at(n_row);
            if (Bj.size() != nnz || Bd.size() != nnz)
                throw std::runtime_error("Array length and nnz do not match.");

            auto Bp_view = make_span<CSX_MAJOR_INDEX>(Bp);
            auto Bj_view = make_mutable_span<CSX_MINOR_INDEX>(Bj);
            auto Bd_view = make_mutable_span<VALUE>(Bd);

            py::gil_scoped_release release;
            return fastercsx::sort_indices(
                ctx->thread_pool().get(),
                n_row,
                nnz,
                Bp_view,
                Bj_view,
                Bd_view);
        },
        value_type,
        csx_major_index_type,
        csx_minor_index_type);
};

void copy_to_dense(
    std::shared_ptr<tiledbsoma::SOMAContext> ctx,
    const uint64_t major_idx_start,
    const int64_t major_idx_end,
    const std::pair<int64_t, int64_t>& shape,
    const std::string& format,
    py::array Bp,
    py::array Bj,
    py::array Bd,
    py::array out) {
    if (format != "csr" && format != "csc")
        throw std::runtime_error("format must be 'csr' or 'csc'");
    const fastercsx::Format cm_format =
        (format == "csr" ? fastercsx::Format::CSR : fastercsx::Format::CSC);

    auto [n_row, n_col] = shape;
    if (n_row < 0 || n_col < 0)
        throw std::runtime_error("n_row and n_col must be >= 0");
    auto n_major = (format == "csr") ? n_row : n_col;

    if (major_idx_start < 0 || major_idx_end > n_major)
        throw std::runtime_error(
            "row_start must be >= 0 and row_end < array shape");
    if (n_major != Bp.size() - 1)
        throw std::runtime_error("n_rows does not match Bp.size()");
    if (!out.writeable())
        throw std::runtime_error("out must be writable");
    if (out.dtype().num() != Bd.dtype().num())
        throw std::runtime_error("out dtype must match Bd dtype");

    // Get dispatch types (TODO: need to throw meaningful errors if missing)
    CsxIndexType csx_major_index_type = lookup_dtype(
        csx_index_type_dispatch, Bp.dtype(), "CSx indptr array");
    CsxIndexType csx_minor_index_type = lookup_dtype(
        csx_index_type_dispatch, Bj.dtype(), "CSx indices array");
    ValueType value_type = lookup_dtype(
        value_type_dispatch, Bd.dtype(), "CSx data array");

    std::visit(
        [&](auto value_type,
            auto csx_major_index_type,
            auto csx_minor_index_type) {
            using CSX_MAJOR_INDEX =
                typename decltype(csx_major_index_type)::type;
            using CSX_MINOR_INDEX =
                typename decltype(csx_minor_index_type)::type;
            using VALUE = typename decltype(value_type)::type;

            auto Bp_view = make_span<CSX_MAJOR_INDEX>(Bp);
            auto Bj_view = make_span<CSX_MINOR_INDEX>(Bj);
            auto Bd_view = make_span<VALUE>(Bd);
            auto out_view = make_mutable_span<VALUE>(out);

            py::gil_scoped_release release;
            return fastercsx::copy_to_dense(
                ctx->thread_pool().get(),
                major_idx_start,
                major_idx_end,
                shape,
                cm_format,
                Bp_view,
                Bj_view,
                Bd_view,
                out_view);
        },
        value_type,
        csx_major_index_type,
        csx_minor_index_type);
};

void count_rows(
    std::shared_ptr<tiledbsoma::SOMAContext> ctx,
    const int64_t n_row,
    const uint64_t nnz,
    py::tuple Ai_,
    py::array Bp) {
    std::vector<py::array> Ai;
    try  // convert py::tuple[py::array] to vector[py::array]
    {
        Ai = to_vector<py::array>(Ai_);
    } catch (const py::cast_error& e) {
        throw std::runtime_error(
            std::string("All COO data must be tuple of ndarray (") + e.what() +
            ")");
    }

    if (Bp.ndim() != 1)
        throw std::runtime_error("All arrays must be of dimension rank 1.");

    if (Bp.size() != (n_row + 1))
        throw std::runtime_error("Pointer array size does not match n_rows.");

    if (!Bp.writeable())
        throw std::runtime_error("Output arrays must be writable.");

    CooIndexType coo_index_type = lookup_dtype(
        coo_index_type_dispatch, Ai[0].dtype(), "COO index (row, col) arrays");
    CsxIndexType csx_major_index_type = lookup_dtype(
        csx_index_type_dispatch, Bp.dtype(), "CSx indptr array");

    std::visit(
        [&](auto csx_major_index_type, auto coo_index_type) {
            using COO_INDEX = typename decltype(coo_index_type)::type;
            using CSX_MAJOR_INDEX =
                typename decltype(csx_major_index_type)::type;

            std::vector<tcb::span<COO_INDEX const>> Ai_views, Aj_views;
            for (size_t i = 0; i < Ai.size(); ++i) {
                Ai_views.push_back(make_span<COO_INDEX>(Ai[i]));
            }
            auto Bp_view = make_mutable_span<CSX_MAJOR_INDEX>(Bp);

            py::gil_scoped_release release;
            fastercsx::count_rows(
                ctx->thread_pool().get(), n_row, nnz, Ai_views, Bp_view);
        },
        csx_major_index_type,
        coo_index_type);
}

void load_fastercsx(py::module& m) {
    py::module fastercsx_m = m.def_submodule("fastercsx", "CSX primitives");

    fastercsx_m.def(
        "compress_coo",
        compress_coo,
        py::arg("ctx").noconvert(),
        py::arg("shape"),
        py::arg("Ai").noconvert(),
        py::arg("Aj").noconvert(),
        py::arg("Ad").noconvert(),
        py::arg("Bp").noconvert(),
        py::arg("Bj").noconvert(),
        py::arg("Bd").noconvert(),
        "Create CSX elements");

    fastercsx_m.def(
        "sort_indices",
        sort_indices,
        py::arg("ctx").noconvert(),
        py::arg("Bp").noconvert(),
        py::arg("Bj").noconvert(),
        py::arg("Bd").noconvert(),
        "Sort minor axis indices and data");

    fastercsx_m.def(
        "copy_to_dense",
        copy_to_dense,
        py::arg("ctx").noconvert(),
        py::arg("major_idx_start"),
        py::arg("major_idx_end"),
        py::arg("shape"),
        py::arg("format"),
        py::arg("Bp").noconvert(),
        py::arg("Bj").noconvert(),
        py::arg("Bd").noconvert(),
        py::arg("out").noconvert(),
        "Copy major axis slice to dense");

    fastercsx_m.def(
        "count_rows",
        count_rows,
        py::arg("ctx").noconvert(),
        py::arg("n_row"),
        py::arg("nnz"),
        py::arg("Ai").noconvert(),
        py::arg("Bp").noconvert(),
        "Count rows");
}

}  // namespace libtiledbsomacpp
