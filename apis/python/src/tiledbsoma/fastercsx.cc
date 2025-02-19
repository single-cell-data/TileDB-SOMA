/**
 * @file   fastercsx.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * Python bindings for CSX conversion primitives.
 */
#include <bit>
#include <variant>

// Define to include extra debugging bindings (e.g., count_rows)
/* #define FASTERCSX__DEBUGGING_HOOKS 1 */

#include <tiledbsoma/utils/fastercsx.h>
#include <span>
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
std::vector<T> to_vector_(const py::tuple& tup) {
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

/**
 * @brief Various helpers which convert from py:array to spans.
 */

template <typename T>
std::span<T const> make_span_(py::array arr) {
    assert(py::isinstance<py::array_t<T>>(arr));
    return std::span<T const>(arr.unchecked<T, 1>().data(0), arr.size());
}

template <typename T>
std::span<T> make_mutable_span_(py::array arr) {
    assert(py::isinstance<py::array_t<T>>(arr));
    return std::span<T>(
        arr.mutable_unchecked<T, 1>().mutable_data(0), arr.size());
}

template <typename T, typename R>
std::span<R const> make_casted_span_(py::array arr) {
    static_assert(sizeof(T) == sizeof(R));
    assert(py::isinstance<py::array_t<T>>(arr));
    std::remove_cv_t<T>* p = (std::remove_cv_t<T>*)arr
                                 .unchecked<std::remove_cv_t<T>, 1>()
                                 .data(0);
    return std::span<R const>(reinterpret_cast<R*>(p), arr.size());
}

template <typename T, typename R>
std::span<R> make_mutable_casted_span_(py::array arr) {
    static_assert(sizeof(T) == sizeof(R));
    assert(py::isinstance<py::array_t<T>>(arr));
    std::remove_cv_t<T>* p = (std::remove_cv_t<T>*)arr
                                 .mutable_unchecked<std::remove_cv_t<T>, 1>()
                                 .data(0);
    return std::span<R>(reinterpret_cast<R*>(p), arr.size());
}

/**
 * @brief Return true if the NP byteorder is native (or equivalent).
 */
bool is_native_byteorder(const char byteorder) {
    if (byteorder == '=')  // native
        return true;
    if (byteorder == '|')  // not-applicable
        return true;
    if constexpr (std::endian::native == std::endian::big)
        return byteorder == '>';  // big
    else
        return byteorder == '<';  // little
}

/**
 * @brief Check enddianness/byteorder is native, and raise exception if not.
 * Necessary becuase we dispatch on dtype().num(), which doesn't confirm
 * byteorder is native.
 */
void check_byteorder(const py::dtype& dtype) {
    if (!is_native_byteorder(dtype.byteorder()))
        throw invalid_argument(
            "All arrays must have native byteorder (endianness).");
}

/*
 * Value/data arrays are cast to an unsigned of the same width as the actual
 * value type. This is solely to reduce the combinatorics of template
 * instantiation, improving final code size and compile time.
 */
template <typename T>
struct remap_value {
    typedef std::make_unsigned_t<T> type;
};

template <>
struct remap_value<float> {
    typedef uint32_t type;
};

template <>
struct remap_value<const float> {
    typedef const uint32_t type;
};

template <>
struct remap_value<double> {
    typedef uint64_t type;
};

template <>
struct remap_value<const double> {
    typedef const uint64_t type;
};

template <class T>
using remap_value_t = typename remap_value<T>::type;

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

static const std::unordered_map<int, CsxIndexType> csx_index_type_dispatch = {
    {npy_api::NPY_INT32_, type_identity<int32_t>{}},
    {npy_api::NPY_INT64_, type_identity<int64_t>{}},
    {npy_api::NPY_UINT16_, type_identity<uint16_t>{}},
    {npy_api::NPY_UINT32_, type_identity<uint32_t>{}}};

static const std::unordered_map<int, CooIndexType> coo_index_type_dispatch = {
    {npy_api::NPY_INT32_, type_identity<int32_t>{}},
    {npy_api::NPY_INT64_, type_identity<int64_t>{}}};

template <typename T>
T lookup_dtype_(
    const std::unordered_map<int, T>& index,
    const py::dtype& dtype,
    const std::string& array_name) {
    check_byteorder(dtype);
    try {
        return index.at(dtype.num());
    } catch (const std::out_of_range& oor) {
        // will bubble up as a ValueError
        std::string name = dtype.attr("name").cast<std::string>();
        throw std::invalid_argument(
            "Unsupported type: " + array_name + " has an unsupported dtype: '" +
            name + "'");
    }
}

/**
 * @brief Perform all checks required to ensure safe access to caller-provided
 * array data for the `coo_compress` function.
 */
void compress_coo_validate_args_(
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
    5. ensure B* are writeable
    6. Ensure each element in A* tuples are same type
    7. Ensure each element in the A* tuples are the same length
    8. byteorder
    etc...

    Not checked:
    1. C-style vs F-style - as all arrays are 1D, there is no need to check.
    */
    if (n_row < 0 || n_col < 0)
        throw std::range_error("n_row and n_col must be >= 0");

    auto n_chunks = Ai.size();
    for (auto& vec : {Ai, Aj, Ad}) {
        if (vec.size() != n_chunks)
            throw std::length_error(
                "All COO array tuples must contain same number of chunks.");
        for (auto& arr : vec) {
            if (arr.ndim() != 1)
                throw std::length_error(
                    "All arrays must be of dimension rank 1.");
            if (arr.dtype().num() != vec[0].dtype().num())
                throw pybind11::type_error(
                    "All chunks of COO arrays must be of same type.");
            check_byteorder(arr.dtype());
        }
    }
    for (uint64_t chunk_idx = 0; chunk_idx < n_chunks; chunk_idx++) {
        if ((Ai[chunk_idx].size() != Aj[chunk_idx].size()) ||
            (Ai[chunk_idx].size() != Ad[chunk_idx].size()))
            throw std::length_error(
                "All COO array tuple elements must be of the same size.");
    }

    if (Bp.ndim() != 1 || Bj.ndim() != 1 || Bd.ndim() != 1)
        throw std::length_error("All arrays must be of dimension rank 1.");

    check_byteorder(Bp.dtype());
    check_byteorder(Bj.dtype());
    check_byteorder(Bd.dtype());

    for (auto& arr : Ad)
        if (arr.dtype().num() != Bd.dtype().num())
            throw pybind11::type_error("All data arrays must be of same type.");

    uint64_t nnz = Bd.size();
    for (auto& vec : {Ai, Aj, Ad}) {
        std::size_t s = std::transform_reduce(
            vec.cbegin(), vec.cend(), 0ul, std::plus<>{}, [](py::array a) {
                return a.size();
            });
        if (s != nnz)
            throw std::length_error(
                "All COO arrays must have same size (nnz).");
    }

    if (static_cast<uint64_t>(Bj.size()) != nnz)
        throw std::length_error(
            "All output data arrays must have same size (nnz).");
    if (Bp.size() != (n_row + 1))
        throw std::length_error("Pointer array size does not match n_rows.");

    if (!Bp.writeable() || !Bj.writeable() || !Bd.writeable())
        throw std::invalid_argument("Output arrays must be writeable.");

    if (Ai.size() > 0) {
        if (!Ai[0].dtype().is(Aj[0].dtype()))
            throw pybind11::type_error(
                "COO index arrays must have same dtype.");
    }
}

/**
 * @brief python binding for coo_compress, which converts COO chunked array (I,
 * J, D) into CSX arrays (P, J, D). The output arrays are NOT sorted on the
 * minor dimension -- see `sort_csx_indices` if that is required.
 */
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
    try {
        // convert py::tuple[py::array] to vector[py::array]
        Ai = to_vector_<py::array>(Ai_);
        Aj = to_vector_<py::array>(Aj_);
        Ad = to_vector_<py::array>(Ad_);
    } catch (const py::cast_error& e) {
        throw pybind11::type_error(
            std::string("All COO data must be tuple of ndarray (") + e.what() +
            ")");
    }

    const auto [n_row, n_col] = shape;
    const auto nnz = Bd.size();

    // Check all user-provided values for constraints, typing, etc.
    compress_coo_validate_args_(n_row, n_col, Ai, Aj, Ad, Bp, Bj, Bd);

    // Lookup the underlying C++ primitive type for tag dispatch.
    CooIndexType coo_index_type = lookup_dtype_(
        coo_index_type_dispatch, Ai[0].dtype(), "COO index (row, col) arrays");
    CsxIndexType csx_major_index_type = lookup_dtype_(
        csx_index_type_dispatch, Bp.dtype(), "CSx indptr array");
    CsxIndexType csx_minor_index_type = lookup_dtype_(
        csx_index_type_dispatch, Bj.dtype(), "CSx indices array");
    ValueType value_type = lookup_dtype_(
        value_type_dispatch, Bd.dtype(), "CSx data array");

    // Dispatch by type
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

            std::vector<std::span<COO_INDEX const>> Ai_views, Aj_views;
            std::vector<std::span<remap_value_t<VALUE> const>> Ad_views;
            for (size_t i = 0; i < Ai.size(); ++i) {
                Ai_views.push_back(make_span_<COO_INDEX>(Ai[i]));
                Aj_views.push_back(make_span_<COO_INDEX>(Aj[i]));
                Ad_views.push_back(
                    make_casted_span_<VALUE, remap_value_t<VALUE>>(Ad[i]));
            }
            auto Bp_view = make_mutable_span_<CSX_MAJOR_INDEX>(Bp);
            auto Bj_view = make_mutable_span_<CSX_MINOR_INDEX>(Bj);
            auto Bd_view =
                make_mutable_casted_span_<VALUE, remap_value_t<VALUE>>(Bd);

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

/**
 * @brief Python binding for sort_csx_indices, which sorts minor dimension of a
 * compressed matrix.
 *
 * Returns false if the matrix contains duplicate coordinates, true if all
 * coordinates are unique.
 */
bool sort_csx_indices(
    std::shared_ptr<tiledbsoma::SOMAContext> ctx,
    py::array Bp,
    py::array Bj,
    py::array Bd) {
    // Error checks first
    //
    if (Bp.ndim() != 1 || Bj.ndim() != 1 || Bd.ndim() != 1)
        throw std::length_error("All arrays must be 1D");
    if (!Bp.writeable() || !Bj.writeable() || !Bd.writeable())
        throw std::invalid_argument("Output arrays must be writeable.");

    check_byteorder(Bp.dtype());
    check_byteorder(Bj.dtype());
    check_byteorder(Bd.dtype());

    // Get dispatch types
    CsxIndexType csx_major_index_type = lookup_dtype_(
        csx_index_type_dispatch, Bp.dtype(), "CSx indptr array");
    CsxIndexType csx_minor_index_type = lookup_dtype_(
        csx_index_type_dispatch, Bj.dtype(), "CSx indices array");
    ValueType value_type = lookup_dtype_(
        value_type_dispatch, Bd.dtype(), "CSx data array");

    // Dispatch by type
    auto no_duplicates = std::visit(
        [&](auto value_type,
            auto csx_major_index_type,
            auto csx_minor_index_type) -> bool {
            using CSX_MAJOR_INDEX =
                typename decltype(csx_major_index_type)::type;
            using CSX_MINOR_INDEX =
                typename decltype(csx_minor_index_type)::type;
            using VALUE = typename decltype(value_type)::type;

            auto n_row = Bp.size() - 1;
            int64_t nnz = py::cast<py::array_t<CSX_MAJOR_INDEX>>(Bp).at(n_row);
            if (Bj.size() != nnz || Bd.size() != nnz)
                throw std::length_error("Array length and nnz do not match.");

            auto Bp_view = make_span_<CSX_MAJOR_INDEX>(Bp);
            auto Bj_view = make_mutable_span_<CSX_MINOR_INDEX>(Bj);
            auto Bd_view =
                make_mutable_casted_span_<VALUE, remap_value_t<VALUE>>(Bd);

            py::gil_scoped_release release;
            return fastercsx::sort_csx_indices(
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

    return no_duplicates;
};

/**
 * @brief Python binding for parallel slice+to_dense operation.
 */
void copy_csx_to_dense(
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
        throw std::invalid_argument("format must be 'csr' or 'csc'");
    const fastercsx::Format cm_format =
        (format == "csr" ? fastercsx::Format::CSR : fastercsx::Format::CSC);

    auto [n_row, n_col] = shape;
    if (n_row < 0 || n_col < 0)
        throw std::length_error("n_row and n_col must be >= 0");
    auto n_major = (format == "csr") ? n_row : n_col;

    if (major_idx_start < 0 || major_idx_end > n_major)
        throw std::range_error(
            "row_start must be >= 0 and row_end < array shape");
    if (n_major != Bp.size() - 1)
        throw std::length_error("n_rows does not match Bp.size()");
    if (!out.writeable())
        throw std::invalid_argument("out must be writeable");
    if (out.dtype().num() != Bd.dtype().num())
        throw pybind11::type_error("out dtype must match Bd dtype");

    // Get dispatch types
    CsxIndexType csx_major_index_type = lookup_dtype_(
        csx_index_type_dispatch, Bp.dtype(), "CSx indptr array");
    CsxIndexType csx_minor_index_type = lookup_dtype_(
        csx_index_type_dispatch, Bj.dtype(), "CSx indices array");
    ValueType value_type = lookup_dtype_(
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

            auto Bp_view = make_span_<CSX_MAJOR_INDEX>(Bp);
            auto Bj_view = make_span_<CSX_MINOR_INDEX>(Bj);
            auto Bd_view = make_casted_span_<VALUE, remap_value_t<VALUE>>(Bd);
            auto out_view =
                make_mutable_casted_span_<VALUE, remap_value_t<VALUE>>(out);

            py::gil_scoped_release release;
            return fastercsx::copy_csx_to_dense(
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

#ifdef FASTERCSX__DEBUGGING_HOOKS
void count_rows(
    std::shared_ptr<tiledbsoma::SOMAContext> ctx,
    const int64_t n_row,
    const uint64_t nnz,
    py::tuple Ai_,
    py::array Bp) {
    std::vector<py::array> Ai;
    try  // convert py::tuple[py::array] to vector[py::array]
    {
        Ai = to_vector_<py::array>(Ai_);
    } catch (const py::cast_error& e) {
        throw pybind11::type_error(
            std::string("All COO data must be tuple of ndarray (") + e.what() +
            ")");
    }

    if (Bp.ndim() != 1)
        throw std::range_error("All arrays must be of dimension rank 1.");

    if (Bp.size() != (n_row + 1))
        throw std::length_error("Pointer array size does not match n_rows.");

    if (!Bp.writeable())
        throw std::invalid_argument("Output arrays must be writeable.");

    CooIndexType coo_index_type = lookup_dtype_(
        coo_index_type_dispatch, Ai[0].dtype(), "COO index (row, col) arrays");
    CsxIndexType csx_major_index_type = lookup_dtype_(
        csx_index_type_dispatch, Bp.dtype(), "CSx indptr array");

    std::visit(
        [&](auto csx_major_index_type, auto coo_index_type) {
            using COO_INDEX = typename decltype(coo_index_type)::type;
            using CSX_MAJOR_INDEX =
                typename decltype(csx_major_index_type)::type;

            std::vector<std::span<COO_INDEX const>> Ai_views, Aj_views;
            for (size_t i = 0; i < Ai.size(); ++i) {
                Ai_views.push_back(make_span_<COO_INDEX>(Ai[i]));
            }
            auto Bp_view = make_mutable_span_<CSX_MAJOR_INDEX>(Bp);

            py::gil_scoped_release release;
            fastercsx::count_rows(
                ctx->thread_pool().get(), n_row, nnz, Ai_views, Bp_view);
        },
        csx_major_index_type,
        coo_index_type);
}
#endif

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
        "sort_csx_indices",
        sort_csx_indices,
        py::arg("ctx").noconvert(),
        py::arg("Bp").noconvert(),
        py::arg("Bj").noconvert(),
        py::arg("Bd").noconvert(),
        "Sort minor axis indices and data");

    fastercsx_m.def(
        "copy_csx_to_dense",
        copy_csx_to_dense,
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

#ifdef FASTERCSX__DEBUGGING_HOOKS
    fastercsx_m.def(
        "count_rows",
        count_rows,
        py::arg("ctx").noconvert(),
        py::arg("n_row"),
        py::arg("nnz"),
        py::arg("Ai").noconvert(),
        py::arg("Bp").noconvert(),
        "Count rows");
#endif
}

}  // namespace libtiledbsomacpp
