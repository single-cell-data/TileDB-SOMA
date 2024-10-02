#include <tiledbsoma/tiledbsoma>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include "common.h"

namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void load_soma_context(py::module&);
void load_soma_object(py::module&);
void load_soma_array(py::module&);
void load_soma_dataframe(py::module&);
void load_soma_point_cloud_dataframe(py::module&);
void load_soma_dense_ndarray(py::module&);
void load_soma_sparse_ndarray(py::module&);
void load_soma_group(py::module&);
void load_soma_collection(py::module&);
void load_query_condition(py::module&);
void load_reindexer(py::module&);

PYBIND11_MODULE(pytiledbsoma, m) {
    py::register_exception<TileDBSOMAError>(m, "SOMAError");

    /* We need to make sure C++ TileDBSOMAError is translated to a
     * correctly-typed Python error
     *
     * We're aware of
     * https://pybind11.readthedocs.io/en/stable/advanced/exceptions.html
     * -- we find empirically that despite this translator, we still
     * find it necessary to do explicit catch-and-rethrow within our
     * pybind11 functions. See also
     * https://github.com/single-cell-data/TileDB-SOMA/pull/2963
     */
    py::register_exception_translator([](std::exception_ptr p) {
        auto tiledb_soma_error = (py::object)py::module::import("tiledbsoma")
                                     .attr("SOMAError");

        try {
            if (p)
                std::rethrow_exception(p);
        } catch (const TileDBSOMAError& e) {
            PyErr_SetString(tiledb_soma_error.ptr(), e.what());
        } catch (py::builtin_exception& e) {
            throw;
        };
    });

    py::enum_<OpenMode>(m, "OpenMode")
        .value("read", OpenMode::read)
        .value("write", OpenMode::write);

    py::enum_<ResultOrder>(m, "ResultOrder")
        .value("automatic", ResultOrder::automatic)
        .value("rowmajor", ResultOrder::rowmajor)
        .value("colmajor", ResultOrder::colmajor);

    py::enum_<URIType>(m, "URIType")
        .value("automatic", URIType::automatic)
        .value("absolute", URIType::absolute)
        .value("relative", URIType::relative);

    m.doc() = "SOMA acceleration library";

    m.def("version", []() { return tiledbsoma::version::as_string(); });

    m.def(
        "config_logging",
        [](const std::string& level, const std::string& logfile) {
            LOG_CONFIG(level, logfile);
        },
        "level"_a,
        "logfile"_a = "");

    m.def("info", &LOG_INFO, "message"_a = "");
    m.def("debug", &LOG_DEBUG, "message"_a = "");

    m.def(
        "tiledbsoma_stats_enable",
        []() { tiledbsoma::stats::enable(); },
        "Enable TileDB internal statistics. Lifecycle: experimental.");
    m.def(
        "tiledbsoma_stats_disable",
        []() { tiledbsoma::stats::disable(); },
        "Disable TileDB internal statistics. Lifecycle: experimental.");
    m.def(
        "tiledbsoma_stats_reset",
        []() { tiledbsoma::stats::reset(); },
        "Reset all TileDB internal statistics to 0. Lifecycle: experimental.");
    m.def(
        "tiledbsoma_stats_dump",
        []() {
            py::print(tiledbsoma::version::as_string());
            std::string stats = tiledbsoma::stats::dump();
            py::print(stats);
        },
        "Print TileDB internal statistics. Lifecycle: experimental.");
    m.def(
        "tiledbsoma_stats_string",
        []() {
            std::string stats = tiledbsoma::stats::dump();
            return stats;
        },
        "Print TileDB internal statistics. Lifecycle: experimental.");

    py::class_<PlatformConfig>(m, "PlatformConfig")
        .def(py::init<>())
        .def_readwrite(
            "dataframe_dim_zstd_level",
            &PlatformConfig::dataframe_dim_zstd_level)
        .def_readwrite(
            "sparse_nd_array_dim_zstd_level",
            &PlatformConfig::sparse_nd_array_dim_zstd_level)
        .def_readwrite(
            "dense_nd_array_dim_zstd_level",
            &PlatformConfig::sparse_nd_array_dim_zstd_level)
        .def_readwrite("write_X_chunked", &PlatformConfig::write_X_chunked)
        .def_readwrite("goal_chunk_nnz", &PlatformConfig::goal_chunk_nnz)
        .def_readwrite("remote_cap_nbytes", &PlatformConfig::remote_cap_nbytes)
        .def_readwrite("capacity", &PlatformConfig::capacity)
        .def_readwrite("offsets_filters", &PlatformConfig::offsets_filters)
        .def_readwrite("validity_filters", &PlatformConfig::validity_filters)
        .def_readwrite("attrs", &PlatformConfig::attrs)
        .def_readwrite("dims", &PlatformConfig::dims)
        .def_readwrite("allows_duplicates", &PlatformConfig::allows_duplicates)
        .def_readwrite("tile_order", &PlatformConfig::tile_order)
        .def_readwrite("cell_order", &PlatformConfig::cell_order)
        .def_readwrite(
            "consolidate_and_vacuum", &PlatformConfig::consolidate_and_vacuum);

    m.def("_update_dataframe_schema", &SOMADataFrame::update_dataframe_schema);

    load_soma_context(m);
    load_soma_object(m);
    load_soma_array(m);
    load_soma_dataframe(m);
    load_soma_dense_ndarray(m);
    load_soma_sparse_ndarray(m);
    load_soma_point_cloud_dataframe(m);
    load_soma_group(m);
    load_soma_collection(m);
    load_query_condition(m);
    load_reindexer(m);
}

};  // namespace libtiledbsomacpp
