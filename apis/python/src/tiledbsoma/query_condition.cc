/**
 * @file   query_condition.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file implements the TileDB-Py query condition.
 */

#include "common.h"

#if TILEDB_VERSION_MAJOR == 2 && TILEDB_VERSION_MINOR >= 2
#if !defined(NDEBUG)
#endif
namespace libtiledbsomacpp {

using namespace std;
using namespace tiledb;
namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

void load_query_condition(py::module& m) {
    py::class_<PyQueryCondition>(m, "PyQueryCondition", py::module_local())
        .def(py::init<py::object>(), py::arg("ctx") = py::none())

        /* TODO surely there's a better way to deal with templated PyBind11
         * functions? but maybe not?
         * https://github.com/pybind/pybind11/issues/1667
         */

        .def(
            "init_string",
            static_cast<void (PyQueryCondition::*)(
                const string&, const string&, tiledb_query_condition_op_t)>(
                &PyQueryCondition::init))

        .def(
            "init_uint64",
            static_cast<void (PyQueryCondition::*)(
                const string&, uint64_t, tiledb_query_condition_op_t)>(
                &PyQueryCondition::init))
        .def(
            "init_int64",
            static_cast<void (PyQueryCondition::*)(
                const string&, int64_t, tiledb_query_condition_op_t)>(
                &PyQueryCondition::init))
        .def(
            "init_uint32",
            static_cast<void (PyQueryCondition::*)(
                const string&, uint32_t, tiledb_query_condition_op_t)>(
                &PyQueryCondition::init))
        .def(
            "init_int32",
            static_cast<void (PyQueryCondition::*)(
                const string&, int32_t, tiledb_query_condition_op_t)>(
                &PyQueryCondition::init))
        .def(
            "init_uint16",
            static_cast<void (PyQueryCondition::*)(
                const string&, uint16_t, tiledb_query_condition_op_t)>(
                &PyQueryCondition::init))
        .def(
            "init_int16",
            static_cast<void (PyQueryCondition::*)(
                const string&, int16_t, tiledb_query_condition_op_t)>(
                &PyQueryCondition::init))
        .def(
            "init_uint8",
            static_cast<void (PyQueryCondition::*)(
                const string&, uint8_t, tiledb_query_condition_op_t)>(
                &PyQueryCondition::init))
        .def(
            "init_int8",
            static_cast<void (PyQueryCondition::*)(
                const string&, int8_t, tiledb_query_condition_op_t)>(
                &PyQueryCondition::init))

        .def(
            "init_float32",
            static_cast<void (PyQueryCondition::*)(
                const string&, float, tiledb_query_condition_op_t)>(
                &PyQueryCondition::init))
        .def(
            "init_float64",
            static_cast<void (PyQueryCondition::*)(
                const string&, double, tiledb_query_condition_op_t)>(
                &PyQueryCondition::init))
        .def(
            "init_bool",
            static_cast<void (PyQueryCondition::*)(
                const string&, bool, tiledb_query_condition_op_t)>(
                &PyQueryCondition::init))

        .def("combine", &PyQueryCondition::combine)

        .def_static(
            "create_string",
            static_cast<PyQueryCondition (*)(
                const std::string&,
                const std::vector<std::string>&,
                tiledb_query_condition_op_t)>(&PyQueryCondition::create))
        .def_static(
            "create_uint64",
            static_cast<PyQueryCondition (*)(
                const std::string&,
                const std::vector<uint64_t>&,
                tiledb_query_condition_op_t)>(&PyQueryCondition::create))
        .def_static(
            "create_int64",
            static_cast<PyQueryCondition (*)(
                const std::string&,
                const std::vector<int64_t>&,
                tiledb_query_condition_op_t)>(&PyQueryCondition::create))
        .def_static(
            "create_uint32",
            static_cast<PyQueryCondition (*)(
                const std::string&,
                const std::vector<uint32_t>&,
                tiledb_query_condition_op_t)>(&PyQueryCondition::create))
        .def_static(
            "create_int32",
            static_cast<PyQueryCondition (*)(
                const std::string&,
                const std::vector<int32_t>&,
                tiledb_query_condition_op_t)>(&PyQueryCondition::create))
        .def_static(
            "create_uint16",
            static_cast<PyQueryCondition (*)(
                const std::string&,
                const std::vector<uint16_t>&,
                tiledb_query_condition_op_t)>(&PyQueryCondition::create))
        .def_static(
            "create_int8",
            static_cast<PyQueryCondition (*)(
                const std::string&,
                const std::vector<int8_t>&,
                tiledb_query_condition_op_t)>(&PyQueryCondition::create))
        .def_static(
            "create_uint16",
            static_cast<PyQueryCondition (*)(
                const std::string&,
                const std::vector<uint16_t>&,
                tiledb_query_condition_op_t)>(&PyQueryCondition::create))
        .def_static(
            "create_int8",
            static_cast<PyQueryCondition (*)(
                const std::string&,
                const std::vector<int8_t>&,
                tiledb_query_condition_op_t)>(&PyQueryCondition::create))
        .def_static(
            "create_float32",
            static_cast<PyQueryCondition (*)(
                const std::string&,
                const std::vector<float>&,
                tiledb_query_condition_op_t)>(&PyQueryCondition::create))
        .def_static(
            "create_float64",
            static_cast<PyQueryCondition (*)(
                const std::string&,
                const std::vector<double>&,
                tiledb_query_condition_op_t)>(&PyQueryCondition::create))

        .def("__capsule__", &PyQueryCondition::__capsule__);

    py::enum_<tiledb_query_condition_op_t>(
        m, "tiledb_query_condition_op_t", py::arithmetic(), py::module_local())
        .value("TILEDB_LT", TILEDB_LT)
        .value("TILEDB_LE", TILEDB_LE)
        .value("TILEDB_GT", TILEDB_GT)
        .value("TILEDB_GE", TILEDB_GE)
        .value("TILEDB_EQ", TILEDB_EQ)
        .value("TILEDB_NE", TILEDB_NE)
        .value("TILEDB_IN", TILEDB_IN)
        .value("TILEDB_NOT_IN", TILEDB_NOT_IN)
        .export_values();

    py::enum_<tiledb_query_condition_combination_op_t>(
        m,
        "tiledb_query_condition_combination_op_t",
        py::arithmetic(),
        py::module_local())
        .value("TILEDB_AND", TILEDB_AND)
        .value("TILEDB_OR", TILEDB_OR)
        .export_values();
}
}  // namespace libtiledbsomacpp

#endif
