/**
 * @file   transformer.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the TransformerPipeline, Transformer, and bindings for
 * derived classes.
 */

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "common.h"

namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

void load_transformers(py::module& m) {
    py::class_<TransformerPipeline>(m, "TransformerPipeline")
        .def(py::init([](py::handle py_batch) {
            ArrowSchema arrow_schema;
            ArrowArray arrow_array;
            uintptr_t arrow_schema_ptr = (uintptr_t)(&arrow_schema);
            uintptr_t arrow_array_ptr = (uintptr_t)(&arrow_array);
            py_batch.attr("_export_to_c")(arrow_array_ptr, arrow_schema_ptr);

            auto array = std::make_unique<ArrowArray>(arrow_array);
            auto schema = std::make_unique<ArrowSchema>(arrow_schema);

            return TransformerPipeline(std::move(array), std::move(schema));
        }))
        .def(
            "transform",
            [](TransformerPipeline& pipeline,
               std::shared_ptr<Transformer> transformer)
                -> TransformerPipeline& {
                return pipeline.transform(transformer);
            })
        .def("asTable", [](TransformerPipeline& pipeline) {
            auto pa = py::module::import("pyarrow");
            auto pa_table_from_arrays = pa.attr("Table").attr("from_arrays");
            auto pa_array_import = pa.attr("Array").attr("_import_from_c");
            auto pa_schema_import = pa.attr("Schema").attr("_import_from_c");

            auto [array, schema] = pipeline.asTable();

            py::list array_list;
            py::list names;

            for (int64_t i = 0; i < schema->n_children; ++i) {
                // Should happen before pyarrow array construction because
                // py::capsule gets ownership of the memory
                names.append(std::string(schema->children[i]->name));

                auto pa_array = pa_array_import(
                    py::capsule(array->children[i]),
                    py::capsule(schema->children[i]));

                array_list.append(pa_array);
            }

            return pa_table_from_arrays(array_list, names);
        });

    py::class_<Transformer, std::shared_ptr<Transformer>>(m, "Transformer");
    py::class_<
        OutlineTransformer,
        Transformer,
        std::shared_ptr<OutlineTransformer>>(m, "OutlineTransformer")
        .def(py::init([](std::string coord_space) {
            auto coordinate_space = SOMACoordinateSpace::from_string(
                coord_space);

            return std::make_shared<OutlineTransformer>(coordinate_space);
        }));
}
}  // namespace libtiledbsomacpp
