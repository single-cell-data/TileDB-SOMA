# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.


from ._registration import (
    ExperimentAmbientLabelMapping,
)
from .ingest import (
    add_matrix_to_collection,
    add_X_layer,
    from_anndata,
    from_h5ad,
    register_anndatas,
    register_h5ads,
    update_matrix,
    update_obs,
    update_var,
)
from .outgest import (
    to_anndata,
    to_h5ad,
)
from .shaping import (
    get_experiment_shapes,
    resize_experiment,
    show_experiment_shapes,
    upgrade_experiment_shapes,
)

__all__ = (
    "ExperimentAmbientLabelMapping",
    "add_X_layer",
    "add_matrix_to_collection",
    "from_anndata",
    "from_h5ad",
    "get_experiment_shapes",
    "register_anndatas",
    "register_h5ads",
    "resize_experiment",
    "show_experiment_shapes",
    "to_anndata",
    "to_h5ad",
    "update_matrix",
    "update_obs",
    "update_var",
    "upgrade_experiment_shapes",
)
