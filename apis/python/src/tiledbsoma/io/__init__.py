# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

import anndata
import scipy
from packaging.version import Version

from ._registration import (
    ExperimentAmbientLabelMapping,
)
from .ingest import (
    add_matrix_to_collection,
    add_X_layer,
    append_obs,
    append_var,
    append_X,
    create_from_matrix,
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

# https://github.com/single-cell-data/TileDB-SOMA/issues/3920
scipy_version = Version(scipy.__version__)
anndata_version = Version(anndata.__version__)
if anndata_version >= Version("0.10.7") and anndata_version < Version("0.11.2"):
    if scipy_version >= Version("1.15.0"):
        raise RuntimeError(
            f"anndata '{anndata.__version__}' is incompatible with '{scipy.__version__}', and more generally, anndata 0.10.7-0.11.1 are incompatible with scipy 0.15. Please upgrade your anndata to >= 0.11.2, or downgrade your scipy to < 0.15."
        )

__all__ = (
    "add_matrix_to_collection",
    "add_X_layer",
    "append_obs",
    "append_var",
    "append_X",
    "create_from_matrix",
    "from_anndata",
    "from_h5ad",
    "register_h5ads",
    "register_anndatas",
    "to_anndata",
    "to_h5ad",
    "update_matrix",
    "update_obs",
    "update_var",
    "upgrade_experiment_shapes",
    "get_experiment_shapes",
    "show_experiment_shapes",
    "resize_experiment",
    "ExperimentAmbientLabelMapping",
)
