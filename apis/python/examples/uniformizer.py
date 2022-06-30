#!/usr/bin/env python

"""
The following shows an example of how to take multiple H5AD files and add them to a single SOMA
Collection.  This is a simple example that assumes all source H5ADs comply with the [cellxgene 2.0
schema](https://pypi.org/project/cellxgene-schema/), although you can modify this to conform with
your own organization's schema.

To keep the example simple, any data outside that schema are discarded.  The example peforms three
steps:

* retains `obs` and `var` annotations specifically defined in the cellxgene 2.0 schema;
* discards transformed `X` values, and retains raw `X` counts;
* creates a trivial `rankit` normalization of `X` values to demonstrate an additional `X` layer.

Examples:

    python uniformizer.py -v atlas add-h5ad data1.h5ad
    python uniformizer.py -v atlas add-h5ad data2.h5ad
    python uniformizer.py -v atlas add-soma soma3 # If you've already ingested into a separate SOMA

Note: This code for populating an atlas is independent of querying an atlas.  See also
[https://github.com/single-cell-data/TileDB-SingleCell](https://github.com/single-cell-data/TileDB-SingleCell)
for query examples.
"""

import argparse
import logging
import os.path
import sys
from typing import Optional

import anndata
import numba as nb
import numpy as np
import scipy.sparse
import scipy.stats
import tiledb

import tiledbsc
import tiledbsc.io

# ================================================================
# MAIN ENTRYPOINT

logger = logging.getLogger("tiledbsc")


# ----------------------------------------------------------------
def main() -> int:
    parser = _create_args_parser()
    args = parser.parse_args()

    if args.verbose:
        tiledbsc.logging.info()

    uniformizer = Uniformizer(args.atlas_uri, args.allow_non_primary_data)

    if "func_name" not in args:
        logging.error("Error: unknown subcommand.")
        parser.print_help()
        return 1

    if args.func_name == "add-h5ad":
        return uniformizer.add_h5ad(args.dataset_id, args.h5ad)
    elif args.func_name == "add-soma":
        return uniformizer.add_soma(args.dataset_id, args.soma)
    else:
        raise Exception(
            f'Internal coding error: handler for "{args.func_name}" not found.'
        )


# ----------------------------------------------------------------
def _create_args_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create a uni-modal single-cell atlas."
    )
    parser.add_argument("atlas_uri", type=str, help="atlas (SOCO) URI (required)")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="verbose logging",
    )

    # An atlas, by definition, contains no duplicate cells. This voids-the-warranty
    # flag lets the user bypass the is-primary-data check -- only for testing/debug purposes.
    parser.add_argument(
        "--allow-non-primary-data",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )

    subparsers = parser.add_subparsers()

    sp = subparsers.add_parser("add-h5ad", help="Add an H5AD to the atlas")
    sp.add_argument("h5ad", type=str, help="H5AD URI")
    sp.add_argument(
        "--dataset-id",
        type=str,
        help="Specify H5AD name/ID. (default: H5AD URI basename)",
    )
    sp.set_defaults(func_name="add-h5ad")

    sp = subparsers.add_parser("add-soma", help="Add a SOMA to the atlas")
    sp.add_argument("soma", type=str, help="SOMA URI")
    sp.add_argument(
        "--dataset-id",
        type=str,
        help="Specify SOMA name/ID. (default: SOMA URI basename)",
    )
    sp.set_defaults(func_name="add-soma")

    return parser


# ================================================================
class Uniformizer:
    """
    Makes an atlas.

    Demo code to make a _simple_ (DEMO ONLY) uni-modal (scRNAseq) atlas from
    H5AD files or SOMAs having certain required columns.

    * Use a SOMACollection for the top-level container (a.k.a. "atlas")
    * Each SOMA in the container contains uniformized data (i.e. they have a common schema)

    IMPORTANT CAVEAT: this implements a simplified (uni-modal) atlas -- both the schema and code
    lack multi-modal support. When added, there will be additional schema complexity, such as nested
    collections.

    Method -- for each H5AD/SOMA added:
    1. Load H5AD/SOMA
    2. Drop all AnnData slots that are not part of the atlas:
      * Retain obs, raw.var, and raw.X
      * Drop all others, eg, drop obsm, obsp, varm, varp, uns, X, etc.
    3. Drop all columns from obs/var DataFrames which are not defined in a list of required columns.
    4. Create a new SOMA
    5. Add rankit layer to X
    6. Add the new SOMA to the SOCO
    """

    ctx: tiledb.Ctx
    atlas_uri: str

    # You can adapt these to match your organization's schema
    OBS_COLUMNS = [
        "assay_ontology_term_id",
        "cell_type_ontology_term_id",
        "development_stage_ontology_term_id",
        "disease_ontology_term_id",
        "ethnicity_ontology_term_id",
        "organism_ontology_term_id",
        "sex_ontology_term_id",
        "tissue_ontology_term_id",
        "assay",
        "cell_type",
        "development_stage",
        "disease",
        "ethnicity",
        "organism",
        "sex",
        "tissue",
    ]
    VAR_COLUMNS = ["feature_biotype", "feature_name", "feature_reference"]

    # ----------------------------------------------------------------
    def __init__(self, atlas_uri: str, allow_non_primary_data: bool = False):
        self.ctx = self._create_tiledb_ctx()
        self.atlas_uri = atlas_uri
        self.allow_non_primary_data = allow_non_primary_data

    # ----------------------------------------------------------------
    def _create_tiledb_ctx(self) -> tiledb.Ctx:
        return tiledb.Ctx(
            {
                "vfs.s3.region": os.environ.get("AWS_DEFAULT_REGION", "us-west-2"),
                "py.init_buffer_bytes": 4 * 1024**3,  # per-column buffer size
            }
        )

    # ----------------------------------------------------------------
    def add_h5ad(self, input_dataset_id: Optional[str], input_h5ad_path) -> int:
        """Add an H5AD to the atlas"""

        soco = self._init_soco()
        soma_name = self._get_soma_name(input_dataset_id, input_h5ad_path)
        if soma_name in soco:
            raise Exception(f"SOMA {soma_name} is already in SOMACollection {soco.uri}")

        logger.info("Loading H5AD")
        ann = anndata.read_h5ad(input_h5ad_path)

        self._clean_and_add(ann, soma_name, soco)
        return 0

    # ----------------------------------------------------------------
    def add_soma(self, input_dataset_id: Optional[str], input_soma_uri: str) -> int:
        """Add a SOMA to the atlas"""

        soco = self._init_soco()
        soma_name = self._get_soma_name(input_dataset_id, input_soma_uri)
        if soma_name in soco:
            raise Exception(f"SOMA {soma_name} is already in SOMACollection {soco.uri}")

        logger.info("Loading SOMA")
        input_soma = tiledbsc.SOMA(input_soma_uri)
        ann = tiledbsc.io.to_anndata(input_soma)

        self._clean_and_add(ann, soma_name, soco)
        return 0

    # ----------------------------------------------------------------
    def _init_soco(self) -> tiledbsc.SOMACollection:
        """
        Makes sure the destination SOMACollection exists for first write.
        """
        soco = tiledbsc.SOMACollection(self.atlas_uri, name="atlas", ctx=self.ctx)
        soco.create_unless_exists()  # Must be done first, to create the parent directory
        if not soco.exists():
            raise Exception(f"Could not create SOCO at {soco.uri}")
        return soco

    # ----------------------------------------------------------------
    def _get_soma_name(self, dataset_id: Optional[str], uri: str) -> str:
        if dataset_id is not None:
            return dataset_id
        basename = os.path.basename(uri.rstrip("/"))
        id, ext = os.path.splitext(basename)
        return id

    # ----------------------------------------------------------------
    def _clean_and_add(
        self,
        ann: anndata.AnnData,
        soma_name: str,
        soco: tiledbsc.SOMACollection,
    ):
        """
        Cleans and uniformizes the data (whether obtained from H5AD or SOMA), writes a new SOMA, adds an
        X/rankit layer, and adds the new SOMA to the SOMACollection.
        """
        logger.info("Cleaning data")
        ann = self._clean_and_uniformize(ann)

        logger.info("Creating rankit")
        X_rankit = _rankit(ann.X)

        logger.info("Saving SOMA")
        soma_uri = f"{self.atlas_uri}/{soma_name}"
        atlas_soma = tiledbsc.SOMA(uri=soma_uri, name=soma_name, ctx=self.ctx)
        tiledbsc.io.from_anndata(atlas_soma, ann)

        logger.info(f"Adding SOMA name {atlas_soma.name} at SOMA URI {atlas_soma.uri}")
        soco.add(atlas_soma)

        # Create rankit X layer and save
        logger.info("Saving rankit layer")
        if "rankit" in atlas_soma.X.keys():
            raise Exception(
                f"rankit layer already exists in the SOMA {atlas_soma.name} {atlas_soma.uri}"
            )
        atlas_soma.X.add_layer_from_matrix_and_dim_values(
            X_rankit, ann.obs.index, ann.var.index, "rankit"
        )

    # ----------------------------------------------------------------
    def _clean_and_uniformize(
        self,
        original_ann: anndata.AnnData,
    ) -> anndata.AnnData:
        """
        Given an AnnData structure (whether from H5AD or SOMA), drops all non-primary,
        non-integrable slots, layers and other unstructured information, and cleans up obs/var
        dataframes.
        """

        # Remove non-primary data.
        if self.allow_non_primary_data:
            ann = original_ann
        else:
            # This one line -- right here -- is the primary reason for using AnnData as an in-memory
            # container even for data read from a SOMA. This is MARVELOUSLY expressive. :)
            ann = original_ann[original_ann.obs.is_primary_data]
            if ann.obs.index.empty:
                raise Exception(
                    "The is_primary_data == True filter left no data remaining for this dataset."
                )

        # Clean up obs by dropping extraneous columns.
        obs = ann.obs[self.OBS_COLUMNS].copy()

        # The raw.var may or may not be present.
        original_var = ann.raw.var if ann.raw is not None else ann.var

        # .raw.var is missing mandatory columns in some datasets.
        # Work-around: populate required columns by joining with .var.
        var_columns_present = [v for v in self.VAR_COLUMNS if v in original_var]
        var_columns_not_present = [
            v for v in self.VAR_COLUMNS if v not in set(var_columns_present)
        ]
        var = original_var[var_columns_present].copy()
        if len(var_columns_not_present) > 0:
            var = var.join(ann.var[var_columns_not_present])

        rawX = ann.raw.X if ann.raw is not None else ann.X
        return anndata.AnnData(obs=obs, var=var, X=rawX)


# ================================================================
# UTILITIES

# ----------------------------------------------------------------
@nb.jit
def _quantiles(n: int) -> np.ndarray:
    """
    Computes quantiles, accelerated using @nb.jit.
    """
    return np.array([np.round((i - 0.5) / n, 5) for i in range(1, n + 1)])


# ----------------------------------------------------------------
def _rankit(Xraw: scipy.sparse.spmatrix) -> scipy.sparse.csr_matrix:
    """
    https://en.wikipedia.org/wiki/Rankit

    Caveat: equal values are ranked in undefined order.
    """
    X = Xraw.tocsr(copy=True)
    indptr = X.indptr
    for row in range(0, indptr.shape[0] - 1):
        data = X.data[indptr[row] : indptr[row + 1]]
        normal_quantiles = scipy.stats.norm.ppf(_quantiles(len(data)))
        rank = np.argsort(data)
        X.data[indptr[row] : indptr[row + 1]][rank] = normal_quantiles

    return X


# ================================================================
if __name__ == "__main__":
    sys.exit(main())
