import os

import pytest

import tiledbsoma

from ._util import ROOT_DATA_DIR


@pytest.mark.parametrize("version", ["1.7.3", "1.12.3", "1.14.5", "1.15.0", "1.15.7"])
@pytest.mark.parametrize(
    "name_and_expected_shape",
    [["pbmc3k_unprocessed", (2700, 13714)], ["pbmc3k_processed", (2638, 1838)]],
)
def test_to_anndata(version, name_and_expected_shape):
    """Checks that experiments written by older versions are still readable,
    in the particular form of doing an outgest."""

    name, expected_shape = name_and_expected_shape
    path = ROOT_DATA_DIR / "soma-experiment-versions-2025-04-04" / version / name
    uri = str(path)
    if not os.path.isdir(uri):
        raise RuntimeError(
            f"Missing '{uri}' directory. Try running `make data` "
            "from the TileDB-SOMA project root directory."
        )

    with tiledbsoma.Experiment.open(uri) as exp:
        adata = tiledbsoma.io.to_anndata(
            exp,
            measurement_name="RNA",
            X_layer_name="data",
        )
        assert adata.shape == expected_shape

        expected_nobs = expected_shape[0]
        expected_nvar = expected_shape[1]

        assert adata.obs.shape[0] == expected_nobs

        assert adata.var.shape[0] == expected_nvar

        assert adata.X.shape[0] == expected_nobs
        assert adata.X.shape[1] == expected_nvar

        if name == "pbmc3k_processed":
            for key in ["X_pca"]:
                assert adata.obsm[key].shape == (expected_nobs, 50)

            for key in ["X_tsne", "X_umap", "X_draw_graph_fr"]:
                assert adata.obsm[key].shape == (expected_nobs, 2)

            for key in ["connectivities", "distances"]:
                assert adata.obsp[key].shape == (expected_nobs, expected_nobs)

            assert adata.varm["PCs"].shape == (expected_nvar, 50)


@pytest.mark.parametrize(
    "version_and_upgraded",
    [
        ["1.7.3", False],
        ["1.12.3", False],
        ["1.14.5", False],
        ["1.15.0", True],
        ["1.15.7", True],
    ],
)
def test_get_experiment_shapes(version_and_upgraded):
    """Checks that experiments written by older versions are still readable,
    in the particular form of doing an outgest."""

    version, upgraded = version_and_upgraded
    name = "pbmc3k_processed"
    path = ROOT_DATA_DIR / "soma-experiment-versions-2025-04-04" / version / name
    uri = str(path)
    if not os.path.isdir(uri):
        raise RuntimeError(
            f"Missing '{uri}' directory. Try running `make data` "
            "from the TileDB-SOMA project root directory."
        )

    # The output includes URIs which vary randomly from one test run to another;
    # mask out the URIs.
    dict_output = tiledbsoma.io.get_experiment_shapes(uri)

    dict_output["obs"]["uri"] = "test"
    dict_output["ms"]["RNA"]["var"]["uri"] = "test"
    dict_output["ms"]["RNA"]["X"]["data"]["uri"] = "test"
    dict_output["ms"]["RNA"]["obsm"]["X_draw_graph_fr"]["uri"] = "test"
    dict_output["ms"]["RNA"]["obsm"]["X_pca"]["uri"] = "test"
    dict_output["ms"]["RNA"]["obsm"]["X_tsne"]["uri"] = "test"
    dict_output["ms"]["RNA"]["obsm"]["X_umap"]["uri"] = "test"
    dict_output["ms"]["RNA"]["obsp"]["connectivities"]["uri"] = "test"
    dict_output["ms"]["RNA"]["obsp"]["distances"]["uri"] = "test"
    dict_output["ms"]["RNA"]["varm"]["PCs"]["uri"] = "test"
    dict_output["ms"]["raw"]["var"]["uri"] = "test"
    dict_output["ms"]["raw"]["X"]["data"]["uri"] = "test"

    if not upgraded:
        expect = {
            "obs": {
                "uri": "test",
                "type": "DataFrame",
                "count": 2638,
                "non_empty_domain": ((0, 2637),),
                "domain": ((0, 2147483646),),
                "maxdomain": ((0, 2147483646),),
                "upgraded": False,
            },
            "ms": {
                "RNA": {
                    "var": {
                        "uri": "test",
                        "type": "DataFrame",
                        "count": 1838,
                        "non_empty_domain": ((0, 1837),),
                        "domain": ((0, 2147483646),),
                        "maxdomain": ((0, 2147483646),),
                        "upgraded": False,
                    },
                    "X": {
                        "data": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (0, 1837)),
                            "shape": (2147483646, 2147483646),
                            "maxshape": (2147483646, 2147483646),
                            "upgraded": False,
                        }
                    },
                    "obsm": {
                        "X_draw_graph_fr": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (0, 1)),
                            "shape": (2147483646, 2147483646),
                            "maxshape": (2147483646, 2147483646),
                            "upgraded": False,
                        },
                        "X_pca": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (0, 49)),
                            "shape": (2147483646, 2147483646),
                            "maxshape": (2147483646, 2147483646),
                            "upgraded": False,
                        },
                        "X_tsne": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (0, 1)),
                            "shape": (2147483646, 2147483646),
                            "maxshape": (2147483646, 2147483646),
                            "upgraded": False,
                        },
                        "X_umap": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (0, 1)),
                            "shape": (2147483646, 2147483646),
                            "maxshape": (2147483646, 2147483646),
                            "upgraded": False,
                        },
                    },
                    "obsp": {
                        "connectivities": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (0, 2637)),
                            "shape": (2147483646, 2147483646),
                            "maxshape": (2147483646, 2147483646),
                            "upgraded": False,
                        },
                        "distances": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (1, 2637)),
                            "shape": (2147483646, 2147483646),
                            "maxshape": (2147483646, 2147483646),
                            "upgraded": False,
                        },
                    },
                    "varm": {
                        "PCs": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 1837), (0, 49)),
                            "shape": (2147483646, 2147483646),
                            "maxshape": (2147483646, 2147483646),
                            "upgraded": False,
                        }
                    },
                },
                "raw": {
                    "var": {
                        "uri": "test",
                        "type": "DataFrame",
                        "count": 13714,
                        "non_empty_domain": ((0, 13713),),
                        "domain": ((0, 2147483646),),
                        "maxdomain": ((0, 2147483646),),
                        "upgraded": False,
                    },
                    "X": {
                        "data": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (0, 13713)),
                            "shape": (2147483646, 2147483646),
                            "maxshape": (2147483646, 2147483646),
                            "upgraded": False,
                        }
                    },
                },
            },
        }
    else:
        expect = {
            "obs": {
                "uri": "test",
                "type": "DataFrame",
                "count": 2638,
                "non_empty_domain": ((0, 2637),),
                "domain": ((0, 2637),),
                "maxdomain": ((0, 9223372036854773758),),
                "upgraded": True,
            },
            "ms": {
                "RNA": {
                    "var": {
                        "uri": "test",
                        "type": "DataFrame",
                        "count": 1838,
                        "non_empty_domain": ((0, 1837),),
                        "domain": ((0, 1837),),
                        "maxdomain": ((0, 9223372036854773968),),
                        "upgraded": True,
                    },
                    "X": {
                        "data": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (0, 1837)),
                            "shape": (2638, 1838),
                            "maxshape": (9223372036854773759, 9223372036854773759),
                            "upgraded": True,
                        }
                    },
                    "obsm": {
                        "X_draw_graph_fr": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (0, 1)),
                            "shape": (2638, 2),
                            "maxshape": (9223372036854773759, 9223372036854773759),
                            "upgraded": True,
                        },
                        "X_pca": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (0, 49)),
                            "shape": (2638, 50),
                            "maxshape": (9223372036854773759, 9223372036854773759),
                            "upgraded": True,
                        },
                        "X_tsne": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (0, 1)),
                            "shape": (2638, 2),
                            "maxshape": (9223372036854773759, 9223372036854773759),
                            "upgraded": True,
                        },
                        "X_umap": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (0, 1)),
                            "shape": (2638, 2),
                            "maxshape": (9223372036854773759, 9223372036854773759),
                            "upgraded": True,
                        },
                    },
                    "obsp": {
                        "connectivities": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (0, 2637)),
                            "shape": (2638, 2638),
                            "maxshape": (9223372036854773759, 9223372036854773759),
                            "upgraded": True,
                        },
                        "distances": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (1, 2637)),
                            "shape": (2638, 2638),
                            "maxshape": (9223372036854773759, 9223372036854773759),
                            "upgraded": True,
                        },
                    },
                    "varm": {
                        "PCs": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 1837), (0, 49)),
                            "shape": (1838, 50),
                            "maxshape": (9223372036854773759, 9223372036854773759),
                            "upgraded": True,
                        }
                    },
                },
                "raw": {
                    "var": {
                        "uri": "test",
                        "type": "DataFrame",
                        "count": 13714,
                        "non_empty_domain": ((0, 13713),),
                        "domain": ((0, 13713),),
                        "maxdomain": ((0, 9223372036854773758),),
                        "upgraded": True,
                    },
                    "X": {
                        "data": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2637), (0, 13713)),
                            "shape": (2638, 13714),
                            "maxshape": (9223372036854773759, 9223372036854773759),
                            "upgraded": True,
                        }
                    },
                },
            },
        }

    assert dict_output == expect
