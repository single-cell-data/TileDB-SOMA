import contextlib
import datetime
import json
import math
import os
import shutil
import struct
import time
from pathlib import Path
from typing import Any, List

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
import somacore
from numpy.testing import assert_array_equal
from pandas.api.types import union_categoricals
from typeguard import suppress_type_checks

import tiledbsoma as soma

from tests._util import raises_no_typeguard

from ._util import ROOT_DATA_DIR


@pytest.fixture
def arrow_schema():
    def _schema():
        return pa.schema(
            [
                pa.field("myint", pa.int64()),
                pa.field("myfloat", pa.float64()),
                pa.field("mystring", pa.string()),
                pa.field("mybool", pa.bool_()),
            ]
        )

    return _schema


def test_dataframe(tmp_path, arrow_schema):
    uri = tmp_path.as_posix()
    # Create
    asch = arrow_schema()
    with pytest.raises(ValueError):
        # requires one or more index columns
        soma.DataFrame.create(uri, schema=asch, index_column_names=[])
    with raises_no_typeguard(TypeError):
        # invalid schema type
        soma.DataFrame.create(uri, schema=asch.to_string(), index_column_names=[])
    with pytest.raises(ValueError):
        # nonexistent indexed column
        soma.DataFrame.create(uri, schema=asch, index_column_names=["bogus"])
    soma.DataFrame.create(
        uri, schema=asch, index_column_names=["myint"], domain=[[0, 99]]
    ).close()

    assert soma.DataFrame.exists(uri)
    assert not soma.Collection.exists(uri)
    assert not soma.SparseNDArray.exists(uri)

    with soma.DataFrame.open(uri) as sdf:
        assert sdf.count == 0
        assert len(sdf) == 0

        assert sorted(sdf.schema.names) == sorted(
            ["myint", "myfloat", "mystring", "soma_joinid", "mybool"]
        )
        assert sorted(sdf.keys()) == sorted(sdf.schema.names)

    pydict = {}
    pydict["soma_joinid"] = [0, 1, 2, 3, 4]
    pydict["myint"] = [10, 20, 30, 40, 50]
    pydict["myfloat"] = [4.1, 5.2, 6.3, 7.4, 8.5]
    pydict["mystring"] = ["apple", "ball", "cat", "dog", "egg"]
    pydict["mybool"] = [True, False, False, True, False]
    rb = pa.Table.from_pydict(pydict)

    with soma.DataFrame.open(uri, "w") as sdf:
        # Write
        for _ in range(3):
            sdf.tiledbsoma_resize_soma_joinid_shape(len(rb))
            sdf.write(rb)
            with raises_no_typeguard(TypeError):
                # non-arrow write
                sdf.write(rb.to_pandas)

    # Array write should fail if array opened in read mode
    with soma.DataFrame.open(uri, "r") as sdf:
        with pytest.raises(soma.SOMAError):
            sdf.write(rb)

    with soma.DataFrame.open(uri) as sdf:
        assert sdf.count == 5
        assert len(sdf) == 5

        assert sdf.tiledbsoma_has_upgraded_domain

        with pytest.raises(AttributeError):
            assert sdf.shape is None

        # soma_joinid is not a dim here
        assert sdf._maybe_soma_joinid_shape is None
        assert sdf._maybe_soma_joinid_maxshape is None

        # Read all
        table = sdf.read().concat()
        assert table.num_rows == 5
        assert table.num_columns == 5
        assert [e.as_py() for e in table["soma_joinid"]] == pydict["soma_joinid"]
        assert [e.as_py() for e in table["myint"]] == pydict["myint"]
        assert [e.as_py() for e in table["myfloat"]] == pydict["myfloat"]
        assert [e.as_py() for e in table["mystring"]] == pydict["mystring"]
        assert [e.as_py() for e in table["mybool"]] == pydict["mybool"]

        # Read ids
        table = sdf.read(coords=[[30, 10]]).concat()
        assert table.num_rows == 2
        assert table.num_columns == 5
        assert sorted([e.as_py() for e in table["soma_joinid"]]) == [0, 2]
        assert sorted([e.as_py() for e in table["myint"]]) == [10, 30]
        assert sorted([e.as_py() for e in table["myfloat"]]) == [4.1, 6.3]
        assert sorted([e.as_py() for e in table["mystring"]]) == ["apple", "cat"]
        assert [e.as_py() for e in table["mybool"]] == [True, False]

    # Open and read with bindings
    with contextlib.closing(
        soma.pytiledbsoma.SOMADataFrame.open(
            uri, soma.pytiledbsoma.OpenMode.read, soma.pytiledbsoma.SOMAContext()
        )
    ) as sdf:
        mq = soma.pytiledbsoma.ManagedQuery(sdf, sdf.context())
        table = mq.next()
        assert table.num_rows == 5
        assert table.num_columns == 5
        assert [e.as_py() for e in table["soma_joinid"]] == pydict["soma_joinid"]
        assert [e.as_py() for e in table["myint"]] == pydict["myint"]
        assert [e.as_py() for e in table["myfloat"]] == pydict["myfloat"]
        assert [e.as_py() for e in table["mystring"]] == pydict["mystring"]
        assert [e.as_py() for e in table["mybool"]] == pydict["mybool"]

    with soma.DataFrame.open(uri) as A:
        cfg = A.schema_config_options()
        assert not cfg.allows_duplicates
        assert json.loads(cfg.dims)["myint"]["filters"] == [
            {"COMPRESSION_LEVEL": 3, "name": "ZSTD"}
        ]
        assert json.loads(cfg.attrs)["myfloat"]["filters"] == [
            {"COMPRESSION_LEVEL": -1, "name": "ZSTD"}
        ]

    with soma.DataFrame.open(uri) as sdf:
        assert sdf.count == 5
        assert len(sdf) == 5

    # Ensure read mode uses clib object
    with soma.DataFrame.open(tmp_path.as_posix(), "r") as A:
        assert isinstance(A._handle._handle, soma.pytiledbsoma.SOMADataFrame)

    # Ensure write mode uses clib object
    with soma.DataFrame.open(tmp_path.as_posix(), "w") as A:
        assert isinstance(A._handle._handle, soma.pytiledbsoma.SOMADataFrame)

    # Ensure it cannot be opened by another type
    with pytest.raises(soma.SOMAError):
        soma.SparseNDArray.open(tmp_path.as_posix())

    with pytest.raises(soma.SOMAError):
        soma.DenseNDArray.open(tmp_path.as_posix())

    with pytest.raises(soma.SOMAError):
        soma.PointCloudDataFrame.open(tmp_path.as_posix())

    with pytest.raises(soma.SOMAError):
        soma.Collection.open(tmp_path.as_posix())

    with pytest.raises(soma.SOMAError):
        soma.Experiment.open(tmp_path.as_posix())

    with pytest.raises(soma.SOMAError):
        soma.Measurement.open(tmp_path.as_posix())

    with pytest.raises(soma.SOMAError):
        soma.Scene.open(tmp_path.as_posix())

    with pytest.raises(soma.SOMAError):
        soma.MultiscaleImage.open(tmp_path.as_posix())


def test_dataframe_reopen(tmp_path, arrow_schema):
    soma.DataFrame.create(
        tmp_path.as_posix(), schema=arrow_schema(), tiledb_timestamp=1
    )

    with soma.DataFrame.open(tmp_path.as_posix(), "r", tiledb_timestamp=1) as sdf1:
        with raises_no_typeguard(ValueError):
            sdf1.reopen("invalid")

        with sdf1.reopen("w", tiledb_timestamp=2) as sdf2:
            with sdf2.reopen("r", tiledb_timestamp=3) as sdf3:
                assert sdf1.mode == "r"
                assert sdf2.mode == "w"
                assert sdf3.mode == "r"
                assert sdf1.tiledb_timestamp_ms == 1
                assert sdf2.tiledb_timestamp_ms == 2
                assert sdf3.tiledb_timestamp_ms == 3

    ts1 = datetime.datetime(2023, 1, 1, 1, 0, tzinfo=datetime.timezone.utc)
    ts2 = datetime.datetime(2024, 1, 1, 1, 0, tzinfo=datetime.timezone.utc)
    with soma.DataFrame.open(tmp_path.as_posix(), "r", tiledb_timestamp=ts1) as sdf1:
        with sdf1.reopen("r", tiledb_timestamp=ts2) as sdf2:
            assert sdf1.mode == "r"
            assert sdf2.mode == "r"
            assert sdf1.tiledb_timestamp == ts1
            assert sdf2.tiledb_timestamp == ts2

    with soma.DataFrame.open(tmp_path.as_posix(), "w") as sdf1:
        with sdf1.reopen("w", tiledb_timestamp=None) as sdf2:
            with sdf1.reopen("w") as sdf3:
                assert sdf1.mode == "w"
                assert sdf2.mode == "w"
                assert sdf3.mode == "w"
                now = datetime.datetime.now(datetime.timezone.utc)
                assert sdf1.tiledb_timestamp <= now
                assert sdf2.tiledb_timestamp <= now
                assert sdf3.tiledb_timestamp <= now


def test_dataframe_with_float_dim(tmp_path, arrow_schema):
    sdf = soma.DataFrame.create(
        tmp_path.as_posix(), schema=arrow_schema(), index_column_names=("myfloat",)
    )
    assert sdf.index_column_names == ("myfloat",)


def test_dataframe_with_enumeration(tmp_path):
    schema = pa.schema(
        [
            pa.field("myint", pa.dictionary(pa.int64(), pa.large_string())),
            pa.field("myfloat", pa.dictionary(pa.int64(), pa.large_string())),
        ]
    )
    enums = {"enmr1": ("a", "bb", "ccc"), "enmr2": ("cat", "dog")}
    with soma.DataFrame.create(
        tmp_path.as_posix(), schema=schema, domain=[[0, 5]]
    ) as sdf:
        data = {}
        data["soma_joinid"] = [0, 1, 2, 3, 4]
        data["myint"] = ["a", "bb", "ccc", "bb", "a"]
        data["myfloat"] = ["cat", "dog", "cat", "cat", "cat"]
        with pytest.raises(soma.SOMAError):
            sdf.write(pa.Table.from_pydict(data))

        data["myint"] = pd.Categorical(["a", "bb", "ccc", "bb", "a"])
        data["myfloat"] = pd.Categorical(["cat", "dog", "cat", "cat", "cat"])
        sdf.write(pa.Table.from_pydict(data))

    with soma.DataFrame.open(tmp_path.as_posix()) as sdf:
        df = sdf.read().concat()
        assert_array_equal(df["myint"].chunk(0).dictionary, enums["enmr1"])
        assert_array_equal(df["myfloat"].chunk(0).dictionary, enums["enmr2"])


# The functionality being tested here doesn't depend on whether the enumerations
# are ordered or not.  (There are no conditional on ordered within the test
# function.) We vary it anyway.
#
# Users should be able to access the schema even when the dataframe is opened in
# write mode.
@pytest.mark.parametrize("ordered", [True, False])
@pytest.mark.parametrize("mode", ["r", "w"])
def test_get_enumeration_values(tmp_path, ordered, mode):
    uri = tmp_path.as_posix()

    schema = pa.schema(
        [
            pa.field("not_an_enum", pa.large_string()),
            pa.field("string_enum", pa.dictionary(pa.int32(), pa.large_string())),
            pa.field("int64_enum", pa.dictionary(pa.int32(), pa.int64())),
            pa.field("float64_enum", pa.dictionary(pa.int16(), pa.float64())),
            pa.field("bool_enum", pa.dictionary(pa.int8(), pa.bool_())),
        ]
    )

    domain = [[0, 7]]

    # Create the dataframe with no levels for any enumerated column
    with soma.DataFrame.create(uri, schema=schema, domain=domain) as sdf:
        pass

    with soma.DataFrame.open(uri, mode) as sdf:
        with pytest.raises(KeyError):
            sdf.get_enumeration_values(["nonesuch"])
        with pytest.raises(KeyError):
            sdf.get_enumeration_values(["not_an_enum"])
        with pytest.raises(KeyError):
            sdf.get_enumeration_values(["string_enum", "not_an_enum"])

        actual = sdf.get_enumeration_values(
            ["string_enum", "int64_enum", "float64_enum", "bool_enum"]
        )
        expect = {
            "string_enum": pa.array([], type=pa.large_string()),
            "int64_enum": pa.array([], type=pa.int64()),
            "float64_enum": pa.array([], type=pa.float64()),
            "bool_enum": pa.array([], type=pa.bool_()),
        }
        assert actual == expect

    # Check with the dataframe closed
    with pytest.raises(soma.SOMAError):
        sdf.get_enumeration_values(["bool_enum"])

    # Write once
    pd_data = {
        "soma_joinid": [0, 1, 2, 3, 4],
        "not_an_enum": ["a", "nn", "zzz", "nn", "a"],
        "string_enum": pd.Categorical(["a", "nn", "zzz", "nn", "a"], ordered=ordered),
        "int64_enum": pd.Categorical(
            [111111111, 99999, 3333333, 111111111, 99999], ordered=ordered
        ),
        # Note: some older versions (I can vouch for pandas 1.5.3 and numpy 1.25.0) do something
        # very sad here:
        #
        # >>> pd.Categorical(np.array([1.5, 0.5, 99.0, 1.5, 99.0], dtype=np.float32))
        # [1.5, 0.5, 99.0, 1.5, 99.0]
        # Categories (3, float64): [0.5, 1.5, 99.0]
        #
        # >>> pd.Categorical(np.array([1.5, 0.5, 99.0, 1.5, 99.0], dtype=np.float64))
        # [1.5, 0.5, 99.0, 1.5, 99.0]
        # Categories (3, float64): [0.5, 1.5, 99.0]
        #
        # i.e. you ask for float32 or float64, you get float64 either way.
        # So we _cannot_ construct a categorical of type float32 _even when we explicitly ask
        # for it. And for as long as this suite must pass on these older versions
        # of pandas/numpy, we must not insist that it do something which it is
        # demonstrably incapable of doing.
        "float64_enum": pd.Categorical(
            np.array([1.5, 0.5, 99.0, 1.5, 99.0], dtype=np.float64), ordered=ordered
        ),
        "bool_enum": pd.Categorical(
            [True, True, True, True, True],
            ordered=ordered,
        ),
    }
    arrow_data = pa.Table.from_pydict(pd_data)

    with soma.DataFrame.open(uri, "w") as sdf:
        sdf.write(arrow_data)
    t2 = int(time.time() * 1000)

    with soma.DataFrame.open(uri, mode) as sdf:
        with pytest.raises(KeyError):
            sdf.get_enumeration_values(["nonesuch"])
        with pytest.raises(KeyError):
            sdf.get_enumeration_values(["not_an_enum"])
        with pytest.raises(KeyError):
            sdf.get_enumeration_values(["string_enum", "not_an_enum"])

        actual = sdf.get_enumeration_values(
            ["string_enum", "int64_enum", "float64_enum", "bool_enum"]
        )
        expect = {
            "string_enum": pa.array(["a", "nn", "zzz"], type=pa.large_string()),
            "int64_enum": pa.array([111111111, 3333333, 99999], type=pa.int64()),
            "float64_enum": pa.array([1.5, 0.5, 99.0], type=pa.float64()),
            "bool_enum": pa.array([False, True, False], type=pa.bool_()),
        }

    # Write again
    pd_data = {
        "soma_joinid": [5, 6, 7],
        "not_an_enum": ["dddd", "nn", "zzz"],
        "string_enum": pd.Categorical(["dddd", "nn", "zzz"], ordered=ordered),
        "int64_enum": pd.Categorical([555555555, 111111111, 99999], ordered=ordered),
        "float64_enum": pd.Categorical(
            np.array([44.25, 0.5, 99.0], dtype=np.float64), ordered=ordered
        ),
        "bool_enum": pd.Categorical([True, False, True], ordered=ordered),
    }
    arrow_data = pa.Table.from_pydict(pd_data)

    with soma.DataFrame.open(uri, "w") as sdf:
        sdf.write(arrow_data)

    with soma.DataFrame.open(uri, mode) as sdf:
        with pytest.raises(KeyError):
            sdf.get_enumeration_values(["nonesuch"])
        with pytest.raises(KeyError):
            sdf.get_enumeration_values(["not_an_enum"])
        with pytest.raises(KeyError):
            sdf.get_enumeration_values(["string_enum", "not_an_enum"])

        actual = sdf.get_enumeration_values(
            ["string_enum", "int64_enum", "float64_enum", "bool_enum"]
        )
        expect = {
            "string_enum": pa.array(["a", "nn", "zzz", "dddd"], type=pa.large_string()),
            "int64_enum": pa.array(
                [99999, 3333333, 111111111, 555555555], type=pa.int64()
            ),
            "float64_enum": pa.array([0.5, 1.5, 99.0, 44.25], type=pa.float64()),
            "bool_enum": pa.array([True, False], type=pa.bool_()),
        }

        assert actual == expect

    # Check that we can read from before the second write
    with soma.DataFrame.open(uri, mode, tiledb_timestamp=t2) as sdf:
        with pytest.raises(KeyError):
            sdf.get_enumeration_values(["nonesuch"])
        with pytest.raises(KeyError):
            sdf.get_enumeration_values(["not_an_enum"])
        with pytest.raises(KeyError):
            sdf.get_enumeration_values(["string_enum", "not_an_enum"])

        actual = sdf.get_enumeration_values(
            ["string_enum", "int64_enum", "float64_enum", "bool_enum"]
        )
        expect = {
            "string_enum": pa.array(["a", "nn", "zzz"], type=pa.large_string()),
            "int64_enum": pa.array([111111111, 3333333, 99999], type=pa.int64()),
            "float64_enum": pa.array([1.5, 0.5, 99.0], type=pa.float64()),
            "bool_enum": pa.array([False, True, False], type=pa.bool_()),
        }


@pytest.mark.parametrize("version", ["1.7.3", "1.12.3", "1.14.5", "1.15.0", "1.15.7"])
@pytest.mark.parametrize("name", ["pbmc3k_unprocessed", "pbmc3k_processed"])
def test_get_enumeration_values_historical(version, name):
    """Checks that experiments written by older versions are still readable,
    in the particular form of doing an outgest."""

    path = ROOT_DATA_DIR / "soma-experiment-versions-2025-04-04" / version / name
    uri = str(path)
    if not os.path.isdir(uri):
        raise RuntimeError(
            f"Missing '{uri}' directory. Try running `make data` "
            "from the TileDB-SOMA project root directory."
        )

    with soma.Experiment.open(uri) as exp:
        if name == "pbmc3k_unprocessed":
            values = exp.obs.get_enumeration_values(
                ["orig.ident", "seurat_annotations"]
            )

            expect = ["pbmc3k"]
            assert values["orig.ident"].to_pylist() == expect

            expect = [
                "B",
                "CD8 T",
                "CD14+ Mono",
                "DC",
                "FCGR3A+ Mono",
                "Memory CD4 T",
                "NA",
                "NK",
                "Naive CD4 T",
                "Platelet",
            ]
            assert values["seurat_annotations"].to_pylist() == expect

        elif name == "pbmc3k_processed":
            values = exp.obs.get_enumeration_values(["louvain"])
            expect = [
                "CD4 T cells",
                "CD14+ Monocytes",
                "B cells",
                "CD8 T cells",
                "NK cells",
                "FCGR3A+ Monocytes",
                "Dendritic cells",
                "Megakaryocytes",
            ]
            assert values["louvain"].to_pylist() == expect


@pytest.mark.parametrize(
    "version",
    ["1.7.3", "1.12.3", "1.14.5", "1.15.0", "1.15.7"],
)
def test_extend_enumeration_values_historical(tmp_path, version):

    original_data_uri = str(
        ROOT_DATA_DIR
        / "soma-experiment-versions-2025-04-04"
        / version
        / "pbmc3k_processed"
    )

    if not os.path.isdir(original_data_uri):
        raise RuntimeError(
            f"Missing '{original_data_uri}' directory. Try running `make data` "
            "from the TileDB-SOMA project root directory."
        )

    # Make a copy of the Experiment as to not write over the data in ROOT_DATA_DIR
    uri = (tmp_path / version).as_posix()
    shutil.copytree(original_data_uri, uri)

    with soma.Experiment.open(uri, "w") as exp:
        with pytest.raises(soma.SOMAError):
            exp.obs.extend_enumeration_values(
                {
                    "louvain": pa.array(
                        [
                            "CD4 T cells",
                            "CD14+ Monocytes",
                            "B cells",
                        ]
                    )
                }
            )

            exp.obs.extend_enumeration_values(
                {
                    "louvain": pa.array(
                        [
                            "CD4 T cells",
                            "CD14+ Monocytes",
                            "B cells",
                        ]
                    )
                },
                deduplicate=True,
            )
    with soma.Experiment.open(uri) as exp:
        values = exp.obs.get_enumeration_values(["louvain"])
        expect = [
            "CD4 T cells",
            "CD14+ Monocytes",
            "B cells",
            "CD8 T cells",
            "NK cells",
            "FCGR3A+ Monocytes",
            "Dendritic cells",
            "Megakaryocytes",
        ]
        assert values["louvain"].to_pylist() == expect

    with soma.Experiment.open(uri, "w") as exp:
        for deduplicate in [True, False]:
            with pytest.raises(ValueError):
                exp.obs.extend_enumeration_values(
                    {"louvain": pa.array(["even more cells", "even more cells"])},
                    deduplicate=deduplicate,
                )

    with soma.Experiment.open(uri, "w") as exp:
        exp.obs.extend_enumeration_values(
            {"louvain": pa.array(["even more cells", "so many more cells"])},
        )

    with soma.Experiment.open(uri) as exp:
        values = exp.obs.get_enumeration_values(["louvain"])
        expect = [
            "CD4 T cells",
            "CD14+ Monocytes",
            "B cells",
            "CD8 T cells",
            "NK cells",
            "FCGR3A+ Monocytes",
            "Dendritic cells",
            "Megakaryocytes",
            "even more cells",
            "so many more cells",
        ]
        assert values["louvain"].to_pylist() == expect


@pytest.mark.parametrize(
    "data_and_expected_levels",
    [
        [[True], [True]],
        [[False], [False]],
        [[True, True], [True]],
        [[True, False], [False, True]],
        [[False, True], [False, True]],
        [[False, False], [False]],
        [[False, True, False, True, False], [False, True]],
    ],
)
def test_bool_enums(tmp_path, data_and_expected_levels):
    uri = tmp_path.as_posix()

    data, expected_levels = data_and_expected_levels
    n = len(data)
    domain = [[0, n - 1]]

    schema = pa.schema(
        [
            pa.field("bool_enum", pa.dictionary(pa.int8(), pa.bool_())),
        ]
    )

    # Create the dataframe with no expected_levels for any enumerated column
    with soma.DataFrame.create(uri, schema=schema, domain=domain) as sdf:
        pass

    pd_data = {
        "soma_joinid": list(range(n)),
        "bool_enum": pd.Categorical(data),
    }
    arrow_data = pa.Table.from_pydict(pd_data)

    with soma.DataFrame.open(uri, "w") as sdf:
        sdf.write(arrow_data)

    with soma.DataFrame.open(uri) as sdf:
        actual = sdf.get_enumeration_values(["bool_enum"])
        expect = {
            "bool_enum": pa.array(expected_levels, type=pa.bool_()),
        }
        assert actual == expect


@pytest.mark.parametrize("extend_not_write", [False, True])
@pytest.mark.parametrize("ordered", [True, False])
def test_extend_enumeration_values(tmp_path, extend_not_write, ordered):
    """Compares the older create+write path against create+extend path"""
    uri = tmp_path.as_posix()
    domain = [[0, 2]]

    schema = pa.schema(
        {
            "soma_joinid": pa.int64(),
            "not_an_enum": pa.large_string(),
            "int64_enum1": pa.dictionary(pa.int8(), pa.int64(), ordered=ordered),
            "int64_enum2": pa.dictionary(pa.int8(), pa.int64(), ordered=ordered),
            "int64_enum3": pa.dictionary(pa.int8(), pa.int64(), ordered=ordered),
            "float64_enum": pa.dictionary(pa.int8(), pa.float64(), ordered=ordered),
            "string_enum1": pa.dictionary(
                pa.int32(), pa.large_string(), ordered=ordered
            ),
            "string_enum2": pa.dictionary(
                pa.int32(), pa.large_string(), ordered=ordered
            ),
            "bool_enum1": pa.dictionary(pa.int32(), pa.bool_(), ordered=ordered),
            "bool_enum2": pa.dictionary(pa.int32(), pa.bool_(), ordered=ordered),
        }
    )
    enum_column_names = [
        name for name in schema.names if pa.types.is_dictionary(schema.field(name).type)
    ]

    pandas_data = {
        "soma_joinid": [0, 1, 2],
        "not_an_enum": ["quick", "brown", "fox"],
        "int64_enum1": pd.Categorical([55555, 55555, 7777777], ordered=ordered),
        "int64_enum2": pd.Categorical(
            np.array([55555, 55555, 7777777]), ordered=ordered
        ),
        "int64_enum3": pd.Categorical(
            np.array([7777777, 55555, 55555], dtype=np.int64),
            ordered=ordered,
        ),
        "float64_enum": pd.Categorical(
            np.array([2.5, 8.875, 2.5], dtype=np.float64), ordered=ordered
        ),
        "string_enum1": pd.Categorical(["hello", "hello", "goodbye"], ordered=ordered),
        "string_enum2": pd.Categorical(
            ["goodbye", "goodbye", "hello"], ordered=ordered
        ),
        "bool_enum1": pd.Categorical([True, True, False], ordered=ordered),
        "bool_enum2": pd.Categorical([False, False, True], ordered=ordered),
    }
    arrow_data = pa.Table.from_pydict(pandas_data)

    expected_levels = {
        "int64_enum1": [55555, 7777777],
        "int64_enum2": [55555, 7777777],
        "int64_enum3": [55555, 7777777],
        "float64_enum": [2.5, 8.875],
        "string_enum1": ["goodbye", "hello"],
        "string_enum2": ["goodbye", "hello"],
        "bool_enum1": [False, True],
        "bool_enum2": [False, True],
    }
    arrow_data = pa.Table.from_pydict(pandas_data)

    # Note: doing the .categories bit is implicitly deduplicating the data
    # before it's sent to tiledbsoma. Other tests, not this one, stress the
    # deduplicate-or-not logic.
    values = {
        name: pa.array(pandas_data[name].categories) for name in enum_column_names
    }

    # Do the extend without write, or write with implicit extend.
    # Note: this function does not test multiple extensions, but another one does.
    with soma.DataFrame.create(
        uri,
        schema=schema,
        domain=domain,
        index_column_names=["soma_joinid"],
    ) as sdf:
        if extend_not_write:
            sdf.extend_enumeration_values(values)
        else:
            sdf.write(arrow_data)

    # The dataframe must be open for write
    with soma.DataFrame.open(uri) as sdf:
        with pytest.raises(soma.SOMAError):
            sdf.extend_enumeration_values({})
    with pytest.raises(soma.SOMAError):
        assert sdf.closed
        sdf.extend_enumeration_values({})

    with soma.DataFrame.open(uri, "w") as sdf:
        # The column must exist
        with pytest.raises(KeyError):
            sdf.extend_enumeration_values({"nonesuch": pa.array(["abc"])})

        # The column must be of enumerated type
        with pytest.raises(KeyError):
            sdf.extend_enumeration_values({"not_an_enum": pa.array(["abc"])})

        # The values provided must be non-dictionary
        xvalues = {name: pa.array(pandas_data[name]) for name in enum_column_names}
        with pytest.raises(
            ValueError,
            match=r"is of dictionary type: pass its dictionary array instead",
        ):
            sdf.extend_enumeration_values(xvalues)

        # The values provided must be Arrow arrays. Our unit tests run with typeguard,
        # but our end users nominally do not -- so we have to ask typeguard to take
        # a breather here so we can test the UX our users will have.
        with suppress_type_checks():
            with pytest.raises(ValueError):
                sdf.extend_enumeration_values({"string_enum1": ["plain", "strings"]})

        # The values provided must all be non-null
        for nvalues in [
            {"string_enum1": pa.array(["greetings", None])},
            {"int64_enum1": pa.array([54321, None])},
        ]:
            with pytest.raises(
                soma.SOMAError,
                match=r"null values are not supported",
            ):
                sdf.extend_enumeration_values(nvalues)

        # Types must match
        type_mismatch_pandas_data = {
            "int64_enum3": pd.Categorical(np.array([4444, 333, 333], dtype=np.float32)),
        }
        type_mismatch_values = {
            name: pa.array(type_mismatch_pandas_data[name].categories)
            for name in type_mismatch_pandas_data.keys()
        }
        with pytest.raises(soma.SOMAError):
            sdf.extend_enumeration_values(type_mismatch_values)

    # Verify
    with soma.DataFrame.open(uri) as sdf:
        table = sdf.read(column_names=enum_column_names).concat()
        for column_name in values.keys():
            field = sdf.schema.field(column_name)
            assert pa.types.is_dictionary(field.type)

            chunk = table.column(column_name).chunk(0)
            assert sorted(chunk.dictionary.to_pylist()) == expected_levels[column_name]
            if extend_not_write:
                assert len(chunk) == 0
            else:
                assert len(chunk) == 3


@pytest.mark.parametrize("deduplicate", [False, True])
@pytest.mark.parametrize("ordered", [False, True])
@pytest.mark.parametrize("extend_not_first_write", [False, True])
def test_extend_enumeration_values_deduplication(
    tmp_path, deduplicate, ordered, extend_not_first_write
):
    uri = tmp_path.as_posix()
    domain = [[0, 9]]

    schema = pa.schema(
        {
            "soma_joinid": pa.int64(),
            "not_an_enum": pa.large_string(),
            "string_enum": pa.dictionary(
                pa.int32(), pa.large_string(), ordered=ordered
            ),
            "int64_enum": pa.dictionary(pa.int8(), pa.int64(), ordered=ordered),
            "float32_enum": pa.dictionary(pa.int8(), pa.float32(), ordered=ordered),
            "bool_enum": pa.dictionary(pa.int32(), pa.bool_(), ordered=ordered),
        }
    )

    with soma.DataFrame.create(
        uri,
        schema=schema,
        domain=domain,
    ) as sdf:

        # Dupes in the inputs are disallowed regardless of the deduplicate flag.
        values_list = [
            {"string_enum": pa.array(["hello", "hello", "goodbye"], type=pa.string())},
            {"int64_enum": pa.array([55555, 55555, 22, 333, 7777777], type=pa.int64())},
            {"float32_enum": pa.array([2.25, 3.75, 2.25], type=pa.float32())},
            {"bool_enum": pa.array([True, True], type=pa.bool_())},
        ]
        for values in values_list:
            # Run separate asserts for each column name. We want to make sure _each_ of them throws.
            with pytest.raises(ValueError):
                sdf.extend_enumeration_values(values, deduplicate=deduplicate)

        values = {"float32_enum": pa.array([2.25, 3.75, 2.25], type=pa.float64())}
        with pytest.raises(soma.SOMAError):
            sdf.extend_enumeration_values(values)

        # Success
        values = {
            "string_enum": pa.array(["hello", "goodbye"], type=pa.large_string()),
            "int64_enum": pa.array([55555, 22, 333, 7777777], type=pa.int64()),
            "float32_enum": pa.array([2.25, 3.75], type=pa.float32()),
            "bool_enum": pa.array([True], type=pa.bool_()),
        }

        if extend_not_first_write:
            sdf.extend_enumeration_values(values)
        else:
            data = pa.Table.from_pydict(
                {
                    "soma_joinid": pa.array([0, 1, 2, 3], type=pa.int64()),
                    "not_an_enum": pa.array(
                        ["the", "quick", "brown", "fox"], type=pa.large_string()
                    ),
                    "string_enum": pa.DictionaryArray.from_arrays(
                        [0, 1, 0, 1], values["string_enum"]
                    ),
                    "int64_enum": pa.DictionaryArray.from_arrays(
                        [0, 2, 3, 1], values["int64_enum"]
                    ),
                    "float32_enum": pa.DictionaryArray.from_arrays(
                        [0, 0, 0, 1], values["float32_enum"]
                    ),
                    "bool_enum": pa.DictionaryArray.from_arrays(
                        [0, 0, 0, 0], values["bool_enum"]
                    ),
                }
            )
            sdf.write(data)

    with soma.DataFrame.open(uri, "w") as sdf:

        # Dupes between the inputs and the existing schema are allowed only
        # if the deduplicate flag is set.
        values_list = [
            {"string_enum": pa.array(["farewell", "goodbye"], type=pa.string())},
            {"int64_enum": pa.array([4444, 7777777], type=pa.int64())},
            {"float32_enum": pa.array([9.00, 2.25], type=pa.float32())},
            {"bool_enum": pa.array([False, True], type=pa.bool_())},
        ]
        for values in values_list:
            # Run separate asserts for each column name. We want to make sure _each_ of them throws.
            if deduplicate:
                sdf.extend_enumeration_values(values, deduplicate=deduplicate)
            else:
                with pytest.raises(soma.SOMAError):
                    sdf.extend_enumeration_values(values, deduplicate=deduplicate)

    if deduplicate:
        expect = {
            "string_enum": pa.array(
                ["hello", "goodbye", "farewell"], type=pa.large_string()
            ),
            "int64_enum": pa.array([55555, 22, 333, 7777777, 4444], type=pa.int64()),
            "float32_enum": pa.array([2.25, 3.75, 9.0], type=pa.float32()),
            "bool_enum": pa.array([True, False], type=pa.bool_()),
        }
    else:
        expect = {
            "string_enum": pa.array(["hello", "goodbye"], type=pa.large_string()),
            "int64_enum": pa.array([55555, 22, 333, 7777777], type=pa.int64()),
            "float32_enum": pa.array([2.25, 3.75], type=pa.float32()),
            "bool_enum": pa.array([True], type=pa.bool_()),
        }

    with soma.DataFrame.open(uri) as sdf:
        actual = sdf.get_enumeration_values(
            ["string_enum", "int64_enum", "float32_enum", "bool_enum"]
        )
        assert actual == expect

    # Different representations of single-precision NaNs.  There are
    # single-precision and double-precision NaNs, of course, but np.nan is
    # single-precision.
    quiet_nan = struct.unpack(">f", b"\x7f\xc0\x00\x00")[0]
    negative_nan = struct.unpack(">f", b"\xff\xc0\x00\x00")[0]

    with soma.DataFrame.open(uri, "w") as sdf:
        values = {
            "float32_enum": pa.array([quiet_nan, quiet_nan], type=pa.float32()),
        }
        with pytest.raises(soma.SOMAError):
            sdf.extend_enumeration_values(values)

        # Core treats these as distinct so we do too
        values = {
            "float32_enum": pa.array([quiet_nan, negative_nan], type=pa.float32()),
        }
        # Implicit check for no throw
        sdf.extend_enumeration_values(values)

    with soma.DataFrame.open(uri) as sdf:
        readback = sdf.get_enumeration_values(["float32_enum"])["float32_enum"]
        if deduplicate:
            assert len(readback) == 5
            assert readback[:3].to_pylist() == [2.25, 3.75, 9.0]
            nans = readback[3:]
        else:
            assert len(readback) == 4
            assert readback[:2].to_pylist() == [2.25, 3.75]
            nans = readback[2:]
        assert len(nans) == 2
        assert struct.pack(">f", nans[0].as_py()) == struct.pack(">f", quiet_nan)
        assert struct.pack(">f", nans[1].as_py()) == struct.pack(">f", negative_nan)


@pytest.mark.parametrize("ordered", [False, True])
def test_extend_enumeration_values_offsets(tmp_path, ordered):
    uri = tmp_path.as_posix()
    domain = [[0, 0]]

    schema = pa.schema(
        {
            "soma_joinid": pa.int64(),
            "not_an_enum": pa.large_string(),
            "string_enum": pa.dictionary(
                pa.int32(), pa.large_string(), ordered=ordered
            ),
            "int64_enum": pa.dictionary(pa.int8(), pa.int64(), ordered=ordered),
            "float32_enum": pa.dictionary(pa.int8(), pa.float32(), ordered=ordered),
            "bool_enum": pa.dictionary(pa.int32(), pa.bool_(), ordered=ordered),
        }
    )

    values = {
        "string_enum": pa.array(
            ["auf", "wieder", "sehen", "Freunde"], type=pa.large_string()
        )[1:3],
        "int64_enum": pa.array([1234, 2345, 3456, 4567], type=pa.int64())[1:3],
        "float32_enum": pa.array([25.0, 26.0, 27.0, 28.0], type=pa.float32())[1:3],
        "bool_enum": pa.array([False, True, True, False], type=pa.bool_())[1:2],
    }

    with soma.DataFrame.create(
        uri,
        schema=schema,
        domain=domain,
    ) as sdf:
        sdf.extend_enumeration_values(values)

    expect = {
        "string_enum": pa.array(["wieder", "sehen"], type=pa.large_string()),
        "int64_enum": pa.array([2345, 3456], type=pa.int64()),
        "float32_enum": pa.array([26.0, 27.0], type=pa.float32()),
        "bool_enum": pa.array([True], type=pa.bool_()),
    }
    with soma.DataFrame.open(uri) as sdf:
        actual = sdf.get_enumeration_values(
            ["string_enum", "int64_enum", "float32_enum", "bool_enum"]
        )
        assert actual == expect


@pytest.mark.parametrize(
    "config",
    [
        # String-valued enums
        {
            "arrow_type": pa.string(),
            "inidx1": [0, 1, 2, 3],
            "inval1": ["red", "yellow", "green", "blue"],
            "inidx2": [0, 1, 2, 3],
            "inval2": ["red", "yellow", "green", "blue"],
            "expidx": [0, 1, 2, 3, 0, 1, 2, 3],
            "expval": ["red", "yellow", "green", "blue"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [0, 1, 2, 3],
            "inval1": ["red", "yellow", "green", "blue"],
            "inidx2": [],
            "inval2": [],
            "expidx": [0, 1, 2, 3],
            "expval": ["red", "yellow", "green", "blue"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [],
            "inval1": [],
            "inidx2": [0, 1, 2, 3],
            "inval2": ["red", "yellow", "green", "blue"],
            "expidx": [0, 1, 2, 3],
            "expval": ["red", "yellow", "green", "blue"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [0, 1, 2, 3],
            "inval1": ["red", "yellow", "green", "blue"],
            "inidx2": [0, 1, 2, 3],
            "inval2": ["orange", "purple", "indigo", "violet"],
            "expidx": [0, 1, 2, 3, 4, 5, 6, 7],
            "expval": [
                "red",
                "yellow",
                "green",
                "blue",
                "orange",
                "purple",
                "indigo",
                "violet",
            ],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [0, 1, 2, 3],
            "inval1": ["red", "yellow", "green", "blue"],
            "inidx2": [3, 2, 1, 0],
            "inval2": ["red", "yellow", "green", "blue"],
            "expidx": [0, 1, 2, 3, 3, 2, 1, 0],
            "expval": ["red", "yellow", "green", "blue"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [0, 1, 2, 3],
            "inval1": ["red", "yellow", "green", "blue"],
            "inidx2": [0, 1, 2, 3],
            "inval2": ["orange", "grey", "green", "blue"],
            "expidx": [0, 1, 2, 3, 4, 5, 2, 3],
            "expval": ["red", "yellow", "green", "blue", "orange", "grey"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [0, 1, 2, 3],
            "inval1": ["red", "yellow", "green", "blue"],
            "inidx2": [0, None, 2, 3],
            "inval2": ["orange", "grey", "green", "blue"],
            "expidx": [0, 1, 2, 3, 4, None, 2, 3],
            "expval": ["red", "yellow", "green", "blue", "orange"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [0, None, 2, 3],
            "inval1": ["red", "yellow", "green", "blue"],
            "inidx2": [0, 1, 2, 3],
            "inval2": ["orange", "grey", "green", "blue"],
            "expidx": [0, None, 1, 2, 3, 4, 1, 2],
            "expval": ["red", "green", "blue", "orange", "grey"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [0, 1, 2, 3],
            "inval1": ["red", "yellow", "green", "blue"],
            "inidx2": [0, None, 0, 0],
            "inval2": ["orange"],
            "expidx": [0, 1, 2, 3, 4, None, 4, 4],
            "expval": ["red", "yellow", "green", "blue", "orange"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [0, 1, 2, 3],
            "inval1": ["red", "yellow", "green", "blue"],
            "inidx2": [None, None, None, None],
            "inval2": ["orange"],
            "expidx": [0, 1, 2, 3, None, None, None, None],
            "expval": ["red", "yellow", "green", "blue"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [0, 1, 2, 3],
            "inval1": ["red", "yellow", "green", "blue"],
            "inidx2": [0, 1, 2, 2],
            "inval2": ["yellow", "green", "blue"],
            "expidx": [0, 1, 2, 3, 1, 2, 3, 3],
            "expval": ["red", "yellow", "green", "blue"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [None, None, None, None],
            "inval1": ["red", "yellow", "green", "blue"],
            "inidx2": [None, None, None, None],
            "inval2": ["yellow", "green", "blue"],
            "expidx": [None, None, None, None, None, None, None, None],
            "expval": [],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [None, None, None, None],
            "inval1": ["orange"],
            "inidx2": [None, None, None, None],
            "inval2": ["orange"],
            "expidx": [None, None, None, None, None, None, None, None],
            "expval": [],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [0, 1, 2, 3],
            "inval1": ["red", "yellow", "green", "blue"],
            "inidx2": [None, None, None, None],
            "inval2": [],
            "expidx": [0, 1, 2, 3, None, None, None, None],
            "expval": ["red", "yellow", "green", "blue"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [None, None, None, None],
            "inval1": [],
            "inidx2": [0, 1, 2, 3],
            "inval2": ["red", "yellow", "green", "blue"],
            "expidx": [None, None, None, None, 0, 1, 2, 3],
            "expval": ["red", "yellow", "green", "blue"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [None, None, None, None],
            "inval1": [],
            "inidx2": [None, None, None, None],
            "inval2": ["orange"],
            "expidx": [None, None, None, None, None, None, None, None],
            "expval": [],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [None, None, None, None],
            "inval1": ["orange"],
            "inidx2": [None, None, None, None],
            "inval2": [],
            "expidx": [None, None, None, None, None, None, None, None],
            "expval": [],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [None, None, None, None],
            "inval1": [],
            "inidx2": [None, None, None, None],
            "inval2": [],
            "expidx": [None, None, None, None, None, None, None, None],
            "expval": [],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [0, None, 2, 3],
            "inval1": ["red", "yellow", "green", "blue"],
            "inidx2": [],
            "inval2": [],
            "expidx": [0, None, 1, 2],
            "expval": ["red", "green", "blue"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [0, 1],
            "inval1": ["green", "red"],
            "inidx2": [0, None, 2, 3],
            "inval2": ["red", "yellow", "green", "blue"],
            "expidx": [0, 1, 1, None, 0, 2],
            "expval": ["green", "red", "blue"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [0, 1, None, 3],
            "inval1": ["RED", "YELLOW", "BLACK", "GREEN"],
            "inidx2": [None, None, 0, 2],
            "inval2": ["YELLOW", "GREY", "BLUE"],
            "expidx": [0, 1, None, 2, None, None, 1, 3],
            "expval": ["RED", "YELLOW", "GREEN", "BLUE"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [1, 0, 1, None, 1, 1, 0],
            "inval1": ["red", "blue"],
            "inidx2": [0, 1, 1, 0, None, None, None, 0, 1],
            "inval2": ["orange", "blue"],
            "expidx": [1, 0, 1, None, 1, 1, 0, 2, 1, 1, 2, None, None, None, 2, 1],
            "expval": ["red", "blue", "orange"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [3, None, 5],
            "inval1": ["purple", "orange", "grey", "red", "blue", "brown"],
            "inidx2": [4, None, 2],
            "inval2": ["lavender", "teal", "salmon", "puce", "mauve"],
            "expidx": [0, None, 1, 3, None, 2],
            "expval": ["red", "brown", "salmon", "mauve"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [0, None, 2, 3],
            "inval1": ["red", "yellow", "green", "blue"],
            "inidx2": [0, 1, None, 3],
            "inval2": ["green", "blue", "orange", "purple"],
            "expidx": [0, None, 1, 2, 1, 2, None, 3],
            "expval": ["red", "green", "blue", "purple"],
        },
        {
            "arrow_type": pa.string(),
            "inidx1": [3, None, 5],
            "inval1": ["purple", "orange", "grey", "red", "blue", "brown"],
            "inidx2": [],
            "inval2": [],
            "expidx": [0, None, 1],
            "expval": ["red", "brown"],
        },
        # Int-valued enums
        {
            "arrow_type": pa.int64(),
            "inidx1": [0, 1, 2, 3],
            "inval1": [111, 222, 333, 444],
            "inidx2": [0, 1, 2, 3],
            "inval2": [111, 222, 333, 444],
            "expidx": [0, 1, 2, 3, 0, 1, 2, 3],
            "expval": [111, 222, 333, 444],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [0, 1, 2, 3],
            "inval1": [111, 222, 333, 444],
            "inidx2": [],
            "inval2": [],
            "expidx": [0, 1, 2, 3],
            "expval": [111, 222, 333, 444],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [],
            "inval1": [],
            "inidx2": [0, 1, 2, 3],
            "inval2": [111, 222, 333, 444],
            "expidx": [0, 1, 2, 3],
            "expval": [111, 222, 333, 444],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [0, 1, 2, 3],
            "inval1": [111, 222, 333, 444],
            "inidx2": [0, 1, 2, 3],
            "inval2": [555, 666, 1111, 1222],
            "expidx": [0, 1, 2, 3, 4, 5, 6, 7],
            "expval": [
                111,
                222,
                333,
                444,
                555,
                666,
                1111,
                1222,
            ],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [0, 1, 2, 3],
            "inval1": [111, 222, 333, 444],
            "inidx2": [3, 2, 1, 0],
            "inval2": [111, 222, 333, 444],
            "expidx": [0, 1, 2, 3, 3, 2, 1, 0],
            "expval": [111, 222, 333, 444],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [0, 1, 2, 3],
            "inval1": [111, 222, 333, 444],
            "inidx2": [0, 1, 2, 3],
            "inval2": [555, 777, 333, 444],
            "expidx": [0, 1, 2, 3, 4, 5, 2, 3],
            "expval": [111, 222, 333, 444, 555, 777],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [0, 1, 2, 3],
            "inval1": [111, 222, 333, 444],
            "inidx2": [0, None, 2, 3],
            "inval2": [555, 777, 333, 444],
            "expidx": [0, 1, 2, 3, 4, None, 2, 3],
            "expval": [111, 222, 333, 444, 555],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [0, None, 2, 3],
            "inval1": [111, 222, 333, 444],
            "inidx2": [0, 1, 2, 3],
            "inval2": [555, 777, 333, 444],
            "expidx": [0, None, 1, 2, 3, 4, 1, 2],
            "expval": [111, 333, 444, 555, 777],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [0, 1, 2, 3],
            "inval1": [111, 222, 333, 444],
            "inidx2": [0, None, 0, 0],
            "inval2": [555],
            "expidx": [0, 1, 2, 3, 4, None, 4, 4],
            "expval": [111, 222, 333, 444, 555],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [0, 1, 2, 3],
            "inval1": [111, 222, 333, 444],
            "inidx2": [None, None, None, None],
            "inval2": [555],
            "expidx": [0, 1, 2, 3, None, None, None, None],
            "expval": [111, 222, 333, 444],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [0, 1, 2, 3],
            "inval1": [111, 222, 333, 444],
            "inidx2": [0, 1, 2, 2],
            "inval2": [222, 333, 444],
            "expidx": [0, 1, 2, 3, 1, 2, 3, 3],
            "expval": [111, 222, 333, 444],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [None, None, None, None],
            "inval1": [111, 222, 333, 444],
            "inidx2": [None, None, None, None],
            "inval2": [222, 333, 444],
            "expidx": [None, None, None, None, None, None, None, None],
            "expval": [],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [None, None, None, None],
            "inval1": [555],
            "inidx2": [None, None, None, None],
            "inval2": [555],
            "expidx": [None, None, None, None, None, None, None, None],
            "expval": [],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [0, 1, 2, 3],
            "inval1": [111, 222, 333, 444],
            "inidx2": [None, None, None, None],
            "inval2": [],
            "expidx": [0, 1, 2, 3, None, None, None, None],
            "expval": [111, 222, 333, 444],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [None, None, None, None],
            "inval1": [],
            "inidx2": [0, 1, 2, 3],
            "inval2": [111, 222, 333, 444],
            "expidx": [None, None, None, None, 0, 1, 2, 3],
            "expval": [111, 222, 333, 444],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [None, None, None, None],
            "inval1": [],
            "inidx2": [None, None, None, None],
            "inval2": [555],
            "expidx": [None, None, None, None, None, None, None, None],
            "expval": [],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [None, None, None, None],
            "inval1": [555],
            "inidx2": [None, None, None, None],
            "inval2": [],
            "expidx": [None, None, None, None, None, None, None, None],
            "expval": [],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [None, None, None, None],
            "inval1": [],
            "inidx2": [None, None, None, None],
            "inval2": [],
            "expidx": [None, None, None, None, None, None, None, None],
            "expval": [],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [0, None, 2, 3],
            "inval1": [111, 222, 333, 444],
            "inidx2": [],
            "inval2": [],
            "expidx": [0, None, 1, 2],
            "expval": [111, 333, 444],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [0, 1],
            "inval1": [333, 111],
            "inidx2": [0, None, 2, 3],
            "inval2": [111, 222, 333, 444],
            "expidx": [0, 1, 1, None, 0, 2],
            "expval": [333, 111, 444],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [3, 1, 0, 1, 3, 0, None, None, 3],
            "inval1": [111, 222, 999, 333],
            "inidx2": [None, None, 0, 2, 2, 0, 2, None],
            "inval2": [222, 777, 444],
            "expidx": [
                2,
                1,
                0,
                1,
                2,
                0,
                None,
                None,
                2,
                None,
                None,
                1,
                3,
                3,
                1,
                3,
                None,
            ],
            "expval": [111, 222, 333, 444],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [1, 0, 1, None, 1, 1, 0],
            "inval1": [111, 444],
            "inidx2": [0, 1, 1, 0, None, None, None, 0, 1],
            "inval2": [555, 444],
            "expidx": [1, 0, 1, None, 1, 1, 0, 2, 1, 1, 2, None, None, None, 2, 1],
            "expval": [111, 444, 555],
        },
        {
            "arrow_type": pa.int64(),
            "inidx1": [3, None, 5],
            "inval1": [666, 555, 777, 111, 444, 888],
            "inidx2": [4, None, 2],
            "inval2": [1333, 1444, 1555, 1666, 1777],
            "expidx": [0, None, 1, 3, None, 2],
            "expval": [111, 888, 1555, 1777],
        },
        # Bool-valued enums
        {
            "arrow_type": pa.bool_(),
            "inidx1": [],
            "inval1": [],
            "inidx2": [],
            "inval2": [],
            "expidx": [],
            "expval": [],
        },
        {
            "arrow_type": pa.bool_(),
            "inidx1": [0, 0, 0],
            "inval1": [True],
            "inidx2": [0, 0],
            "inval2": [True],
            "expidx": [0, 0, 0, 0, 0],
            "expval": [True],
        },
        {
            "arrow_type": pa.bool_(),
            "inidx1": [0, 0, 0],
            "inval1": [True],
            "inidx2": [0, 0],
            "inval2": [False],
            "expidx": [0, 0, 0, 1, 1],
            "expval": [True, False],
        },
        {
            "arrow_type": pa.bool_(),
            "inidx1": [0, None, 0],
            "inval1": [True],
            "inidx2": [None, 0],
            "inval2": [False],
            "expidx": [0, None, 0, None, 1],
            "expval": [True, False],
        },
        {
            "arrow_type": pa.bool_(),
            "inidx1": [0, None, 0],
            "inval1": [True, False],
            "inidx2": [None, 0],
            "inval2": [False],
            "expidx": [0, None, 0, None, 1],
            "expval": [True, False],
        },
        {
            "arrow_type": pa.bool_(),
            "inidx1": [0, None, 0],
            "inval1": [True, False],
            "inidx2": [None, 1],
            "inval2": [True, False],
            "expidx": [0, None, 0, None, 1],
            "expval": [True, False],
        },
        {
            "arrow_type": pa.bool_(),
            "inidx1": [0, None, 1],
            "inval1": [True, False],
            "inidx2": [None, 1, 0],
            "inval2": [True, False],
            "expidx": [0, None, 1, None, 1, 0],
            "expval": [True, False],
        },
    ],
)
@pytest.mark.parametrize("ordered", [True, False])
def test_extend_enumeration_null_indices(tmp_path, config, ordered):
    uri = tmp_path.as_posix()

    arrow_type = config["arrow_type"]
    inidx1 = config["inidx1"]
    inval1 = config["inval1"]
    inidx2 = config["inidx2"]
    inval2 = config["inval2"]
    expidx = config["expidx"]
    expval = config["expval"]

    n1 = len(inidx1)
    n2 = len(inidx2)

    domain = [[0, 99]]

    schema = pa.schema(
        {
            "soma_joinid": pa.int64(),
            "enum_test": pa.dictionary(pa.int32(), arrow_type, ordered=ordered),
        }
    )

    data1 = pa.Table.from_pydict(
        {
            "soma_joinid": pa.array(list(range(n1)), type=pa.int64()),
            "enum_test": pa.DictionaryArray.from_arrays(
                pa.array(inidx1, type=pa.int32()),
                pa.array(inval1, type=arrow_type),
            ),
        }
    )

    with soma.DataFrame.create(uri, schema=schema, domain=domain) as sdf:
        sdf.write(data1)

    data2 = pa.Table.from_pydict(
        {
            "soma_joinid": pa.array(list(range(n1, n1 + n2)), type=pa.int64()),
            "enum_test": pa.DictionaryArray.from_arrays(
                pa.array(inidx2, type=pa.int32()),
                pa.array(inval2, type=arrow_type),
            ),
        }
    )

    with soma.DataFrame.open(uri, "w") as sdf:
        sdf.write(data2)

    incol1 = data1["enum_test"].to_pylist()
    incol2 = data2["enum_test"].to_pylist()
    expcol = incol1 + incol2

    with soma.DataFrame.open(uri, "r") as sdf:
        table = sdf.read().concat()
        column = table["enum_test"].combine_chunks()
        outcol = column.to_pylist()
        outidx = column.indices.to_pylist()
        outval = column.dictionary.to_pylist()
        getval = sdf.get_enumeration_values(["enum_test"])["enum_test"].to_pylist()

        assert outcol == expcol
        assert outidx == expidx
        assert outval == expval
        assert getval == expval


def test_extend_enumeration_empty(tmp_path):
    uri = tmp_path.as_posix()
    schema = pa.schema(
        {
            "soma_joinid": pa.int64(),
            "string_enum": pa.dictionary(pa.int32(), pa.large_string()),
        }
    )
    domain = [[0, 7]]

    data1 = pa.Table.from_pydict(
        {
            "soma_joinid": pa.array([0, 1, 2, 3], type=pa.int64()),
            "string_enum": pa.DictionaryArray.from_arrays(
                pa.array([0, 1, 2, 3], type=pa.int32()),
                pa.array(["red", "yellow", None, "blue"], type=pa.large_string()),
            ),
        }
    )

    with soma.DataFrame.create(uri, schema=schema, domain=domain) as sdf:
        with pytest.raises(soma.SOMAError):
            sdf.write(data1)

    data1 = pa.Table.from_pydict(
        {
            "soma_joinid": pa.array([0, 1, 2, 3], type=pa.int64()),
            "string_enum": pa.DictionaryArray.from_arrays(
                pa.array([0, 1, 2, 3], type=pa.int32()),
                pa.array(["red", "yellow", "green", "blue"], type=pa.large_string()),
            ),
        }
    )

    with soma.DataFrame.open(uri, "w") as sdf:
        sdf.write(data1)

    with soma.DataFrame.open(uri, "r") as sdf:
        assert sdf.get_enumeration_values(["string_enum"])[
            "string_enum"
        ].to_pylist() == [
            "red",
            "yellow",
            "green",
            "blue",
        ]

    with soma.DataFrame.open(uri, "w") as sdf:
        sdf.extend_enumeration_values(
            {"string_enum": pa.array([], type=pa.large_string())}
        )


@pytest.fixture
def simple_data_frame(tmp_path):
    """
    A pytest fixture which creates a simple DataFrame for use in tests below.
    """
    schema = pa.schema(
        [
            pa.field("index", pa.int64(), False),
            pa.field("soma_joinid", pa.int64()),
            pa.field("A", pa.int64()),
            pa.field("B", pa.float64()),
            pa.field("C", pa.large_string()),
        ]
    )
    index_column_names = ["index"]
    with soma.DataFrame.create(
        tmp_path.as_posix(),
        schema=schema,
        index_column_names=index_column_names,
        domain=[[0, 3]],
    ) as sdf:
        data = {
            "index": [0, 1, 2, 3],
            "soma_joinid": [10, 11, 12, 13],
            "A": [10, 11, 12, 13],
            "B": [100.1, 200.2, 300.3, 400.4],
            "C": ["this", "is", "a", "test"],
        }
        n_data = len(data["index"])
        rb = pa.Table.from_pydict(data)
        sdf.write(rb)
    sdf = soma.DataFrame.open(tmp_path.as_posix())

    return (sdf, n_data)


@pytest.mark.parametrize(
    "ids",
    [
        [None],
        [0],
        [[1, 3]],
    ],
)
@pytest.mark.parametrize(
    "col_names",
    [
        ["A"],
        ["B"],
        ["A", "B"],
        ["index"],
        ["index", "A", "B", "C"],
        ["soma_joinid"],
        ["soma_joinid", "A"],
        None,
    ],
)
def test_DataFrame_read_column_names(simple_data_frame, ids, col_names):
    sdf, n_data = simple_data_frame

    expected_schema = (
        sdf.schema
        if col_names is None
        else pa.schema([sdf.schema.field(name) for name in col_names])
    )

    def _check_tbl(tbl, col_names, ids):
        assert tbl.num_columns == len(sdf.schema if col_names is None else col_names)

        if ids is None:
            assert tbl.num_rows == n_data
        elif ids[0] is None:
            assert tbl.num_rows == n_data
        elif isinstance(ids[0], int):
            assert tbl.num_rows == 1
        else:
            assert tbl.num_rows == len(ids[0])

        assert tbl.schema == expected_schema

    # TileDB ASCII -> Arrow large_string
    _check_tbl(sdf.read(ids, column_names=col_names).concat(), col_names, ids)
    _check_tbl(sdf.read(column_names=col_names).concat(), col_names, None)

    # pa.Table.from_pandas infers nullability from the data if a schema is not
    # provided. If there are no null values in the data, then the data is marked
    # as non-nullable. To ensure nullability is preserved, explicitly pass in
    # the schema

    # TileDB ASCII -> Pandas string -> Arrow string (not large_string)
    _check_tbl(
        pa.Table.from_pandas(
            pd.concat(
                [tbl.to_pandas() for tbl in sdf.read(ids, column_names=col_names)],
            ),
            schema=expected_schema,
        ),
        col_names,
        ids,
    )
    _check_tbl(
        pa.Table.from_pandas(
            sdf.read(column_names=col_names).concat().to_pandas(),
            schema=expected_schema,
        ),
        col_names,
        None,
    )


def test_empty_dataframe(tmp_path):
    soma.DataFrame.create(
        (tmp_path / "A").as_posix(),
        schema=pa.schema([("a", pa.int32())]),
        index_column_names=["a"],
    ).close()
    with soma.DataFrame.open((tmp_path / "A").as_posix()) as a:
        # Must not throw
        assert len(next(a.read())) == 0
        assert len(a.read().concat()) == 0
        assert len(next(a.read()).to_pandas()) == 0
        assert len(a.read().concat().to_pandas()) == 0
        assert isinstance(a.read().concat().to_pandas(), pd.DataFrame)

    with pytest.raises(ValueError):
        # illegal column name
        soma.DataFrame.create(
            (tmp_path / "B").as_posix(),
            schema=pa.schema([("a", pa.int32()), ("soma_bogus", pa.int32())]),
            index_column_names=["a"],
        )


def test_columns(tmp_path):
    """
    1. soma_joinid is int64
    2. soma_joinid will be added by default, if missing in call to create
    3. soma_joinid is explicit in keys/schema
    4. No other soma_ ids allowed
    """

    A = soma.DataFrame.create(
        (tmp_path / "A").as_posix(),
        schema=pa.schema([("a", pa.int32())]),
        index_column_names=["a"],
    )
    assert sorted(A.keys()) == sorted(["a", "soma_joinid"])
    assert A.schema.field("soma_joinid").type == pa.int64()

    with pytest.raises(ValueError):
        soma.DataFrame.create(
            (tmp_path / "B").as_posix(),
            schema=pa.schema([("a", pa.int32()), ("soma_joinid", pa.float32())]),
            index_column_names=["a"],
        )

    D = soma.DataFrame.create(
        (tmp_path / "D").as_posix(),
        schema=pa.schema([("a", pa.int32()), ("soma_joinid", pa.int64())]),
        index_column_names=["a"],
    )
    assert sorted(D.keys()) == sorted(["a", "soma_joinid"])
    assert D.schema.field("soma_joinid").type == pa.int64()

    with pytest.raises(ValueError):
        soma.DataFrame.create(
            (tmp_path / "E").as_posix(),
            schema=pa.schema(
                [("a", pa.int32()), ("soma_is_a_reserved_prefix", pa.bool_())]
            ),
            index_column_names=["a"],
        )


@pytest.fixture
def make_dataframe(request):
    index_type, domain = request.param

    index = {
        pa.string(): ["A", "B", "C"],
        pa.large_string(): ["A", "B", "C"],
        pa.binary(): [b"A", b"B", b"C"],
        pa.large_binary(): [b"A", b"B", b"C"],
        **{
            t: np.arange(3, dtype=t.to_pandas_dtype())
            for t in (
                pa.int8(),
                pa.uint8(),
                pa.int16(),
                pa.uint16(),
                pa.int32(),
                pa.uint32(),
                pa.int64(),
                pa.uint64(),
                pa.float32(),
                pa.float64(),
            )
        },
    }[index_type]

    df = pd.DataFrame(
        data={
            "index": index,
            "soma_joinid": np.arange(3, dtype=np.int64),
            "ascii": ["aa", "bbb", "cccccc"],
            "float32": np.array([0.0, 1.1, 2.2], np.float32),
        }
    )
    return [pa.Table.from_pandas(df), domain]


@pytest.mark.parametrize(
    "make_dataframe",
    [
        [pa.float32(), [[-1000, 1000]]],
        [pa.float64(), [[-1000, 1000]]],
        [pa.int8(), [[-100, 100]]],
        [pa.uint8(), [[0, 100]]],
        [pa.int16(), [[-1000, 1000]]],
        [pa.uint16(), [[0, 1000]]],
        [pa.int32(), [[-1000, 1000]]],
        [pa.uint32(), [[0, 1000]]],
        [pa.int64(), [[-1000, 1000]]],
        [pa.uint64(), [[0, 1000]]],
        [pa.string(), [None]],
        [pa.large_string(), [None]],
        [pa.binary(), [None]],
        [pa.large_binary(), [None]],
    ],
    indirect=True,
)
def test_index_types(tmp_path, make_dataframe):
    """Verify that the index columns can be of various types"""
    sdf = soma.DataFrame.create(
        tmp_path.as_posix(),
        schema=make_dataframe[0].schema,
        index_column_names=["index"],
        domain=make_dataframe[1],
    )
    sdf.write(make_dataframe[0])


def make_multiply_indexed_dataframe(
    tmp_path, index_column_names: list[str], domain: List[Any]
):
    """
    Creates a variably-indexed DataFrame for use in tests below.
    """
    schema = pa.schema(
        [
            # Note: exhaustive type-coverage is in a separate test case.
            # (As of this writing: test_dataframe_column_indexing.py)
            ("0_thru_5", pa.int64()),
            ("strings_aaa", pa.string()),
            ("zero_one", pa.int64()),
            ("thousands", pa.int64()),
            ("both_signs", pa.int64()),
            ("soma_joinid", pa.int64()),
            ("A", pa.int64()),
        ]
    )

    sdf = soma.DataFrame.create(
        uri=tmp_path.as_posix(),
        schema=schema,
        index_column_names=index_column_names,
        domain=domain,
    )

    data: dict[str, list] = {
        "0_thru_5": [0, 1, 2, 3, 4, 5],
        "strings_aaa": ["aaa", "aaa", "bbb", "bbb", "ccc", "ccc"],
        "zero_one": [0, 1, 0, 1, 0, 1],
        "thousands": [1000, 2000, 1000, 1000, 1000, 1000],
        "both_signs": [-1, -2, -3, 1, 2, 3],
        "soma_joinid": [10, 11, 12, 13, 14, 15],
        "A": [10, 11, 12, 13, 14, 15],
    }

    n_data = len(data["0_thru_5"])
    sdf.write(pa.Table.from_pandas(pd.DataFrame(data=data)))

    return (schema, sdf, n_data)


@pytest.mark.parametrize(
    "io",
    [
        {
            "name": "1D indexing slot is None",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [None],
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        {
            "name": "1D indexing slot is int",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [0],
            "A": [10],
            "throws": None,
        },
        {
            "name": "1D no results for 100",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [100],
            "A": [],
            "throws": None,
        },
        {
            "name": "1D no results for -100",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [-100],
            "A": [],
            "throws": None,
        },
        {
            "name": "1D indexing slot is list",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [[1, 3]],
            "A": [11, 13],
            "throws": None,
        },
        {
            "name": "1D no results for -100, 100",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [[-100, 100]],
            "A": [],
            "throws": None,
        },
        {
            "name": "1D empty list returns empty results",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [[]],
            "A": [],
            "throws": None,
        },
        {
            "name": "1D indexing slot is tuple",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [(1, 3)],
            "A": [11, 13],
            "throws": None,
        },
        {
            "name": "1D indexing slot is range",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [range(1, 3)],
            "A": [11, 12],
            "throws": None,
        },
        {
            "name": "1D indexing slot is pa.ChunkedArray",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [pa.chunked_array(pa.array([1, 3]))],
            "A": [11, 13],
            "throws": None,
        },
        {
            "name": "1D indexing slot is pa.Array",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [pa.array([1, 3])],
            "A": [11, 13],
            "throws": None,
        },
        # 1D: indexing slot is np.ndarray
        {
            "name": "1D indexing slot is np.ndarray",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [np.asarray([1, 3])],
            "A": [11, 13],
            "throws": None,
        },
        {
            "name": "1D indexing by 2D np.ndarray",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [
                np.asarray([[1, 3], [2, 4]])
            ],  # Error since 2D array in the slot
            "A": [11, 13],
            "throws": ValueError,
        },
        {
            "name": "1D indexing by slice(None)",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [
                slice(None)
            ],  # Indexing slot is none-slice i.e. `[:]` which is like None
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        {
            "name": "1D indexing by empty coords",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [],
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        {
            "name": "1D indexing by 1:3",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [slice(1, 3)],  # Indexing slot is double-ended slice
            "A": [11, 12, 13],
            "throws": None,
        },
        {
            "name": "1D indexing by [:3]",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [slice(None, 3)],  # Half-slice
            "A": [10, 11, 12, 13],
            "throws": None,
        },
        {
            "name": "1D indexing by [2:]",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [slice(2, None)],  # Half-slice
            "A": [12, 13, 14, 15],
            "throws": None,
        },
        {
            "name": "1D indexing with negatives",
            "index_column_names": ["both_signs"],
            "domain": [[-10, 10]],
            "coords": [slice(-2, 1)],
            "A": [11, 10, 13],
            "throws": None,
        },
        {
            "name": "1D indexing by ['bbb':'c']",
            "index_column_names": ["strings_aaa", "zero_one"],
            "domain": [None, [0, 1]],
            "coords": [slice("bbb", "c")],
            "A": [12, 13],
            "throws": None,
        },
        {
            "name": "1D indexing by ['ccc':]",
            "index_column_names": ["strings_aaa", "zero_one"],
            "domain": [None, [0, 1]],
            "coords": [slice("ccc", None)],
            "A": [14, 15],
            "throws": None,
        },
        {
            "name": "1D indexing by [:'bbd']",
            "index_column_names": ["strings_aaa", "zero_one"],
            "domain": [None, [0, 1]],
            "coords": [slice("bbd")],
            "A": [10, 11, 12, 13],
            "throws": None,
        },
        {
            "name": "1D indexing with one partition",
            "index_column_names": ["0_thru_5"],
            "domain": [[0, 8]],
            "coords": [slice(2, None)],
            "partitions": somacore.IOfN(0, 1),
            "A": [12, 13, 14, 15],
            "throws": None,
        },
        {
            "name": "partitioned reads unimplemented",
            "index_column_names": ["0_thru_5"],
            "domain": [[0, 8]],
            "coords": [],
            "partitions": somacore.IOfN(1, 2),
            "A": None,
            "throws": ValueError,
        },
        {
            "name": "steps forbidden",
            "index_column_names": ["0_thru_5"],
            "domain": [[0, 8]],
            "coords": [slice(1, 5, 2)],
            "A": None,
            "throws": ValueError,
        },
        {
            "name": "slice must overlap domain (negative)",
            "index_column_names": ["soma_joinid"],
            "domain": [[0, 59]],
            "coords": [slice(-2, -1)],
            "A": None,
            "throws": ValueError,
        },
        {
            "name": "backwards slice",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [slice(1, 0)],
            "A": None,
            "throws": ValueError,
        },
        {
            "name": "too many columns",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [(1,), (2,)],
            "A": None,
            "throws": ValueError,
        },
        {
            "name": "wrong coords type",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": "bogus",
            "A": None,
            "throws": TypeError,
        },
        {
            "name": "bad index type dict",
            "index_column_names": ["0_thru_5"],
            "domain": [[-1000, 1000]],
            "coords": [{"bogus": True}],
            "A": None,
            # Disable Typeguard while asserting this error, otherwise a typeguard.TypeCheckError is
            # raised (though that's not what would happen in production)
            "throws": (TypeError, False),
        },
        {
            "name": "bad index type bool",
            "index_column_names": ["strings_aaa", "zero_one"],
            "domain": [None, [0, 1]],
            "coords": [[True], slice(None)],
            "A": None,
            "throws": TypeError,
        },
        {
            "name": "2D index empty",
            "index_column_names": ["strings_aaa", "zero_one"],
            "domain": [None, [0, 1]],
            "coords": (),
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        {
            "name": "2D index None",
            "index_column_names": ["strings_aaa", "zero_one"],
            "domain": [None, [0, 1]],
            "coords": [None, None],
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        {
            "name": "2D index 0, 0",
            "index_column_names": ["0_thru_5", "zero_one"],
            "domain": [[-1000, 1000], [0, 1]],
            "coords": [0, 0],
            "A": [10],
            "throws": None,
        },
        {
            "name": "2D index str, int",
            "index_column_names": ["strings_aaa", "zero_one"],
            "domain": [None, [0, 1]],
            "coords": [["aaa"], 0],
            "A": [10],
            "throws": None,
        },
        {
            "name": "2D index str, not sequence[str]",
            "index_column_names": ["strings_aaa", "zero_one"],
            "domain": [None, [0, 1]],
            "coords": ["aaa", 0],
            "A": [10],
            "throws": None,
        },
        {
            "name": "2D index List[str]",
            "index_column_names": ["strings_aaa", "zero_one"],
            "domain": [None, [0, 1]],
            "coords": [["aaa", "ccc"], None],
            "A": [10, 11, 14, 15],
            "throws": None,
        },
        {
            "name": "3D index List[str]",
            "index_column_names": ["strings_aaa", "zero_one", "thousands"],
            "domain": [None, [0, 1], [0, 9999]],
            "coords": [["aaa", "ccc"], None, None],
            "A": [10, 11, 14, 15],
            "throws": None,
        },
        {
            "name": "3D index mixed",
            "index_column_names": ["strings_aaa", "zero_one", "thousands"],
            "domain": [None, [0, 1], [0, 9999]],
            "coords": [("aaa", "ccc"), None, np.asarray([2000, 9999])],
            "A": [11],
            "throws": None,
        },
        {
            "name": "value filter good",
            "index_column_names": ["0_thru_5", "strings_aaa"],
            "domain": [[-1000, 1000], None],
            "coords": [None, ("ccc", "zzz")],
            "value_filter": "soma_joinid > 13",
            "A": [14, 15],
        },
        {
            "name": "value filter bad",
            "index_column_names": ["0_thru_5", "strings_aaa"],
            "domain": [[-1000, 1000], None],
            "coords": [None, ("bbb", "zzz")],
            "value_filter": "quick brown fox",
            "A": None,
            "throws": soma.SOMAError,
        },
    ],
    ids=lambda d: d.get("name"),
)
def test_read_indexing(tmp_path, io):
    """Test various ways of indexing on read"""

    schema, sdf, n_data = make_multiply_indexed_dataframe(
        tmp_path, io["index_column_names"], io["domain"]
    )
    with soma.DataFrame.open(uri=sdf.uri) as sdf:
        assert list(sdf.index_column_names) == io["index_column_names"]

        read_kwargs = {"column_names": ["A"]}
        read_kwargs.update(
            {k: io[k] for k in ("coords", "partitions", "value_filter") if k in io}
        )

        # `throws` can be `type[Exception]`, or `(type[Exception], bool)` indicating explicitly
        # whether Typeguard should be enabled during the `with raises` check.
        throws = io.get("throws", None)
        if throws:
            if isinstance(throws, tuple) and not throws[1]:
                # Disable Typeguard, verify actual runtime error type (avoid
                # `typeguard.TypeCheckError` short-circuit)
                throws = throws[0]
                throws_ctx = raises_no_typeguard
            else:
                throws_ctx = pytest.raises
        else:
            throws_ctx = None

        if throws_ctx:
            with throws_ctx(throws):
                next(sdf.read(**read_kwargs))
        else:
            table = next(sdf.read(**read_kwargs))
            assert table["A"].to_pylist() == io["A"]

        if throws_ctx:
            with throws_ctx(throws):
                next(sdf.read(**read_kwargs)).to_pandas()
        else:
            table = next(sdf.read(**read_kwargs)).to_pandas()
            assert table["A"].to_list() == io["A"]


def test_write_categorical_types(tmp_path):
    """
    Verify that write path accepts categoricals
    """
    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            (
                "string-ordered",
                pa.dictionary(pa.int8(), pa.large_string(), ordered=True),
            ),
            ("string-unordered", pa.dictionary(pa.int8(), pa.large_string())),
            ("string-compat", pa.large_string()),
            ("int-ordered", pa.dictionary(pa.int8(), pa.int64(), ordered=True)),
            ("int-unordered", pa.dictionary(pa.int8(), pa.int64())),
            ("int-compat", pa.int64()),
            ("bool-ordered", pa.dictionary(pa.int8(), pa.bool_(), ordered=True)),
            ("bool-unordered", pa.dictionary(pa.int8(), pa.bool_())),
            ("bool-compat", pa.bool_()),
        ]
    )
    with soma.DataFrame.create(
        tmp_path.as_posix(),
        schema=schema,
        index_column_names=["soma_joinid"],
        domain=[[0, 3]],
    ) as sdf:
        df = pd.DataFrame(
            data={
                "soma_joinid": [0, 1, 2, 3],
                "string-ordered": pd.Categorical(
                    ["a", "b", "a", "b"], ordered=True, categories=["b", "a"]
                ),
                "string-unordered": pd.Categorical(
                    ["a", "b", "a", "b"], ordered=False, categories=["b", "a"]
                ),
                "string-compat": pd.Categorical(
                    ["a", "b", "a", "b"], ordered=False, categories=["a", "b"]
                ),
                "int-ordered": pd.Categorical(
                    [777777777, 888888888, 777777777, 888888888],
                    ordered=True,
                    categories=[888888888, 777777777],
                ),
                "int-unordered": pd.Categorical(
                    [777777777, 888888888, 777777777, 888888888],
                    ordered=False,
                    categories=[888888888, 777777777],
                ),
                "int-compat": pd.Categorical(
                    [777777777, 888888888, 777777777, 888888888],
                    ordered=False,
                    categories=[777777777, 888888888],
                ),
                "bool-ordered": pd.Categorical(
                    [True, False, True, False],
                    ordered=True,
                    categories=[True, False],
                ),
                "bool-unordered": pd.Categorical(
                    [True, False, True, False],
                    ordered=False,
                    categories=[True, False],
                ),
                "bool-compat": pd.Categorical(
                    [True, False, True, False],
                    ordered=False,
                    categories=[True, False],
                ),
            }
        )
        sdf.write(pa.Table.from_pandas(df))

    with soma.DataFrame.open(tmp_path.as_posix()) as sdf:
        assert (df == sdf.read().concat().to_pandas()).all().all()


# def test_write_categorical_dims(tmp_path):
#     """
#     Categories are not supported as dims. Here we test our handling of what we
#     do when we are given them as input.
#     """
#     schema = pa.schema(
#         [
#             ("soma_joinid", pa.int64()),
#             ("string", pa.dictionary(pa.int8(), pa.large_string())),
#         ]
#     )
#     with soma.DataFrame.create(
#         tmp_path.as_posix(),
#         schema=schema,
#         index_column_names=["soma_joinid"],
#     ) as sdf:
#         df = pd.DataFrame(
#             data={
#                 "soma_joinid": pd.Categorical([0, 1, 2, 3], categories=[0, 1, 2, 3]),
#                 "string": pd.Categorical(["a", "b", "a", "b"], categories=["b", "a"]),
#             }
#         )
#         sdf.write(pa.Table.from_pandas(df))

#     with soma.DataFrame.open(tmp_path.as_posix()) as sdf:
#         assert (df == sdf.read().concat().to_pandas()).all().all()


@pytest.mark.parametrize("index_type", [pa.int8(), pa.int16(), pa.int32(), pa.int64()])
def test_write_categorical_dim_extend(tmp_path, index_type):
    """
    Introduce new categorical values in each subsequent write.
    """
    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            ("string", pa.dictionary(pa.int8(), pa.large_string())),
        ]
    )

    df1 = pd.DataFrame(
        data={
            "soma_joinid": [0, 1, 2, 3],
            "string": pd.Categorical(["a", "b", "a", "b"], categories=["b", "a"]),
        }
    )

    df2 = pd.DataFrame(
        data={
            "soma_joinid": [4, 5],
            "string": pd.Categorical(["c", "b"], categories=["b", "c"]),
        }
    )

    with soma.DataFrame.create(
        tmp_path.as_posix(),
        schema=schema,
        index_column_names=["soma_joinid"],
        domain=[[0, 5]],
    ) as sdf:
        table = pa.Table.from_pandas(df1)
        dtype = pa.dictionary(index_type, pa.string())
        set_index_type = table.set_column(
            1,
            pa.field("string", dtype),
            pa.array(["a", "b", "a", "b"], dtype),
        )
        sdf.write(set_index_type)

    with soma.DataFrame.open(tmp_path.as_posix(), "w") as sdf:
        sdf.write(pa.Table.from_pandas(df2))

    # https://stackoverflow.com/questions/45639350/retaining-categorical-dtype-upon-dataframe-concatenation
    uc = union_categoricals([df1.string, df2.string])
    df1.string = pd.Categorical(df1.string, categories=uc.categories)
    df2.string = pd.Categorical(df2.string, categories=uc.categories)
    expected_df = pd.concat((df1, df2), ignore_index=True)

    with soma.DataFrame.open(tmp_path.as_posix()) as sdf:
        data = sdf.read().concat()
        assert expected_df.compare(data.to_pandas()).empty


def test_result_order(tmp_path):
    # cf. https://docs.tiledb.com/main/background/key-concepts-and-data-format#data-layout
    schema = pa.schema(
        [
            ("row", pa.int64()),
            ("col", pa.int64()),
            ("soma_joinid", pa.int64()),
        ]
    )
    with soma.DataFrame.create(
        uri=tmp_path.as_posix(),
        schema=schema,
        index_column_names=["row", "col"],
        domain=[[0, 15], [0, 15]],
    ) as sdf:
        data = {
            "row": [0] * 4 + [1] * 4 + [2] * 4 + [3] * 4,
            "col": [0, 1, 2, 3] * 4,
            "soma_joinid": list(range(16)),
        }
        sdf.write(pa.Table.from_pydict(data))

    with soma.DataFrame.open(tmp_path.as_posix()) as sdf:
        table = sdf.read(result_order="row-major").concat().to_pandas()
        assert table["soma_joinid"].to_list() == list(range(16))

        table = sdf.read(result_order="column-major").concat().to_pandas()
        assert table["soma_joinid"].to_list() == [
            0,
            4,
            8,
            12,
            1,
            5,
            9,
            13,
            2,
            6,
            10,
            14,
            3,
            7,
            11,
            15,
        ]

        with raises_no_typeguard(ValueError):
            next(sdf.read(result_order="bogus"))


@pytest.mark.parametrize(
    "create_options,expected_schema_fields",
    (
        (
            {"allows_duplicates": True},
            {
                "validity_filters": [{"COMPRESSION_LEVEL": -1, "name": "RLE"}],
                "allows_duplicates": True,
            },
        ),
        (
            {"allows_duplicates": False},
            {
                "validity_filters": [{"COMPRESSION_LEVEL": -1, "name": "RLE"}],
                "allows_duplicates": False,
            },
        ),
        (
            {"validity_filters": ["NoOpFilter"], "allows_duplicates": False},
            {
                "validity_filters": [{"name": "NOOP"}],
                "allows_duplicates": False,
            },
        ),
    ),
)
def test_create_platform_config_overrides(
    tmp_path, create_options, expected_schema_fields
):
    uri = tmp_path.as_posix()
    soma.DataFrame.create(
        uri,
        schema=pa.schema([pa.field("colA", pa.string())]),
        platform_config={"tiledb": {"create": {**create_options}}},
    ).close()

    with soma.DataFrame.open(tmp_path.as_posix()) as A:
        cfg = A.schema_config_options()
        assert expected_schema_fields["validity_filters"] == json.loads(
            cfg.validity_filters
        )
        assert expected_schema_fields["allows_duplicates"] == cfg.allows_duplicates


@pytest.mark.parametrize("allows_duplicates", [False, True])
@pytest.mark.parametrize("consolidate", [False, True])
def test_timestamped_ops(tmp_path, allows_duplicates, consolidate):
    uri = tmp_path.as_posix()

    platform_config = {"tiledb": {"create": {"allows_duplicates": allows_duplicates}}}
    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            ("float", pa.float64()),
            ("string", pa.large_string()),
        ]
    )

    start = datetime.datetime(2021, 3, 10, 19, 0, tzinfo=datetime.timezone.utc)
    with soma.DataFrame.create(
        uri,
        schema=schema,
        index_column_names=["soma_joinid"],
        domain=[[0, 1]],
        tiledb_timestamp=start,
        platform_config=platform_config,
    ) as sidf:
        data = {
            "soma_joinid": [0],
            "float": [100.1],
            "string": ["apple"],
        }
        sidf.write(pa.Table.from_pydict(data))
        assert sidf.tiledb_timestamp_ms == 1615402800000
        assert sidf.tiledb_timestamp.isoformat() == "2021-03-10T19:00:00+00:00"

    end = start + datetime.timedelta(minutes=3, seconds=25)
    with soma.DataFrame.open(uri=uri, mode="w", tiledb_timestamp=end) as sidf:
        data = {
            "soma_joinid": [0, 1],
            "float": [200.2, 300.3],
            "string": ["ball", "cat"],
        }

        # Without consolidate:
        # * There are two fragments:
        #   o One with tiledb.fragment.FragmentInfoList[i].timestamp_range = (10, 10)
        #   o One with tiledb.fragment.FragmentInfoList[i].timestamp_range = (20, 20)
        # With consolidate:
        # * There is one fragment:
        #   o One with tiledb.fragment.FragmentInfoList[i].timestamp_range = (10, 20)
        sidf.write(
            pa.Table.from_pydict(data),
            soma.TileDBWriteOptions(consolidate_and_vacuum=consolidate),
        )

    # read without timestamp (i.e., after final write) & see final image
    with soma.DataFrame.open(uri) as sidf:
        table = sidf.read().concat()
        if allows_duplicates:
            assert sorted(list(x.as_py() for x in table["soma_joinid"])) == [0, 0, 1]
            assert sorted(list(x.as_py() for x in table["float"])) == [
                100.1,
                200.2,
                300.3,
            ]
            assert sorted(list(x.as_py() for x in table["string"])) == [
                "apple",
                "ball",
                "cat",
            ]
            assert sidf.count == 3
        else:
            assert list(x.as_py() for x in table["soma_joinid"]) == [0, 1]
            assert list(x.as_py() for x in table["float"]) == [200.2, 300.3]
            assert list(x.as_py() for x in table["string"]) == ["ball", "cat"]
            assert sidf.count == 2

    middle = 1615402887987
    # read at t=15 & see only the first write
    with soma.DataFrame.open(tmp_path.as_posix(), tiledb_timestamp=middle) as sidf:
        tab = sidf.read().concat()
        assert list(x.as_py() for x in tab["soma_joinid"]) == [0]
        assert list(x.as_py() for x in tab["float"]) == [100.1]
        assert list(x.as_py() for x in tab["string"]) == ["apple"]
        assert sidf.tiledb_timestamp_ms == 1615402887987
        assert sidf.tiledb_timestamp.isoformat() == "2021-03-10T19:01:27.987000+00:00"


def test_extend_enumerations(tmp_path):
    written_df = pd.DataFrame(
        {
            "soma_joinid": pd.Series([0, 1, 2, 3, 4, 5], dtype=np.int64),
            "str": pd.Series(["A", "B", "A", "B", "B", "B"], dtype="category"),
            "byte": pd.Series([b"A", b"B", b"A", b"B", b"B", b"B"], dtype="category"),
            "bool": pd.Series(
                [True, False, True, False, False, False], dtype="category"
            ),
            "int64": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.int64), dtype="category"
            ),
            "uint64": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.uint64), dtype="category"
            ),
            "int32": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.int32), dtype="category"
            ),
            "uint32": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.uint32), dtype="category"
            ),
            "int16": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.int16), dtype="category"
            ),
            "uint16": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.uint16), dtype="category"
            ),
            "int8": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.int8), dtype="category"
            ),
            "uint8": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.uint8), dtype="category"
            ),
            "float32": pd.Series(
                np.array([0, 1.1, 2.1, 0, 1.1, 2.1], dtype=np.float32), dtype="category"
            ),
            "float64": pd.Series(
                np.array([0, 1.1, 2.1, 0, 1.1, 2.1], dtype=np.float64), dtype="category"
            ),
            "float64_w_non_finite": pd.Series(
                np.array([0, 1.1, 2.1, 0, np.inf, -np.inf], dtype=np.float64),
                dtype="category",
            ),
            "str_ordered": pd.Series(
                pd.Categorical(
                    ["A", "B", "A", "B", "B", "B"],
                    categories=["B", "A", "C"],
                    ordered=True,
                ),
            ),
            "int64_ordered": pd.Series(
                pd.Categorical(
                    [1, 2, 3, 3, 2, 1],
                    categories=np.array([3, 2, 1], dtype=np.int64),
                    ordered=True,
                ),
            ),
            "uint64_ordered": pd.Series(
                pd.Categorical(
                    [1, 2, 3, 3, 2, 1],
                    categories=np.array([3, 2, 1], dtype=np.uint64),
                    ordered=True,
                ),
            ),
            "float64_ordered": pd.Series(
                pd.Categorical(
                    [0, 1.1, 2.1, 0, 1.1, 2.1],
                    categories=np.array([1.1, 0, 2.1], dtype=np.float64),
                    ordered=True,
                ),
            ),
        },
    )

    schema = pa.Schema.from_pandas(written_df, preserve_index=False)

    with soma.DataFrame.create(
        str(tmp_path), schema=schema, domain=[[0, 9]]
    ) as soma_dataframe:
        tbl = pa.Table.from_pandas(written_df, preserve_index=False)
        soma_dataframe.write(tbl)

    with soma.open(str(tmp_path)) as soma_dataframe:
        readback_df = soma_dataframe.read().concat().to_pandas()
        for c in readback_df:
            assert readback_df[c].dtype == written_df[c].dtype
            if readback_df[c].dtype == "category":
                assert (
                    readback_df[c].cat.categories.dtype
                    == written_df[c].cat.categories.dtype
                )
            assert (readback_df[c] == written_df[c]).all()


def test_multiple_writes_with_str_enums(tmp_path):
    uri = tmp_path.as_posix()

    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            (
                "obs",
                pa.dictionary(
                    index_type=pa.int8(), value_type=pa.string(), ordered=False
                ),
            ),
        ]
    )
    soma.DataFrame.create(uri, schema=schema, domain=[[0, 7]]).close()

    df1 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([0, 1, 2], dtype=np.int64),
            "obs": pd.Series(["A", "B", "A"], dtype="category"),
        }
    )
    tbl = pa.Table.from_pandas(df1, preserve_index=False)
    with soma.open(uri, mode="w") as A:
        A.write(tbl)

    df2 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([3, 4, 5], dtype=np.int64),
            "obs": pd.Series(["B", "C", "B"], dtype="category"),
        }
    )
    tbl = pa.Table.from_pandas(df2, preserve_index=False)
    with soma.open(uri, mode="w") as A:
        A.write(tbl)

    with soma.open(uri) as A:
        df = A.read().concat().to_pandas()

    # https://stackoverflow.com/questions/45639350/retaining-categorical-dtype-upon-dataframe-concatenation
    uc = union_categoricals([df1.obs, df2.obs])
    df1.obs = pd.Categorical(df1.obs, categories=uc.categories)
    df2.obs = pd.Categorical(df2.obs, categories=uc.categories)
    expected_df = pd.concat((df1, df2), ignore_index=True)

    assert df.equals(expected_df)

    # No new enumerations introduced but we still need to remap
    # the indexes
    df3 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([6, 7], dtype=np.int64),
            "obs": pd.Series(["C", "C"], dtype="category"),
        }
    )
    tbl = pa.Table.from_pandas(df3, preserve_index=False)
    with soma.open(uri, mode="w") as A:
        A.write(tbl)

    with soma.open(uri) as A:
        df = A.read().concat().to_pandas()

    uc = union_categoricals([df1.obs, df2.obs, df3.obs])
    df1.obs = pd.Categorical(df1.obs, categories=uc.categories)
    df2.obs = pd.Categorical(df2.obs, categories=uc.categories)
    df3.obs = pd.Categorical(df3.obs, categories=uc.categories)
    expected_df = pd.concat((df1, df2, df3), ignore_index=True)

    assert df.equals(expected_df)


def test_multiple_writes_with_int_enums(tmp_path):
    uri = tmp_path.as_posix()

    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            (
                "obs",
                pa.dictionary(
                    index_type=pa.int8(), value_type=pa.int64(), ordered=False
                ),
            ),
        ]
    )
    soma.DataFrame.create(uri, schema=schema, domain=[[0, 9]]).close()

    df1 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([0, 1, 2], dtype=np.int64),
            "obs": pd.Series([1, 2, 1], dtype="category"),
        }
    )
    tbl = pa.Table.from_pandas(df1, preserve_index=False)
    with soma.open(uri, mode="w") as A:
        A.write(tbl)

    df2 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([3, 4, 5], dtype=np.int64),
            "obs": pd.Series([2, 3, 2], dtype="category"),
        }
    )
    tbl = pa.Table.from_pandas(df2, preserve_index=False)
    with soma.open(uri, mode="w") as A:
        A.write(tbl)

    with soma.open(uri) as A:
        df = A.read().concat().to_pandas()

    # https://stackoverflow.com/questions/45639350/retaining-categorical-dtype-upon-dataframe-concatenation
    uc = union_categoricals([df1.obs, df2.obs])
    df1.obs = pd.Categorical(df1.obs, categories=uc.categories)
    df2.obs = pd.Categorical(df2.obs, categories=uc.categories)
    expected_df = pd.concat((df1, df2), ignore_index=True)

    assert df.equals(expected_df)

    # No new enumerations introduced but we still need to remap
    # the indexes
    df3 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([6, 7], dtype=np.int64),
            "obs": pd.Series([3, 3], dtype="category"),
        }
    )
    tbl = pa.Table.from_pandas(df3, preserve_index=False)
    with soma.open(uri, mode="w") as A:
        A.write(tbl)

    with soma.open(uri) as A:
        df = A.read().concat().to_pandas()

    uc = union_categoricals([df1.obs, df2.obs, df3.obs])
    df1.obs = pd.Categorical(df1.obs, categories=uc.categories)
    df2.obs = pd.Categorical(df2.obs, categories=uc.categories)
    df3.obs = pd.Categorical(df3.obs, categories=uc.categories)
    expected_df = pd.concat((df1, df2, df3), ignore_index=True)

    assert df.equals(expected_df)


def test_multichunk(tmp_path):
    uri = tmp_path.as_posix()

    # --- three dataframes, all with identical schema
    df_0 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([0, 1, 2, 3], dtype=np.int64),
            "obs": pd.Series(["A", "B", "A", "B"], dtype="str"),
        }
    )
    df_1 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([4, 5, 6, 7], dtype=np.int64),
            "obs": pd.Series(["A", "A", "B", "B"], dtype="str"),
        }
    )
    df_2 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([8, 9, 10, 11], dtype=np.int64),
            "obs": pd.Series(["B", "C", "B", "C"], dtype="str"),
        }
    )
    expected_df = pd.concat([df_0, df_1, df_2], ignore_index=True)

    soma.DataFrame.create(
        uri,
        schema=pa.Schema.from_pandas(df_0, preserve_index=False),
        domain=[[0, 11]],
    ).close()

    with soma.open(uri, mode="w") as A:
        # one-chunk table
        A.write(pa.Table.from_pandas(df_0, preserve_index=False))

    with soma.open(uri, mode="w") as A:
        # two-chunk table
        A.write(
            pa.concat_tables(
                [
                    pa.Table.from_pandas(df_1, preserve_index=False),
                    pa.Table.from_pandas(df_2, preserve_index=False),
                ]
            )
        )

    with soma.open(uri) as A:
        df = A.read().concat().to_pandas()

    assert df.equals(expected_df)


def test_multichunk_with_enums(tmp_path):
    uri = tmp_path.as_posix()

    # --- three dataframes, all with identical schema
    df_0 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([0, 1, 2, 3], dtype=np.int64),
            "obs": pd.Series(["A", "B", "A", "B"], dtype="category"),
        }
    )
    df_1 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([4, 5, 6, 7], dtype=np.int64),
            "obs": pd.Series(["A", "A", "B", "B"], dtype="category"),
        }
    )
    df_2 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([8, 9, 10, 11], dtype=np.int64),
            "obs": pd.Series(["B", "C", "B", "C"], dtype="category"),
        }
    )
    expected_df = pd.concat([df_0, df_1, df_2], ignore_index=True)

    soma.DataFrame.create(
        uri,
        schema=pa.Schema.from_pandas(df_0, preserve_index=False),
        domain=[[0, 11]],
    ).close()

    with soma.open(uri, mode="w") as A:
        # one-chunk table
        A.write(pa.Table.from_pandas(df_0, preserve_index=False))

    with soma.open(uri, mode="w") as A:
        # two-chunk table
        A.write(
            pa.concat_tables(
                [
                    pa.Table.from_pandas(df_1, preserve_index=False),
                    pa.Table.from_pandas(df_2, preserve_index=False),
                ]
            )
        )

    with soma.open(uri) as A:
        df = A.read().concat().to_pandas()

    # https://stackoverflow.com/questions/45639350/retaining-categorical-dtype-upon-dataframe-concatenation
    uc = union_categoricals([df_0.obs, df_1.obs, df_2.obs])
    df_0.obs = pd.Categorical(df_0.obs, categories=uc.categories)
    df_1.obs = pd.Categorical(df_1.obs, categories=uc.categories)
    df_2.obs = pd.Categorical(df_2.obs, categories=uc.categories)
    expected_df = pd.concat((df_0, df_1, df_2), ignore_index=True)

    assert df.equals(expected_df)


def test_enum_extend_past_numerical_limit(tmp_path):
    uri = tmp_path.as_posix()

    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            (
                "obs",
                pa.dictionary(
                    index_type=pa.int8(), value_type=pa.large_string(), ordered=False
                ),
            ),
        ]
    )
    soma.DataFrame.create(uri, schema=schema, domain=[[0, 999]]).close()

    n_elem = 132
    n_cats = 127
    df1 = pd.DataFrame(
        {
            "soma_joinid": pd.Series(np.arange(n_elem), dtype=np.int64),
            "obs": pd.Series(
                [f"enum_{i % n_cats}" for i in range(n_elem)], dtype="category"
            ),
        }
    )

    # use max number of possible categories
    tbl = pa.Table.from_pandas(df1, preserve_index=False)
    with soma.open(uri, mode="w") as A:
        A.write(tbl)

    more_elem = 4
    df2 = pd.DataFrame(
        {
            "soma_joinid": pd.Series(
                np.arange(n_elem, n_elem + more_elem), dtype=np.int64
            ),
            "obs": pd.Series(["TEST"] * more_elem, dtype="category"),
        }
    )

    # cannot add additional categories as already maxed out earlier
    tbl = pa.Table.from_pandas(df2, preserve_index=False)
    with pytest.raises(soma.SOMAError):
        with soma.open(uri, mode="w") as A:
            A.write(tbl)


def test_write_str_empty_ned(tmp_path):
    tmp_path.as_posix()


def test_enum_schema_report(tmp_path):
    uri = tmp_path.as_posix()

    pandas_df = pd.DataFrame(
        {
            "soma_joinid": pd.Series([0, 1, 2, 3, 4, 5], dtype=np.int64),
            "int_cat": pd.Series([10, 20, 10, 20, 20, 20], dtype="category"),
            "int": pd.Series([10, 20, 10, 20, 20, 20]),
            "str_cat": pd.Series(["A", "B", "A", "B", "B", "B"], dtype="category"),
            "str": pd.Series(["A", "B", "A", "B", "B", "B"]),
            "byte_cat": pd.Series(
                [b"A", b"B", b"A", b"B", b"B", b"B"], dtype="category"
            ),
            "byte": pd.Series([b"A", b"B", b"A", b"B", b"B", b"B"]),
        },
    )

    arrow_schema = pa.Schema.from_pandas(pandas_df, preserve_index=False)

    with soma.DataFrame.create(uri, schema=arrow_schema, domain=[[0, 5]]) as sdf:
        arrow_table = pa.Table.from_pandas(pandas_df, preserve_index=False)
        sdf.write(arrow_table)

    # Verify SOMA Arrow schema
    with soma.open(uri) as sdf:
        f = sdf.schema.field("int_cat")
        assert f.type.index_type == pa.int8()
        assert f.type.value_type == pa.int64()

        f = sdf.schema.field("str_cat")
        assert f.type.index_type == pa.int8()
        assert f.type.value_type == pa.string()

        f = sdf.schema.field("byte_cat")
        assert f.type.index_type == pa.int8()
        assert f.type.value_type == pa.binary()


def test_nullable(tmp_path):
    uri = tmp_path.as_posix()

    # Arrow fields are nullable by default.  They can be explicitly set nullable
    # or non-nullable via the nullable kwarg to pa.field.  Also, they can be
    # explicitly set nullable via metadata. The latter, if present, overrides
    # the former.
    asch = pa.schema(
        [
            pa.field("int", pa.int32()),
            pa.field("bool", pa.bool_()),
            pa.field("ord", pa.dictionary(pa.int64(), pa.string())),
            pa.field("no-meta-flag-unspecified", pa.int32()),
            pa.field("no-meta-flag-true", pa.int32(), nullable=True),
            pa.field("no-meta-flag-false", pa.int32(), nullable=False),
            pa.field("yes-meta-flag-unspecified", pa.int32()),
            pa.field("yes-meta-flag-true", pa.int32(), nullable=True),
            pa.field("yes-meta-flag-false", pa.int32(), nullable=False),
        ],
        metadata={
            "int": "nullable",
            "bool": "nullable",
            "ord": "nullable",
            "yes-meta-flag-unspecified": "nullable",
            "yes-meta-flag-true": "nullable",
            "yes-meta-flag-false": "nullable",
        },
    )

    pydict = {}
    pydict["soma_joinid"] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    pydict["int"] = [1, 2, 3, 4, 5, 6, None, 8, None, None]
    pydict["bool"] = [True, True, True, False, True, False, None, False, None, None]
    pydict["ord"] = pd.Categorical(
        ["g1", "g2", "g3", None, "g2", "g3", "g1", None, "g3", "g1"]
    )
    pydict["no-meta-flag-unspecified"] = [1, 2, 3, 4, 5, 6, None, 8, None, None]
    pydict["no-meta-flag-true"] = [1, 2, 3, 4, 5, 6, None, 8, None, None]
    pydict["no-meta-flag-false"] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    pydict["yes-meta-flag-unspecified"] = [1, 2, 3, 4, 5, 6, None, 8, None, None]
    pydict["yes-meta-flag-true"] = [1, 2, 3, 4, 5, 6, None, 8, None, None]
    pydict["yes-meta-flag-false"] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    data = pa.Table.from_pydict(pydict)

    with soma.DataFrame.create(uri, schema=asch, domain=[[0, 9]]) as sdf:
        sdf.write(data)

    with soma.DataFrame.open(uri, "r") as sdf:
        df = sdf.read().concat().to_pandas()
        assert df.compare(data.to_pandas()).empty


def test_only_evolve_schema_when_enmr_is_extended(tmp_path):
    uri = tmp_path.as_posix()

    schema = pa.schema(
        [
            pa.field("myint", pa.dictionary(pa.int64(), pa.large_string())),
            pa.field("myfloat", pa.large_string()),
        ]
    )

    # +1 creating the schema
    # +1 evolving the schema
    with soma.DataFrame.create(uri, schema=schema, domain=[[0, 4]]) as sdf:
        data = {}
        data["soma_joinid"] = [0, 1, 2, 3, 4]
        data["myint"] = pd.Categorical(["a", "bb", "ccc", "bb", "a"])
        data["myfloat"] = ["cat", "dog", "cat", "cat", "cat"]
        sdf.write(pa.Table.from_pydict(data))

    # +1 evolving the schema
    with soma.DataFrame.open(uri, "w") as sdf:
        data = {}
        data["soma_joinid"] = [0, 1, 2, 3, 4]
        data["myint"] = pd.Categorical(["a", "bb", "ccc", "d", "a"])
        data["myfloat"] = ["cat", "dog", "cat", "cat", "cat"]
        sdf.write(pa.Table.from_pydict(data))

    # +0 no changes to enumeration values
    with soma.DataFrame.open(uri, "w") as sdf:
        data = {}
        data["soma_joinid"] = [0, 1, 2, 3, 4]
        data["myint"] = pd.Categorical(["a", "bb", "ccc", "d", "a"])
        data["myfloat"] = ["cat", "dog", "cat", "cat", "cat"]
        sdf.write(pa.Table.from_pydict(data))

    # +0 no changes enumeration values
    with soma.DataFrame.open(uri, "w") as sdf:
        data = {}
        data["soma_joinid"] = [0, 1, 2, 3, 4]
        data["myint"] = pd.Categorical(["a", "bb", "ccc", "d", "d"])
        data["myfloat"] = ["cat", "dog", "cat", "cat", "cat"]
        sdf.write(pa.Table.from_pydict(data))

    # total 3 fragment files

    # subtract 1 for the __schema/__enumerations directory;
    # only looking at fragment files
    assert len(list((Path(uri) / "__schema").iterdir())) - 1 == 3


def test_fix_update_dataframe_with_var_strings(tmp_path):
    uri = tmp_path.as_posix()

    tbl = pa.table(
        {
            "soma_joinid": pa.array([0, 1, 2, 3], pa.int64()),
            "mystring": pa.array(["a", "bb", "ccc", "dddd"], pa.large_utf8()),
            "myint": pa.array([33, 44, 55, 66], pa.int32()),
            "myfloat": pa.array([4.5, 5.5, 6.5, 7.5], pa.float32()),
        }
    )

    with soma.DataFrame.create(uri, schema=tbl.schema, domain=[[0, 3]]) as sdf:
        sdf.write(tbl)

    with soma.DataFrame.open(uri, "r") as sdf:
        updated_sdf = sdf.read().concat().to_pandas()
    updated_sdf["newattr"] = np.array(["a", "b", "c", "d"])

    with soma.DataFrame.open(uri, "w") as sdf:
        soma.io.ingest._update_dataframe(
            sdf,
            updated_sdf,
            "testing",
            platform_config=None,
            context=None,
            default_index_name="mystring",
        )

    with soma.DataFrame.open(uri, "r") as sdf:
        results = sdf.read().concat().to_pandas()
        assert results.equals(updated_sdf)


def test_presence_matrix(tmp_path):
    uri = tmp_path.as_uri()

    # Cerate the dataframe
    soma_df = soma.DataFrame.create(
        uri,
        schema=pa.schema(
            [
                ("soma_joinid", pa.int64()),
                ("scene_id", pa.string()),
                ("data", pa.bool_()),
            ]
        ),
        domain=((0, 99), ("", "")),
        index_column_names=("soma_joinid", "scene_id"),
    )

    # Create datda to write
    joinid_data = pa.array(np.arange(0, 100, 5))
    scene_id_data = 10 * ["scene1"] + 10 * ["scene2"]
    df = pd.DataFrame(
        {
            "soma_joinid": joinid_data,
            "scene_id": scene_id_data,
            "data": 20 * [True],
        }
    )
    arrow_table = pa.Table.from_pandas(df)
    soma_df.write(arrow_table)

    soma_df.close()

    with soma.DataFrame.open(uri) as soma_df:
        actual = soma_df.read().concat().to_pandas()

    assert actual.equals(df)


def test_bounds_on_somajoinid_domain(tmp_path):
    uri = tmp_path.as_posix()

    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            ("mystring", pa.string()),
            ("myint", pa.int32()),
            ("myfloat", pa.float32()),
        ]
    )

    with pytest.raises(ValueError):
        soma.DataFrame.create(
            uri,
            schema=schema,
            domain=[[0, -1]],
        )

    with pytest.raises(ValueError):
        soma.DataFrame.create(
            uri,
            schema=schema,
            domain=[[-1, 2]],
        )

    soma.DataFrame.create(
        uri,
        schema=schema,
        domain=[[2, 99]],
    )

    assert soma.DataFrame.exists(uri)


def test_pass_configs(tmp_path, arrow_schema):
    uri = tmp_path.as_posix()

    with soma.DataFrame.create(uri, schema=arrow_schema()) as sdf:
        pydict = {}
        pydict["soma_joinid"] = [0, 1, 2, 3, 4]
        pydict["myint"] = [10, 20, 30, 40, 50]
        pydict["myfloat"] = [4.1, 5.2, 6.3, 7.4, 8.5]
        pydict["mystring"] = ["apple", "ball", "cat", "dog", "egg"]
        pydict["mybool"] = [True, False, False, True, False]
        rb = pa.Table.from_pydict(pydict)
        sdf.tiledbsoma_resize_soma_joinid_shape(len(rb))
        sdf.write(rb)

    # Pass a custom config to open
    with soma.DataFrame.open(
        uri,
        "r",
        context=soma.SOMATileDBContext(
            {"sm.mem.total_budget": "0", "sm.io_concurrency_level": "0"}
        ),
    ) as sdf:

        # This errors out as 0 is not a valid value to set the total memory
        # budget or number of threads
        with pytest.raises(soma.SOMAError):
            next(sdf.read())

        # This still errors out because read still sees that the number of
        # threads is 0 and therefore invalid
        with pytest.raises(soma.SOMAError):
            next(sdf.read(platform_config={"sm.mem.total_budget": "10000"}))

        # With correct values, this reads without issue
        next(
            sdf.read(
                platform_config={
                    "sm.mem.total_budget": "10000",
                    "sm.io_concurrency_level": "1",
                }
            )
        )


def test_arrow_table_sliced_writer(tmp_path):
    """Tests writes of sliced Arrow tables, with fixed-length and variable-length attributes"""
    uri = tmp_path.as_posix()
    num_rows = 50

    schema = pa.schema(
        [
            ("myint", pa.int32()),
            ("mystring", pa.large_string()),
            ("mybool", pa.bool_()),
            ("myenumint", pa.dictionary(pa.int64(), pa.int32())),
            ("myenumstr", pa.dictionary(pa.int64(), pa.large_string())),
            ("myenumbool", pa.dictionary(pa.int64(), pa.bool_())),
        ]
    )

    pydict = {
        "soma_joinid": list(range(num_rows)),
        "myint": np.random.randint(10, 100, size=num_rows),
        "mystring": [f"s_{np.random.randint(1, 100000):08d}" for _ in range(num_rows)],
        "mybool": np.random.choice([False, True], size=num_rows),
        "myenumint": pd.Categorical(
            np.random.choice([1, 2, 3], size=num_rows, replace=True)
        ),
        "myenumstr": pd.Categorical(
            np.random.choice(["a", "bb", "ccc"], size=num_rows, replace=True)
        ),
        "myenumbool": pd.Categorical(
            np.random.choice([False, True], size=num_rows, replace=True)
        ),
    }

    pydict["myenumint"] = pa.DictionaryArray.from_arrays(
        pa.array(pydict["myenumint"].codes, type=pa.int32()),
        pa.array([1, 2, 3], type=pa.int32()),
    )

    pydict["myenumstr"] = pa.DictionaryArray.from_arrays(
        pa.array(pydict["myenumstr"].codes, type=pa.int32()),
        pa.array(["a", "bb", "ccc"], type=pa.large_string()),
    )

    pydict["myenumbool"] = pa.DictionaryArray.from_arrays(
        pa.array(pydict["myenumbool"].codes, type=pa.int32()),
        pa.array([False, True], type=pa.bool_()),
    )

    table = pa.Table.from_pydict(pydict)

    domain = [[0, len(table) - 1]]

    with soma.DataFrame.create(uri, schema=schema, domain=domain) as sdf:
        sdf.write(table[:])

    with soma.DataFrame.open(uri) as sdf:
        pdf = sdf.read().concat()

        assert_array_equal(pdf["myint"], pydict["myint"])
        assert_array_equal(pdf["mystring"], pydict["mystring"])
        assert_array_equal(pdf["mybool"], pydict["mybool"])

        assert_array_equal(pdf["myenumint"], pydict["myenumint"])
        assert_array_equal(pdf["myenumstr"], pydict["myenumstr"])
        assert_array_equal(pdf["myenumbool"], pydict["myenumbool"])

    with soma.DataFrame.open(uri, mode="w") as sdf:
        mid = num_rows // 2
        sdf.write(table[:mid])
        sdf.write(table[mid:])

    with soma.DataFrame.open(uri) as sdf:
        pdf = sdf.read().concat()

        assert_array_equal(pdf["myint"], pydict["myint"])
        assert_array_equal(pdf["mystring"], pydict["mystring"])
        assert_array_equal(pdf["mybool"], pydict["mybool"])

        assert_array_equal(pdf["myenumint"], pydict["myenumint"])
        assert_array_equal(pdf["myenumstr"], pydict["myenumstr"])
        assert_array_equal(pdf["myenumbool"], pydict["myenumbool"])


def test_arrow_table_validity_with_slicing(tmp_path):
    uri = tmp_path.as_posix()
    num_rows = 10
    domain = ((0, np.iinfo(np.int64).max - 2050),)

    schema = pa.schema(
        [
            ("myint", pa.int32()),
            ("mystring", pa.large_string()),
            ("mybool", pa.bool_()),
            ("mydatetime", pa.timestamp("s")),
            ("myenum", pa.dictionary(pa.int64(), pa.large_string())),
        ]
    )

    soma.DataFrame.create(uri, schema=schema, domain=domain)

    pydict = {}
    pydict["soma_joinid"] = [None, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    pydict["myint"] = [1, 2, 3, 4, 5, 6, None, 8, None, None]
    pydict["mystring"] = ["g1", "g2", "g3", None, "g2", "g3", "g1", None, "g3", "g1"]
    pydict["mybool"] = [True, True, True, False, True, False, None, False, None, None]
    pydict["mydatetime"] = [
        np.datetime64("NaT", "s"),
        np.datetime64(1, "s"),
        np.datetime64(2, "s"),
        np.datetime64("NaT", "s"),
        np.datetime64(4, "s"),
        np.datetime64(5, "s"),
        np.datetime64(6, "s"),
        np.datetime64(7, "s"),
        np.datetime64("NaT", "s"),
        np.datetime64(9, "s"),
    ]
    pydict["myenum"] = pd.Categorical(
        ["g1", "g2", "g3", None, "g2", "g3", "g1", None, "g3", "g1"]
    )
    table = pa.Table.from_pydict(pydict)

    # As of version 1.15.6 we were throwing in this case. However, we found
    # a compatibility issue with pyarrow versions below 17. Thus this is
    # now non-fatal.
    # with soma.DataFrame.open(uri, "w") as A:
    #    with raises_no_typeguard(soma.SOMAError):
    #        # soma_joinid cannot be nullable
    #        A.write(table)

    pydict["soma_joinid"] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    table = pa.Table.from_pydict(pydict)

    with soma.DataFrame.open(uri, "w") as A:
        A.write(table)

    with soma.DataFrame.open(uri) as A:
        pdf = A.read().concat()
        assert_array_equal(pdf["myint"], table["myint"])
        assert_array_equal(pdf["mystring"], table["mystring"])
        assert_array_equal(pdf["mybool"], table["mybool"])
        assert_array_equal(pdf["mydatetime"], table["mydatetime"])
        assert_array_equal(pdf["myenum"], table["myenum"])

    with soma.DataFrame.open(uri, "w") as A:
        mid = num_rows // 2
        A.write(table[:mid])
        A.write(table[mid:])

    with soma.DataFrame.open(uri) as A:
        pdf = A.read().concat()
        assert_array_equal(pdf["myint"], table["myint"])
        assert_array_equal(pdf["mystring"], table["mystring"])
        assert_array_equal(pdf["mybool"], table["mybool"])
        assert_array_equal(pdf["mydatetime"], table["mydatetime"])
        assert_array_equal(pdf["myenum"], table["myenum"])


def test_enum_regression_62887(tmp_path):
    uri = tmp_path.as_posix()

    schema = pa.schema(
        [
            pa.field("soma_joinid", pa.int64(), nullable=False),
            pa.field("A", pa.dictionary(pa.int8(), pa.int8())),
        ]
    )

    tbl = pa.Table.from_pydict(
        {
            "soma_joinid": pa.chunked_array([[0, 1, 2, 3, 4, 5, 6, 7], [8, 9]]),
            "A": pa.chunked_array(
                [
                    pa.DictionaryArray.from_arrays(
                        indices=pa.array([0, 0, 0, 0, 0, 0, 0, 0], type=pa.int8()),
                        dictionary=pa.array(
                            [0, 1, 2, 3, 4, 5, 6, 7, 8], type=pa.int8()
                        ),
                    ),
                    pa.DictionaryArray.from_arrays(
                        indices=pa.array([0, 0], type=pa.int8()),
                        dictionary=pa.array(
                            [0, 1, 2, 3, 4, 5, 6, 7, 8], type=pa.int8()
                        ),
                    ),
                ]
            ),
        }
    )

    with soma.DataFrame.create(
        uri, schema=schema, index_column_names=["soma_joinid"], domain=[(0, 10000000)]
    ) as A:
        A.write(tbl)

    with soma.open(uri) as A:
        assert_array_equal(A.read().concat()["A"], tbl["A"])


def test_enum_handling_category_of_nan_62449(tmp_path):
    uri = tmp_path.as_posix()

    # Different representations of single-precision NaNs
    quiet_nan = struct.unpack(">f", b"\x7f\xc0\x00\x00")[0]
    negative_nan = struct.unpack(">f", b"\xff\xc0\x00\x00")[0]
    signaling_nan = struct.unpack(">f", b"\x7f\x80\x00\x01")[0]

    def nan_check(expected_nan, dict_vals):
        return any(
            math.isnan(val)
            and struct.pack(">f", val) == struct.pack(">f", expected_nan)
            for val in dict_vals
        )

    schema = pa.schema(
        [
            pa.field("soma_joinid", pa.int64(), nullable=False),
            pa.field("A", pa.dictionary(pa.int32(), pa.float32())),
        ]
    )

    # Ensure that unique NaN values are respected as different dictionary values
    expected_data1 = pa.Table.from_pydict(
        {
            "soma_joinid": [0, 1, 2, 3],
            "A": pa.DictionaryArray.from_arrays(
                indices=pa.array([0, 1, 2, 0], type=pa.int32()),
                dictionary=pa.array(
                    [negative_nan, quiet_nan, signaling_nan], type=pa.float32()
                ),
            ),
        }
    )

    with soma.DataFrame.create(
        uri, schema=schema, index_column_names=["soma_joinid"], domain=[(0, 5)]
    ) as A:
        A.write(expected_data1)

    with soma.open(uri) as A:
        actual_data = A.read().concat()["A"]
        actual_dict_vals = actual_data.chunk(0).dictionary.to_pylist()
        assert len(actual_dict_vals) == 3
        assert nan_check(quiet_nan, actual_dict_vals)
        assert nan_check(negative_nan, actual_dict_vals)
        assert nan_check(signaling_nan, actual_dict_vals)
        assert_array_equal(expected_data1["A"], actual_data)

    # Ensure that the dictionary indexes get shifted correctly when appending
    # to the dataframe
    expected_data2 = pa.Table.from_pydict(
        {
            "soma_joinid": [4, 5],
            "A": pa.DictionaryArray.from_arrays(
                indices=pa.array([1, 0], type=pa.int32()),
                dictionary=pa.array([quiet_nan, signaling_nan], type=pa.float32()),
            ),
        }
    )

    with soma.open(uri, mode="w") as A:
        A.write(expected_data2)

    with soma.open(uri) as A:
        actual_data = A.read().concat()["A"]
        actual_dict_vals = actual_data.chunk(0).dictionary.to_pylist()
        assert len(actual_dict_vals) == 3
        assert nan_check(quiet_nan, actual_dict_vals)
        assert nan_check(negative_nan, actual_dict_vals)
        assert nan_check(signaling_nan, actual_dict_vals)
        assert_array_equal(expected_data1["A"], actual_data[:4])
        assert_array_equal(expected_data2["A"], actual_data[4:])


def test_return_datetime_type_for_domain_and_maxdomain_62887(tmp_path):
    uri = tmp_path.as_posix()

    schema = pa.schema([("A", pa.timestamp("s"))])
    index_column_names = ("soma_joinid", "A")
    domain = (
        (0, 1000),
        (np.datetime64("2000-01-01T00:00:00"), np.datetime64("2025-01-04T09:16:49")),
    )

    soma.DataFrame.create(
        uri, schema=schema, index_column_names=index_column_names, domain=domain
    )

    with soma.DataFrame.open(uri) as A:
        assert A.index_column_names == index_column_names
        assert A.schema.field("A").type == pa.timestamp("s")

        assert A.domain[1] == (
            pa.scalar(domain[1][0], pa.timestamp("s")),
            pa.scalar(domain[1][1], pa.timestamp("s")),
        )
        assert A.maxdomain[1] == (
            pa.scalar(-9223372036854775807, pa.timestamp("s")),
            pa.scalar(9223372036853775807, pa.timestamp("s")),
        )


@pytest.mark.parametrize(
    "pa_type,tile",
    (
        (pa.int32(), "1"),
        (pa.float32(), "1"),
        (pa.timestamp("s"), "1"),
        (pa.large_string(), ""),
        (pa.large_binary(), ""),
    ),
)
def test_extents(tmp_path, pa_type, tile):
    uri = tmp_path.as_posix()
    asch = pa.schema([pa.field("dim", pa_type)])
    soma.DataFrame.create(uri, schema=asch, index_column_names=["dim"])

    with soma.DataFrame.open(tmp_path.as_posix()) as A:
        dim_info = json.loads(A.schema_config_options().dims)
        assert dim_info["dim"]["tile"] == tile
