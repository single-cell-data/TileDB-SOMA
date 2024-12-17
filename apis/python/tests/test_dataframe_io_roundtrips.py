# The tests in this file verify issues where an ingest/outgest "round trip" modifies an AnnData's
# "obs" or "var" DataFrames. See https://github.com/single-cell-data/TileDB-SOMA/issues/2829 for more info.

from __future__ import annotations

import json
from copy import deepcopy
from dataclasses import dataclass
from os.path import join
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import pyarrow as pa
from anndata import AnnData
from pandas._testing import assert_frame_equal
from scipy.sparse import csr_matrix

from tiledbsoma import SOMA_JOINID, DataFrame, Experiment
from tiledbsoma.io._common import _DATAFRAME_ORIGINAL_INDEX_NAME_JSON
from tiledbsoma.io._registration import AxisIDMapping
from tiledbsoma.io.ingest import IngestionParams, _write_dataframe, from_anndata
from tiledbsoma.io.outgest import _read_dataframe, to_anndata

from tests._util import assert_adata_equal, make_pd_df
from tests.parametrize_cases import parametrize_cases


@dataclass
class RoundTrip:
    # Test-case name
    id: str
    # DataFrame to ingest
    original_df: pd.DataFrame
    # Expected DataFrame after outgest
    outgested_df: pd.DataFrame
    # Columns that should be persisted; all are expected to be of type "large string" (for the purposes of these test
    # cases); the required `soma_joinid` (with type `int64 not null`) is excluded from this list, but always verified to
    # also exist.
    persisted_column_names: List[str]
    # Expected "original index metadata" (attached to the persisted SOMA DataFrame)
    persisted_metadata: str | None = None
    # Argument passed to `_write_dataframe` on ingest (default here matches `from_anndata`'s "obs" path)
    ingest_id_column_name: str | None = "obs_id"


# fmt: off
# These cases verify issues with ingest/outgest where an AnnData's "obs" or "var" DataFrame is not round-tripped
# correctly. See https://github.com/single-cell-data/TileDB-SOMA/issues/2829 for more info.
ROUND_TRIPS = [
    RoundTrip(
        '1. `df.index` named "index"',
        make_pd_df("index=xx,yy,zz", col0="aa,bb,cc", col1="AA,BB,CC"),
        # ⇒ index name is lost
        make_pd_df(      "xx,yy,zz", col0="aa,bb,cc", col1="AA,BB,CC"),
        [ "obs_id", "col0", "col1", ],
    ),
    RoundTrip(
        '2. DataFrame has a column named `obs_id`',
        make_pd_df("xx,yy,zz", col0="AA,BB,CC", obs_id="aa,bb,cc"),
        # ⇒ `obs_id` column becomes `df.index`, loses name
        # ⇒ Original `df.index` dropped
        make_pd_df("aa,bb,cc", col0="AA,BB,CC"),
        [ "col0", "obs_id", ],
    ),
    RoundTrip(
        '3. DataFrame has a column named "index", and `df.index` is unnamed',
        make_pd_df("xx,yy,zz", col0="aa,bb,cc", index="AA,BB,CC"),
        # ⇒ "index" column is promoted to `df.index` (unnamed)
        # ⇒ Original (unnamed) index becomes a column named `level_0`
        make_pd_df("AA,BB,CC", level_0="xx,yy,zz", col0="aa,bb,cc"),
        [ "level_0", "col0", "obs_id", ],
    ),
    RoundTrip(
        '4. DataFrame has a column named "index", and `df.index` is named `id_column_name` (default `obs_id`)',
        make_pd_df("obs_id=xx,yy,zz", col0="aa,bb,cc", index="AA,BB,CC"),
        # ⇒ "index" column is dropped
        make_pd_df("obs_id=xx,yy,zz", col0="aa,bb,cc"),
        [ "obs_id", "col0", ],
        "obs_id",
    ),
    RoundTrip(
        '5. DataFrame has a column named "index" and df.index has another name',
        make_pd_df("idx=xx,yy,zz", col0="aa,bb,cc", index="AA,BB,CC"),
        # ⇒ "index" column renamed to `obs_id`
        make_pd_df("idx=xx,yy,zz", col0="aa,bb,cc", obs_id="AA,BB,CC"),
        [ "idx", "col0", "obs_id", ],
        "idx",
    ),
    RoundTrip(
        '6. DataFrame has columns named "index" and `id_column_name` (default: obs_id), and `df.index` is unnamed:',
        make_pd_df("xx,yy,zz", obs_id="aa,bb,cc", index="AA,BB,CC"),
        # ⇒ unnamed index → column named `level_0`
        # ⇒ `obs_id` column → unnamed index
        # ⇒ "index" column dropped
        make_pd_df("aa,bb,cc", level_0="xx,yy,zz"),
        [ "level_0", "obs_id", ],
    ),
    RoundTrip(
        '7. DataFrame has columns named "index" and `id_column_name` (default: obs_id), and `df.index` has a name:',
        make_pd_df("idx=xx,yy,zz", obs_id="aa,bb,cc", index="AA,BB,CC"),
        # ⇒ "index" column is dropped
        make_pd_df("idx=xx,yy,zz", obs_id="aa,bb,cc"),
        [ "idx", "obs_id", ],
        "idx",
    ),
]
# fmt: on


def verify_metadata(
    sdf: DataFrame, persisted_column_names: List[str], persisted_metadata: str | None
):
    # Verify column names and types
    schema = sdf.schema
    assert schema.names == [SOMA_JOINID, *persisted_column_names]
    [soma_joinid_type, *string_col_types] = schema.types
    assert soma_joinid_type == pa.int64() and schema.field(0).nullable is False
    for string_col_type in string_col_types:
        assert string_col_type == pa.large_string()

    # Verify "original index metadata"
    actual_index_metadata = json.loads(
        sdf.metadata[_DATAFRAME_ORIGINAL_INDEX_NAME_JSON]
    )
    assert actual_index_metadata == persisted_metadata


@parametrize_cases(ROUND_TRIPS)
def test_adata_io_roundtrips(
    tmp_path: Path,
    original_df: pd.DataFrame,
    persisted_column_names: List[str],
    persisted_metadata: str | None,
    ingest_id_column_name: str | None,
    outgested_df: pd.DataFrame,
):
    """Given an `original_df`, set it as the `obs` DataFrame of an AnnData, ingest it, outgest it back, and compare the
    original and final DataFrames. Also verify the persisted column names and "original index metadata."

    `ingest_id_column_name` and `outgest_default_index_name`
    """
    uri = str(tmp_path)
    n_obs = len(original_df)
    var = pd.DataFrame({"var1": [1, 2, 3], "var2": ["a", "b", "c"]})  # unused
    n_var = len(var)
    X = csr_matrix(np.array([0] * n_obs * n_var).reshape(n_obs, n_var))  # unused
    adata0 = AnnData(X=X, obs=original_df, var=var)
    ingested_uri = from_anndata(uri, adata0, "meas", obs_id_name=ingest_id_column_name)
    assert ingested_uri == uri

    # Verify column names and types
    obs_uri = join(uri, "obs")
    obs = DataFrame.open(obs_uri)
    verify_metadata(obs, persisted_column_names, persisted_metadata)

    # Verify outgested pd.DataFrame
    with Experiment.open(ingested_uri) as exp:
        adata1 = to_anndata(exp, "meas")
        outgested_obs = adata1.obs

    assert_frame_equal(outgested_obs, outgested_df)

    expected = deepcopy(adata0)
    # Patch in the expected outgested DataFrame (which in these test cases is known to differ from
    # what was ingested).
    expected.obs = outgested_df
    assert_adata_equal(expected, adata1)


@parametrize_cases(ROUND_TRIPS)
def test_df_io_roundtrips(
    tmp_path: Path,
    original_df: pd.DataFrame,
    persisted_column_names: List[str],
    persisted_metadata: str | None,
    ingest_id_column_name: str | None,
    outgested_df: pd.DataFrame,
):
    uri = str(tmp_path)
    _write_dataframe(
        uri,
        original_df,
        id_column_name=ingest_id_column_name,
        axis_mapping=AxisIDMapping(data=tuple(range(len(original_df)))),
        ingestion_params=IngestionParams("write", None),
    ).close()

    sdf = DataFrame.open(uri)
    verify_metadata(sdf, persisted_column_names, persisted_metadata)

    # Verify outgested pd.DataFrame
    actual_outgested_df = _read_dataframe(
        sdf,
        default_index_name=None,  # corresponds to `to_anndata`'s default `obs_id_name`
        fallback_index_name="obs_id",  # corresponds to `to_anndata`'s default behavior for "obs"
    )
    assert_frame_equal(actual_outgested_df, outgested_df)
