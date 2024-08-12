import json
from copy import deepcopy
from dataclasses import asdict, dataclass, fields
from inspect import getfullargspec
from os.path import join
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from anndata import AnnData
from pandas._testing import assert_frame_equal
from scipy.sparse import csr_matrix

from tiledbsoma import SOMA_JOINID, DataFrame, Experiment
from tiledbsoma.io._common import _DATAFRAME_ORIGINAL_INDEX_NAME_JSON
from tiledbsoma.io._registration import AxisIDMapping
from tiledbsoma.io.ingest import IngestionParams, _write_dataframe, from_anndata
from tiledbsoma.io.outgest import _read_dataframe, to_anndata

from tests._util import assert_adata_equal


def parse_col(col_str: str) -> Tuple[Optional[str], List[str]]:
    """Parse a "column string" of the form `val1,val2,...` or `name=val1,val2,...`."""
    pcs = col_str.split("=")
    if len(pcs) == 1:
        return None, col_str.split(",")
    elif len(pcs) == 2:
        name, vals_str = pcs
        vals = vals_str.split(",")
        return name, vals
    else:
        raise ValueError(f"Invalid column string: {col_str}")


def make_df(index_str: str, **cols) -> pd.DataFrame:
    """DataFrame construction helper, for tests.

    - index and columns are provided as strings of the form `name=val1,val2,...`.
    - `name=` is optional for the initial (`index_str`) arg.
    """
    cols = dict([(col, parse_col(col_str)[1]) for col, col_str in cols.items()])
    index = None
    index_name = None
    if index_str:
        index_name, index = parse_col(index_str)
    df = pd.DataFrame(cols, index=index)
    df.index.name = index_name
    return df


@dataclass
class RoundTrip:
    # Test-case name
    name: str
    # DataFrame to ingest
    original_df: pd.DataFrame
    # Expected DataFrame after outgest
    outgested_df: pd.DataFrame
    # Columns that should be persisted; all are expected to be of type "large string" (for the purposes of these test
    # cases); the required `soma_joinid` (with type `int64 not null`) is excluded from this list, but always verified to
    # also exist.
    persisted_column_names: List[str]
    # Expected "original index metadata" (attached to the persisted SOMA DataFrame)
    persisted_metadata: Optional[str] = None
    # Argument passed to `_write_dataframe` on ingest (default here matches `from_anndata`'s "obs" path)
    ingest_id_column_name: Optional[str] = "obs_id"
    # Argument passed to `_read_dataframe` on outgest (default here matches `to_anndata`'s "obs_id" path)
    outgest_default_index_name: Optional[str] = None
    # Argument passed to `_read_dataframe` on outgest (default here matches `to_anndata`'s "obs_id" path)
    outgest_fallback_index_name: Optional[str] = "obs_id"


def parametrize_roundtrips(roundtrips: List[RoundTrip]):
    def wrapper(fn):
        # Test-case IDs
        ids = [rt.name for rt in roundtrips]
        # Convert `RoundTrip`s to "values" arrays, filtered and reordered to match kwargs expected by the wrapped
        # function
        fields_names = [f.name for f in fields(RoundTrip)]
        spec = getfullargspec(fn)
        names = [arg for arg in spec.args if arg in fields_names]
        values = [
            {name: rt_dict[name] for name in names}.values()
            for rt_dict in [asdict(rt) for rt in roundtrips]
        ]
        # Delegate to PyTest `parametrize`
        return pytest.mark.parametrize(
            names,  # arg names
            values,  # arg value lists
            ids=ids,  # test-case names
        )(fn)

    return wrapper


# fmt: off
ROUND_TRIPS = [
    RoundTrip(
        '1. `df.index` named "index"',
        make_df("index=xx,yy,zz", col0="aa,bb,cc", col1="AA,BB,CC"),
        # ⇒ index name is lost
        make_df(      "xx,yy,zz", col0="aa,bb,cc", col1="AA,BB,CC"),
        [ "obs_id", "col0", "col1", ],
    ),
    RoundTrip(
        '2. DataFrame has a column named `obs_id`',
        make_df("xx,yy,zz", col0="AA,BB,CC", obs_id="aa,bb,cc"),
        # ⇒ `obs_id` column becomes `df.index`, loses name
        # ⇒ Original `df.index` dropped
        make_df("aa,bb,cc", col0="AA,BB,CC"),
        [ "col0", "obs_id", ],
    ),
    RoundTrip(
        '3. DataFrame has a column named "index", and `df.index` is unnamed',
        make_df("xx,yy,zz", col0="aa,bb,cc", index="AA,BB,CC"),
        # ⇒ "index" column is promoted to `df.index` (unnamed)
        # ⇒ Original (unnamed) index becomes a column named `level_0`
        make_df("AA,BB,CC", level_0="xx,yy,zz", col0="aa,bb,cc"),
        [ "level_0", "col0", "obs_id", ],
    ),
    RoundTrip(
        '4. DataFrame has a column named "index", and `df.index` is named `id_column_name` (default `obs_id`)',
        make_df("obs_id=xx,yy,zz", col0="aa,bb,cc", index="AA,BB,CC"),
        # ⇒ "index" column is dropped
        make_df("obs_id=xx,yy,zz", col0="aa,bb,cc"),
        [ "obs_id", "col0", ],
        "obs_id",
    ),
    RoundTrip(
        '5. DataFrame has a column named "index" and df.index has another name',
        make_df("idx=xx,yy,zz", col0="aa,bb,cc", index="AA,BB,CC"),
        # ⇒ "index" column renamed to `obs_id`
        make_df("idx=xx,yy,zz", col0="aa,bb,cc", obs_id="AA,BB,CC"),
        [ "idx", "col0", "obs_id", ],
        "idx",
    ),
    RoundTrip(
        '6. DataFrame has columns named "index" and `id_column_name` (default: obs_id), and `df.index` is unnamed:',
        make_df("xx,yy,zz", obs_id="aa,bb,cc", index="AA,BB,CC"),
        # ⇒ unnamed index → column named `level_0`
        # ⇒ `obs_id` column → unnamed index
        # ⇒ "index" column dropped
        make_df("aa,bb,cc", level_0="xx,yy,zz"),
        [ "level_0", "obs_id", ],
    ),
    RoundTrip(
        '7. DataFrame has columns named "index" and `id_column_name` (default: obs_id), and `df.index` has a name:',
        make_df("idx=xx,yy,zz", obs_id="aa,bb,cc", index="AA,BB,CC"),
        # ⇒ "index" column is dropped
        make_df("idx=xx,yy,zz", obs_id="aa,bb,cc"),
        [ "idx", "obs_id", ],
        "idx",
    ),
]
# fmt: on


def verify_metadata(
    sdf: DataFrame, persisted_column_names: List[str], persisted_metadata: Optional[str]
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


@parametrize_roundtrips(ROUND_TRIPS)
def test_adata_io_roundtrips(
    tmp_path: Path,
    original_df: pd.DataFrame,
    persisted_column_names: List[str],
    persisted_metadata: Optional[str],
    ingest_id_column_name: Optional[str],
    outgest_default_index_name: Optional[str],
    outgested_df: pd.DataFrame,
):
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
        adata1 = to_anndata(exp, "meas", obs_id_name=outgest_default_index_name)
        outgested_obs = adata1.obs

    assert_frame_equal(outgested_obs, outgested_df)

    expected = deepcopy(adata0)
    # Patch in the expected outgested DataFrame (which in these test cases is known to differ from
    # what was ingested).
    expected.obs = outgested_df
    assert_adata_equal(expected, adata1)


@parametrize_roundtrips(ROUND_TRIPS)
def test_df_io_roundtrips(
    tmp_path: Path,
    original_df: pd.DataFrame,
    persisted_column_names: List[str],
    persisted_metadata: Optional[str],
    ingest_id_column_name: Optional[str],
    outgest_default_index_name: Optional[str],
    outgest_fallback_index_name: Optional[str],
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
        default_index_name=outgest_default_index_name,
        fallback_index_name=outgest_fallback_index_name,
    )
    assert_frame_equal(actual_outgested_df, outgested_df)
