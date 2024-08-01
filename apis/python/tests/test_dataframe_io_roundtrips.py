import json
from dataclasses import asdict, dataclass, fields
from os.path import join
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from anndata import AnnData
from pandas._testing import assert_frame_equal

from tiledbsoma import SOMA_JOINID, DataFrame, Experiment
from tiledbsoma.io._common import _DATAFRAME_ORIGINAL_INDEX_NAME_JSON
from tiledbsoma.io.ingest import from_anndata
from tiledbsoma.io.outgest import to_anndata


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
    # Argument passed to `_extract_pdf` on outgest (default here matches `to_anndata`'s "obs_id" path)
    outgest_default_index_name: Optional[str] = None


def parametrize_roundtrips(*roundtrips: RoundTrip):
    def wrapper(fn):
        names = [f.name for f in fields(RoundTrip)[1:]]
        values = [list(asdict(rt).values()) for rt in roundtrips]
        ids, values = zip(*([(vs[0], vs[1:]) for vs in values]))
        return pytest.mark.parametrize(
            names,
            values,
            ids=ids,
        )(fn)

    return wrapper


# fmt: off
@parametrize_roundtrips(
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
)
# fmt: on
def test_io_roundtrips(
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
    X = np.array([0] * n_obs * n_var).reshape(n_obs, n_var)  # unused
    adata0 = AnnData(X=X, obs=original_df, var=var)
    ingested_uri = from_anndata(uri, adata0, "meas", obs_id_name=ingest_id_column_name)
    assert ingested_uri == uri

    # Verify column names and types
    obs_uri = join(uri, "obs")
    obs = DataFrame.open(obs_uri)
    schema = obs.schema
    assert schema.names == [SOMA_JOINID, *persisted_column_names]
    [soma_joinid_type, *string_col_types] = schema.types
    assert soma_joinid_type == pa.int64() and schema.field(0).nullable is False
    for string_col_type in string_col_types:
        assert string_col_type == pa.large_string()

    # Verify "original index metadata"
    actual_index_metadata = json.loads(
        obs.metadata[_DATAFRAME_ORIGINAL_INDEX_NAME_JSON]
    )
    assert actual_index_metadata == persisted_metadata

    # Verify outgested pd.DataFrame
    with Experiment.open(ingested_uri) as exp:
        adata1 = to_anndata(exp, "meas", obs_id_name=outgest_default_index_name)
        outgested_obs = adata1.obs
    assert_frame_equal(outgested_obs, outgested_df)
