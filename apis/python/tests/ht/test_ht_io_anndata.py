"""
Hypothesis-based tests of some.io.{to_from}_anndata.
"""

from __future__ import annotations

import math
import operator
import string
import sys
from collections import OrderedDict
from collections.abc import Hashable, Mapping, Sequence
from datetime import timedelta
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Literal,
    TypeVar,
    get_args,
)

import anndata as ad
import deepdiff
import hypothesis as ht
import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.sparse as sp
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra import numpy as ht_np
from hypothesis.extra import pandas as ht_pd
from numpy.testing import assert_allclose
from packaging.version import Version
from pandas.testing import assert_frame_equal
from typeguard import suppress_type_checks

import tiledbsoma
import tiledbsoma.io

from tests.ht._ht_test_config import HT_TEST_CONFIG
from tests.ht._ht_util import posix_filename

if TYPE_CHECKING:
    Ex_co = TypeVar("Ex_co", covariant=True, default=Any)
else:
    Ex_co = TypeVar("Ex_co", covariant=True)

T = TypeVar("T")


def df_col_dtypes() -> st.SearchStrategy[npt.DTypeLike]:
    """obs/var column types supported by TileDB-SOMA _and_ hypothesis.extras.pandas.data_frames."""
    return st.one_of(
        ht_np.byte_string_dtypes(endianness="="),  # |S
        ht_np.unicode_string_dtypes(endianness="="),  # |U
        ht_np.boolean_dtypes(),
        ht_np.integer_dtypes(endianness="="),
        ht_np.unsigned_integer_dtypes(endianness="="),
        ht_np.floating_dtypes(sizes=(32, 64), endianness="="),
    )


def df_index_dtypes() -> st.SearchStrategy[npt.DTypeLike]:
    """Can only include types supported by both TileDB-SOMA _and_ hypothesis.extras.pandas.indexes."""
    return st.one_of(
        ht_np.integer_dtypes(endianness="="),
        ht_np.unsigned_integer_dtypes(endianness="="),
        ht_np.floating_dtypes(sizes=(32, 64), endianness="="),
        ht_np.unicode_string_dtypes(endianness="="),  # generates object dtype
    )


def elements_and_dtypes(dtype: np.dtype, is_index: bool) -> tuple[np.dtype, st.SearchStrategy[Any]]:
    """Given dtype, return elements and dtype strategy"""
    elements = None  # default
    if dtype.kind == "m":
        elements = ht_np.from_dtype(
            np.dtype("timedelta64[ns]"),
            min_value=-10_000_000,
            max_value=10_000_000,
            allow_nan=not is_index,
        )
    elif dtype.kind == "M":
        elements = ht_np.from_dtype(
            np.dtype("datetime64[ns]"),
            min_value=-10_000_000,
            max_value=10_000_000,
            allow_nan=not is_index,
        )
    elif dtype.kind == "f" and is_index:
        elements = ht_np.from_dtype(dtype, allow_nan=False)
    elif dtype.kind == "U":
        # disallow surrogate codepoints
        if HT_TEST_CONFIG["sc-63404_workaround"]:
            elements = ht_np.from_dtype(
                dtype,
                alphabet=st.characters(exclude_categories=["C"], exclude_characters=["\x00"]),
                min_size=1,
            )
        else:
            elements = ht_np.from_dtype(
                dtype,
                alphabet=st.characters(exclude_categories=["C"], exclude_characters=["\x00"]),
            )
    elif dtype.kind in ("S", "a"):
        max_size = dtype.itemsize or None
        if HT_TEST_CONFIG["sc-63404_workaround"]:
            elements = st.binary(min_size=1, max_size=max_size).filter(lambda b: b"\0" not in b)
        else:
            elements = st.binary(max_size=max_size).filter(lambda b: b"\0" not in b)

    return elements, dtype


@st.composite
def dataframe_indexes(draw: st.DrawFn, size: int, name: str) -> st.SearchStrategy[pd.Index]:
    """Strategy that returns an index-building strategy"""
    if draw(st.booleans()):
        return ht_pd.range_indexes(min_size=size, max_size=size, name=st.just(f"{name}_index"))

    elements, dtype = elements_and_dtypes(draw(df_index_dtypes()), is_index=True)
    return ht_pd.indexes(
        min_size=size,
        max_size=size,
        name=st.just(f"{name}_index"),
        dtype=dtype,
        elements=elements,
        unique=True,
    )


@st.composite
def dataframe_columns(draw: st.DrawFn, name: str) -> ht_pd.column:
    elements, dtype = elements_and_dtypes(draw(df_col_dtypes()), is_index=False)
    return ht_pd.column(name=name, dtype=dtype, elements=elements)


@st.composite
def dataframes(draw: st.DrawFn, size: int, name: str, n_cols: int | None = None) -> pd.DataFrame:
    """
    Strategy to generate dataframe compatible with obs/var requirements
    (must be made of types/values supported by both TileDB-SOMA _and_ AnnData).
    """

    index = draw(dataframe_indexes(size=size, name=name))
    n_cols = draw(st.integers(min_value=1, max_value=10)) if n_cols is None else n_cols
    columns = [draw(dataframe_columns(name=f"{name}_{c}")) for c in range(n_cols)]
    df = draw(ht_pd.data_frames(columns=columns, index=index))

    # Hypothesis.extras.pandas does not support generating categorical or (pandas) string
    # types, so do those as post-processing steps.
    for c in columns:
        if c.dtype.kind in ("S", "U", "b", "i", "u") and df[c.name].dtype != "category" and draw(st.booleans()):
            df[c.name] = df[c.name].astype("category")

        # AnnData does not (yet) support Pandas string-type. See https://github.com/scverse/anndata/issues/1571
        # if (
        #     c.dtype.kind == "U"
        #     and df[c.name].dtype != "category"
        #     and draw(st.booleans())
        # ):
        #     df[c.name] = df[c.name].astype("string")

    return df


def keys() -> st.SearchStrategy[str]:
    # if sc-63410_workaround, we can't handle anything that isn't a legal posix path name.
    if HT_TEST_CONFIG["sc-63410_workaround"]:
        stgy = posix_filename()
    else:
        stgy = st.text(string.printable, min_size=1, max_size=20)
    return stgy.filter(lambda s: not s.startswith("soma_"))


# AnnData <= 0.10 does not support scipy sparse_array
if Version(ad.__version__) >= Version("0.11.0"):
    MatrixFormats = Literal["csr_matrix", "csc_matrix", "csr_array", "csc_array", "ndarray", "ma"]
else:
    MatrixFormats = Literal["csr_matrix", "csc_matrix", "ndarray", "ma"]


@st.composite
def matrixes(
    draw: st.DrawFn,
    shape: tuple[int, int] | st.SearchStrategy[tuple[int, ...]],
    formats: MatrixFormats | None = None,
) -> sp.csr_matrix | sp.csc_matrix | np.ndarray | np.ma.MaskedArray:
    """Random 2D array in a variety of formats supported by AnnData."""

    if isinstance(shape, st.SearchStrategy):
        shape = draw(shape)

    rng = np.random.default_rng(seed=draw(st.integers(min_value=0)))
    dtype = draw(
        st.one_of(
            ht_np.floating_dtypes(sizes=(32, 64), endianness="="),  # float16 not supported
            ht_np.integer_dtypes(endianness="="),
            ht_np.unsigned_integer_dtypes(endianness="="),
            # ht_np.boolean_dtypes(),  # bools not currently supported by TileDB-SOMA
        ),
    )

    size = math.prod(shape)
    if dtype == np.bool_:

        def data_sampler(size):
            return rng.integers(0, high=1, size=size, dtype=dtype, endpoint=True)

    elif np.issubdtype(dtype, np.integer):

        def data_sampler(size):
            return rng.integers(
                np.iinfo(dtype).min,
                high=np.iinfo(dtype).max,
                size=size,
                dtype=dtype,
                endpoint=False,
            )

    elif np.issubdtype(dtype, np.floating):
        data_sampler = rng.uniform
    else:
        raise TypeError("Unsupported matrix type.")

    vals = data_sampler(size=size).astype(dtype, copy=False).reshape(shape)

    # TODO: AnnData accepts a pd.Dataframe as X value

    formats = formats or get_args(MatrixFormats)
    format = draw(st.sampled_from(formats))
    if format == "ndarray":
        return vals

    # the remainder are sparse, so create a mask
    density = 0.2
    masked = rng.choice(size, size=int(size * (1 - density)), replace=False)  # coordinates
    mask = np.ones((size,), dtype=np.bool_)
    mask[masked] = 0
    mask = mask.reshape(vals.shape)

    if format == "ma":
        return np.ma.masked_array(vals, mask=mask)
    if format in ("csr_matrix", "csr_array"):
        vals.reshape((size,))[masked] = 0
        return sp.csr_matrix(vals) if format == "csr_matrix" else sp.csr_array(vals)
    if format in ("csc_matrix", "csc_array"):
        vals.reshape((size,))[masked] = 0
        return sp.csc_matrix(vals) if format == "csc_matrix" else sp.csc_array(vals)
    raise ValueError(f"oops, bad code - missing format path for {format}")


@st.composite
def matrix_shapes(
    draw: st.DrawFn,
    prelude: tuple[int, ...],
    min_dims: int = 0,
    max_dims: int | None = None,
    min_side: int = 1,
    max_side: int | None = None,
) -> tuple[int, int]:
    if max_dims is None:
        max_dims = min(min_dims + 2, 5)
    if max_side is None:
        max_side = min_side + 5

    min_dims -= len(prelude)
    max_dims -= len(prelude)
    if min_dims <= 0:
        return prelude

    postlude = draw(st.lists(st.integers(min_side, max_side), min_size=min_dims, max_size=max_dims))
    return tuple(list(prelude) + postlude)


@st.composite
def map_of_matrixes(
    draw: st.DrawFn,
    shape_prelude: tuple[int, ...],
    formats: MatrixFormats | None = None,
) -> dict[str, np.ndarray | sp.sparray | sp.spmatrix | np.ma.MaskedArray]:
    def key_unique_by(i: str) -> str:
        return i[0].lower() if HT_TEST_CONFIG["sc-63402_workaround"] else i[0]

    return draw(
        st.one_of(
            dictionaries_unique_by(
                keys=keys(),
                values=matrixes(
                    shape=matrix_shapes(
                        prelude=shape_prelude,
                        min_dims=2,
                        max_dims=2 if HT_TEST_CONFIG["sc-63409_workaround"] else None,
                    ),
                    formats=formats,
                ),
                max_size=3,
                unique_by=key_unique_by,
            ),
            st.none(),
        ),
    )


@st.composite
def monomorphic_list(
    draw: st.DrawFn,
    element_type: Sequence[st.SearchStrategy[Any]],
    min_size: int = 0,
    max_size: int = 10,
) -> list[Any]:
    """List of elements of a single type. from_anndata saving of `uns` does not support save anything polymorphic."""
    elmt_type = draw(st.sampled_from(element_type))
    return draw(st.lists(elmt_type, min_size=min_size, max_size=max_size))


def dictionaries_unique_by(
    keys: st.SearchStrategy[Ex_co],
    values: st.SearchStrategy[T],
    *,
    dict_class: type = dict,
    min_size: int = 0,
    max_size: int | None = None,
    unique_by: Callable[[Ex_co], Hashable] | tuple[Callable[[Ex_co], Hashable], ...] | None = None,
) -> st.SearchStrategy[dict[Ex_co, T]]:
    """
    Identical to hypothesis.strategies.dictionaries, except it allows user-configurable
    `unique_by` param, AND has less error checking.
    """
    if max_size == 0:
        return st.fixed_dictionaries(dict_class())

    if unique_by is None:
        unique_by = operator.itemgetter(0)

    return st.lists(
        st.tuples(keys, values),
        min_size=min_size,
        max_size=max_size,
        unique_by=unique_by,
    ).map(dict_class)


@st.composite
def unses(draw: st.DrawFn) -> dict[str, Any]:
    # the types AnnData.write will allow
    np_dtypes = st.one_of(
        ht_np.boolean_dtypes(),
        ht_np.integer_dtypes(endianness="=" if HT_TEST_CONFIG["sc-63459_workaround"] else "?"),
        ht_np.unsigned_integer_dtypes(endianness="=" if HT_TEST_CONFIG["sc-63459_workaround"] else "?"),
        ht_np.floating_dtypes(
            sizes=(32, 64),
            endianness="=" if HT_TEST_CONFIG["sc-63459_workaround"] else "?",
        ),
    )
    value_types = (
        # st.none(),  # both anndata.write and tiledbsoma.io.from_anndata drop any elements with a value of None
        st.booleans(),
        st.floats(),
        st.integers(min_value=-(2**63), max_value=2**63 - 1),  # outside this range unsupported by TileDB metadata
        st.text(string.printable),
        # | st.binary()  # Currently unsupported by Anndata write
    )
    base = st.one_of(
        *value_types,
        # from_anndata saves "list" as a TileDB array, which must be monomorphic
        monomorphic_list(value_types, min_size=1),
        ht_np.arrays(
            dtype=np_dtypes,
            shape=ht_np.array_shapes(max_dims=2 if HT_TEST_CONFIG["sc-63409_workaround"] else None),
        ),
    )

    def key_unique_by(i: str) -> str:
        return i[0].lower() if HT_TEST_CONFIG["sc-63402_workaround"] else i[0]

    return draw(
        st.one_of(
            st.none(),
            dictionaries_unique_by(
                keys=keys(),
                values=st.recursive(
                    base,
                    lambda children: dictionaries_unique_by(keys=keys(), values=children, unique_by=key_unique_by),
                    max_leaves=10,
                ),
                unique_by=key_unique_by,
            ),
        ),
    )


@st.composite
def anndatas(draw: st.DrawFn) -> ad.AnnData:
    """
    Generate a non-empty (i.e., minimum shape (1,1) anndata)
    """
    n_obs = draw(st.integers(min_value=1, max_value=100))
    n_vars = draw(st.integers(min_value=1, max_value=100))

    obs = draw(dataframes(size=n_obs, name="obs"))
    var = draw(dataframes(size=n_vars, name="var"))
    # AnnData likes a string index
    obs.index = obs.index.astype(str)
    var.index = var.index.astype(str)

    X = draw(st.one_of(matrixes(shape=(n_obs, n_vars))), st.none())
    layers = draw(map_of_matrixes(shape_prelude=(n_obs, n_vars)))
    obsm = draw(map_of_matrixes(shape_prelude=(n_obs,), formats=("ndarray",)))
    varm = draw(map_of_matrixes(shape_prelude=(n_vars,), formats=("ndarray",)))
    obsp = draw(map_of_matrixes(shape_prelude=(n_obs, n_obs), formats=("ndarray",)))
    varp = draw(map_of_matrixes(shape_prelude=(n_vars, n_vars), formats=("ndarray",)))

    uns = draw(unses())

    # TODO: this approach does not accurately represent real-world use cases, where
    # X is often a transformation of X.raw, and both X and var may be a slice
    # of raw.X/raw.var on the obs axis.
    raw = {"X": X, "var": var}

    return ad.AnnData(
        X=X,
        obs=obs,
        var=var,
        uns=uns,
        obsp=obsp,
        obsm=obsm,
        layers=layers,
        varm=varm,
        varp=varp,
        raw=raw,
    )


def assert_matrix_equal(m1: np.ndarray | sp.sparray | None, m2: np.ndarray | sp.sparray | None) -> None:
    if m1 is None and m2 is None:
        return

    assert m1 is not None and m2 is not None

    if sp.issparse(m1):
        m1 = m1.toarray()
    if sp.issparse(m2):
        m2 = m2.toarray()

    assert_allclose(m1, m2)


def assert_deep_eq(a: Any, b: Any, **kwargs: Any) -> None:
    diff = deepdiff.DeepDiff(a, b, **kwargs)
    assert diff == {}, repr(diff.to_dict())


def assert_frame_equal_strict(f1: pd.DataFrame, f2: pd.DataFrame) -> None:
    assert_frame_equal(f1, f2, check_index_type=True, check_column_type=True, check_exact=True)


def assert_map_of_matrix_equal(m1: Mapping[str, np.ndarray], m2: dict[str, np.ndarray]) -> None:
    assert m1.keys() == m2.keys()
    for k in m1:
        assert_matrix_equal(m1[k], m2[k])


def ndarray_equal(a: np.ndarray, b: np.ndarray) -> bool:
    """See https://github.com/numpy/numpy/issues/16377

    numpy.array_equal w/ equal_nan=True will fail on some types e.g., strig.
    We only need it for floats, so specialize.
    """
    if np.issubdtype(a.dtype, np.floating) or np.issubdtype(b.dtype, np.floating):
        return np.allclose(a, b, equal_nan=True)
    return np.array_equal(a, b)


def assert_uns_equal(src_adata: ad.AnnData, read_adata: ad.Anndata) -> None:
    """
    src_adata is original in-mem AnnData (oracle)
    read_adata is the read-back round-tripped AnnData created by tiledbsoma.io
    This order is important as there are certain transformations expected from
    source to read-back (e.g., list->ndarray).
    """

    src_uns, read_uns = src_adata.uns, read_adata.uns

    # AnnData will generate an OrderedDict in some situations
    if isinstance(src_uns, OrderedDict):
        src_uns = dict(src_uns)
    if isinstance(read_uns, OrderedDict):
        read_uns = dict(read_uns)

    diff = deepdiff.DeepDiff(src_uns, read_uns, ignore_nan_inequality=True)

    # Ignore expected differences
    for key in list(diff.get("type_changes", ())):
        chng = diff["type_changes"][key]

        # All `list` objects are saved as SOMA array, and read back as ndarray. Ignore if values are "eq"-ish
        if (
            chng["old_type"] is list
            and isinstance(chng["new_value"], np.ndarray)
            and ndarray_equal(np.asarray(chng["old_value"]), chng["new_value"])
        ):
            del diff["type_changes"][key]
            continue

        # bools incorrectly read back as uint8
        if HT_TEST_CONFIG["sc-63447_workaround"] and (
            chng["old_type"] in (np.bool_, list)
            and chng["new_type"] == np.uint8
            and ndarray_equal(np.asarray(chng["old_value"]).astype(np.uint8), chng["new_value"])
        ):
            del diff["type_changes"][key]
            continue

    if diff.get("type_changes", None) == {}:
        del diff["type_changes"]


def assert_anndata_equal(src_adata: ad.AnnData, read_adata: ad.AnnData) -> None:
    assert (
        src_adata.shape == read_adata.shape
        and src_adata.n_obs == read_adata.n_obs
        and src_adata.n_vars == read_adata.n_vars
    ), f"AnnData shape is not eq. Got {read_adata.shape}, expected {src_adata.shape}"
    assert src_adata.X is None or sp.issparse(read_adata.X)

    assert_frame_equal_strict(src_adata.obs, read_adata.obs)
    assert_frame_equal_strict(src_adata.var, read_adata.var)
    assert_matrix_equal(src_adata.X, read_adata.X)

    assert_map_of_matrix_equal(src_adata.layers, read_adata.layers)
    assert_map_of_matrix_equal(src_adata.obsm, read_adata.obsm)
    assert_map_of_matrix_equal(src_adata.obsp, read_adata.obsp)
    assert_map_of_matrix_equal(src_adata.varm, read_adata.varm)
    assert_map_of_matrix_equal(src_adata.varp, read_adata.varp)

    # sc-63461 - io.to_anndata does not load raw.
    #
    # if src_adata.raw or read_adata.raw:
    #     assert_frame_equal_strict(src_adata.raw.var, read_adata.raw.var)
    #     assert_matrix_equal(src_adata.raw.X, read_adata.raw.X)

    assert_uns_equal(src_adata, read_adata)


@settings(
    suppress_health_check=(ht.HealthCheck.function_scoped_fixture,),
    deadline=timedelta(milliseconds=2500),
)
@given(
    data=st.data(),
    adata=anndatas(),
    measurement_name=posix_filename(),
    context=st.one_of(st.from_type(tiledbsoma.SOMAContext), st.none()),
)
def test_roundtrip_from_anndata_to_anndata(
    data: st.DataFn,
    adata: ad.AnnData,
    measurement_name: str,
    context: tiledbsoma.SOMAContext | None,
    tmp_path_factory,  # fixture
) -> None:
    """
    Round-trip an AnnData.

    Of note:
    * to_anndata does not read "raw". So it can't round-trip this part of the AnnData
    * AnnData X can be many types, but when read back will always be a sparse array
    """

    test_path = tmp_path_factory.mktemp("anndata-")
    experiment_uri = (test_path / "soma").as_posix()

    # save the anndata for test failure debugging convenince and to ensure that the
    # generated AnnData comply with the anndata HDF5 writer.
    try:
        adata.write(test_path / "adata.h5ad")
    except Exception as e:
        print("Unable to save AnnData!", e)
        print(adata)
        assert False, repr(e)

    # Pick X and raw.X layer names that are NOT already used by the layers mapping
    # Case insensitive on MacOS.
    if sys.platform == "darwin":
        used_layer_names = set(lname.upper() for lname in adata.layers)
        X_layer_name = data.draw(posix_filename().filter(lambda x: x.upper() not in used_layer_names))
        raw_X_layer_name = data.draw(
            posix_filename().filter(lambda x: x.upper() not in used_layer_names and x.upper() != X_layer_name.upper())
        )
    else:
        X_layer_name = data.draw(posix_filename().filter(lambda x: x not in adata.layers))
        raw_X_layer_name = data.draw(posix_filename().filter(lambda x: x not in adata.layers and x != X_layer_name))

    with suppress_type_checks():
        tiledbsoma.io.from_anndata(
            experiment_uri,
            adata,
            measurement_name=measurement_name,
            context=context,
            X_layer_name=X_layer_name,
            raw_X_layer_name=raw_X_layer_name,
        )

        with tiledbsoma.Experiment.open(experiment_uri, context=context) as E:
            read_adata = tiledbsoma.io.to_anndata(
                E,
                measurement_name=measurement_name,
                X_layer_name=X_layer_name if adata.X is not None else None,
                extra_X_layer_names=list(adata.layers.keys()),
            )

            # TODO: io.to_anndata does not load raw. Do a manual verification
            # of the arrays created by io.from_anndata.
            #
            # NB: from_anndata allows user-specified naming of the X layer of raw,
            # but NOT the actual collection name. Hard-wire for now. See sc-63483
            if adata.raw is not None:
                assert_frame_equal_strict(
                    adata.raw.var.reset_index(),
                    E.ms["raw"]["var"].read().concat().to_pandas().drop(columns="soma_joinid"),
                )
                if adata.raw.X is not None:
                    tbl = E.ms["raw"]["X"][raw_X_layer_name].read().tables().concat()
                    data, i, j = (
                        tbl["soma_data"].to_numpy(),
                        tbl["soma_dim_0"].to_numpy(),
                        tbl["soma_dim_1"].to_numpy(),
                    )
                    X = sp.csr_matrix(
                        (data, (i, j)),
                        dtype=tbl.schema.field("soma_data").type.to_pandas_dtype(),
                        shape=E.ms["raw"]["X"][raw_X_layer_name].shape,
                    )
                    assert_matrix_equal(adata.raw.X, X)

    assert_anndata_equal(adata, read_adata)
