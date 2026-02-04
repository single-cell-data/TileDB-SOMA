"""Hypothesis tests for SOMADataFrame."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Generic, TypeVar, Union

import attrs
import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from hypothesis import strategies as st
from hypothesis.extra import numpy as ht_np
from hypothesis.stateful import initialize, invariant, precondition, rule
from more_itertools import pairwise
from packaging.version import Version
from somacore.options import OpenMode

import tiledbsoma as soma

from tests.ht._array_state_machine import SOMAArrayStateMachine
from tests.ht._ht_test_config import HT_TEST_CONFIG
from tests.ht._ht_util import (
    arrow_array,
    arrow_schema,
    dataframe_datatype,
    from_datatype,
    pad_array,
    schemas_equal,
    splitss,
    tables_equal,
)
from tests.ht._ledger import ArrowTableLedgerEntry, Ledger, get_entries

# Only a subset of Arrow types are allowed as an indexed column (TileDB dimension)
DataFrameIndexTypes = [
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
    pa.string(),
    pa.large_string(),
    pa.timestamp("ns"),
]
if not HT_TEST_CONFIG["sc-62236_workaround"]:
    DataFrameIndexTypes += [
        pa.binary(),
        pa.large_binary(),
    ]

if Version(pd.__version__) >= Version("2.0.0"):
    DataFrameIndexTypes += [
        pa.timestamp("s"),
        pa.timestamp("ms"),
        pa.timestamp("us"),
    ]

AxisDomain = Union[tuple[Any, Any], list[Any], None]
Domain = Sequence[AxisDomain]


T = TypeVar("T")


@attrs.define(kw_only=True, frozen=True)
class EnumerationMetadata(Generic[T]):
    type: pa.DictionaryType
    max_categories: int = attrs.field(init=False)
    categories: tuple[T, ...] = attrs.field(factory=tuple)

    def __attrs_post_init__(self):
        # we are frozen, so use __setattr__ to bypass.
        max_categories = np.iinfo(self.type.index_type.to_pandas_dtype()).max

        # catch the corner case where the cardinality of the value type
        # is smaller than the index type.
        if pa.types.is_integer(self.type.value_type):
            max_categories = min(max_categories, np.iinfo(self.type.value_type.to_pandas_dtype()).max)
        object.__setattr__(self, "max_categories", max_categories)

    @property
    def ordered(self) -> bool:
        return self.type.ordered != 0

    @property
    def index_type(self) -> pa.DataType:
        return self.type.index_type

    @property
    def value_type(self) -> pa.DataType:
        return self.type.value_type

    @property
    def num_categories(self) -> int:
        return len(self.categories)

    def extend_categories(self, additional_categories: Sequence[T]) -> EnumerationMetadata[T]:
        return attrs.evolve(self, categories=tuple(list(self.categories) + list(additional_categories)))


@st.composite
def dataframe_schema(
    draw: st.DrawFn,
) -> tuple[tuple[str], pa.Schema, dict[str, EnumerationMetadata[Any]]]:
    """Strategy will generate a legal DataFrame schema and accompanying index names.

    Will comply with SOMA/TileDB conventions:
    * index columns must not be nullable
    * schema field order must start with index colum names and be in same order
    * must contain a `soma_joinid` column
    * must have at least two columns, one indexed, one not indexed
    """

    # initial schema draw
    schema = draw(
        arrow_schema(
            required_fields=(pa.field("soma_joinid", pa.int64(), nullable=False),),
            unique_field_names=True,
            elements=dataframe_datatype(),
        ),
    )
    assert len(schema) > 1

    # randomly choose index columns
    if draw(st.booleans()):
        # common choice; treat as such
        index_column_names = ("soma_joinid",)
    else:
        # find candidate fields to be indexed, and select random subset
        candidate_index_fields = [
            f.name
            for f in schema
            if f.type in DataFrameIndexTypes
            and not f.name.startswith(".")  # Arrow compute functions choke on table columns beginning with '.'
        ]
        assert len(candidate_index_fields) > 0  # at least one index must exist
        n_indices = draw(st.integers(min_value=1, max_value=min(len(candidate_index_fields), len(schema) - 1)))
        rng = np.random.default_rng(seed=draw(st.integers(min_value=0)))
        index_column_names = tuple(
            candidate_index_fields[i] for i in rng.choice(len(candidate_index_fields), size=n_indices, replace=False)
        )
        # TileDB dimensions may not be nullable, so just rewrite those we have selected
        for name in index_column_names:
            idx = schema.get_field_index(name)
            schema = schema.set(idx, schema.field(idx).with_nullable(False))

    # reorder schema to match index_column_names for ease of read eq tests
    reordered_fields = [schema.field(name) for name in index_column_names] + [
        f for f in schema if f.name not in index_column_names
    ]
    schema = pa.schema(reordered_fields)

    # define enumerations metadata
    enumeration_metadata: dict[str, EnumerationMetadata[Any]] = {}
    for field_idx, field in enumerate(schema):
        if field.name in index_column_names:
            continue
        if not pa.types.is_dictionary(field.type):
            continue
        enumeration_metadata[field.name] = EnumerationMetadata(type=field.type)

    assert len(schema) > 1
    assert len(index_column_names) > 0
    assert len(index_column_names) < len(schema)

    return index_column_names, schema, enumeration_metadata


def default_max_domain(datatype: pa.DataType) -> AxisDomain:
    """Return the accepted default for the domain of a given Arrow DataType.

    NB:
    * there are bugs that prescribe some values (noted inline), e.g. sc-61331
    """
    if datatype in [pa.string(), pa.large_string(), pa.binary(), pa.large_binary()]:
        return ("", "")
    if pa.types.is_floating(datatype):
        dtype = datatype.to_pandas_dtype()
        return (np.finfo(dtype).min, np.finfo(dtype).max)
    if pa.types.is_integer(datatype):
        dtype = datatype.to_pandas_dtype()
        md = (
            np.iinfo(dtype).min,
            np.iinfo(dtype).max - 2,  # sc-61331 - can't use entire range(!).
        )
        # Also, sc-61334, which has different limit for create() than change_domain().
        # Seemingly only affects int64.
        if dtype in [np.int16, np.int32, np.int64, np.uint16, np.uint32, np.uint64]:
            md = (md[0], md[1] - 2048)

        return md
    if pa.types.is_timestamp(datatype):
        return (
            pa.scalar(
                -(2**63) + 1,
                type=pa.timestamp(datatype.unit),
            ),  # NB: -2**63 is NaT, per NEP-7, and indices can't be nullable
            pa.scalar(
                2**63 - 1_000_001,
                type=pa.timestamp(datatype.unit),
            ),  # sc-61331: 1_000_001 appears to be a weird buggy magic number?
        )

    raise ValueError("Unsupported type.")


@st.composite
def dataframe_domain(
    draw: st.DrawFn,
    *,
    schema: pa.Schema,
    index_column_names: Sequence[str],
    max_domain: Domain | None = None,
    current_domain: Domain | None = None,
    apply_defaults: bool = False,
) -> Domain:
    """Strategy to generate DataFrame domains.

    If current_domain specified, will never shrink. Will not exceed max_domain.

    NB:
    * domain can't be set for string or binary index columns - use None or ('','').
    * domain can only expand.
    * all other domain values must be native python types, not pyarrow.Scalar
    """
    if max_domain is None:
        max_domain = tuple(default_max_domain(schema.field(n).type) for n in index_column_names)
    assert len(index_column_names) == len(max_domain)
    new_domain = []
    for field_index, field_name in enumerate(index_column_names):
        field = schema.field(field_name)
        if pa.types.is_primitive(field.type):
            zero = (
                pa.scalar(0, type=field.type)
                if pa.types.is_timestamp(field.type)
                else pa.scalar(0, type=field.type).as_py()
            )
            max_lower, max_upper = max_domain[field_index]
            if field_name == "soma_joinid":
                max_lower = max(0, max_lower)  # per SOMA spec
            current_lower, current_upper = current_domain[field_index] if current_domain is not None else (zero, zero)
            lower = (
                draw(
                    from_datatype(
                        field.type,
                        min_value=max_lower,
                        max_value=current_lower,
                        allow_nan=False,
                    ),
                )
                if current_lower is None or draw(st.booleans())
                else current_lower
            )
            upper = (
                draw(
                    from_datatype(
                        field.type,
                        min_value=current_upper,
                        max_value=max_upper,
                        allow_nan=False,
                    ),
                )
                if current_upper is None or draw(st.booleans())
                else current_upper
            )

            if pa.types.is_timestamp(field.type):
                # TypeError: '<=' not supported between instances of
                # 'pyarrow.lib.TimestampScalar' and 'pyarrow.lib.TimestampScalar'
                # so compare as int here
                assert int(lower.value) <= int(upper.value)
                assert int(max_lower.value) <= int(lower.value) <= int(current_lower.value)
                assert int(max_upper.value) >= int(upper.value) >= int(current_upper.value)
            else:
                lower = lower.as_py() if isinstance(lower, pa.Scalar) else lower
                upper = upper.as_py() if isinstance(upper, pa.Scalar) else upper
                assert lower <= upper
                assert max_lower <= lower <= current_lower
                assert max_upper >= upper >= current_upper

            new_domain.append((lower, upper))

        else:
            # no idea what this is, so specify default
            new_domain.append(None)

    assert len(new_domain) == len(index_column_names)
    return tuple(new_domain)


@st.composite
def column_values(
    draw: st.DrawFn,
    type: pa.DataType,
    size: int,
    is_index: bool,
    unique: bool,
    domain: tuple[int, int] | tuple[None, None],
    is_dict_value: bool = False,  # only used for bug workarounds
) -> pa.Array:
    min_value, max_value = domain

    if pa.types.is_timestamp(type):
        # don't generate NaT. ht_np.from_dtype doesn't obey min/max value
        # params, so draw ints, and then convert. NEB-7 says NaT is -2**63.
        min_value = -(2**63) + 1 if min_value is None else max(-(2**63) + 1, int(min_value.cast("int64").as_py()))
        max_value = 2**63 - 1 if max_value is None else min(2**63 - 1, int(max_value.cast("int64").as_py()))
        dtype = np.dtype(type.to_pandas_dtype())
        elements = st.builds(
            dtype.type,
            st.integers(min_value=min_value, max_value=max_value),
            st.just(type.unit),
        )
        return draw(arrow_array(type, size, elements=elements, unique=unique, padding=False))

    if pa.types.is_floating(type) and (HT_TEST_CONFIG["sc-61506_workaround"] or HT_TEST_CONFIG["sc-62449_workaround"]):
        dtype = np.dtype(type.to_pandas_dtype())
        elements = ht_np.from_dtype(dtype, min_value=min_value, max_value=max_value)
        if HT_TEST_CONFIG["sc-61506_workaround"]:
            # Array dimensions do not de-dup -0. and 0. as the same. Disable any generation
            # of negative zero until this is resolved. NB: ledger de-dup treats them a equivalent
            # per IEEE 754 semantics.
            elements = elements.filter(lambda x: not (x == 0 and np.signbit(x)))

        if HT_TEST_CONFIG["sc-62449_workaround"] and is_dict_value:
            # NaN as categorical values fails to evolve the enum correctly, do disable NaN values
            # ONLY when generating cat values.
            elements = elements.filter(lambda x: not np.isnan(x))

        return draw(arrow_array(type, size, elements=elements, unique=unique, padding=False))

    if pa.types.is_primitive(type):
        dtype = np.dtype(type.to_pandas_dtype())
        elements = ht_np.from_dtype(dtype, min_value=min_value, max_value=max_value)
        return draw(arrow_array(type, size, elements=elements, unique=unique, padding=False))

    if type in [pa.binary(), pa.large_binary()]:
        if HT_TEST_CONFIG["sc-62447_workaround"]:
            return draw(
                arrow_array(
                    np.dtype(bytes),
                    size,
                    elements=st.binary(min_size=1).filter(lambda b: b"\x00" not in b),
                    unique=unique,
                    padding=False,
                ),
            )
        return draw(arrow_array(np.dtype(bytes), size, unique=unique, padding=False))

    if type in [pa.string(), pa.large_string()]:
        dtype = np.dtype(str)
        if is_index:
            # TileDB string index columns are restricted to "7 bit ASCII". These tests use
            # Pandas indexing, which are foobared on anything containing a null.
            # So in practice, use [1,127]
            if HT_TEST_CONFIG["sc-62265_workaround"]:
                elements = st.text(alphabet=st.characters(codec="ascii", min_codepoint=1, max_codepoint=126))
            else:
                elements = st.text(alphabet=st.characters(codec="ascii", min_codepoint=1))
        else:
            # Disallow surrogate codepoints.  Arrow doesn't implement them in the
            # encoder/decoder, and will throw if they are present.
            if is_dict_value and HT_TEST_CONFIG["sc-62447_workaround"]:
                # don't allow empty string due to sc-62447
                elements = st.text(alphabet=st.characters(exclude_categories=["C"]), min_size=1)
            else:
                elements = st.text(alphabet=st.characters(exclude_categories=["C"]))

        return draw(arrow_array(dtype, size, elements=elements, unique=unique, padding=False))

    assert False, f"Unknown type: no arrow_table strategy for this type {type}"


def setdiff(a: set[Any], b: set[Any]) -> set[Any]:
    """Set diff (a-b) with nan equivalence."""

    def wo_nan(s):
        return {v for v in s if v == v}  # v!=v means Nan

    a_wo_nan, b_wo_nan = wo_nan(a), wo_nan(b)
    if a_wo_nan != a and b_wo_nan != b:
        # both had a Nan, so diff the wo_nan sets
        return a_wo_nan - b_wo_nan
    if a_wo_nan == a and b_wo_nan == b:
        # neither had a NaN, so diff the original sets
        return a - b
    if a_wo_nan != a and b_wo_nan == b:
        # a had a NaN, b did not, diff the wo sets and add a NaN.
        # this handles the case where set a had multiple (different)
        # NaNs
        return (a_wo_nan - b_wo_nan) | {np.nan}
    # b had a NaN, a did not, so just diff the wo sets
    return a_wo_nan - b_wo_nan


@st.composite
def arrow_table2(
    draw: st.DrawFn,
    schema: pa.Schema,
    index_column_names: Sequence[str],
    enumeration_metadata: dict[str, EnumerationMetadata[Any]],
    domain: Domain,
    *,
    min_size: int = 0,
) -> tuple[pa.Table, dict[str, EnumerationMetadata[Any]]]:
    index_domains = {k: v if v is not None else (None, None) for k, v in zip(index_column_names, domain)}
    is_unique = {f.name: (f.name in index_domains or f.name == "soma_joinid") for f in schema}

    # First, decide if we have any dictionary/categoricals that we want to extend
    for field in schema:
        field_name = field.name
        if pa.types.is_dictionary(field.type):
            assert field_name not in index_domains
            enmr = enumeration_metadata[field_name]
            assert enmr.type == field.type

            # extend enum categories if it is len == 0 or draw says to do it
            if enmr.num_categories == 0 or ((enmr.num_categories < enmr.max_categories) and draw(st.booleans())):
                MAX_CATEGORIES = 129
                new_cat_count = enmr.num_categories + draw(
                    st.integers(
                        min_value=1,
                        max_value=min(MAX_CATEGORIES, enmr.max_categories - enmr.num_categories),
                    ),
                )
                assert new_cat_count <= enmr.max_categories

                # draw until we have sufficient unique values
                while enmr.num_categories < new_cat_count:
                    new_cats = draw(
                        column_values(
                            field.type.value_type,
                            new_cat_count - enmr.num_categories,
                            is_index=False,
                            unique=True,
                            domain=(None, None),
                            is_dict_value=True,
                        ),
                    )
                    new_unique_cats = setdiff(set(new_cats.to_pylist()), set(enmr.categories))
                    enmr = enmr.extend_categories(new_unique_cats)

                enumeration_metadata[field_name] = enmr

    # Second, calculate size of table based upon uniqueness requirements
    def get_max_size() -> int:
        """max_size is mininimum of:
        * domain range of all int/uint/ts index domains
        * number of categories for any column with a unique draw
        """
        max_size = 1024  # default max
        for f in schema:
            if not is_unique[f.name]:
                continue
            if f.name in index_domains:
                d = index_domains[f.name]
                if pa.types.is_integer(f.type):
                    max_size = min(max_size, d[1] - d[0] + 1)
                elif pa.types.is_floating(f.type):
                    with np.errstate(over="ignore"):
                        max_size = int(
                            min(
                                max_size,
                                (d[1] - d[0]) / np.finfo(f.type.to_pandas_dtype()).tiny + 1,
                            ),
                        )
                elif pa.types.is_timestamp(f.type):
                    delta = int(d[1].cast("int64").as_py()) - int(d[0].cast("int64").as_py())
                    assert delta >= 0
                    max_size = min(max_size, delta + 1)
            elif pa.types.is_dictionary(f.type):
                max_size = min(max_size, enumeration_metadata[f.name].num_categories)

        return max_size

    size = draw(st.integers(min_value=min_size, max_value=get_max_size()))

    # Third, draw table columns
    columns = {}
    for field in schema:
        field_name = field.name
        is_index = field_name in index_domains

        if pa.types.is_dictionary(field.type):
            assert not is_index
            enmr = enumeration_metadata[field_name]

            dictionary = pa.array(enmr.categories, type=field.type.value_type)
            indices = draw(
                ht_np.arrays(
                    dtype=field.type.index_type.to_pandas_dtype(),
                    shape=(size,),
                    unique=is_unique[field_name],
                    elements=st.integers(min_value=0, max_value=enmr.num_categories - 1),
                ),
            )
            columns[field_name] = pa.DictionaryArray.from_arrays(indices, dictionary, ordered=field.type.ordered)
        else:
            domain = index_domains.get(field_name, (None, None))
            if field_name == "soma_joinid" and domain == (None, None):
                domain = (0, 2**56 - 1)
            columns[field_name] = draw(column_values(field.type, size, is_index, is_unique[field_name], domain))

    assert all(len(columns[k]) == size for k in columns)
    tbl = pa.Table.from_pydict(columns, schema)
    assert tbl.schema == schema

    # split, sometimes
    if len(tbl) > 3 and draw(st.booleans()):
        n_splits = draw(st.integers(min_value=0, max_value=max(0, len(tbl) // 10)))
        if n_splits > 0:
            split_points = draw(splitss(n_splits=n_splits, max_value=len(tbl)))
            split_points = [0, *split_points, len(tbl)]
            tbl = pa.concat_tables([tbl[start:end] for start, end in pairwise(split_points)])

    # pad, sometimes
    if draw(st.booleans()):
        batches = tbl.to_batches()
        batch_to_pad = draw(st.integers(min_value=0, max_value=len(batches) - 1))
        batch_arrays = [pad_array(arr, draw) for arr in batches[batch_to_pad].columns]
        batches[batch_to_pad] = pa.RecordBatch.from_arrays(batch_arrays, schema=tbl.schema)
        tbl = pa.Table.from_batches(batches)

    return tbl, enumeration_metadata


class SOMADataFrameStateMachine(SOMAArrayStateMachine):
    def __init__(self) -> None:
        super().__init__()

    @initialize(data=st.data(), dataframe_schema=dataframe_schema())
    def setup(
        self,
        data: st.DataObject,
        dataframe_schema: tuple[Sequence[str], pa.Schema, dict[str, EnumerationMetadata[Any]]],
    ) -> None:
        # Schema in total includes:  arrow schema, index column names and current enumerations for
        # any dictionary columns. These must be evolved as a unit, as there are lots of cross-
        # dependencies (e.g, some types may not be an index column).
        self.index_column_names, self.schema, self.enumeration_metadata = dataframe_schema
        self.domain = data.draw(  # TODO XXX: should be a ledger
            dataframe_domain(schema=self.schema, index_column_names=self.index_column_names),
        )
        super().setup(
            soma.DataFrame.create(
                self.uri,
                schema=self.schema,
                domain=self.domain,
                index_column_names=self.index_column_names,
                context=self.context,
                tiledb_timestamp=None,  # TODO: no time-travel for now
            ),
        )
        self.domain = self.A.domain
        assert not self.A.closed
        assert self.A.mode == "w"
        self.data_ledger = Ledger[ArrowTableLedgerEntry](
            initial_entry=ArrowTableLedgerEntry(
                data=self.schema.empty_table(),
                timestamp_ms=self.A.tiledb_timestamp_ms,
                name="initial entry",
                index_columns=self.index_column_names,
            ),
            allows_duplicates=False,
        )

    def _array_exists(self, uri: str, context: soma.SOMAContext, tiledb_timestamp: int | None) -> bool:
        return soma.DataFrame.exists(uri, context=context, tiledb_timestamp=tiledb_timestamp)

    def _array_open(self, *, mode: OpenMode, tiledb_timestamp: int | None = None) -> None:
        self.A = soma.DataFrame.open(self.uri, mode=mode, context=self.context, tiledb_timestamp=tiledb_timestamp)

    ##
    # --- schema
    ##

    @precondition(lambda self: not self.closed)
    @invariant()
    def check_schema(self) -> None:
        assert isinstance(self.A, soma.DataFrame)
        assert self.A.soma_type == "SOMADataFrame"
        assert schemas_equal(
            self.schema,
            self.A.schema,
            ignore_field_order=True,
            large_type_equivalence=True,
        )
        assert sorted(self.schema.names) == sorted(self.A.keys())
        assert self.index_column_names == self.A.index_column_names

    ##
    # --- domain
    ##

    @precondition(lambda self: not self.closed)
    @invariant()
    def check_domain(self) -> None:
        domain = []
        for iname, idomain in zip(self.index_column_names, self.domain):
            if idomain is not None:
                domain.append(idomain)
            else:
                type = self.schema.field(iname).type
                if type in [
                    pa.string(),
                    pa.large_string(),
                    pa.binary(),
                    pa.large_binary(),
                ]:
                    domain.append(("", ""))
                elif pa.type.is_primitive(type):
                    domain.append((0, 0))
                elif pa.types.is_timestamp(type):
                    zero_ts = pa.scalar(0, type=type)
                    domain.append((zero_ts, zero_ts))
                else:
                    domain.append(None)

        assert self.A.domain == tuple(
            domain,
        ), f"Unexpected domain in {self.A}: had {self.A.domain}, expected {self.domain}"

    @rule(data=st.data())
    def expand_domain(self, data: st.DataObject) -> None:
        assert self.index_column_names == self.A.index_column_names
        new_domain = data.draw(
            dataframe_domain(
                schema=self.schema,
                index_column_names=self.index_column_names,
                current_domain=self.domain,
                max_domain=self.A.maxdomain,
                apply_defaults=True,
            ),
        )

        # Always re-open at latest. Without this, there is a good chance we will end up with a
        # schema fragment with the same timestamp. Concrete case: a write has been done to an
        # array that has a dict field, which triggers an automatic schema evolution.
        if not self.closed:
            self._close()
        self._open(mode="w")
        assert self.mode == "w"

        self.A.change_domain(new_domain)
        self.domain = new_domain  # TODO XXX should be a ledger
        self._close()  # domain is committed upon close
        assert self.A.closed

    ##
    # --- data
    ##

    @precondition(lambda self: not self.closed and self.mode == "r")
    @invariant()
    def check_read_all(self) -> None:
        timestamp_ms = self.A.tiledb_timestamp_ms
        sort_order = [(name, "ascending") for name in self.index_column_names]
        expected = self.data_ledger.read(timestamp_ms=timestamp_ms).to_table().sort_by(sort_order)
        found = self.A.read().concat().sort_by(sort_order)
        assert tables_equal(found, expected, equal_nan=True), f"{found}\n is not equal to {expected}"

    @precondition(lambda self: not self.closed and self.mode == "r")
    @invariant()
    def check_count(self) -> None:
        expected = len(self.data_ledger.read(timestamp_ms=self.A.tiledb_timestamp_ms).to_table())
        assert expected == self.A.count, "count mismatch"
        assert self.A.count <= self.A._handle.fragment_cell_count(), "count vs fragment_cell_count inconsistency"

    @precondition(lambda self: not self.closed and self.mode == "w")
    @precondition(
        lambda self: self.A.tiledb_timestamp_ms not in self.data_ledger.timestamps,
    )  # only one write per timestamp until sc-61223 (is FIXED) and sc-61226 are fixed
    @rule(data=st.data())
    def write(self, data: st.DataObject) -> None:
        df_tbl, self.enumeration_metadata = data.draw(
            arrow_table2(
                self.schema,
                self.index_column_names,
                self.enumeration_metadata,
                self.domain,
                min_size=1,
            ),
        )
        fragments_before_write = get_entries(f"{self.uri}/__fragments")
        self.A.write(df_tbl)
        new_fragments = set(get_entries(f"{self.uri}/__fragments")) - set(fragments_before_write)
        assert len(new_fragments) == (1 if not HT_TEST_CONFIG["sc-61462_workaround"] else df_tbl[0].num_chunks)
        self.data_ledger.write(
            ArrowTableLedgerEntry(
                timestamp_ms=self.A.tiledb_timestamp_ms,
                name=new_fragments.pop(),
                data=df_tbl,
                index_columns=self.index_column_names,
            ),
        )


TestSOMADataFrame = pytest.mark.usefixtures("setup_fixtures")(SOMADataFrameStateMachine.TestCase)
