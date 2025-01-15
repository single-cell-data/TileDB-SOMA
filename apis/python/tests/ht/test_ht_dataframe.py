"""Hypothesis tests for SOMADataFrame."""

from __future__ import annotations

from typing import Any, Sequence, Union

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from hypothesis import strategies as st
from hypothesis.extra import numpy as ht_np
from hypothesis.extra import pandas as ht_pd
from hypothesis.stateful import initialize, invariant, precondition, rule
from more_itertools import pairwise
from packaging.version import Version

import tiledbsoma as soma

from tests.ht._array_state_machine import SOMAArrayStateMachine
from tests.ht._ht_test_config import HT_TEST_CONFIG
from tests.ht._ht_util import (
    arrow_schema,
    df_to_table,
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
    pa.binary(),
    pa.large_binary(),
    pa.string(),
    pa.large_string(),
    pa.timestamp("ns"),
]
if Version(pd.__version__) >= Version("2.0.0"):
    DataFrameIndexTypes += [
        pa.timestamp("s"),
        pa.timestamp("ms"),
        pa.timestamp("us"),
    ]

AxisDomain = Union[None, tuple[Any, Any], list[Any]]
Domain = Sequence[AxisDomain]


@st.composite
def dataframe_schema(draw: st.DrawFn) -> tuple[Sequence[str], pa.Schema]:
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
            elements=st.from_type(pa.DataType).filter(
                lambda t: (
                    pa.types.is_primitive(t)
                    and not (pa.types.is_timestamp(t) and t.tz is not None)
                    and not pa.types.is_time(t)
                    and not pa.types.is_date(t)
                    and t
                    not in [
                        pa.float16(),
                    ]
                )
            ),
        )
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
            and not f.name.startswith(
                "."
            )  # Arrow compute functions choke on table columns beginning with '.'
        ]
        assert len(candidate_index_fields) > 0  # at least one index must exist
        n_indices = draw(
            st.integers(
                min_value=1, max_value=min(len(candidate_index_fields), len(schema) - 1)
            )
        )
        rng = np.random.default_rng(seed=draw(st.integers(min_value=0)))
        index_column_names = tuple(
            candidate_index_fields[i]
            for i in rng.choice(
                len(candidate_index_fields), size=n_indices, replace=False
            )
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

    assert len(schema) > 1
    assert len(index_column_names) > 0
    assert len(index_column_names) < len(schema)

    return index_column_names, schema


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
        # return Numpy! See sc-61328 and sc-61329
        return (
            np.datetime64(
                -(2**63) + 1, datatype.unit
            ),  # NB: -2**63 is NaT, per NEP-7, and indices can't be nullable
            np.datetime64(
                2**63 - 1_000_001, datatype.unit
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
) -> Domain:
    """Strategy to generate DataFrame domains.

    If current_domain specified, will never shrink. Will not exceed max_domain.

    NB:
    * domain can't be set for string or binary index columns - use None or ('','').
    * domain can only expand.
    * timestamp64 domain must be specified as a numpy.datetime64 (see sc-61328 and sc-61329)
    * all other domain values must be native python types, not pyarrow.Scalar
    """
    if max_domain is None:
        max_domain = tuple(
            default_max_domain(schema.field(n).type) for n in index_column_names
        )
    assert len(index_column_names) == len(max_domain)
    new_domain = []
    for field_index, field_name in enumerate(index_column_names):
        field = schema.field(field_name)
        if not pa.types.is_primitive(field.type):
            new_domain.append(None)  # i.e., noop, use default
        else:
            zero = (
                np.datetime64(0, field.type.unit)
                if pa.types.is_timestamp(field.type)
                else pa.scalar(0, type=field.type).as_py()
            )
            max_lower, max_upper = max_domain[field_index]
            if field_name == "soma_joinid":
                max_lower = max(0, max_lower)  # per SOMA spec
            current_lower, current_upper = (
                current_domain[field_index]
                if current_domain is not None
                else (zero, zero)
            )
            lower = (
                draw(
                    from_datatype(
                        field.type,
                        min_value=max_lower,
                        max_value=current_lower,
                        allow_nan=False,
                    )
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
                    )
                )
                if current_upper is None or draw(st.booleans())
                else current_upper
            )

            # timestamp64 columns only accept np.datetime64 for domain (see sc-61328 and sc-61329)
            # In addition, pa.TimestampScalar overflows in a variety of situations, so don't use it
            # (e.g., `pa.scalar(-161650356352888167,type=pa.timestamp('s')).as_py()` )
            if pa.types.is_timestamp(field.type):
                lower = (
                    np.datetime64(lower.value, field.type.unit)
                    if isinstance(lower, pa.TimestampScalar)
                    else lower
                )
                upper = (
                    np.datetime64(upper.value, field.type.unit)
                    if isinstance(upper, pa.TimestampScalar)
                    else upper
                )
            else:
                lower = lower.as_py() if isinstance(lower, pa.Scalar) else lower
                upper = upper.as_py() if isinstance(upper, pa.Scalar) else upper

            assert lower <= upper
            assert max_lower <= lower <= current_lower
            assert max_upper >= upper >= current_upper
            new_domain.append((lower, upper))

    assert len(new_domain) == len(index_column_names)
    return tuple(new_domain)


@st.composite
def arrow_table(
    draw: st.DrawFn,
    schema: pa.Schema,
    index_column_names: Sequence[str],
    domain: Domain,
    *,
    min_size: int | None = None,
) -> pa.Table:
    """Strategy to generate Arrow Tables which:
    * match the schema
    * have unique values in the index columns
    * have values within the domain for the index columns
    """
    index_domains = {k: v for k, v in zip(index_column_names, domain)}
    columns = []
    for field in schema:
        name = field.name
        dtype = np.dtype(field.type.to_pandas_dtype())
        unique = name in index_column_names or name == "soma_joinid"
        elements = None

        min_value, max_value = index_domains.get(name, (None, None))
        assert name in index_domains or (min_value is None and max_value is None)

        if pa.types.is_timestamp(field.type):
            # don't generate NaT. ht_np.from_dtype doesn't obey min/max value
            # params, so draw ints, and then convert. NEB-7 says NaT is -2**63.
            min_value = (
                -(2**63) + 1
                if min_value is None
                else max(-(2**63) + 1, int(min_value.astype(np.int64)))
            )
            max_value = (
                2**63 - 1
                if max_value is None
                else min(2**63 - 1, int(max_value.astype(np.int64)))
            )
            elements = st.builds(
                dtype.type,
                st.integers(min_value=min_value, max_value=max_value),
                st.just(field.type.unit),
            )

        elif pa.types.is_primitive(field.type):
            elements = ht_np.from_dtype(dtype, min_value=min_value, max_value=max_value)
            # Array dimensions do not de-dup -0. and 0. as the same. Disable any generation
            # of negative zero until this is resolved. NB: ledger de-dup treats them a equivalent
            # per IEEE 754 semantics.
            if HT_TEST_CONFIG["sc-61506_workaround"] and pa.types.is_floating(
                field.type
            ):
                elements = elements.filter(lambda x: not (x == 0 and np.signbit(x)))

        # else, use default

        columns.append(
            ht_pd.column(name=name, dtype=dtype, unique=unique, elements=elements)
        )

    df = draw(
        ht_pd.data_frames(columns=columns, index=ht_pd.range_indexes(min_size=min_size))
    )
    assert min_size is None or len(df) >= min_size
    tbl = df_to_table(df, schema=schema)
    assert schemas_equal(schema, tbl.schema)
    if len(tbl) == 0:
        return tbl

    # split, sometimes
    if (
        len(tbl) > 3
        and draw(st.booleans())
        and not HT_TEST_CONFIG["sc-61239_workaround"]
    ):
        n_splits = draw(st.integers(min_value=0, max_value=max(0, len(tbl) // 10)))
        if n_splits > 0:
            split_points = draw(splitss(n_splits=n_splits, max_value=len(tbl)))
            split_points = [0] + split_points + [len(tbl)]
            tbl = pa.concat_tables([tbl[st:sp] for st, sp in pairwise(split_points)])

    # pad, sometimes
    if draw(st.booleans()) and not HT_TEST_CONFIG["sc-61239_workaround"]:
        batches = tbl.to_batches()
        batch_to_pad = draw(st.integers(min_value=0, max_value=len(batches) - 1))
        batch_arrays = [
            pad_array(arr.to_numpy(zero_copy_only=(arr.type != pa.bool_())), draw)
            for arr in batches[batch_to_pad].columns
        ]
        batches[batch_to_pad] = pa.RecordBatch.from_arrays(
            batch_arrays, schema=tbl.schema
        )
        tbl = pa.Table.from_batches(batches)

    return tbl


class SOMADataFrameStateMachine(SOMAArrayStateMachine):

    def __init__(self) -> None:
        super().__init__()

    @initialize(data=st.data(), index_cols_and_schema=dataframe_schema())
    def setup(
        self,
        data: st.DataObject,
        index_cols_and_schema: tuple[Sequence[str], pa.Schema],
    ) -> None:
        self.index_column_names, self.schema = index_cols_and_schema
        self.domain = data.draw(  # TODO XXX: should be a ledger
            dataframe_domain(
                schema=self.schema, index_column_names=self.index_column_names
            )
        )
        super().setup(
            soma.DataFrame.create(
                self.uri,
                schema=self.schema,
                domain=self.domain,
                index_column_names=self.index_column_names,
                context=self.context,
                tiledb_timestamp=None,  # TODO: no time-travel for now
            )
        )
        self.domain = self.A.domain
        assert not self.A.closed
        assert self.A.mode == "w"
        assert schemas_equal(self.schema, self.A.schema, ignore_field_order=True)

        self.data_ledger = Ledger[ArrowTableLedgerEntry](
            initial_entry=ArrowTableLedgerEntry(
                data=self.schema.empty_table(),
                timestamp_ms=self.A.tiledb_timestamp_ms,
                name="initial entry",
                index_columns=self.index_column_names,
            ),
            allows_duplicates=False,
        )

    def _array_exists(
        uri: str, context: soma.SOMATileDBContext, tiledb_timestamp: int | None
    ) -> bool:
        return soma.DataFrame.exists(
            uri, context=context, tiledb_timestamp=tiledb_timestamp
        )

    def _array_open(self, *, mode: str, tiledb_timestamp: int | None = None) -> None:
        self.A = soma.DataFrame.open(
            self.uri, mode=mode, context=self.context, tiledb_timestamp=tiledb_timestamp
        )

    ##
    ## --- schema
    ##

    @precondition(lambda self: not self.closed)
    @invariant()
    def check_schema(self) -> None:
        assert isinstance(self.A, soma.DataFrame)
        assert self.A.soma_type == "SOMADataFrame"
        assert schemas_equal(self.schema, self.A.schema, ignore_field_order=True)
        assert sorted(self.schema.names) == sorted(self.A.keys())
        assert self.index_column_names == self.A.index_column_names

    ##
    ## --- domain
    ##

    @precondition(lambda self: not self.closed)
    @invariant()
    def check_domain(self) -> None:
        assert (
            self.A.domain == self.domain
        ), f"Unexpected domain in {self.A}: had {self.A.domain}, expected {self.domain}"

    @precondition(lambda self: self.closed or self.mode == "w")
    @rule(data=st.data())
    def expand_domain(self, data: st.DataObject) -> None:
        assert self.index_column_names == self.A.index_column_names
        new_domain = data.draw(
            dataframe_domain(
                schema=self.schema,
                index_column_names=self.index_column_names,
                current_domain=self.domain,
                max_domain=self.A.maxdomain,
            )
        )
        if self.closed:
            self._open(mode="w")
        assert self.mode == "w"

        self.A.change_domain(new_domain)
        self.domain = new_domain  # TODO XXX should be a ledger
        self._close()  # domain is committed upon close

    ##
    ## --- data
    ##

    @precondition(lambda self: not self.closed and self.mode == "r")
    @invariant()
    def check_read_all(self) -> None:
        timestamp_ms = self.A.tiledb_timestamp_ms
        sort_order = [(name, "ascending") for name in self.index_column_names]
        expected = (
            self.data_ledger.read(timestamp_ms=timestamp_ms)
            .to_table()
            .sort_by(sort_order)
        )
        found = self.A.read().concat().sort_by(sort_order)
        assert tables_equal(
            found, expected, equal_nan=True
        ), f"{found}\n is not equal to {expected}"

    @precondition(lambda self: not self.closed and self.mode == "r")
    @invariant()
    def check_count(self) -> None:
        expected = len(
            self.data_ledger.read(timestamp_ms=self.A.tiledb_timestamp_ms).to_table()
        )
        assert expected == self.A.count, "count mismatch"

    @precondition(lambda self: not self.closed and self.mode == "w")
    @precondition(
        lambda self: self.A.tiledb_timestamp_ms not in self.data_ledger.timestamps
    )  # only one write per timestamp until sc-61223 and sc-61226 are fixed
    @rule(data=st.data())
    def write(self, data: st.DataObject) -> None:
        df_tbl = data.draw(
            arrow_table(self.schema, self.index_column_names, self.domain, min_size=1)
        )
        fragments_before_write = get_entries(f"{self.uri}/__fragments")
        self.A.write(df_tbl)
        new_fragments = set(get_entries(f"{self.uri}/__fragments")) - set(
            fragments_before_write
        )
        assert len(new_fragments) == (
            1 if not HT_TEST_CONFIG["sc-61462_workaround"] else df_tbl[0].num_chunks
        )
        self.data_ledger.write(
            ArrowTableLedgerEntry(
                timestamp_ms=self.A.tiledb_timestamp_ms,
                name=new_fragments.pop(),
                data=df_tbl,
                index_columns=self.index_column_names,
            )
        )


TestSOMADataFrame = pytest.mark.usefixtures("make_tmp_dir")(
    SOMADataFrameStateMachine.TestCase
)
