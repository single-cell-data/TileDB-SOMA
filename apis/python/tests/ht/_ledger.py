"""Ledger/log to model fragment/schema/metadata log entries."""

from __future__ import annotations

import pathlib
import re
from abc import ABCMeta, abstractmethod
from typing import Generic, Sequence, TypeVar

import numpy as np
import pandas as pd
import pyarrow as pa

from tests.ht._ht_util import schemas_equal


def get_entries(path: str | pathlib.Path) -> set[str]:
    """Get log entry names from directory, and return in canonical order.

    This is used to determine, by inspection, the names that TileDB Core
    assigns to log entries, such as write fragments, metadata tiles, etc.
    """
    dir = pathlib.Path(path)
    children = [p.relative_to(dir).as_posix() for p in dir.iterdir()]
    entries = [c for c in children if re.match(r"__[0-9]+_[0-9]+_[0-9a-fA-F]+", c)]
    entries.sort()
    return entries


LedgerEntryDataType = TypeVar("LedgerEntryDataType")


class LedgerEntry(Generic[LedgerEntryDataType], metaclass=ABCMeta):
    """An abstract consistent unit of written data, such as a fragment."""

    def __init__(self, timestamp_ms: int, name: str, data: LedgerEntryDataType) -> None:
        self.timestamp_ms: int = timestamp_ms
        self.name = name
        self.data = data

    @abstractmethod
    def consolidate_with(
        self, other: LedgerEntry[LedgerEntryDataType], allow_duplicates: bool
    ) -> LedgerEntry[LedgerEntryDataType]:
        pass


class Ledger(Generic[LedgerEntryDataType]):
    def __init__(
        self,
        initial_entry: LedgerEntry[LedgerEntryDataType],
        *,
        allows_duplicates: bool = False,
    ) -> None:
        self.entries: list[LedgerEntry[LedgerEntryDataType]] = [initial_entry]
        self.initial_entry = (
            initial_entry  # XXX: do we need this or can we use entries[0]?
        )
        self.allows_duplicates = allows_duplicates

        # multiple fragments with same timestamp are trouble. Just disallow for now
        # or we will have unstable tests.  See sc-61223.  When this is fixed, and there
        # is a stable read order, we could in principle reproduce that order (assuming
        # we know the UUID of each fragment, which requires sc-61226)
        self.timestamps = set()

    def __repr__(self) -> str:
        return (
            f"Ledger(n_entries={len(self.entries)}, "
            + f"allows_duplicates={self.allows_duplicates}):\n"
            + "\n".join(repr(f) for f in self.entries)
            + "\n"
        )

    def read(self, timestamp_ms: int) -> LedgerEntry[LedgerEntryDataType]:
        """Return a single ledger entry representing all writes <= timestamp"""
        assert len(self.entries) > 0
        entries_to_consolidate = sorted(
            filter(lambda f: f.timestamp_ms <= timestamp_ms, self.entries),
            key=lambda f: (f.timestamp_ms, f.name),
        )
        consolidated_result = entries_to_consolidate[0]
        for entry in entries_to_consolidate[1:]:
            consolidated_result = consolidated_result.consolidate_with(
                entry, self.allows_duplicates
            )
        return consolidated_result

    def write(self, entry: LedgerEntry[LedgerEntryDataType]) -> None:
        """Write new entry to the ledger."""
        assert entry.timestamp_ms >= 0
        assert type(entry) is type(self.initial_entry)

        if entry.timestamp_ms in self.timestamps:
            raise ValueError("Timestamp already written - may lead to unstable test.")

        self.entries.append(entry)
        self.timestamps.add(entry.timestamp_ms)


class ArrowTableLedgerEntry(LedgerEntry[pa.Table]):
    """Ledger entry based upon an Arrow Table."""

    def __init__(
        self, timestamp_ms: int, name: str, data: pa.Table, index_columns: Sequence[str]
    ) -> None:
        super().__init__(timestamp_ms, name, data)
        self.index_columns: list[str] = list(index_columns)

    def __repr__(self) -> str:
        return f"ArrowTableLedgerEntry(timestamp_ms={self.timestamp_ms}, index_columns={self.index_columns}):\n{self.data}"

    def consolidate_with(
        self, other: ArrowTableLedgerEntry, allow_duplicates: bool
    ) -> ArrowTableLedgerEntry:

        assert (self.timestamp_ms, self.name) < (other.timestamp_ms, other.name)
        assert schemas_equal(self.data.schema, other.data.schema)
        assert self.index_columns == other.index_columns

        earliest, latest = self, other

        if allow_duplicates:
            combined_table = pa.concat_tables((earliest, latest))
        else:
            if len(earliest.data) == 0:
                combined_table = latest.data
            elif len(latest.data) == 0:
                combined_table = earliest.data
            else:
                schema = self.data.schema
                latest_indexed = latest.to_pandas().set_index(self.index_columns)
                earliest_indexed = earliest.to_pandas().set_index(self.index_columns)

                # Table.from_pandas attempts to infer nulled values (e.g., NaN->null, NaT->null).
                # We do not want this behavior, so explicitly override it with `from_pandas=False`
                combined_table = pa.Table.from_pydict(
                    {
                        k: pa.array(v, from_pandas=False)
                        for k, v in combine_first(latest_indexed, earliest_indexed)
                        .reset_index()
                        .items()
                    },
                    schema=schema,
                )

        return ArrowTableLedgerEntry(
            timestamp_ms=latest.timestamp_ms,
            data=combined_table,
            name="consolidated",
            index_columns=self.index_columns,
        )

    def to_pandas(self) -> pd.DataFrame:
        return self.data.to_pandas(ignore_metadata=True)

    def to_table(self) -> pa.Table:
        return self.data


class ArrowTensorLedgerEntry(LedgerEntry[pa.Tensor]):
    """Ledger entry based upon an Arrow Tensor."""

    def __init__(
        self,
        timestamp_ms: int,
        name: str,
        data: pa.Tensor,
    ) -> None:
        super().__init__(timestamp_ms, name, data)

    def __repr__(self) -> str:
        return f"ArrowTensorLedgerEntry(timestamp_ms={self.timestamp_ms}):\n{self.data}"

    def consolidate_with(
        self, other: ArrowTensorLedgerEntry, allow_duplicates: bool
    ) -> ArrowTensorLedgerEntry:
        assert not allow_duplicates, "Unsupported"
        assert (self.timestamp_ms, self.name) < (other.timestamp_ms, other.name)
        return other

    def to_tensor(self) -> pa.Tensor:
        return self.data

    def to_numpy(self) -> np.ndarray:
        return self.data.to_numpy()


def combine_first(first: pd.DataFrame, second: pd.DataFrame) -> pd.DataFrame:
    """Combine dataframes - similar to pandas.DataFrame.combine_first,
    except fixes pandas#60128 and ignores NA values (they are copied as is).

    NB: the two dataframes MUST have the same structure, and we aren't
    too careful about checking for that.
    """

    assert first.columns.equals(second.columns)
    assert first.dtypes.equals(second.dtypes)
    assert first.index.nlevels == second.index.nlevels

    new_index = first.index.union(second.index)
    new_data = {}
    for col in first.columns:
        first_series = first[col]
        second_series = second[col]

        keep_second_index = second.index.difference(first.index)
        keep_first_index = first.index

        first_series = first_series.reindex(keep_first_index, copy=False)
        second_series = second_series.reindex(keep_second_index, copy=False)

        if first_series.dtype.kind == "M" and second_series.dtype.kind == "M":
            second_series = pd.to_datetime(second_series)

        combined_series = pd.concat([first_series, second_series])
        combined_series = combined_series.reindex(new_index, copy=False)

        new_data[col] = combined_series

    return pd.DataFrame(new_data, index=new_index)
