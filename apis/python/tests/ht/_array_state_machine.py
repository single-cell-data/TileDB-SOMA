"""Hypothesis rule-based statemachine ABC for SOMAArray.

Intended to be specialized for SparseNDArray, et al.
"""

from __future__ import annotations

import math
import re
from abc import abstractmethod
from typing import Any, Literal, Protocol, Union

import numpy as np
import pyarrow as pa
from hypothesis import strategies as st
from hypothesis.stateful import RuleBasedStateMachine, invariant, precondition, rule
from somacore.options import OpenMode
from typing_extensions import TypeAlias

import tiledbsoma as soma

from tests.ht._ledger import Ledger, PyDictLedgerEntry

SOMAArray: TypeAlias = Union[soma.DataFrame, soma.SparseNDArray, soma.DenseNDArray]


class SOMAArrayStateMachine(RuleBasedStateMachine):
    """Abstract base class for a soma array Hypothesis state machine"""

    def __init__(self) -> None:
        super().__init__()
        self.context = self.TestCase.soma_tiledb_context
        self.closed: bool = True
        self.mode: Literal["r", "w"] | None = None
        self.A: SOMAArray | None = None
        self.uri = self.TestCase.tmp_path_factory.mktemp(f"{self.__class__.__name__}-").as_posix()

    def setup(self, A: SOMAArray) -> None:
        assert isinstance(A, (soma.DataFrame, soma.SparseNDArray, soma.DenseNDArray))
        assert A.mode == "w" and not A.closed
        self.A = A
        self.create_timestamp_ms = self.A.tiledb_timestamp_ms
        self.closed = self.A.closed
        self.mode = self.A.mode
        self.metadata_ledger = Ledger[PyDictLedgerEntry](
            initial_entry=PyDictLedgerEntry(
                data=dict(self.A.metadata),
                timestamp_ms=self.A.tiledb_timestamp_ms,
                name="initial entry",
            ),
            allows_duplicates=False,
        )
        self.pending_metadata: dict[str, Any] | None = None

    def teardown(self) -> None:
        if self.A is not None:
            if not self.closed:
                self.A.close()
            self.A = None

        super().teardown()

    @property
    def is_initialized(self) -> bool:
        return self.A is not None

    @abstractmethod
    def _array_exists(self, uri: str, context: soma.SOMAContext, tiledb_timestamp: int | None) -> bool:
        pass

    @abstractmethod
    def _array_open(self, mode: OpenMode, tiledb_timestamp: int | None = None) -> None:
        pass

    def _open(self, *, mode: OpenMode, tiledb_timestamp: int | None = None) -> None:
        assert self.A.closed
        tiledb_timestamp = None  # TODO/XXX: no time-travel for now. FIXME
        self._array_open(mode=mode, tiledb_timestamp=tiledb_timestamp)
        assert self.A is not None
        self.closed = False
        self.mode = mode

    def _close(self) -> None:
        assert not self.A.closed
        if self.pending_metadata is not None:
            self.metadata_ledger.write(PyDictLedgerEntry(self.A.tiledb_timestamp_ms, "", self.pending_metadata))
            self.pending_metadata = None

        self.A.close()
        self.closed = True
        self.mode = None

    ##
    # ---- Open/close state
    ##

    @precondition(lambda self: self.is_initialized)
    @invariant()
    def check_exists(self) -> None:
        assert self._array_exists(self.uri, self.context, None)

    @precondition(lambda self: self.is_initialized)
    @invariant()
    def check_mode(self) -> None:
        assert self.closed or self.mode == self.A.mode

    @precondition(lambda self: self.is_initialized)
    @invariant()
    def check_closed(self) -> None:
        assert self.closed == self.A.closed

    @precondition(lambda self: not self.closed)
    @rule()
    def close(self) -> None:
        self._close()

    @precondition(lambda self: self.closed)
    @rule(mode=st.sampled_from(["r", "w", "d"]))
    def open(self, mode: OpenMode) -> None:
        # TODO: time travel
        self._open(mode=mode)

    @precondition(lambda self: self.is_initialized)
    @rule(mode=st.sampled_from(["r", "w", "d"]))
    def reopen(self, mode: OpenMode) -> None:
        self.A.reopen(
            mode,
            tiledb_timestamp=None,  # no time-travel for now
        )
        assert self.A.mode == mode and not self.A.closed
        self.closed = False
        self.mode = mode

    ##
    # --- metadata
    ##
    METADATA_KEY_ALPHABET = st.characters(codec="utf-8", exclude_characters=["\x00"])
    METADATA_KEYS = st.text(min_size=0, max_size=4096, alphabet=METADATA_KEY_ALPHABET).filter(
        lambda k: not k.startswith("soma_"),
    )
    METADATA_VALUE_ALPHABET = st.characters(codec="utf-8", exclude_characters=["\x00"])
    METADATA_VALUES = st.one_of(
        st.text(alphabet=METADATA_VALUE_ALPHABET, min_size=0)
        | st.integers(min_value=np.iinfo(np.int64).min, max_value=np.iinfo(np.int64).max)
        | st.floats(),
    )
    IGNORE_KEYS = re.compile(r"^soma_.*$")

    @classmethod
    def filter_metadata(cls, d: dict[str, Any]) -> dict[str, Any]:
        """Apply the "ignore" regex to dict keys, returning the filtered dict."""
        return {k: v for k, v in d.items() if not cls.IGNORE_KEYS.match(k)}

    @precondition(lambda self: not self.closed)
    @invariant()
    def check_metadata(self) -> None:
        array_metadata = self.filter_metadata(dict(self.A.metadata))
        expected_metadata = self.filter_metadata(
            (
                self.metadata_ledger.read(timestamp_ms=self.A.tiledb_timestamp_ms).to_dict()
                if self.pending_metadata is None
                else self.pending_metadata
            ),
        )
        assert set(array_metadata.keys()) == set(expected_metadata.keys())
        for k in array_metadata:
            if isinstance(array_metadata[k], float) and math.isnan(array_metadata[k]):
                assert math.isnan(expected_metadata[k])
                continue
            assert array_metadata[k] == expected_metadata[k]

    @precondition(lambda self: not self.closed and self.mode == "w" and len(self.A.metadata) < 100)
    @rule(k=METADATA_KEYS, v=METADATA_VALUES)
    def set_metadata(self, k: str, v: str | int | float) -> None:
        self.A.metadata[k] = v
        if self.pending_metadata is None:
            self.pending_metadata = self.metadata_ledger.read(self.A.tiledb_timestamp_ms).to_dict()
        self.pending_metadata[k] = v

    @precondition(lambda self: not self.closed and self.mode == "w" and len(self.filter_metadata(self.A.metadata)))
    @precondition(lambda self: not self.closed)
    @rule(data=st.data())
    def del_metadata(self, data: st.DataObject) -> None:
        if self.pending_metadata is None:
            self.pending_metadata = self.metadata_ledger.read(self.A.tiledb_timestamp_ms).to_dict()

        k = data.draw(st.sampled_from(sorted(list(self.filter_metadata(self.pending_metadata).keys()))))
        del self.A.metadata[k]
        del self.pending_metadata[k]


class ShapesFactory(Protocol):
    """Factory for a strategy returning ndarray shape."""

    def __call__(
        self,
        *,
        min_shape: tuple[int, ...] | None = None,
        max_shape: tuple[int, ...] | None = None,
    ) -> st.SearchStrategy[tuple[int | None, ...]]: ...


class SOMANDArrayStateMachine(SOMAArrayStateMachine):
    """Abstract base class for NDArray Hypothesis state machine."""

    def __init__(self, shapes_factory: ShapesFactory) -> None:
        super().__init__()
        self.shapes_factory = shapes_factory

    def setup(self, type: pa.DataType, shape: tuple[int, ...], array) -> None:
        super().setup(array)
        self.type = type
        self.schema = pa.schema(
            [pa.field(f"soma_dim_{n}", pa.int64(), nullable=False) for n in range(len(shape))]
            + [pa.field("soma_data", self.type, nullable=False)],
        )
        assert all((shape[i] or 1) == self.A.shape[i] for i in range(len(shape)))
        assert self.schema == self.A.schema
        self.shape = tuple((shape[i] or 1) for i in range(len(shape)))  # XXX TODO: shape should be a ledger

    ##
    # --- schema
    ##

    @precondition(lambda self: not self.closed)
    @invariant()
    def check_schema(self) -> None:
        schema = self.A.schema
        assert len(schema.types) == len(self.shape) + 1
        assert schema.field("soma_data").type == self.type
        for idx in range(len(self.shape)):
            assert schema.names[idx] == f"soma_dim_{idx}"
            assert schema.types[idx] == pa.int64()
            assert schema.field(f"soma_dim_{idx}").type == pa.int64()
        assert self.A.schema == self.schema

    ##
    # --- shape
    ##

    @precondition(lambda self: not self.closed)
    @invariant()
    def check_shape(self) -> None:
        assert hasattr(self.A, "shape")  # sc-61123
        assert self.A.shape == tuple((s or 1) for s in self.shape), (
            f"Unexpected shape in {self.A}: had {self.A.shape}, expected {self.shape}"
        )
        assert self.A.ndim == len(self.shape)

    @precondition(lambda self: self.closed or self.mode == "w")
    @rule(data=st.data())
    def expand_shape(self, data: st.DataObject) -> None:
        # Always re-open at latest. Without this, there is a good chance we will end up with a
        # schema fragment with the same timestamp. Concrete case: a write has been done to an
        # array that has a dict field, which triggers an automatic schema evolution.
        if not self.closed:
            self._close()
        self._open(mode="w")
        assert self.mode == "w"

        new_shape = data.draw(self.shapes_factory(min_shape=self.shape, max_shape=self.A.maxshape))
        self.A.resize(new_shape)
        self.shape = new_shape
        self._close()  # resize is committed upon close
