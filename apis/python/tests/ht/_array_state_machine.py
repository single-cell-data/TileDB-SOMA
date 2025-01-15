"""Hypothesis rule-based statemachine ABC for SOMAArray.

Intended to be specialized for SparseNDArray, et al.
"""

from __future__ import annotations

import re
from abc import abstractmethod
from typing import Any, Literal, Protocol, Union

import numpy as np
import pyarrow as pa
from hypothesis import strategies as st
from hypothesis.stateful import RuleBasedStateMachine, invariant, precondition, rule
from typing_extensions import TypeAlias

import tiledbsoma as soma

from tests.ht._ht_test_config import HT_TEST_CONFIG

SOMAArray: TypeAlias = Union[soma.DataFrame, soma.SparseNDArray, soma.DenseNDArray]


class SOMAArrayStateMachine(RuleBasedStateMachine):
    """Abstract base class for a soma array Hypothesis state machine"""

    def __init__(self) -> None:
        super().__init__()
        self.context = soma.SOMATileDBContext()
        self.closed: bool = True
        self.mode: Literal["r", "w"] | None = None
        self.A: SOMAArray | None = None
        self.uri = self.TestCase.tmp_path_factory.mktemp(
            f"{self.__class__.__name__}-"
        ).as_posix()
        self.metadata: dict[str, Any] = (
            {}
        )  # XXX TODO: should be a ledger to allow for time travel
        self.initial_metadata_keys: set[str] = set()

    def setup(self, A: SOMAArray) -> None:
        assert isinstance(A, (soma.DataFrame, soma.SparseNDArray, soma.DenseNDArray))
        assert A.mode == "w" and not A.closed
        self.A = A
        self.create_timestamp_ms = self.A.tiledb_timestamp_ms
        self.closed = self.A.closed
        self.mode = self.A.mode
        self.metadata = dict(self.A.metadata)
        self.initial_metadata_keys = set(self.metadata)

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
    def _array_exists(
        uri: str, context: soma.SOMATileDBContext, tiledb_timestamp: int | None
    ) -> bool:
        pass

    @abstractmethod
    def _array_open(self, mode: str) -> None:
        pass

    def _open(self, *, mode: str, tiledb_timestamp: int | None = None) -> None:
        assert self.A.closed
        tiledb_timestamp = None  # TODO/XXX: no time-travel for now. FIXME
        self._array_open(mode=mode, tiledb_timestamp=tiledb_timestamp)
        assert self.A is not None
        self.closed = False
        self.mode = mode

    def _close(self) -> None:
        assert not self.A.closed
        self.A.close()
        self.closed = True
        self.mode = None

    @abstractmethod
    def _reopen(self, mode: str) -> None:
        pass

    ##
    ## ---- Open/close state
    ##

    @precondition(lambda self: self.is_initialized)
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
    @rule(mode=st.sampled_from(["r", "w"]))
    def open(self, mode: str) -> None:
        # TODO: time travel
        self._open(mode=mode)

    @precondition(
        lambda self: not HT_TEST_CONFIG["sc-61123_workaround"]
    )  # TODO: this entire rule disabled until sc-61123 fixed.
    @precondition(lambda self: not self.closed)
    @precondition(
        lambda self: not HT_TEST_CONFIG["sc-61118_workaround"] or self.mode != "w"
    )  # TODO - fails due to loss of metadata on reopen from w->r. See sc-61118. Remove when fixed.
    @rule(mode=st.sampled_from(["r", "w"]))
    def reopen(self, mode: str) -> None:
        assert not self.A.closed
        assert not self.closed
        assert self.mode is not None
        self.A = self.A.reopen(
            mode,
            tiledb_timestamp=None,  # no time-travel for now
        )
        self.mode = mode
        assert self.A.mode == mode and not self.A.closed

    ##
    ## --- metadata
    ##
    # TODO: sc-61092 causes SOMA to fail on writing a metadata value with a non-ASCII codepoint.
    # TODO: due to sc-61093, zero length bytes and strings are mishandled (not written correctly). Remove the `min_size` when fixed.
    # TODO: due to sc-61094, strings containing a zero code point also fail.

    METADATA_KEY_ALPHABET = (
        st.characters(codec="utf-8", exclude_characters=["\x00"])
        if HT_TEST_CONFIG["sc-61094_workaround"]
        else st.characters(codec="utf-8")
    )
    METADATA_KEYS = st.text(min_size=1, max_size=4096, alphabet=METADATA_KEY_ALPHABET)

    METADATA_VALUE_ALPHABET = (
        st.characters(codec="ascii", exclude_characters=["\x00"])
        if (
            HT_TEST_CONFIG["sc-61092_workaround"]
            or HT_TEST_CONFIG["sc-61094_workaround"]
        )
        else st.characters(codepoint="utf-8")
    )
    METADATA_VALUES = st.one_of(
        st.text(
            alphabet=METADATA_VALUE_ALPHABET,
            min_size=1 if HT_TEST_CONFIG["sc-61093_workaround"] else 0,
        )
        | st.integers(
            min_value=np.iinfo(np.int64).min, max_value=np.iinfo(np.int64).max
        )
        | st.floats(
            allow_nan=False
        )  # FIXME: disabled NaNs make assertions easier (they are supported)
    )

    IGNORE_KEYS = re.compile(r"^soma_dim_[0-9]+_domain_(upper|lower)$")

    @precondition(lambda self: not self.closed)
    @invariant()
    def check_metadata(self) -> None:
        # Prior to tiledbsoma 1.16, the "used domain" keys were still included. Ignore them.
        # TODO: we could generalize this by removing _all_ keys that are reserved soma_* keys.
        array_metadata = {
            k: v for k, v in self.A.metadata.items() if not self.IGNORE_KEYS.match(k)
        }
        assert array_metadata == self.metadata

    @precondition(
        lambda self: not self.closed and self.mode == "w" and len(self.metadata) < 100
    )
    @rule(k=METADATA_KEYS, v=METADATA_VALUES)
    def set_metadata(self, k: str, v: str | int | float) -> None:
        self.metadata[k] = v
        self.A.metadata[k] = v

    @precondition(
        lambda self: not self.closed
        and self.mode == "w"
        and len(self.metadata) > len(self.initial_metadata_keys)
    )
    @precondition(lambda self: not self.closed)
    @rule(data=st.data())
    def del_metadata(self, data: st.DataObject) -> None:
        k = data.draw(
            st.sampled_from(
                [
                    kn
                    for kn in self.metadata.keys()
                    if kn not in self.initial_metadata_keys
                ]
            )
        )
        del self.metadata[k]
        del self.A.metadata[k]


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

    def setup(self, type, shape, array) -> None:
        super().setup(array)
        self.type = type
        self.schema = pa.schema(
            [
                pa.field(f"soma_dim_{n}", pa.int64(), nullable=False)
                for n in range(len(shape))
            ]
            + [pa.field("soma_data", self.type, nullable=False)]
        )
        assert all((shape[i] or 1) == self.A.shape[i] for i in range(len(shape)))
        assert self.schema == self.A.schema
        self.shape = tuple(
            (shape[i] or 1) for i in range(len(shape))
        )  # XXX TODO: shape should be a ledger

    ##
    ## --- schema
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
    ## --- shape
    ##

    @precondition(lambda self: not self.closed)
    @invariant()
    def check_shape(self) -> None:
        assert hasattr(self.A, "shape")  # sc-61123
        assert self.A.shape == tuple(
            (s or 1) for s in self.shape
        ), f"Unexpected shape in {self.A}: had {self.A.shape}, expected {self.shape}"
        assert self.A.ndim == len(self.shape)

    @precondition(lambda self: self.closed or self.mode == "w")
    @rule(data=st.data())
    def expand_shape(self, data: st.DataObject) -> None:
        if self.closed:
            self._open(mode="w")
        assert self.mode == "w"
        new_shape = data.draw(
            self.shapes_factory(min_shape=self.shape, max_shape=self.A.maxshape)
        )
        self.A.resize(new_shape)
        self.shape = new_shape
        self._close()  # resize is committed upon close
