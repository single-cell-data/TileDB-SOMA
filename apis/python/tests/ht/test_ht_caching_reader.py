import io

from hypothesis import strategies as st
from hypothesis.stateful import (
    RuleBasedStateMachine,
    initialize,
    invariant,
    precondition,
    rule,
)

from tiledbsoma.io._caching_reader import CachingReader

from .._util import TESTDATA


class CachingReaderStateMachine(RuleBasedStateMachine):

    def __init__(self) -> None:
        super().__init__()
        self._fname = TESTDATA / "pbmc3k.h5ad"

    @initialize()
    def _setup(self):
        self._closed = True

        # the two file handles that are compared
        self.raw = None
        self.cached = None

    def teardown(self) -> None:
        if not self._closed:
            self.raw.close()
            self.cached.close()
            self.raw = None
            self.cache = None

        super().teardown()

    @precondition(lambda self: self.raw is not None and self._closed)
    @invariant()
    def check_closed_state(self) -> None:
        """Invariants that should be true when we have a file handle (in any state)"""
        assert self.raw is not None and self.cached is not None
        assert self.raw.closed == self.cached.closed

    @precondition(lambda self: self.raw is not None and not self._closed)
    @invariant()
    def check_open_state(self) -> None:
        """Invariants that should be true if we have an open handle"""
        assert self.raw is not None and not self._closed
        assert self.raw.closed == self.cached.closed
        assert self.raw.readable() == self.cached.readable()
        assert self.raw.tell() == self.cached.tell()

    @precondition(lambda self: self._closed)
    @rule()
    def open(self) -> None:
        self._closed = False
        self.raw = open(self._fname, "rb")
        self.cached = CachingReader(open(self._fname, "rb"))

    @precondition(lambda self: not self._closed)
    @rule()
    def close(self) -> None:
        self._closed = True
        self.raw.close()
        self.cached.close()
        assert self.raw.closed == self.cached.closed

    @precondition(lambda self: not self._closed)
    @rule(size=st.integers(min_value=-1, max_value=8192))
    def read(self, size: int) -> None:
        raw_buf = self.raw.read(size)
        cached_buf = self.cached.read(size)

        assert raw_buf == cached_buf
        assert self.raw.tell() == self.cached.tell()

    @precondition(lambda self: not self._closed)
    @rule(size=st.integers(min_value=0, max_value=8192))
    def readinto(self, size: int) -> None:
        raw_buf = bytearray(size)
        cached_buf = bytearray(size)

        raw_nbytes = self.raw.readinto(raw_buf)
        cached_nbytes = self.cached.readinto(cached_buf)

        assert raw_nbytes == cached_nbytes
        assert raw_buf == cached_buf
        assert self.raw.tell() == self.cached.tell()

    @precondition(lambda self: not self._closed)
    @rule(
        offset=st.integers(min_value=-8192, max_value=8192),
        whence=st.sampled_from([io.SEEK_SET, io.SEEK_CUR, io.SEEK_END]),
    )
    def seek(self, offset: int, whence: int) -> None:
        raw_offset: int | None = None
        cached_offset: int | None = None
        try:
            raw_err = False
            raw_offset = self.raw.seek(offset, whence)
        except Exception:
            raw_err = True

        try:
            cached_err = False
            cached_offset = self.cached.seek(offset, whence)
        except Exception:
            cached_err = True

        assert raw_offset == cached_offset
        assert self.raw.tell() == self.cached.tell()
        assert raw_err == cached_err


TestCachingReader = CachingReaderStateMachine.TestCase
