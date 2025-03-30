from __future__ import annotations

import ctypes
import io
import threading
from collections import OrderedDict
from typing import Any, cast

import attrs
import pyarrow as pa
from typing_extensions import Buffer


@attrs.define
class CacheStats:
    block_idx: int
    hit: int = 0
    miss: int = 0


class CachingReader:
    """Buffering reader which maintains a LRU cache of previously read data, wrapping
    and presenting a file-like interface. Typical use is similar to ``io.BufferedReader``:

        f = CachingReader(open("a file", mode="rb"), memory_budget=64*1024**2, cache_block_size=1024**2)

    The LRU cache is user-defined by:
        * memory budget - total cache size in bytes
        * block_size - the unit of cached data. Also the read size

    Instances present as a file-like object, and all public methods conform to the normal Python
    ``io`` package semantics. Most, but not all methods required for io.RawIOBase (read-only) are
    implemented.

    Class assumes that the underlying file is not changing concurrent with read.

    Primary use case is wrapping a VFS FileBuffer object, to improve performance of read
    operations by doing fewer, larger reads.

    The implemtation stores cached data as PyArrow arrays, primarily to benefit from the
    easy, minimal-copy concat and slice operations upon underlying buffers.
    """

    def __init__(
        self,
        file: Any,  # file-like object. Unfortunately, Python lacks a good typing signature for this concept.
        *,
        memory_budget: int = 64 * 1024**2,
        cache_block_size: int = 1024**2,
    ):
        if not file.readable():
            raise io.UnsupportedOperation("file must be readable")
        if memory_budget < cache_block_size:
            raise ValueError("memory_budget must be >= cache_block_size")

        self._file = file
        self._file_length = file.seek(0, io.SEEK_END)
        file.seek(0)
        self._pos = 0

        self._cache_block_size = cache_block_size
        self._max_cache_blocks = max(1, memory_budget // cache_block_size)
        _n_blocks = (self._file_length + cache_block_size - 1) // cache_block_size

        self._cache_lock = threading.Lock()
        self._cache: OrderedDict[int, pa.UInt8Array] = OrderedDict()
        self._cache_stats: list[CacheStats] = [
            CacheStats(block_idx) for block_idx in range(_n_blocks)
        ]

    def _read_block(self, block_idx: int) -> pa.UInt8Array:
        nbytes = min(
            self._cache_block_size,
            self._file_length - block_idx * self._cache_block_size,
        )
        assert nbytes > 0
        buffer = pa.allocate_buffer(nbytes)
        ctypes.memset(buffer.address, 0, len(buffer))  # better safe than sorry
        self._file.seek(block_idx * self._cache_block_size)
        bytes_read = self._file.readinto(memoryview(buffer))
        assert nbytes == bytes_read == len(buffer)
        a = pa.UInt8Array.from_buffers(pa.uint8(), len(buffer), [None, buffer])
        return a

    def _load_cache(self, start: int, end: int) -> list[pa.UInt8Array]:
        end = min(end, self._file_length)
        start_block = start // self._cache_block_size
        end_block = (end + self._cache_block_size - 1) // self._cache_block_size
        with self._cache_lock:
            missing_blocks = [
                i for i in range(start_block, end_block) if i not in self._cache
            ]
            for block_idx in missing_blocks:
                self._cache_stats[block_idx].miss += 1
                self._cache[block_idx] = self._read_block(block_idx)

            requested_blocks = [
                self._cache[block_idx] for block_idx in range(start_block, end_block)
            ]

            self._mark_and_sweep_blocks(start_block, end_block)

        return requested_blocks

    def _mark_and_sweep_blocks(self, start: int, stop: int) -> None:
        for block_idx in range(start, stop):
            self._cache.move_to_end(block_idx)
            self._cache_stats[block_idx].hit += 1

        for i in range(max(0, len(self._cache) - self._max_cache_blocks)):
            self._cache.popitem(last=False)

    def _reset_cache(self) -> None:
        with self._cache_lock:
            self._cache.clear()

    def read(self, size: int = -1) -> bytes:
        if size is None:
            size = -1  # type: ignore[unreachable]
        if size < 0:
            size = self._file_length - self._pos
        if size == 0:
            return b""

        blocks = self._load_cache(self._pos, self._pos + size)

        start = self._pos % self._cache_block_size
        end = start + size
        arr = pa.chunked_array(blocks)[start:end].combine_chunks()  # NB: copy
        assert arr.offset == 0
        b = arr.buffers()[1].to_pybytes()  # NB: copy
        self._pos += len(b)
        return cast(bytes, b)

    def readinto(self, buf: Buffer | memoryview) -> int | None:
        if not isinstance(buf, memoryview):
            buf = memoryview(buf)
        if buf.nbytes == 0:
            return 0
        buf = buf.cast("B")

        size = buf.nbytes
        blocks = self._load_cache(self._pos, self._pos + size)
        start = self._pos % self._cache_block_size
        end = start + size
        carr = pa.chunked_array(blocks)[start:end]

        dst = ctypes.addressof(ctypes.c_char.from_buffer(buf))
        dstidx = 0
        for c in carr.chunks:
            src = c.buffers()[1].address + c.offset
            count = len(c)
            ctypes.memmove(dst + dstidx, src, count)
            dstidx += count

        assert dstidx == len(carr)
        self._pos += dstidx
        return dstidx

    def close(self) -> None:
        self._reset_cache()
        self._file.close()

    def tell(self) -> int:
        return self._pos

    def seek(self, offset: int, whence: int = 0) -> int:
        if whence not in [io.SEEK_SET, io.SEEK_CUR, io.SEEK_END]:
            raise ValueError("Invalid whence value")

        if whence == io.SEEK_END:
            offset = self._file_length + offset
        elif whence == io.SEEK_CUR:
            offset += self._pos
        elif whence == io.SEEK_SET:
            offset = offset

        if offset < 0:
            raise OSError("seek() returned invalid position")

        self._pos = offset
        return offset

    @property
    def closed(self) -> bool:
        return bool(self._file.closed)

    def readable(self) -> bool:
        return bool(self._file.readable())

    def seekable(self) -> bool:
        return bool(self._file.seekable())
