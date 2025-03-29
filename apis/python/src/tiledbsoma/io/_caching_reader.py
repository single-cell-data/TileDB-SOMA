from __future__ import annotations

import ctypes
import io
import threading
from collections import OrderedDict
from typing import Any, cast

import pyarrow as pa
from typing_extensions import Buffer


class CachingReader:
    """Buffering reader which will maintain an LRU cache of previously read data, wrapping
    and presenting a file-like interface.

    Much like ``io.BufferedReader``, wraps file-like object. Class assumes that the underlying
    file is not changing concurrent with read. Any read to the file-like object will be a
    minimum of "cache block" size bytes.

    Primary use case is wrapping a VFS FileBuffer object, to improve performance of read
    operations by doing fewer, larger reads.

    Cached read bytes are stored as PyArrow arrays, primarily to benefit from the
    easy, minimal-copy concat and slice operations upon underlying buffers.
    """

    def __init__(
        self,
        file: Any,  # file-like object. Unfortunately, Python lacks a good typing signature for this concept.
        *,
        memory_budget: int = 128 * 1024**2,
        cache_block_size: int = 1024**2,
    ):
        if not file.readable():
            raise io.UnsupportedOperation("file must be readable")
        if memory_budget < cache_block_size:
            raise ValueError("memory_budget must be >= cache_block_size")

        self._file = file
        self._file_length = file.seek(0, io.SEEK_END)
        file.seek(0)
        self._offset = 0

        self._cache_block_size = cache_block_size
        self._max_cache_blocks = max(1, memory_budget // cache_block_size)
        self._n_blocks = (self._file_length + cache_block_size - 1) // cache_block_size

        self._cache_lock = threading.Lock()
        self._cached_blocks: list[bytes | None] = [None] * self._n_blocks
        self._cache_lru: OrderedDict[int, None] = OrderedDict()

    def _read_block(self, block: int) -> pa.UInt8Array:
        nbytes = min(
            self._cache_block_size,
            self._file_length - block * self._cache_block_size,
        )
        assert nbytes > 0
        buffer = pa.allocate_buffer(nbytes)
        ctypes.memset(buffer.address, 0, len(buffer))  # better safe than sorry
        self._file.seek(block * self._cache_block_size)
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
                i
                for i in range(start_block, end_block)
                if self._cached_blocks[i] is None
            ]
            for block in missing_blocks:
                self._cached_blocks[block] = self._read_block(block)

            requested_blocks = [
                self._cached_blocks[i] for i in range(start_block, end_block)
            ]

            self._mark_and_sweep_blocks(start_block, end_block)

        assert all(b is not None for b in requested_blocks)
        return cast(list[pa.UInt8Array], requested_blocks)

    def _mark_and_sweep_blocks(self, start: int, stop: int) -> None:
        for block_idx in range(start, stop):
            if block_idx in self._cache_lru:
                self._cache_lru.move_to_end(block_idx)

            elif len(self._cache_lru) < self._max_cache_blocks:
                self._cache_lru[block_idx] = None

            else:
                remove_idx = self._cache_lru.popitem(last=False)
                self._cache_lru[block_idx] = None
                self._cached_blocks[remove_idx[0]] = None

    def _reset_cache(self) -> None:
        with self._cache_lock:
            self._cached_blocks = [None] * self._n_blocks
            self._cache_lru.clear()

    def read(self, size: int = -1) -> bytes:
        if size == -1:
            size = self._file_length - self._offset
        blocks = self._load_cache(self._offset, self._offset + size)

        start = self._offset % self._cache_block_size
        end = start + size
        arr = pa.chunked_array(blocks)[start:end].combine_chunks()  # NB: copy
        assert arr.offset == 0
        b = arr.buffers()[1].to_pybytes()  # NB: copy
        self._offset += size
        return cast(bytes, b)

    def readinto(self, buf: Buffer) -> int:
        buf = memoryview(buf)
        size = len(buf)
        blocks = self._load_cache(self._offset, self._offset + size)
        start = self._offset % self._cache_block_size
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
        self._offset += dstidx
        return dstidx

    def close(self) -> None:
        self._reset_cache()
        self._file.close()

    def tell(self) -> int:
        return self._offset

    def seek(self, offset: int, whence: int = 0) -> int:
        if whence not in [io.SEEK_SET, io.SEEK_CUR, io.SEEK_END]:
            raise ValueError("Invalid whence value")

        if whence == io.SEEK_END:
            offset = self._file_length + offset
        elif whence == io.SEEK_CUR:
            offset += self._offset
        elif whence == io.SEEK_SET:
            offset = offset

        if offset < 0:
            raise OSError("seek() returned invalid position")

        self._offset = offset
        return offset
