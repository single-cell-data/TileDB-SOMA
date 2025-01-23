# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licenced under the MIT License.

from __future__ import annotations

from concurrent import futures
from typing import Iterator, TypeVar

_T = TypeVar("_T")


class EagerIterator(Iterator[_T]):
    def __init__(
        self,
        iterator: Iterator[_T],
        pool: futures.Executor | None = None,
    ):
        super().__init__()
        self.iterator = iterator
        self._pool = pool or futures.ThreadPoolExecutor()
        self._own_pool = pool is None
        self._preload_future = self._pool.submit(self.iterator.__next__)

    def __next__(self) -> _T:
        stopped = False
        try:
            if self._preload_future.cancel():
                # If `.cancel` returns True, cancellation was successful.
                # The self.iterator.__next__ call has not yet been started,
                # and will never be started, so we can compute next ourselves.
                # This prevents deadlocks if the thread pool is too small
                # and we can never create a preload thread.
                return next(self.iterator)
            # `.cancel` returned false, so the preload is already running.
            # Just wait for it.
            return self._preload_future.result()
        except StopIteration:
            self._cleanup()
            stopped = True
            raise
        finally:
            if not stopped:
                # If we have more to do, go for the next thing.
                self._preload_future = self._pool.submit(self.iterator.__next__)

    def _cleanup(self) -> None:
        if self._own_pool:
            self._pool.shutdown()

    def __del__(self) -> None:
        # Ensure the threadpool is cleaned up in the case where the
        # iterator is not exhausted. For more information on __del__:
        # https://docs.python.org/3/reference/datamodel.html#object.__del__
        self._cleanup()
        super_del = getattr(super(), "__del__", lambda: None)
        super_del()
