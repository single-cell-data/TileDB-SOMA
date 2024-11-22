import threading
import unittest
from concurrent import futures
from unittest import mock

from tiledbsoma._eager_iter import EagerIterator


class EagerIterTest(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.kiddie_pool = futures.ThreadPoolExecutor(1)
        """Tiny thread pool for testing."""
        self.verify_pool = futures.ThreadPoolExecutor(1)
        """Separate thread pool so verification is not blocked."""

    def tearDown(self):
        self.verify_pool.shutdown(wait=False)
        self.kiddie_pool.shutdown(wait=False)
        super().tearDown()

    def test_thread_starvation(self):
        sem = threading.Semaphore()
        try:
            # Monopolize the threadpool.
            sem.acquire()
            self.kiddie_pool.submit(sem.acquire)
            eager = EagerIterator(iter("abc"), pool=self.kiddie_pool)
            got_a = self.verify_pool.submit(lambda: next(eager))
            self.assertEqual("a", got_a.result(0.1))
            got_b = self.verify_pool.submit(lambda: next(eager))
            self.assertEqual("b", got_b.result(0.1))
            got_c = self.verify_pool.submit(lambda: next(eager))
            self.assertEqual("c", got_c.result(0.1))
            with self.assertRaises(StopIteration):
                self.verify_pool.submit(lambda: next(eager)).result(0.1)
        finally:
            sem.release()

    def test_nesting(self):
        inner = EagerIterator(iter("abc"), pool=self.kiddie_pool)
        outer = EagerIterator(inner, pool=self.kiddie_pool)
        self.assertEqual(
            "a, b, c", self.verify_pool.submit(", ".join, outer).result(0.1)
        )

    def test_exceptions(self):
        flaky = mock.MagicMock()
        flaky.__next__.side_effect = [1, 2, ValueError(), 3, 4]

        eager_flaky = EagerIterator(flaky, pool=self.kiddie_pool)
        got_1 = self.verify_pool.submit(lambda: next(eager_flaky))
        self.assertEqual(1, got_1.result(0.1))
        got_2 = self.verify_pool.submit(lambda: next(eager_flaky))
        self.assertEqual(2, got_2.result(0.1))
        with self.assertRaises(ValueError):
            self.verify_pool.submit(lambda: next(eager_flaky)).result(0.1)
        got_3 = self.verify_pool.submit(lambda: next(eager_flaky))
        self.assertEqual(3, got_3.result(0.1))
        got_4 = self.verify_pool.submit(lambda: next(eager_flaky))
        self.assertEqual(4, got_4.result(0.1))
        for _ in range(5):
            with self.assertRaises(StopIteration):
                self.verify_pool.submit(lambda: next(eager_flaky)).result(0.1)
