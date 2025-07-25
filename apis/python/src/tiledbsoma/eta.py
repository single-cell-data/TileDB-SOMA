# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Estimated time helpers."""

import numpy as np


class Tracker:
    """Computes estimated time to completion for chunked writes."""

    chunk_percents: list[float]
    cumulative_seconds: list[float]

    def __init__(self) -> None:
        """Initialize."""
        self.chunk_percents = []
        self.cumulative_seconds = []

    def ingest_and_predict(self, chunk_percent: float, chunk_seconds: float) -> str:
        """Updates from most recent chunk percent-done and chunk
        completion-seconds, then does a linear regression
        on all chunks done so far and estimates time to completion.

        Args:
            chunk_percent:
                A percent done like 6.1 or 10.3.
            chunk_seconds:
                Number of seconds it took to do the current chunk operation.
        """
        self._ingest(chunk_percent, chunk_seconds)
        eta_seconds = self._predict()
        return self._format_seconds(eta_seconds)

    def _ingest(self, chunk_percent: float, chunk_seconds: float) -> None:
        """Takes the current percent done like 10.3 and current chunk seconds
        like 58.4 and grows an array of percent-dones and cumulative seconds.
        This means self.chunk_percents is a list of all the chunk_percent
        arguments from calling _ingest, while each self.chunk_seconds slot
        is the sum of all previous chunk_seconds arguments from calling _ingest.
        """
        if len(self.chunk_percents) == 0:
            self.chunk_percents = [chunk_percent]
            self.cumulative_seconds = [chunk_seconds]
        else:
            self.chunk_percents.append(chunk_percent)
            self.cumulative_seconds.append(self.cumulative_seconds[-1] + chunk_seconds)

    def _predict(self) -> float:
        """Does a linear regression on all chunks done so far and estimates t
        ime to completion.  Returns ETA seconds as a number.
        """
        # Linear regression where x is cumulative seconds and y is percent done.
        x = np.array(self.cumulative_seconds)
        y = np.array(self.chunk_percents)
        A = np.vstack([x, np.ones(len(x))]).T
        m, b = np.linalg.lstsq(A, y, rcond=None)[0]
        # Solve for x where y == 100.
        done_cumu_seconds = (100.0 - b) / m if m > 0 else 0

        return float(done_cumu_seconds) - self.cumulative_seconds[-1]

    def _format_seconds(self, seconds: float) -> str:
        """Formats the ETA seconds as a compact, human-readable string."""
        if seconds >= 86400:
            return "%.2f days" % (seconds / 86400)
        if seconds >= 3600:
            return "%.2f hours" % (seconds / 3600)
        if seconds >= 60:
            return "%.2f minutes" % (seconds / 60)
        return "%.2f seconds" % (seconds)

    def __str__(self) -> str:
        """User-friendly string."""
        return str(self.chunk_percents) + " " + str(self.cumulative_seconds)

    def __repr__(self) -> str:
        """User-friendly string."""
        return self.__str__()
