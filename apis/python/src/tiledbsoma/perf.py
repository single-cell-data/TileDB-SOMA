# ----------------------------------------------------------------
# in some foo class:
#     s = self.startfunc({'foo': 'bar'}) # uri=self.uri implicit in there
#     self.endfunc(s)

import json
import time
from typing import Any, Dict, List, Optional


class Entry:
    """
    A simple container for a single labeled elapsed-time stat.
    """

    t1: float
    t2: Optional[float] = None

    def __init__(self, **kwargs: Any):
        self.t1 = time.time()
        self.t2 = None
        self.state = {}
        for k, v in kwargs.items():
            self.state[k] = v

    def finish(self) -> None:
        self.t2 = time.time()

    def as_line(self) -> str:
        line = ",".join([f"{k}={v}" for k, v in self.state.items()])
        if self.t2 is None:
            return line

        dt = "%.6f" % (self.t2 - self.t1)
        if len(self.state.items()) > 0:
            return f"seconds={dt},{line}"
        else:
            return f"seconds={dt}"

    def as_dict(self) -> Dict[str, Any]:
        retval = {}
        if self.t1 is not None and self.t2 is not None:
            retval["seconds"] = "%.6f" % (self.t2 - self.t1)
        for k, v in self.state.items():
            retval[k] = v
        return retval


# ----------------------------------------------------------------
class Tracker:
    """
    A simple container for a multiple labeled elapsed-time stats.
    """

    _entries: List[Entry]

    def __init__(self) -> None:
        self._entries = []

    def reset(self) -> None:
        self._entries = []

    def track(self, entry: Entry) -> None:
        self._entries.append(entry)

    def as_json(self) -> str:
        return json.dumps([entry.as_dict() for entry in self._entries], indent=2)

    def show_as_lines(self) -> None:
        """
        Displays the stats in key1=value1,key2=value2,... format which is splendid
        for when the key-list varies from one row to the next.
        """
        for entry in self._entries:
            print(entry.as_line())

    def show(self) -> None:
        """
        Displays the stats in pretty-printed tabular format, starting a new pretty-printed batch
        whenever the key-list changes.
        """
        homogeneous_batch = []
        last_header: List[str] = []
        printed_a_batch = False
        for entry in self._entries:
            row = entry.as_dict()
            # Build up a list of rows all having the same key-list, then format them
            # using column alignment.
            current_header = list(row.keys())
            if (last_header == []) or (current_header == last_header):
                homogeneous_batch.append(row)
            else:
                if printed_a_batch:
                    print()
                self._show_tabular_batch(homogeneous_batch, last_header)
                printed_a_batch = True
                homogeneous_batch = []
            last_header = current_header
        if len(homogeneous_batch) > 0:
            if printed_a_batch:
                print()
            self._show_tabular_batch(homogeneous_batch, last_header)

    def _show_tabular_batch(
        self, batch: List[Dict[str, Any]], header: List[str]
    ) -> None:
        """
        Helper method for `show`. Prints a list of perf entries all having the
        same key-list, using column alignment.
        """
        # Find the widths of the keys in the header line.
        assert len(batch) > 0
        max_widths = {}
        for key in header:
            max_widths[key] = len(key)

        # Find the widths of the values in the data lines, taking the max length
        # down columns.
        for row in batch:
            for key in header:
                value = str(row[key])
                max_widths[key] = max(max_widths[key], len(value))

        aligned_header = " ".join(["%-*s" % (max_widths[key], key) for key in header])
        print(aligned_header)

        for row in batch:
            aligned_data = " ".join(
                ["%-*s" % (max_widths[key], value) for key, value in row.items()]
            )
            print(aligned_data)


# ----------------------------------------------------------------
# Canonical singleton, although others can be instantiated if desired
tracker = Tracker()
