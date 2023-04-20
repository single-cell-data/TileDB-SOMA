import builtins
import cProfile
import pstats
import sys
import tracemalloc
from functools import wraps
from io import StringIO
from timeit import timeit
from typing import Any, Callable, Dict, Optional

import psutil
from line_profiler import LineProfiler
from typing_extensions import TypeVarTuple


class Measurement:
    """A class representing measurements and their labels"""

    _data: Dict[str, Any] = {}

    def __init__(self) -> None:
        self._data = {}

    @classmethod
    def read_single_stats(self, name: str) -> Optional[Any]:
        return self._data.get(name)

    @classmethod
    def update_single_stats(self, name: str, value: Any) -> None:
        self._data[name] = value

    @classmethod
    def update_stats_single_name(self, name: str, values) -> None:  # type: ignore
        for idx, stat in enumerate(values):
            self.update_single_stats(f"{name}{idx}", stat)

    @classmethod
    def __str__(self) -> str:
        result = ""
        for name, value in self._data.items():
            result += f"{name} {value}\n"
        return result


measurements: Measurement = Measurement()
Ts = TypeVarTuple("Xs")  # type: ignore


def stdout_string_wrapper(func) -> str:  # type: ignore
    """Function to convert and std output to string"""
    tmp = sys.stdout
    my_result = StringIO()
    sys.stdout = my_result
    func()
    sys.stdout = tmp
    return my_result.getvalue()


def profiler(
    param: Optional[str] = None,
) -> Callable[[Callable[[], Callable[[*Ts], None]]], Callable[[], None]]:  # type: ignore
    """The parameterized decorator for memory and performance profiling"""

    def measure(func: Callable[[], Callable[[*Ts], None]]) -> Callable[[], None]:  # type: ignore
        @wraps(func)
        def wrapper(*args, **kwargs) -> None:  # type: ignore
            if param == "timeit":
                print("timeit")
                time = timeit(lambda: func(*args, **kwargs))
                measurements.update_single_stats("timeit", time)

            # Collecting tracemalloc stats
            # each stat is stored with tracemalloc prefix as key"""

            if param == "tracemalloc":
                print("tracemalloc")
                tracemalloc.start()
                snapshot1 = tracemalloc.take_snapshot()
                func(*args, **kwargs)
                snapshot2 = tracemalloc.take_snapshot()
                top_stats = snapshot2.compare_to(snapshot1, "lineno")
                assert top_stats
                measurements.update_stats_single_name("tracemalloc", top_stats[:100])

            # Collecting cprofile stats
            # each stat is stored with cprofile prefix as key"""

            if param == "cprofile":
                cprofiler = cProfile.Profile()
                cprofiler.enable()
                func(*args, **kwargs)
                cprofiler.disable()
                stats = pstats.Stats(cprofiler).sort_stats("ncalls")
                print("cprofile output")
                stats_strings = stdout_string_wrapper(stats.print_stats)
                lines = stats_strings.splitlines()
                measurements.update_stats_single_name("cprofile", lines)
                print(stats_strings)

            # Collecting line_profiler stats
            # each stat is stored with line_profiler prefix as key"""

            if param == "line_profiler":
                prof = LineProfiler()
                builtins.__dict__["profile"] = prof
                func(*args, **kwargs)
                print("line_profiler:")
                prof.print_stats
                stats_strings = stdout_string_wrapper(prof.print_stats)
                lines = stats_strings.splitlines()
                assert lines
                measurements.update_stats_single_name("line_profiler", lines)

            # Collecting psutil stats
            # each stat is stored with cpu_times, virtual_memory, swap_memory prefix as key"""

            if param == "mem":
                func(*args, **kwargs)
                assert (
                    psutil.cpu_times()
                    and psutil.virtual_memory()
                    and psutil.swap_memory()
                )
                measurements.update_stats_single_name(
                    "cpu_times", [str(x) for x in psutil.cpu_times()]
                )
                measurements.update_stats_single_name(
                    "virtual_memory", [str(x) for x in psutil.virtual_memory()]
                )
                measurements.update_stats_single_name(
                    "swap_memory", [str(x) for x in psutil.swap_memory()]
                )

        return wrapper

    return measure
