from functools import wraps
from timeit import timeit
import sys
import tracemalloc
import cProfile, pstats
from line_profiler import LineProfiler
from io import StringIO
import builtins
import psutil

"""A class representing measurements and their labels"""
class Measurement:
    _data = {}

    def __init__(self):
        self._data = {}

    @classmethod
    def readSingleStats(self, name):
        return self._data.get(name)

    @classmethod
    def updateSingleStats(self, name, value):
        self._data[name] = value

    @classmethod
    def updateStatsSingleName(self, name, values):
        for idx, stat in enumerate(values):
            self.updateSingleStats(f"{name}{idx}", stat)

    @classmethod
    def __str__(self):
        result = ""
        for name, value in self._data.items():
            result += f"{name} {value}\n"
        return result


measurements = Measurement()

"""Function to convert and std output to string"""


def stdout_string_wrapper(func):
    tmp = sys.stdout
    my_result = StringIO()
    sys.stdout = my_result
    func()
    sys.stdout = tmp
    return my_result.getvalue()


"""The parameterized decorator for memory and performance profiling"""


def profiler(param=None):
    def measure(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if param == "timeit":
                print("timeit")
                time = timeit(
                    lambda: func(*args, **kwargs)
                )
                measurements.updateSingleStats("timeit", time)

            """Collecting tracemalloc stats
            each stat is stored with tracemalloc prefix as key"""

            if param == "tacemalloc":
                print("tacemalloc")
                tracemalloc.start()
                snapshot1 = tracemalloc.take_snapshot()
                func(*args, **kwargs)
                snapshot2 = tracemalloc.take_snapshot()
                top_stats = snapshot2.compare_to(snapshot1, 'lineno')
                measurements.updateStatsSingleName("tracemalloc", top_stats[:100])

            """Collecting cprofile stats
            each stat is stored with cprofile prefix as key"""

            if param == "cprofile":
                cprofiler = cProfile.Profile()
                cprofiler.enable()
                func(*args, **kwargs)
                cprofiler.disable()
                stats = pstats.Stats(cprofiler).sort_stats('ncalls')
                print("cprofile output")
                stats_strings = stdout_string_wrapper(stats.print_stats)
                lines = stats_strings.splitlines()
                measurements.updateStatsSingleName("cprofile", lines)
                print(stats_strings)

            """Collecting line_profiler stats
            each stat is stored with line_profiler prefix as key"""

            if param == "line_profiler":
                prof = LineProfiler()
                builtins.__dict__['profile'] = prof
                func(*args, **kwargs)
                print("line_profiler:")
                stats_strings = stdout_string_wrapper(prof.print_stats)
                lines = stats_strings.splitlines()
                measurements.updateStatsSingleName("line_profiler", lines)

            """Collecting psutil stats
            each stat is stored with cpu_times, virtual_memory, swap_memory prefix as key"""

            if param == "mem":
                func(*args, **kwargs)
                measurements.updateStatsSingleName("cpu_times", [str(x) for x in psutil.cpu_times()])
                measurements.updateStatsSingleName("virtual_memory", [str(x) for x in psutil.virtual_memory()])
                measurements.updateStatsSingleName("swap_memory", [str(x) for x in psutil.swap_memory()])

        return wrapper

    return measure
