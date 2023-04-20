from tiledbsoma.profiler import measurements, profiler

"""A group of basic tests for the memory and performance profiler decorator"""


@profiler(param="tracemalloc")
def test_tracemalloc():
    my_list = []
    for i in range(100000):
        my_list += [0]
    return 1


@profiler(param="mem")
def test_mem():
    my_list = []
    for i in range(100000):
        my_list += [0]
    return 1


@profiler(param="timeit")
def test_timeit():
    [x for x in range(10)]
    return 1


def inner():
    [x for x in range(10000)]
    inner2()


def inner2():
    [x for x in range(100)]


@profiler(param="cprofile")
def test_cprofile():
    [x for x in range(10000)]
    inner()
    return 1


@profiler(param="line_profiler")
def test_line_profile():
    [x for x in range(10000)]
    inner()
    return 1


def test_profiles():
    test_tracemalloc()
    test_timeit()
    test_cprofile()
    test_line_profile()
    test_mem()
    print(str(measurements))
