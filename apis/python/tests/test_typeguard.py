import typeguard

from tiledbsoma import SparseNDArray


def example_fn(n: int):
    """Intentionally broken type hints, "control" for typeguard test."""
    assert n == "abc"


def test_typeguard1():
    """Verify typeguard objects when type-hints are broken, at call-site below or in called function."""
    example_fn("abc")


@typeguard.typechecked
def test_typeguard2():
    """Verify typeguard objects when type-hints are broken, at call-site below or in called function."""
    example_fn("abc")


def test_typeguard_method():
    SparseNDArray.example_method("abc")
