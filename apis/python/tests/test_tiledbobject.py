import pytest

import tiledbsoma as soma


# Checking that objects _do_ exist is already done (thoroughly) in other
# tests. Here we check the negative cases.
#
# Specifically, we're checking these don't raise, and all return False.


@pytest.mark.parametrize(
    "uri",
    ["/nonesuch/no/nope/never/ever", "foo://bar", "s3://@@@@ILLEGAL@@@@"],
)
@pytest.mark.parametrize(
    "somaclass",
    [
        soma.DataFrame,
        soma.SparseNDArray,
        soma.DenseNDArray,
        soma.Collection,
        soma.Measurement,
        soma.Experiment,
    ],
)
def test_tiledbobject_exists(uri, somaclass):
    assert not somaclass.exists(uri)
