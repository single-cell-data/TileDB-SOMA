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
def test_tiledbobject_exists_nonexistent_path(uri, somaclass):
    assert not somaclass.exists(uri)


@pytest.mark.parametrize("uri", [b"/dev/null", 123.45, ["path"], {}])
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
def test_tiledbobject_exists_invalid_uri_type(uri, somaclass):
    with pytest.raises(TypeError):
        somaclass.exists(uri)
