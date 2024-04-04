import pyarrow as pa
import pytest

import tiledbsoma as soma

from tests._util import raises_no_typeguard

# Checking that objects _do_ exist is already done (thoroughly) in other tests. Here
# we primarily focus on the negative cases.


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
    # We're implicitly checking these don't raise, and explicitly checking they all
    # return False.
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
    with raises_no_typeguard(TypeError):
        somaclass.exists(uri)


def _make_object(name, uri):
    if name == "dataframe":
        pydict = {
            "soma_joinid": pa.array([0, 1, 2, 3, 4], pa.int64()),
            "string": ["apple", "ball", "cat", "dog", "egg"],
        }
        arrow_df = pa.Table.from_pydict(pydict)
        return (soma.DataFrame, soma.DataFrame.create(uri, schema=arrow_df.schema))

    elif name == "sparsendarray":
        return (
            soma.SparseNDArray,
            soma.SparseNDArray.create(uri, type=pa.int64(), shape=(10, 20)),
        )

    elif name == "densendarray":
        return (
            soma.DenseNDArray,
            soma.DenseNDArray.create(uri, type=pa.int64(), shape=(10, 20)),
        )

    elif name == "collection":
        return (soma.Collection, soma.Collection.create(uri))

    elif name == "measurement":
        return (soma.Measurement, soma.Measurement.create(uri))

    elif name == "experiment":
        return (soma.Experiment, soma.Experiment.create(uri))

    else:
        raise Exception("interal unit-test coding error")


@pytest.mark.parametrize(
    "name1",
    [
        "dataframe",
        "sparsendarray",
        "densendarray",
        "collection",
        "measurement",
        "experiment",
    ],
)
@pytest.mark.parametrize(
    "name2",
    [
        "dataframe",
        "sparsendarray",
        "densendarray",
        "collection",
        "measurement",
        "experiment",
    ],
)
def test_tiledbobject_exists_cross_types(tmp_path, name1, name2):
    if name1 == name2:
        uri1 = (tmp_path / name1).as_posix()
        (cls1, obj1) = _make_object(name1, uri1)
        assert cls1.exists(uri1)

    else:
        uri1 = (tmp_path / name1).as_posix()
        uri2 = (tmp_path / name2).as_posix()

        (cls1, obj1) = _make_object(name1, uri1)
        (cls2, obj2) = _make_object(name2, uri2)

        assert not cls1.exists(uri2)
        assert not cls2.exists(uri1)
