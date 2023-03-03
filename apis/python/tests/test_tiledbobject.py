import tiledbsoma as soma

# Checking that objects _do_ exist is already done (thoroughly) in other
# tests. Here we check the negative cases.
#
# Specifically, we're checking these don't raise, and all return False.


def test_dataframe_exists():
    assert not soma.DataFrame.exists("/nonesuch/no/nope/never/ever")
    assert not soma.DataFrame.exists("foo://bar")
    assert not soma.DataFrame.exists("s3://@@@@ILLEGAL@@@@")


def test_sparse_nd_array_exists():
    assert not soma.SparseNDArray.exists("/nonesuch/no/nope/never/ever")
    assert not soma.SparseNDArray.exists("foo://bar")
    assert not soma.SparseNDArray.exists("s3://@@@@ILLEGAL@@@@")


def test_dense_nd_array_exists():
    assert not soma.DenseNDArray.exists("/nonesuch/no/nope/never/ever")
    assert not soma.DenseNDArray.exists("foo://bar")
    assert not soma.DenseNDArray.exists("s3://@@@@ILLEGAL@@@@")


def test_collection_exists():
    assert not soma.Collection.exists("/nonesuch/no/nope/never/ever")
    assert not soma.Collection.exists("foo://bar")
    assert not soma.Collection.exists("s3://@@@@ILLEGAL@@@@")


def test_measurement_exists():
    assert not soma.Measurement.exists("/nonesuch/no/nope/never/ever")
    assert not soma.Measurement.exists("foo://bar")
    assert not soma.Measurement.exists("s3://@@@@ILLEGAL@@@@")


def test_experiment_exists():
    assert not soma.Experiment.exists("/nonesuch/no/nope/never/ever")
    assert not soma.Experiment.exists("foo://bar")
    assert not soma.Experiment.exists("s3://@@@@ILLEGAL@@@@")
