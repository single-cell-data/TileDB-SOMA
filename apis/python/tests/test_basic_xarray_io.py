from urllib.parse import urljoin

import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma

soma_xarray = pytest.importorskip("tiledbsoma.experimental._xarray_backend")
xr = pytest.importorskip("xarray")


@pytest.fixture(scope="module")
def sample_1d_dense_array(tmp_path_factory):
    """Create a sample DenseNDArray with data and return open for reading."""
    baseuri = tmp_path_factory.mktemp("basic_xarray_io").as_uri()
    array_uri = urljoin(baseuri, "sample_1d_dense_array")
    array = soma.DenseNDArray.create(array_uri, type=pa.uint8(), shape=(8,))
    data = pa.Tensor.from_numpy(np.arange(8, dtype=np.uint8))
    array.write((None,), data)
    array.close()
    array = soma.DenseNDArray.open(array_uri)
    return array


@pytest.fixture(scope="module")
def sample_3d_dense_array(tmp_path_factory):
    """Create a sample DenseNDArray with data and return open for reading."""
    baseuri = tmp_path_factory.mktemp("basic_xarray_io").as_uri()
    array_uri = urljoin(baseuri, "sample_3d_dense_array")
    array = soma.DenseNDArray.create(array_uri, type=pa.uint8(), shape=(8, 2, 4))
    data = pa.Tensor.from_numpy(np.reshape(np.arange(64, dtype=np.uint8), (8, 2, 4)))
    array.write((None,), data)
    array.close()
    array = soma.DenseNDArray.open(array_uri)
    return array


class Test1DSampleDenseNDArrayWrapper:

    @pytest.fixture(scope="class")
    def soma_variable(self, sample_1d_dense_array):
        array_wrapper = soma_xarray.DenseNDArrayWrapper(sample_1d_dense_array)
        assert isinstance(array_wrapper, xr.backends.BackendArray)
        return xr.Variable(
            ("xdim",), xr.core.indexing.LazilyIndexedArray(array_wrapper)
        )

    @pytest.fixture(scope="class")
    def numpy_variable(self):
        data = np.arange(8, dtype=np.uint8)
        return xr.Variable(("xdim",), data)

    def test_variable_basics(self, soma_variable):
        assert soma_variable.shape == (8,)
        assert soma_variable.dims == ("xdim",)
        assert isinstance(soma_variable._data, xr.core.indexing.LazilyIndexedArray)

    def test_getitem_all(self, soma_variable, numpy_variable):
        expected = numpy_variable[...]

        actual = soma_variable[...]
        assert actual.shape == (8,)
        np.testing.assert_equal(actual.data, expected.data)

        actual = soma_variable[:]
        assert actual.shape == (8,)
        np.testing.assert_equal(actual.data, expected.data)

    @pytest.mark.parametrize(
        "key",
        [
            (2,),
            (-1,),
            (slice(1, 2),),
            (slice(None, None),),
            (slice(1, -1),),
            (slice(0, 4, 2),),
            ([0, 3],),
            ([0, -3],),
        ],
    )
    def test_getitem(self, soma_variable, numpy_variable, key):
        actual = soma_variable[*key]
        expected = numpy_variable[*key]

        assert actual.dims == expected.dims
        assert actual.shape == expected.shape

        np.testing.assert_equal(actual.data, expected.data)


class TestSample3DDenseNDArrayWrapper:

    @pytest.fixture(scope="class")
    def soma_variable(self, sample_3d_dense_array):
        array_wrapper = soma_xarray.DenseNDArrayWrapper(sample_3d_dense_array)
        assert isinstance(array_wrapper, xr.backends.BackendArray)
        return xr.Variable(
            ("xdim", "ydim", "zdim"), xr.core.indexing.LazilyIndexedArray(array_wrapper)
        )

    @pytest.fixture(scope="class")
    def numpy_variable(self):
        data = np.reshape(np.arange(64, dtype=np.uint8), (8, 2, 4))
        return xr.Variable(("xdim", "ydim", "zdim"), data)

    def test_variable_basics(self, soma_variable):
        assert soma_variable.shape == (8, 2, 4)
        assert soma_variable.dims == ("xdim", "ydim", "zdim")
        assert isinstance(soma_variable._data, xr.core.indexing.LazilyIndexedArray)

    def test_getitem_all(self, soma_variable, numpy_variable):
        expected = numpy_variable[...]

        actual = soma_variable[...]
        assert actual.shape == (8, 2, 4)
        np.testing.assert_equal(actual.data, expected.data)

        actual = soma_variable[:]
        assert actual.shape == (8, 2, 4)
        np.testing.assert_equal(actual.data, expected.data)

        actual = soma_variable[:, :, :]
        assert actual.shape == (8, 2, 4)
        np.testing.assert_equal(actual.data, expected.data)

    @pytest.mark.parametrize(
        "key",
        [
            (2, 1, 3),
            (0, 0, slice(1, 2)),
            (2, 1, slice(None, None)),
            (0, slice(None, None), 1),
            ([1, 3], 1, 1),
            ([1, 3], 1, [0, 2]),
            (-1, -1, -1),
            (0, 1, slice(0, 4, 2)),
            ([0, 3, 4], [0, 1, 0], [1, 3, 2]),
            (slice(0, 5, 2), [0, 1, 0], [1, 3, 2]),
        ],
    )
    def test_getitem(self, soma_variable, numpy_variable, key):
        actual = soma_variable[*key]
        expected = numpy_variable[*key]

        assert actual.dims == expected.dims
        assert actual.shape == expected.shape

        np.testing.assert_equal(actual.data, expected.data)


def test_1d_xarray_dataset(sample_1d_dense_array):
    datastore = soma_xarray.DenseNDArrayDatastore(
        sample_1d_dense_array.uri, attrs={"test1": 1, "test2": 2}
    )
    ds = xr.Dataset.load_store(datastore)
    assert isinstance(ds, xr.Dataset)
    assert len(ds.attrs) == 2
    assert ds.attrs == {"test1": 1, "test2": 2}

    assert len(ds.data_vars) == 1
    assert "soma_data" in ds

    data = ds["soma_data"]
    assert isinstance(data, xr.DataArray)
    assert data.dims == ("soma_dim_0",)
    assert len(data.attrs) == 0

    expected = np.arange(8, dtype=np.uint8)
    np.testing.assert_equal(data.data, expected)


def test_3d_xarray_dataset(sample_3d_dense_array):
    datastore = soma_xarray.DenseNDArrayDatastore(sample_3d_dense_array.uri)
    ds = xr.Dataset.load_store(datastore)
    assert isinstance(ds, xr.Dataset)
    assert len(ds.attrs) == 0

    assert len(ds.data_vars) == 1
    assert "soma_data" in ds

    data = ds["soma_data"]
    assert isinstance(data, xr.DataArray)
    assert data.dims == ("soma_dim_0", "soma_dim_1", "soma_dim_2")
    assert len(data.attrs) == 0

    expected = np.reshape(np.arange(64, dtype=np.uint8), (8, 2, 4))
    np.testing.assert_equal(data.data, expected)
