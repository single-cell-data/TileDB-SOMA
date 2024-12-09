from urllib.parse import urljoin

import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma

soma_xarray = pytest.importorskip("tiledbsoma.io.spatial._xarray_backend")
xr = pytest.importorskip("xarray")


class TestDenseNDDataArray1D:

    @pytest.fixture(scope="class")
    def xr_soma_data_array(self, tmp_path_factory):
        baseuri = tmp_path_factory.mktemp("basic_xarray_io").as_uri()
        array_uri = urljoin(baseuri, "sample_1d_dense_array")
        platform_config = {"tiledb": {"create": {"dims": {"soma_dim_0": {"tile": 4}}}}}
        with soma.DenseNDArray.create(
            array_uri, type=pa.uint8(), shape=(8,), platform_config=platform_config
        ) as array:
            data = pa.Tensor.from_numpy(np.arange(8, dtype=np.uint8))
            array.write((None,), data)
        return soma_xarray.dense_nd_array_to_data_array(
            array_uri,
            dim_names=("xdim",),
            attrs={"meta1": 1, "meta2": "abc"},
        )

    @pytest.fixture(scope="class")
    def xr_numpy_data_array(self):
        data = np.arange(8, dtype=np.uint8)
        return xr.DataArray(data, dims=("xdim",))

    def test_data_array_basics(self, xr_soma_data_array):
        assert isinstance(xr_soma_data_array, xr.DataArray)
        assert xr_soma_data_array.shape == (8,)
        assert xr_soma_data_array.dims == ("xdim",)
        assert xr_soma_data_array.chunks == ((4, 4),)
        assert xr_soma_data_array.attrs == {"meta1": 1, "meta2": "abc"}

    def test_getitem_all(self, xr_soma_data_array, xr_numpy_data_array):
        expected = xr_numpy_data_array[...]
        actual = xr_soma_data_array[...]
        assert actual.shape == (8,)
        np.testing.assert_equal(actual.data.compute(), expected.data)

        actual = xr_soma_data_array[:]
        assert actual.shape == (8,)
        np.testing.assert_equal(actual.data.compute(), expected.data)

    @pytest.mark.parametrize(
        "key",
        [
            (2,),
            (-1,),
            (slice(1, 2),),
            (slice(None, None),),
            (slice(1, -1),),
            ([0, 3],),
            ([0, -3],),
        ],
    )
    def test_getitem(self, xr_soma_data_array, xr_numpy_data_array, key):
        actual = xr_soma_data_array[key]
        expected = xr_numpy_data_array[key]

        assert actual.dims == expected.dims
        assert actual.shape == expected.shape

        np.testing.assert_equal(actual.data.compute(), expected.data)

    # TEMPORARY PENDING https://github.com/single-cell-data/TileDB-SOMA/issues/3398
    # @pytest.mark.parametrize(
    #    "key",
    #    [
    #        (slice(0, 4, 2),),
    #    ],
    # )
    # def test_getitem_with_steps(self, xr_soma_data_array, key):
    #    with pytest.raises(NotImplementedError):
    #        xr_soma_data_array[key].data.compute()


class TestDenseNDDataArray3D:

    @pytest.fixture(scope="class")
    def xr_soma_data_array(self, tmp_path_factory):
        baseuri = tmp_path_factory.mktemp("basic_xarray_io").as_uri()
        array_uri = urljoin(baseuri, "sample_3d_dense_array")
        with soma.DenseNDArray.create(
            array_uri, type=pa.uint8(), shape=(8, 2, 4)
        ) as array:
            data = pa.Tensor.from_numpy(
                np.reshape(np.arange(64, dtype=np.uint8), (8, 2, 4))
            )
            array.write((None,), data)
        return soma_xarray.dense_nd_array_to_data_array(
            array_uri, dim_names=("xdim", "ydim", "zdim")
        )

    @pytest.fixture(scope="class")
    def xr_numpy_data_array(self):
        data = np.reshape(np.arange(64, dtype=np.uint8), (8, 2, 4))
        return xr.DataArray(data, dims=("xdim", "ydim", "zdim"))

    def test_data_array_basics(self, xr_soma_data_array):
        assert isinstance(xr_soma_data_array, xr.DataArray)
        assert xr_soma_data_array.shape == (8, 2, 4)
        assert xr_soma_data_array.chunks == ((8,), (2,), (4,))
        assert xr_soma_data_array.dims == ("xdim", "ydim", "zdim")
        assert len(xr_soma_data_array.attrs) == 0

    def test_getitem_all(self, xr_soma_data_array, xr_numpy_data_array):
        expected = xr_numpy_data_array[...]

        actual = xr_soma_data_array[...]
        assert actual.shape == (8, 2, 4)
        np.testing.assert_equal(actual.data, expected.data)

        actual = xr_soma_data_array[:]
        assert actual.shape == (8, 2, 4)
        np.testing.assert_equal(actual.data, expected.data)

        actual = xr_soma_data_array[:, :, :]
        assert actual.shape == (8, 2, 4)
        np.testing.assert_equal(actual.data, expected.data)

    @pytest.mark.parametrize(
        "key",
        [
            (2, 1, 3),
            (0, 0, slice(1, 2)),
            (2, 1, slice(None)),
            (0, slice(None, None), 1),
            ([1, 3], 1, 1),
            ([1, 3], 1, [0, 2]),
            (-1, -1, -1),
            ([0, 3, 4], [0, 1, 0], [1, 3, 2]),
            (..., slice(0, 3)),
            (slice(1, 2), ...),
            (slice(0, 5, 2), [0, 1, 0], [1, 3, 2]),
        ],
    )
    def test_getitem(self, xr_soma_data_array, xr_numpy_data_array, key):
        actual = xr_soma_data_array[key]
        expected = xr_numpy_data_array[key]

        assert actual.dims == expected.dims
        assert actual.shape == expected.shape

        np.testing.assert_equal(actual.data.compute(), expected.data)

    # TEMPORARY PENDING https://github.com/single-cell-data/TileDB-SOMA/issues/3398
    # @pytest.mark.parametrize(
    #    "key",
    #    [
    #        (0, 1, slice(0, 4, 2)),
    #    ],
    # )
    # def test_getitem_with_steps(self, xr_soma_data_array, key):
    #    with pytest.raises(NotImplementedError):
    #        xr_soma_data_array[key].data.compute()
