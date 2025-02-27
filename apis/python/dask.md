# `to_anndata` Dask notes

## Example runs

### Scanpy HVG

#### 1.5M heart cells, 100k cells per Dask chunk

##### 4x4x4 (Dask procs x Dask threads x Tile threads)

Top allocations ([hvg/heart/100k_4x4x4.json]):
<!-- `top-allocs.py -m out/hvg/heart/100k_4x4x4.json` -->
<div>

- [tiledbsoma/\_read\_iters.py](http://github.com/single-cell-data/TileDB-SOMA/blob/a0348ff5/apis/python/src/tiledbsoma/_read_iters.py#L586): 37 GiB
- [cloudpickle/cloudpickle.py](http://github.com/cloudpipe/cloudpickle/blob/v3.1.1/cloudpickle/cloudpickle.py#L1303): 5.1 GiB
- `<stack trace unavailable>`: 2.75 GiB
- [dask/tokenize.py](http://github.com/dask/dask/blob/2024.11.2/dask/tokenize.py#L242): 2.35 GiB
- [msgpack/\_\_init\_\_.py](http://github.com/msgpack/msgpack-python/blob/v1.1.0/msgpack/__init__.py#L36): 1.32 GiB
</div>


<details><summary>Video</summary>

https://github.com/user-attachments/assets/355c4151-56a6-4642-b744-2861a5acd42d
</details>

##### 2x2x2

Top allocations ([hvg/heart/100k_2x2x2.json]):
<!-- `top-allocs.py -m out/hvg/heart/100k_2x2x2.json` -->
<div>

- [tiledbsoma/\_read\_iters.py](http://github.com/single-cell-data/TileDB-SOMA/blob/a0348ff5/apis/python/src/tiledbsoma/_read_iters.py#L586): 37 GiB
- [cloudpickle/cloudpickle.py](http://github.com/cloudpipe/cloudpickle/blob/v3.1.1/cloudpickle/cloudpickle.py#L1303): 5.18 GiB
- [dask/tokenize.py](http://github.com/dask/dask/blob/2024.11.2/dask/tokenize.py#L242): 2.35 GiB
- `<stack trace unavailable>`: 2.11 GiB
- [msgpack/\_\_init\_\_.py](http://github.com/msgpack/msgpack-python/blob/v1.1.0/msgpack/__init__.py#L36): 1.38 GiB
</div>


<details><summary>Video</summary>

https://github.com/user-attachments/assets/ba090066-9811-41b8-b538-c5efaff56ef1
</details>

#### 3.9M lung cells, 100k cells per Dask chunk

##### 4x4x4

Top allocations ([hvg/lung/100k_4x4x4.json]):

<!-- `top-allocs.py -m out/hvg/lung/100k_4x4x4.json` -->
<div>

- [tiledbsoma/\_read\_iters.py](http://github.com/single-cell-data/TileDB-SOMA/blob/a0348ff5/apis/python/src/tiledbsoma/_read_iters.py#L586): 37 GiB
- [cloudpickle/cloudpickle.py](http://github.com/cloudpipe/cloudpickle/blob/v3.1.1/cloudpickle/cloudpickle.py#L1303): 25.2 GiB
- [msgpack/\_\_init\_\_.py](http://github.com/msgpack/msgpack-python/blob/v1.1.0/msgpack/__init__.py#L36): 6.32 GiB
- [dask/tokenize.py](http://github.com/dask/dask/blob/2024.11.2/dask/tokenize.py#L242): 5.44 GiB
- `<stack trace unavailable>`: 4.72 GiB
</div>

[Video](https://drive.google.com/file/d/1zjU29bxrvF52RzmyysKaYAxEEAhCk7Br/view)

##### 2x2x2

Top allocations ([hvg/lung/100k_2x2x2.json]):

<!-- `top-allocs.py -m out/hvg/lung/100k_2x2x2.json` -->
<div>

- [tiledbsoma/\_read\_iters.py](http://github.com/single-cell-data/TileDB-SOMA/blob/a0348ff5/apis/python/src/tiledbsoma/_read_iters.py#L586): 37 GiB
- [cloudpickle/cloudpickle.py](http://github.com/cloudpipe/cloudpickle/blob/v3.1.1/cloudpickle/cloudpickle.py#L1303): 26 GiB
- [msgpack/\_\_init\_\_.py](http://github.com/msgpack/msgpack-python/blob/v1.1.0/msgpack/__init__.py#L36): 6.67 GiB
- [dask/tokenize.py](http://github.com/dask/dask/blob/2024.11.2/dask/tokenize.py#L242): 5.44 GiB
- `<stack trace unavailable>`: 3.73 GiB
</div>


[Video](https://drive.google.com/file/d/1uVtV_iubQAAN2AnaPO50gr2sA7c_7TDE/view)

[hvg/heart/100k_4x4x4.json]: https://github.com/single-cell-data/TileDB-SOMA/blob/rw/dsk2/apis/python/out/hvg/heart/100k_4x4x4.json
[hvg/heart/100k_2x2x2.json]: https://github.com/single-cell-data/TileDB-SOMA/blob/rw/dsk2/apis/python/out/hvg/heart/100k_2x2x2.json
[hvg/lung/100k_4x4x4.json]: https://github.com/single-cell-data/TileDB-SOMA/blob/rw/dsk2/apis/python/out/hvg/lung/100k_4x4x4.json
[hvg/lung/100k_2x2x2.json]: https://github.com/single-cell-data/TileDB-SOMA/blob/rw/dsk2/apis/python/out/hvg/lung/100k_2x2x2.json
