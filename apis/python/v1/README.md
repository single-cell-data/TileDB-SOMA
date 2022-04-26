# Status

* This is a prototype version of the ingestor.
* It's known to work with a particular raw 10X data file: `../../data/10x-pbmc-multiome-v1.0/subset_100_100.h5ad`.
* It is not known support cellxgene data files.

# Running the v1 code

Note that `anndata_to_tiledb.py` has an `if __name__ == "__main__":` allowing us to run it from
the command line, with input HDF5 file as `argv[1]` and output TileDB directory name as `argv[2]`.

```
$ pwd
/Users/johnkerl/git/single-cell-data/TileDB-SingleCell/util

$ python anndata_to_tiledb.py ../../data/10x-pbmc-multiome-v1.0/subset_100_100.h5ad subset_100_100

$ ls -l subset_100_100/
total 0
drwxr-xr-x  7 johnkerl  staff  224 Apr 25 10:33 X
drwxr-xr-x  2 johnkerl  staff   64 Apr 25 10:33 __group
drwxr-xr-x  2 johnkerl  staff   64 Apr 25 10:33 __meta
-rw-r--r--  1 johnkerl  staff    0 Apr 25 10:33 __tiledb_group.tdb
drwxr-xr-x  7 johnkerl  staff  224 Apr 25 10:33 obs
drwxr-xr-x  7 johnkerl  staff  224 Apr 25 10:33 var

>>> import tiledb

>>> print(tiledb.open("subset_100_100/X").schema)
ArraySchema(
  domain=Domain(*[
    Dim(name='obs', domain=(None, None), tile=None, dtype='|S0', var=True),
    Dim(name='var', domain=(None, None), tile=None, dtype='|S0', var=True),
  ]),
  attrs=[
    Attr(name='data', dtype='float32', var=False, nullable=False),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=True,
  allows_duplicates=False,
)
```

Discussion:

* `../../data/10x-pbmc-multiome-v1.0/subset_100_100.h5ad subset_100_100` is raw off-the-hardware data
* In the HDF5 file we read, there is `X`, `obs`, `var`
* In the TileDB group we write, there is `X` as an array, `obs`, and `var`
