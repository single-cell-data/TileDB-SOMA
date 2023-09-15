soco.tgz was created by running the following from the repo base directory (.. from here):

```
apis/python/devtools/ingestor --debug --soco -o test/soco -n \
  data/pbmc3k_processed.h5ad data/10x-pbmc-multiome-v1.0/subset_100_100.h5ad
```

while versions were at

```
>>> import tiledbsoma
>>> tiledbsoma.show_package_versions()
tiledbsoma.__version__        1.2.7
TileDB-Py tiledb.version()    (0, 21, 6)
TileDB core version           2.15.4
libtiledbsoma version()       libtiledb=2.15.2
python version                3.10.6.final.0
OS version                    Linux 5.19.0-1025-aws
```
