# What is SOMA?

SOMA -- for _stack of matrices, annotated_ -- is a unified data model and API for single-cell data.

If you know about `obs`, `var`, and `X`, you'll recognize what you're seeing.

The data model and API -- here as implemented using the TileDB storage engine -- allow you do persist and share single-cell data

* at scale
* with auditable data-sharing
* using the same storage across multiple high-level languages (currently Python and R)
* allowing interop with multiple tools including AnnData, Scanpy, Seurat, and Bioconductor.

See also [the schema specification](https://github.com/single-cell-data/SOMA/blob/main/README.md).

# Examples overview

In these example we will offer how-to's on the [TileDB SingleCell Python package](https://github.com/single-cell-data/TileDB-SingleCell/tree/main/apis/python):

* How to install the software
* How to retrieve sample H5AD inputs from a public S3 bucket
* How to ingest these into SOMA storage on local disk, S3, and TileDB Cloud for increasing levels of scalability and shareability
* How to slice and query SOMA and SOMA-collection objects in new and empowering ways

See also for the R package:

* [R repo](https://github.com/TileDB-Inc/tiledbsc)
* [R docs](https://tiledb-inc.github.io/tiledbsc)

# Notebook

Examples with screenshots and copy/pasteable reusable samples are shown here. As well, you can use
the [public TileDB Cloud notebook](https://cloud.tiledb.com/notebooks/details/johnkerl-tiledb/d3d7ff44-dc65-4cd9-b574-98312c4cbdbd/preview) from which most of these screenshots are taken.
