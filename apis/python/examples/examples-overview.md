# What is SOMA?

SOMA -- for _stack of matrices, annotated_ -- is a unified data model and API for single-cell data.

If you know about `obs`, `var`, and `X`, you'll recognize what you're seeing.

The [data model and API](https://github.com/single-cell-data) -- here as implemented using the [TileDB storage engine](https://tiledb.com) -- allow you to persist, investigate, and share annotated 2D matrices, commonly used in single-cell biology.

Features:

* flexible, extensible, and open-source API
* supports access to persistent, cloud-resident annotated 2D matrix datasets
* enables use within popular data science environments (e.g., R, Python), using the tools of that environment (e.g., Python Pandas integration), with the same storage regardless of language
* allows interop with multiple tools including AnnData, Scanpy, Seurat, and Bioconductor
* cloud-native TileDB arrays allow you to slice straight from remote storage
* reduces costs and processing time by utilizing cost-efficient object storage services like S3
* enables out-of-core access to data aggregations much larger than single-host main memory
* enables distributed computation over datasets

# Examples overview

In these example we will offer how-to's on the [TileDB SingleCell Python package](https://github.com/single-cell-data/TileDB-SingleCell/tree/main/apis/python):

* How to install the software
* How to retrieve sample H5AD inputs from a public S3 bucket
* How to ingest these into SOMA storage on local disk, S3, and TileDB Cloud for increasing levels of scalability and shareability
* How to slice and query SOMA and SOMA-collection objects in new and empowering ways

# Notebook

Examples with screenshots and copy/pasteable reusable samples are shown here. As well, you can use
the [public TileDB Cloud notebook](https://cloud.tiledb.com/notebooks/details/johnkerl-tiledb/33c4fe81-d15f-43cd-a588-5c277cf70cb6/preview) from which most of these screenshots are taken.
