# SOMA C++ API Notes
This document describes the SOMA C++ API implementation of the [SOMA API specification](https://github.com/single-cell-data/SOMA/blob/spec-revision/brainstorming.md), including:
* Notes for for SOMA High-Level Language API implementations (Python and R), aka. **SOMA HLL APIs**.
* SOMA API specification questions, ambiguities, and issues.

## SOMA API Specification Types

> The foundational types are:
> 
> * SOMACollection - a string-keyed container (key-value map) of other SOMA data types, e.g., SOMADataFrame, SOMADataMatrix and SOMACollection.
> * SOMADataFrame and SOMAIndexedDataFrame - a multi-column table -- essentially a dataframe, available in types which support offset indexing/slicing, or indexing based upon user-defined columns.
> * SOMADenseNdArray and SOMASparseNdArray- an offset addressed (zero-based), single-type N-D array, available in either sparse or dense instantiations
> 
> The composed types are:
> 
> * SOMAExperiment - a specialization and extension of SOMACollection, codifying a set of naming and indexing conventions to represent annotated, 2-D matrix of observations across multiple sets of variables.
> 

The `SOMACollection` and `SOMAExperiment` SOMA API types are implemented with TileDB Groups. The C++ API will implement a `SOMAGroup` class that provides functional operations for these **SOMA API Group** types.

The `SOMADataFrame`, `SOMAIndexedDataFrame`, `SOMADenseNdArray`, and `SOMASparseNdArray` SOMA API types are implemented with TileDB Arrays. The C++ API will implement `SOMAReader` and `SOMAWriter` classes that provide functional operations for these **SOMA API Array** types.

## SOMAReader (SOMA C++ API)
The `SOMAReader` class implements read operations for the **SOMA API Array** types:

### SOMADataFrame (SOMA API Spec)
```
read(
    slices=[row_slices, ...],
    column_names=[`string`, ...]|all,
    batch_size,
    partitions,
    result_order,
    value_filter,
    platform_config,
) -> delayed iterator over Arrow.RecordBatch
```

### SOMAIndexedDataFrame (SOMA API Spec)
```
read(
    ids=[[id,...]|all, ...],
    column_names=[`string`, ...]|all,
    batch_size,
    partitions,
    result_order,
    value_filter,
    platform_config,
) -> delayed iterator over Arrow.RecordBatch
```

### SOMADenseNdArray (SOMA API Spec)
```
read(
    [slice, ...],
    batch_size,
    partitions,
    result_order,
    batch_format,
    platform_config,
) -> delayed iterator over ReadResult
```

### SOMASparseNdArray (SOMA API Spec)
```
read(
    [slice, ...],
    batch_size,
    partitions,
    result_order,
    batch_format,
    platform_config,
) -> delayed iterator over ReadResult
```

### SOMA API Spec `read` parameters

A combined set of `read` parameters for all **SOMA API Spec** array types are implemented by the `SOMAReader` class. The **SOMA HLL APIs** will expose the appropriate parameters for each **SOMA API Array** type.

* `ids` - the rows to read. Defaults to 'all'. Coordinates for each dimension may be specified by value, a value range (slice) or a list of both.
* `slices` - the rows to read. Defaults to 'all'. Rows are addressable with a row offset (uint), a row-offset range (slice) or a list of both.
* `column_names` - the named columns to read and return. Defaults to all. The pseudo-column soma_rowid may be included in this list.
* `batch_size` - a SOMABatchSize, indicating the size of each "batch" returned by the read iterator. Defaults to auto.
* `partitions` - an optional SOMAReadPartitions to partition read operations.
* `result_order` - order of read results. If dataframe is indexed, can be one of row-major, col-major or unordered. If dataframe is non-indexed, can be one of rowid-ordered or unordered.
* `value_filter` - an optional value filter to apply to the results. Defaults to no filter. The soma_rowid pseudo-column can not be filtered.
* `platform_config` - optional storage-engine specific configuration

> **Question**: Can we allow `slices` to be passed as an `Arrow.Array` or `Arrow.ChunkedArray`? This would avoid copying `obs` or `var` query results (already in `Arrow` format) to a list.

> **Question**: How will the `read` memory budget be set? Maybe with `platform_config`? Setting the `batch_size` will require dividing the memory budget between TileDB (`sm.mem.total_budget`) and TileDB-SOMA buffers for the selected `column_names`.

> **Note**: The **SOMA Python API** `read` iterator currently returns an `PyArrow.Table`, which is a collection of `Arrow.RecordBatch` objects.

> **Note**: The `PyArrow.Table` returned by the **SOMA C++ API** `read` iterator must be fully processed/analyzed or copied (for example, `to_pandas`) before reading the next batch from the iterator. This is required because the buffers backing the Arrow results are reused by the next read by the iterator. This approach helps control memory usage, but prevents the concatenation of `PyArrow.Table`s from different read batches.

> **Note**: Assuming `partitions` will partition reads by "rows to read", read partitioning divides the reads based on `ids` or `slices` or `all`. For each read, the **SOMA HLL APIs** must set the `ids` or `slices` **SOMA C++ API** parameter *once* (including for the default `all`) and provide the partitioning information (`IofN`).

> **Note**: Assuming `partitions` will partition reads by "rows to read", read partitioning divides the reads based on `ids` or `slices` or `all`. For each read, the **SOMA HLL APIs** must set the `ids` or `slices` **SOMA C++ API** parameter *once* (including for the default `all`) and provide the partitioning information (`IofN`).

> **Note**: Initially, partitioning is only supported along one dimension.


```python
reader = soma.SOMAReader(
    uri, 
    name='unnamed',
    column_names=[], # empty list is all columns
    query_condition=None,
    schema=None,
    batch_size='auto',
    result_order='auto',
    platform_config={},
    )
```
