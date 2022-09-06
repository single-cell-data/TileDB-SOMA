TileDB-SC Python API Reference
==============================

**SOMA --- for stack of matrices, annotated --- is a unified data model and API for single-cell data.**

If you know about ``obs``, ``var``, and ``X``, you'll recognize what you're seeing.

The `data model and API <https://github.com/single-cell-data>`_ --- here as implemented using the `TileDB storage engine <https://tiledb.com>`_ --- allow you to persist, investigate, and share annotated 2D matrices, commonly used in single-cell biology.

Features:

* Flexible, extensible, and open-source API
* Supports access to persistent, cloud-resident annotated 2D matrix datasets
* Enables use within popular data science environments (e.g., R, Python), using the tools of that environment (e.g., Python Pandas integration), with the same storage regardless of language
* Allows interoperability with multiple tools including AnnData, Scanpy, Seurat, and Bioconductor
* Cloud-native TileDB arrays allow you to slice straight from remote storage
* Reduces costs and processing time by utilizing cost-efficient object storage services like S3
* Enables out-of-core access to data aggregations much larger than single-host main memory
* Enables distributed computation over datasets

Modules
-------

Typical usage of the Python interface to TileDB-SC will use the top-level module ``tiledbsc``, e.g.

.. code-block:: python

   import tiledbsc

There is also a submodule ``io`` which contains logic for importing data from ``AnnData`` to SOMA structure, and exporting back to ``AnnData``.

.. code-block:: python

   import tiledbsc.io

SOMA
----

.. autoclass:: tiledbsc.SOMA
   :members:

SOMACollection
--------------

.. autoclass:: tiledbsc.SOMACollection
   :members:

SOMASlice
---------

.. autoclass:: tiledbsc.SOMASlice
   :members:

I/O functions
-------------

.. autofunction:: tiledbsc.io.from_h5ad
.. autofunction:: tiledbsc.io.from_anndata
.. autofunction:: tiledbsc.io.to_h5ad
.. autofunction:: tiledbsc.io.to_anndata

Options
-------

.. autoclass:: tiledbsc.SOMAOptions
   :members:
.. autoclass:: tiledbsc.logging
   :members:

SOMA-element classes
--------------------

.. autoclass:: tiledbsc.AssayMatrixGroup
   :members:
.. autoclass:: tiledbsc.AssayMatrix
   :members:
.. autoclass:: tiledbsc.AnnotationDataFrame
   :members:
.. autoclass:: tiledbsc.AnnotationMatrixGroup
   :members:
.. autoclass:: tiledbsc.AnnotationMatrix
   :members:
.. autoclass:: tiledbsc.AnnotationPairwiseMatrixGroup
   :members:
.. autoclass:: tiledbsc.RawGroup
   :members:
.. autoclass:: tiledbsc.UnsGroup
   :members:
.. autoclass:: tiledbsc.UnsArray
   :members:

Implementation-level classes
----------------------------

.. autoclass:: tiledbsc.TileDBArray
   :members:
.. autoclass:: tiledbsc.TileDBGroup
   :members:
.. autoclass:: tiledbsc.TileDBObject
   :members:
.. autoclass:: tiledbsc.util
   :members:
.. autoclass:: tiledbsc.util_ann
   :members:
.. autoclass:: tiledbsc.util_tiledb
   :members:
