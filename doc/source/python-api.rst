TileDB-SOMA Python API Reference
================================

**SOMA is a unified data model and API for single-cell data.**

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

Typical usage of the Python interface to TileDB-SOMA will use the top-level module ``tiledbsoma``, e.g.

.. code-block:: python

   import tiledbsoma

There is also a submodule ``io`` which contains logic for importing data from ``AnnData`` to SOMA structure, and exporting back to ``AnnData``.

.. code-block:: python

   import tiledbsoma.io


Composed Types
--------------

Collection
^^^^^^^^^^^^^^
.. autoclass:: tiledbsoma.Collection
   :members:

Experiment
^^^^^^^^^^^^^^
.. autoclass:: tiledbsoma.Experiment
   :members:

Measurement
^^^^^^^^^^^^^^^
.. autoclass:: tiledbsoma.Measurement
   :members:

Foundational Types
------------------

DataFrame
^^^^^^^^^^^^^
.. autoclass:: tiledbsoma.DataFrame
   :members:

SparseNDArray
^^^^^^^^^^^^^^^^^
.. autoclass:: tiledbsoma.SparseNDArray
   :members:

DenseNDArray
^^^^^^^^^^^^^^^^
.. autoclass:: tiledbsoma.DenseNDArray
   :members:

I/O functions
-------------

from_h5ad
^^^^^^^^^
.. autofunction:: tiledbsoma.io.from_h5ad

from_anndata
^^^^^^^^^^^^
.. autofunction:: tiledbsoma.io.from_anndata

to_h5ad
^^^^^^^
.. autofunction:: tiledbsoma.io.to_h5ad

to_anndata
^^^^^^^^^^
.. autofunction:: tiledbsoma.io.to_anndata

Options
-------

PlatformConfig
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: tiledbsoma.PlatformConfig
   :members:

logging
^^^^^^^
.. automodule:: tiledbsoma.logging
   :members:
