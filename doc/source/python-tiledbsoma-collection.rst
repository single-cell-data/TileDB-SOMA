tiledbsoma.Collection
=====================

.. currentmodule:: tiledbsoma

.. autoclass:: Collection

   .. automethod:: __init__

   .. rubric:: Methods

   .. autosummary::
      :toctree: _generated

      ~Collection.__init__
      ~Collection.exists
      ~Collection.create
      ~Collection.open
      ~Collection.reopen
      ~Collection.close
      ~Collection.verify_open_for_writing

      ~Collection.items
      ~Collection.keys
      ~Collection.values
      ~Collection.members
      ~Collection.get
      ~Collection.set
      ~Collection.setdefault
      ~Collection.update
      ~Collection.clear
      ~Collection.pop
      ~Collection.popitem

      ~Collection.add_new_collection
      ~Collection.add_new_dataframe
      ~Collection.add_new_dense_ndarray
      ~Collection.add_new_sparse_ndarray

   .. rubric:: Attributes

   .. autosummary::
      :toctree: _generated

      ~Collection.uri
      ~Collection.closed
      ~Collection.context
      ~Collection.metadata

      ~Collection.mode
      ~Collection.soma_type
      ~Collection.tiledb_timestamp
      ~Collection.tiledb_timestamp_ms

