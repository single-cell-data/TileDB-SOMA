tiledbsoma.DataFrame
====================

.. currentmodule:: tiledbsoma

.. autoclass:: DataFrame

   .. automethod:: __init__

   .. rubric:: Methods

   .. autosummary::
      :toctree: _generated

      ~DataFrame.__init__
      ~DataFrame.exists
      ~DataFrame.create
      ~DataFrame.open
      ~DataFrame.reopen
      ~DataFrame.close
      ~DataFrame.read
      ~DataFrame.write
      ~DataFrame.verify_open_for_writing

      ~DataFrame.keys

      ~DataFrame.tiledbsoma_upgrade_domain
      ~DataFrame.change_domain
      ~DataFrame.tiledbsoma_resize_soma_joinid_shape
      ~DataFrame.tiledbsoma_upgrade_soma_joinid_shape
      ~DataFrame.non_empty_domain

      ~DataFrame.config_options_from_schema

   .. rubric:: Attributes

   .. autosummary::
      :toctree: _generated

      ~DataFrame.uri
      ~DataFrame.soma_type
      ~DataFrame.schema
      ~DataFrame.index_column_names

      ~DataFrame.count
      ~DataFrame.shape
      ~DataFrame.maxshape
      ~DataFrame.domain
      ~DataFrame.maxdomain
      ~DataFrame.tiledbsoma_has_upgraded_domain

      ~DataFrame.mode
      ~DataFrame.closed

      ~DataFrame.context
      ~DataFrame.tiledb_timestamp
      ~DataFrame.tiledb_timestamp_ms

      ~DataFrame.metadata

