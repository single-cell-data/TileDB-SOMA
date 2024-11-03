tiledbsoma.Scene
================

.. currentmodule:: tiledbsoma

.. autoclass:: Scene

   .. automethod:: __init__

   .. rubric:: Methods

   .. autosummary::
      :toctree: _generated

      ~Scene.__init__
      ~Scene.exists
      ~Scene.create
      ~Scene.open
      ~Scene.reopen
      ~Scene.close
      ~Scene.verify_open_for_writing

      ~Scene.items
      ~Scene.keys
      ~Scene.values
      ~Scene.members
      ~Scene.get
      ~Scene.set
      ~Scene.setdefault
      ~Scene.update
      ~Scene.clear
      ~Scene.pop
      ~Scene.popitem

      ~Scene.add_new_collection
      ~Scene.add_new_dataframe
      ~Scene.add_new_dense_ndarray
      ~Scene.add_new_geometry_dataframe
      ~Scene.add_new_multiscale_image
      ~Scene.add_new_point_cloud_dataframe
      ~Scene.add_new_sparse_ndarray

      ~Scene.get_transform_from_geometry_dataframe
      ~Scene.get_transform_from_multiscale_image
      ~Scene.get_transform_from_point_cloud_dataframe
      ~Scene.get_transform_to_geometry_dataframe
      ~Scene.get_transform_to_multiscale_image
      ~Scene.get_transform_to_point_cloud_dataframe
      ~Scene.set_transform_to_geometry_dataframe
      ~Scene.set_transform_to_multiscale_image
      ~Scene.set_transform_to_point_cloud_dataframe

   .. rubric:: Attributes

   .. autosummary::
      :toctree: _generated

      ~Scene.uri
      ~Scene.context
      ~Scene.metadata

      ~Scene.mode
      ~Scene.closed
      ~Scene.soma_type
      ~Scene.tiledb_timestamp
      ~Scene.tiledb_timestamp_ms

      ~Scene.coordinate_space
      ~Scene.img
      ~Scene.obsl
      ~Scene.varl

