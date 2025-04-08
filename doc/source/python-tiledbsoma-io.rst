The tiledbsoma.io module
========================

.. currentmodule:: tiledbsoma.io

.. automodule:: tiledbsoma.io


Functions
---------

.. rubric:: Data conversion to/from TileDB-SOMA

.. autosummary::
    :toctree: _autosummary/
    :nosignatures:

    tiledbsoma.io.from_h5ad
    tiledbsoma.io.from_anndata
    tiledbsoma.io.to_h5ad
    tiledbsoma.io.to_anndata

.. rubric:: Updating values within a TileDB-SOMA Experiment

.. autosummary::
    :toctree: _autosummary/
    :nosignatures:

    tiledbsoma.io.add_X_layer
    tiledbsoma.io.add_matrix_to_collection
    tiledbsoma.io.create_from_matrix
    tiledbsoma.io.update_obs
    tiledbsoma.io.update_var
    tiledbsoma.io.update_matrix

.. rubric:: Growing a TileDB-SOMA Experiment

.. autosummary::
    :toctree: _autosummary/
    :nosignatures:

    tiledbsoma.io.register_anndatas
    tiledbsoma.io.register_h5ads
    tiledbsoma.io.append_X
    tiledbsoma.io.append_obs
    tiledbsoma.io.append_var
    tiledbsoma.io.show_experiment_shapes
    tiledbsoma.io.upgrade_experiment_shapes
    tiledbsoma.io.resize_experiment


Classes
-------

.. currentmodule:: tiledbsoma.io

.. autoclass:: ExperimentAmbientLabelMapping

   .. rubric:: Methods

   .. autosummary::
      :toctree: _generated

      ~ExperimentAmbientLabelMapping.prepare_experiment
      ~ExperimentAmbientLabelMapping.subset_for_anndata
      ~ExperimentAmbientLabelMapping.subset_for_h5ad
