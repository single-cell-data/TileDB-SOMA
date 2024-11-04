The tiledbsoma module
=====================

.. currentmodule:: tiledbsoma

.. automodule:: tiledbsoma

Classes
-------

.. autosummary::
    :toctree: _autosummary/
    :nosignatures:

    tiledbsoma.Collection
    tiledbsoma.Experiment
    tiledbsoma.Measurement

    tiledbsoma.DataFrame
    tiledbsoma.SparseNDArray
    tiledbsoma.SparseNDArrayRead
    tiledbsoma.DenseNDArray

    tiledbsoma.Axis
    tiledbsoma.CoordinateSpace
    tiledbsoma.MultiscaleImage
    tiledbsoma.PointCloudDataFrame
    tiledbsoma.Scene

    tiledbsoma.ResultOrder

    tiledbsoma.AxisColumnNames
    tiledbsoma.AxisQuery
    tiledbsoma.ExperimentAxisQuery

    tiledbsoma.SOMATileDBContext
    tiledbsoma.TileDBCreateOptions
    tiledbsoma.TileDBWriteOptions

    tiledbsoma.IntIndexer
    tiledbsoma.tiledbsoma_build_index

Exceptions
----------

.. autosummary::
    :toctree: _autosummary/
    :nosignatures:

    tiledbsoma.SOMAError
    tiledbsoma.DoesNotExistError
    tiledbsoma.AlreadyExistsError
    tiledbsoma.NotCreateableError

Functions
---------

.. autosummary::
    :toctree: _autosummary/
    :nosignatures:

    tiledbsoma.open

    tiledbsoma.show_package_versions

    tiledbsoma.get_SOMA_version
    tiledbsoma.get_implementation
    tiledbsoma.get_implementation_version
    tiledbsoma.get_storage_engine

    tiledbsoma.tiledbsoma_stats_enable
    tiledbsoma.tiledbsoma_stats_disable
    tiledbsoma.tiledbsoma_stats_reset
    tiledbsoma.tiledbsoma_stats_dump
    tiledbsoma.tiledbsoma_stats_as_py
    tiledbsoma.tiledbsoma_stats_json
