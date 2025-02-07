## ================================================================
### UDFs
# def remote_obs_schema(exp_uri):
#    import tiledbsoma
#    exp = tiledbsoma.Experiment.open(exp_uri)
#    return exp.obs.schema
# import tiledb.cloud
# import tiledb.cloud.udf
# tiledb.cloud.udf.exec(
#    remote_obs_schema,
#    soma_pbmc3k_uri,
# )
#
# def remote_query(exp_uri):
#    import tiledbsoma
#    exp = tiledbsoma.Experiment.open(exp_uri)
#
#    query = tiledbsoma.ExperimentAxisQuery(
#        experiment=exp,
#        measurement_name="RNA",
#        obs_query=tiledbsoma.AxisQuery(
#            value_filter="n_genes_by_counts > 1000",
#        ),
#        var_query=tiledbsoma.AxisQuery(
#            value_filter="n_cells_by_counts > 100",
#        ),
#    )
#
#    return (query.n_obs, query.n_vars)
# tiledb.cloud.udf.exec(
#    remote_query,    soma_pbmc3k_uri,
# )
#
## ================================================================
## Collection-mapper test
# from tiledb.cloud.taskgraphs import client_executor as executor
# soco_uri = 'tiledb://TileDB-Inc/stack-small-soco-staging'
# res = tiledb.cloud.udf.exec(
#    'TileDB-Inc/soma_experiment_collection_mapper',
#    soco_uri=soco_uri,
#    measurement_name="RNA",
#    X_layer_name="data",
#    # callback = lambda x: x.obs.shape,
#    # callback = lambda x: x,
#    callback = lambda adata: [adata.obs.shape, adata.var.shape, adata.X.shape],
#    # callback = lambda adata: adata.var,
#    args_dict={},
#    reducer = lambda x: x,
#    obs_attrs = ['obs_id', 'cell_type', 'is_primary_data'],
#    var_attrs = ['var_id', 'means'],
# )
# dag = executor.LocalExecutor(res, namespace = "TileDB-Inc")
# dag.visualize()
##%%time
# dag.execute()
# dag.wait()
# dag.node("output").result()

# * Show, upgrade, resize
