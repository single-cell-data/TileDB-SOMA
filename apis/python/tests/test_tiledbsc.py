import anndata

import tiledbsc as api


def test_group():
    # mapping between layer name and (array, attribute)
    normalized_X_layers = {
        "counts": ("counts_array", "counts"),
        "normed": ("normed_array", "normed"),
    }

    # compose a sc_group
    obs_array_name = "obs_array_name"
    normalized_group = api.Group(
        obs=obs_array_name,
        obs_index="obs_dim name",
        var=obs_array_name,
        var_index="var_dim_name",
        X=normalized_X_layers,
    )


# def test_create():
#     # Compose sc_group from a TileDB group containing required arrays.
#     # Validates that the triple <obs, var, X> are all consistent &
#     # compatible.
#     #
#     # Example validation/constraints:
#     # 	* indexing dims exist and are of same type, etc
#     # 	* all arrays are in the group
#     # 	* all X layer defn's are valid
#     # etc.

#     # mapping between layer name and (array, attribute)
#     array_name = "array_name"
#     normalized_X_layers = {
#         "counts": (array_name, "counts"),
#         "normed": (array_name, "normed"),
#     }

#     # compose a sc_group
#     obs_array_name = "obs_array_name"
#     normalized_group = api.Group(
#         obs=obs_array_name,
#         obs_index="obs_dim name",
#         var=obs_array_name,
#         var_index="var_dim_name",
#         X=normalized_X_layers,
#     )

#     unfiltered_group = api.Group(
#         obs=obs_array_name,
#         obs_index="obs_dim name",
#         var=obs_array_name,
#         var_index="var_dim_name",
#         X={"counts": (array_name, "counts")},
#     )

#     # compose a sc_super_group
#     super_group_uri = "test_create"
#     api.create(
#         super_group_uri,
#         groups={"normalized": normalized_group, "unfiltered": unfiltered_group},
#         default_group="normalized",
#     )


# def test_import_export():
#     # Create from in-memory AnnData/ScanPy. Suggest a similar approach
#     # for R ecosystem (Seurat & SingleCellExperiment)

#     adata = anndata.read_h5ad("file.h5ad")
#     api.create(<group_uri>, anndata=adata)

#     # To AnnData in-memory object.
#     # Similar approach for R (Seurat and SingleCellExperiment)
#     # Requires further elaboration on conventions - some initial
#     # thoughts in a later section.
#     adata = sc_group.to_anndata()

# def test_shape():
#     # load from URI.  If <uri> is a supergroup, will return the default
#     # sc_group in the sc_super_group.  Also works as a Python context
#     # manager.
#     # This open function, similar to the functionality supported in
#     # TileDB, can take as arguments timestamps to support versioning
#     # and time traveling.
#     sc_group = api.open(<uri>)

#     # shape
#     sc_group.shape  # shared across obs, var and X

#     # closes all arrays in group
#     sc_group.close()

# def test_query():
#     with api.open(<uri>) as sc_group:
#     sc_group_slice = sc_group[
#     <obs_filter_condition>,
#     <var_filter_condition>
#     ]

#     # or, to configure query (eg, incremental results, constrain result
#     # to partial attributes, etc).
#     with api.open(<uri>) as sc_group:
#         sc_group_slice = sc_group.query(
#             return_incomplete=True,
#             <option>=<value>, ...
#         )[<obs_filter_condition>, <var_filter_condition>]

#         for partial_result in sc_group_slice:
#             ...

#     # export slice as in-mem AnnData object
#     adata = sc_group_slice.to_anndata()

#     # export obs/var query as dataframe or dict-like
#     sc_group_slice.obs.to_dict() # return multi-dict
#     sc_group_slice.obs.to_df() # return df indexed by var label
#     sc_group_slice.var.to_dict() # return multi-dict
#     sc_group_slice.var.to_df() # return df indexed by var label

#     # export X or aux arrays
#     # dense
#     sc_group_slice.X["layer"].to_array() # dense ndarray
#     sc_group_slice.aux['arrayname'].to_array() # dense ndarray
#     # sparse - dataframe, with 2D index (obs_label, var_label)
#     # or dict-like (COO): {
#     #    "obs_labels": 1d_ndarray,
#     #    "var_labels": 1d_ndarray,
#     #    "layer_data": 1d_ndarray
#     # }
#     sc_group_slice.X["layer"].to_df()
#     sc_group_slice.X["layer"].to_dict()
#     sc_group_slice.aux['arrayname'].to_df()
#     sc_group_slice.aux['arrayname'].to_dict()

# def test_subarray_access():
#     # Accessing sc_group subarray: X, obs, var.
#     # Returns TileDB array
#     sc_group.X["layer"]  # return opened array
#     sc_group.var  # return opened array
#     sc_group.obs  # return opened array
#     # Other auxiliary arrays associated with the sc_group
#     sc_group.aux['array name']  # return opened array
