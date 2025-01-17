"""Config settings for all Hypothesis tests. Primarily used to toggle bug work-arounds, etc.
"""

HT_TEST_CONFIG = {
    #
    # Defect work-arounds, while awaiting a fix
    #
    # data corruption due to incorrect Arrow array offset handling
    # See also sc-62104
    "sc-61239_workaround": True,
    # creating array with timestamp==0 fails in 1.15 (regression)
    "sc-61054_workaround": True,
    # Tables returned by SparseNDArray.read have incorrect nullability in schema fields
    "sc-61222_workaround": True,
    # SparseNDArray.read returns table with type==int64 when array schema has type==timestamp[us]
    "sc-61227_workaround": True,
    # reopen return mangled object
    "sc-61123_workaround": True,
    # reopen w->r loses all metadata modifications
    "sc-61118_workaround": True,
    # dataframe column names of \x00 silently mutated to empty Python string
    "sc-61291_workaround": True,
    # DataFrame.write creates 1+ fragments (one per table chunk)
    "sc-61462_workaround": True,
    # Core does not correctly de-dup 0. and -0. on float dimensions
    "sc-61506_workaround": True,
    # DenseNDArray can't read timestamps
    "sc-61743_workaround": True,
    # Read of new array returns incorrect info
    "sc-61676_workaround": True,
    #
    # Enable/disable partially implemented features
    #
    "allow_nullable": False,
}
