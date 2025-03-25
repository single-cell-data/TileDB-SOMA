"""Config settings for all Hypothesis tests. Primarily used to toggle bug workarounds, etc."""

HT_TEST_CONFIG = {
    #
    # Defect workarounds, while awaiting fixes
    #
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
    # index columns of type binary/large_binary are reported as large_string
    "sc-62236_workaround": True,
    # string index values starting with 0x7F barf
    "sc-62265_workaround": True,
    # dictionary of timestamps is not working
    "sc-62364_workaround": True,
    # string categories of value '' fail in write
    "sc-62447_workaround": True,
    # float categoricals fail with NaN
    "sc-62449_workaround": True,
    # enum handling fails with an enum of ''
    "sc-63404_workaround": True,
    # from_anndata fails on any high-D obms/varm/obsp/varp/uns
    "sc-63409_workaround": True,
    # Non-posix path names used as obsm/varm/obsp/varp keys fails
    "sc-63410_workaround": True,
    # soma.io.to_anndata incorrectly reads bool_ arrays in uns
    "sc-63447_workaround": True,
    # tiledbsoma.io.from_anndata ignores byte order
    "sc-63459_workaround": True,
    # path generated from key in collections
    "sc-63402_workaround": True,
    #
    # Enable/disable partially implemented features
    #
    "allow_nullable": False,
}
