import tiledbsoma

# This test is extracted from tutorial_soma_objects.ipynb notebook
# To test the profile run: `python main.py python tests/objects.py` in the profile folder
experiment = tiledbsoma.open("../apis/python/notebooks/data/dense/pbmc3k")

obs = experiment.obs
print(f"obs schema {obs.schema}")
p1 = obs.read().concat().to_pandas()
p2 = obs.read((slice(0, 10),), column_names=["obs_id", "n_genes"]).concat().to_pandas()
p3 = obs.read((slice(None),), value_filter="n_genes > 1500").concat().to_pandas()

print(p1)
print(p2)
print(p3)

raw = experiment.ms["raw"]
X = experiment["ms"]["RNA"].X["data"]

print(X.schema)
print(X.read().to_numpy())
sliced_X = X.read((slice(0, 9),)).to_numpy()
print(sliced_X)
print(sliced_X.shape)

var = experiment.ms["RNA"].var
idx = var.read(value_filter="var_id == 'ICOSLG'").concat()["soma_joinid"].to_numpy()

print(X.read((None, int(idx[0]))).to_numpy())

experiment = tiledbsoma.open("../apis/python/notebooks/data/sparse/pbmc3k")
X = experiment.ms["RNA"].X["data"]

print(X.schema)
print(X.shape)
print(X.nnz)
print(X.read)

print("Extract tensor!")

tensor = X.read().coos().concat()
tensor.to_scipy()

print("SlicedX!")

sliced_X = X.read((slice(0, 9),)).coos().concat().to_scipy()
print(sliced_X)
print(sliced_X.nonzero()[0])
