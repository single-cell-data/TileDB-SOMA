import tiledb

ref = tiledb.open("ref")
print("ref schema: ")
print(ref.schema)

tgt = tiledb.open("tgt")
print("tgt schema: ")
print(tgt.schema)
