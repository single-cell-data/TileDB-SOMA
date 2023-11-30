import tiledb

def foo(uri):
    print("AAA300", uri)
    tiledb.open(uri)
    print("AAA301")
