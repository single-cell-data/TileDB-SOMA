# Build `tiledbsoma` from source

Run the [../scripts/setup](../scripts/setup) script to execute these steps.

## Setup Python Environment (optional)

For example, create and activate a python venv:
```
python -m venv test/tiledbsoma
source test/tiledbsoma/bin/activate
```

## Clean

Remove old build files:
```
rm -rf build dist
rm -f apis/python/src/tiledbsoma/libtiledb.*
rm -f apis/python/src/tiledbsoma/libtiledbsoma.*
```

## Build

```
cd apis/python
pip install -e .
```

## Check

```
python -c "import tiledbsoma.libtiledbsoma; print(tiledbsoma.libtiledbsoma.version())"
```

## Test

Install test data:
```
cd ../..
mkdir -p test
./apis/python/tools/ingestor \
  --soco \
  -o test/soco \
  -n \
  data/pbmc3k_processed.h5ad \
  data/10x-pbmc-multiome-v1.0/subset_100_100.h5ad
```

Run pytests:
```
pytest -v --durations=0 libtiledbsoma
```

## Notes

### C++ Prerequisites

A build system with `libfmt` and `libspdlog` installed may conflict with the required versions for `tiledbsoma`. If that is the case, this configuration is known-good on Ubuntu 22.04:

```
sudo apt install libfmt-dev libspdlog-dev
```

Check:

```
$ dpkg -l | egrep "lib(spdlog|fmt)" | cut -c-80
ii  libfmt-dev:amd64                  8.1.1+ds1-2                             am
ii  libfmt8:amd64                     8.1.1+ds1-2                             am
ii  libspdlog-dev:amd64               1:1.9.2+ds-0.2                          am
ii  libspdlog1:amd64                  1:1.9.2+ds-0.2                          am
```
