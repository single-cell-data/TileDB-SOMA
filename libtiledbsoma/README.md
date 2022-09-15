# Build `tiledbsoma` from source

## Setup Python Environment (optional)

For example, create and activate activate a python venv:
```
python -m venv test/tiledbsoma-venv
source test/tiledbsoma-venv/bin/activate
```

## Install prerequisites:
```
python -m pip install pybind11 pytest
```

## Build

```
cd apis/python
python -m pip install .
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
