# Prerequisites

Known-good on Ubuntu 22.04:

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

# Setup

This populates test data:

```
../scripts/setup
```

# Build

```
../scripts/bld
../scripts/test
```

Without `../scripts/bld`:

```
cd libtiledbsc
rm -rvf build
mkdir build
cd build
cmake ..
export pybind11_DIR=$(python -m pybind11 --cmakedir)
make -j 16
```

Then in Python `import tiledb.libtiledb`.
