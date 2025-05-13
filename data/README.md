# TileDB-SOMA Data

This folder contains data for use in tests and examples.

## How to add new data

### Check new data and upload to GitHub for hosting

Temporarily add a copy of your data to this directory in the way you will eventually want to use it. Create your new tests or examples and check that the data serves your purpose.

Once you have verified the data works as intended, add it to a [TileDB-SOMA-Test-Data](https://github.com/single-cell-data/TileDB-SOMA-Test-Data) release as described in that project's [README](https://github.com/single-cell-data/TileDB-SOMA-Test-Data/blob/main/README.md).

### Update `make data` to include new dataset

The Makefile `data` target calls [scripts/prepare-test-data.sh](../scripts/prepare-test-data.sh). You will need to update the script to include a section to prepare your new data. Your update must do the following:

1. Create a new file or folder in the `data/` directory (this directory) that includes the desired data if the data has not already been prepared.
1. Gracefully skip preparing data that has already been added.

If you are including data from TileDB-SOMA-Test-Data, the easiest way to access the data is to download all necessary to a new directory using "wget". Extract any data that is compressed.

### Update `make clean_data` to remove new dataset

The Makefile `clean_data` target calls [scripts/clean-test-data.sh](../scripts/clean-test-data). You will need to update the script to include a section that removes your new data.

### Checks

Before committing your changes to `scripts/prepare-test-data.sh` verify the following:

1. The command `make data` fully prepares your data so that if can be easily accessed in your new test and/or example.
1. The command `make clean_data` fully removes your data.
1. The data is being ignored by Git.

It is a good idea to include you new test or example in the same PR that adds the new data to help ensure you can check it has been successfully prepared.
