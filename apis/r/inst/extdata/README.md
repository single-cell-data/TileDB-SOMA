How to recreate these data files:

* `soma-exp-pbmc-small.tar.gz`
  * `data-raw/create-soma-exp-pbmc-small.R`
* `soma-dataframe-pbmc3k-processed-obs.tar.gz`
  * `data-raw/create-soma-dataframe-pbmc3k-processed-obs.R`

Note for unit tests: if I unit-test case compares the length of `list_datasets()` to the length of `dir(example_data_dir())`, then that unit-test case should subtract 1 to account for this `README.md` file itself.
