* Most things are configured using GitHub Actions at `../../.github/workflows`
* Pre-push suggestions:
  * `black . tools/[a-z]*`
  * `isort . tools/[a-z]*`
  * `flake8 . tools/[a-z]*`
  * `python -m pytest tests`
* PyPI:
  * https://pypi.org/project/tiledbsc/
  * `tiledbinc` is an owner
  * See https://github.com/single-cell-data/TileDB-SingleCell/pull/178 for setups that were done
