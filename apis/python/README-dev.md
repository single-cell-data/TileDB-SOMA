* Most things are configured using GitHub Actions at `../../.github/workflows`
* Pre-push suggestions:
  * `black . tools/[a-z]*`
  * `isort . tools/[a-z]*`
  * `flake8 . tools/[a-z]*`
  * `mypy .`
  * `python -m pytest tests`
* And/or you can follow https://github.com/single-cell-data/TileDB-SOMA/pull/193:
  * `pip install pre-commit && pre-commit install`
  * Then on every `git commit` all configured linters will be run locally on the changed files, and the commit will be blocked if there is any error.
* PyPI:
  * https://pypi.org/project/tiledbsoma/
  * `tiledbinc` is an owner
  * See https://github.com/single-cell-data/TileDB-SOMA/pull/178 for setups that were done
