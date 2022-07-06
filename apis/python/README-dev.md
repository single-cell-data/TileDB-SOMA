* Most things are configured using GitHub Actions at `../../.github/workflows`
* Pre-push suggestions:
  * `black . tools/[a-z]*`
  * `isort . tools/[a-z]*`
  * `flake8 . tools/[a-z]*`
  * `mypy .`
  * `python -m pytest tests`
* And/or you can follow https://github.com/single-cell-data/TileDB-SingleCell/pull/193:
  * `pip install pre-commit && pre-commit install`
  * Then on every `git commit` all configured linters will be run locally on the changed files, and the commit will be blocked if there is any error.
* PyPI:
  * https://pypi.org/project/tiledbsc/
  * `tiledbinc` is an owner
  * See https://github.com/single-cell-data/TileDB-SingleCell/pull/178 for setups that were done
* Public-bucket status quo:
  * `s3://tiledb-singlecell-data`

```
aws-login
cd ~/tiledb-singlecell-data
aws --profile prod-admin s3 sync . s3://tiledb-singlecell-data
```

  * `s3://tiledb-singlecell-docs`

```
aws-login
cd ~/git/single-cell-data/TileDB-SingleCell/apis/python
sh mkmd.sh
cd ~/git/single-cell-data/TileDB-SingleCell
quarto render
aws --profile prod-admin s3 sync docs s3://tiledb-singlecell-docs/docs
reload and double-check https://tiledb-singlecell-docs.s3.amazonaws.com/docs/overview.html
```
