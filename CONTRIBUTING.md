# Contributing to TileDB-SOMA

Thanks for your interest in TileDB-SOMA. The notes below give some pointers for filing issues and bug reports, or contributing to the code.

## Contribution Checklist
- Reporting a bug? Please include the following information
  - operating system and version (windows, linux, macos, etc.)
  - TileDB and TileDB-SOMA version (for example, the output of `conda list` or `git status`).
  - if possible, a minimal working example demonstrating the bug or issue (along with any data to re-create, when feasible)
- Please paste code blocks with triple backquotes (```) so that github will format it nicely. See [GitHub's guide on Markdown](https://guides.github.com/features/mastering-markdown) for more formatting tricks.

## Contributing Code
*By contributing code to TileDB-SOMA, you are agreeing to release it under the [MIT License](https://github.com/TileDB-Inc/TileDB/tree/dev/LICENSE).*

### Contribution Workflow

- [Please follow these instructions to build from source](https://docs.tiledb.com/developer/tiledbsoma/installation)
- Make changes locally, then rebuild as appropriate for the level of changes (e.g.: `make` for `libtilebsc` or `python setup.py develop` for `apis/python`).
- Make sure to run `make check`, or `pytest` to verify changes against tests (add new tests where applicable).
- Please submit [pull requests](https://help.github.com/en/desktop/contributing-to-projects/creating-a-pull-request) against the default [`main` branch of TileDB-SOMA](https://github.com/TileDB-Inc/TileDB-SOMA/tree/master).
- If you edit the Python files, please run the pre-commit hooks

  ```sh
  python -m venv ./pre-commit
  source ./pre-commit/bin/activate
  python -m pip -v install pre-commit
  pre-commit run -a -v
  deactivate
  ```
