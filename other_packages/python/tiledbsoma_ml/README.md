
# tiledbsoma_ml

A Python package containing ML tools for use with `tiledbsoma`.

## Description

The package currently contains a prototype PyTorch `IterableDataset` for use with the
[`torch.utils.data.DataLoader`](https://pytorch.org/docs/stable/data.html#torch.utils.data.DataLoader)
API.

## Getting Started

### Installing

Install using your favorite package installer.  For exapmle, with pip:

> pip install tiledbsoma-ml

Developers may install editable, from source, in the usual manner:

> pip install -e .

### Documentation

TBD

## Builds

This is a pure Python package. To build a wheel, ensure you have the `build` package installed, and then:

> python -m build .

## Version History

See the [CHANGELOG.md](CHANGELOG.md) file.

## License

This project is licensed under the MIT License.

## Acknowledgements

The SOMA team is grateful to the Chan Zuckerberg Initiative Foundation [CELLxGENE Census](https://cellxgene.cziscience.com)
team for their initial contribution.
