# SOMA Hypothesis-based tests

This folder contains Hypothesis-based tests and supporting code. All will run within the standard pytest
framework and will run in the course of normal pytest execution.

## Configuration

The default configuration is suitable for use in CI, i.e., run fairly quickly. Please do not
change this behavior.

In the course of development, it is often useful to more exhaustively search for test cases.
A Hypothesis profile has been defined for this case called `expensive`. You can run the tests in this
mode:

> pytest tests/ --hypothesis-profile=expensive

## For More Information

See the [Hypothesis documentation](https://hypothesis.readthedocs.io/)
