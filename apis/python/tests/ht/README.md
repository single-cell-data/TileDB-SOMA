# SOMA Hypothesis-based tests

This folder contains Hypothesis-based tests and supporting code. All will run within the standard pytest
framework and will run in the course of normal pytest execution (including during CI).

## Recreating test failures

Property-based tests generate unique test invocations each time they are run. As such, test failures in CI
may not immediately reproduce if you re-run a test. The test log files should contain enough
information to reproduce the error, including the random number seed used for test selection, and the
so-called "blob" that allows Hypothesis to recreate a particular test.

The [Hypothesis documentation] has information on reproducing test failures.

## Exhaustive testing

The more test cases Hypothesis generates, the more likely it is to find a bug. This is at odds with
the need for our CI pipeline to complete in a "reasonable" amount of time.

The default configuration is suitable for use in CI, i.e., all tests will complete fairly quickly. Please do not
change this behavior.

In the course of development, it is often useful to more exhaustively search for test cases.
A Hypothesis profile has been defined for this case called `expensive`. You can run the tests in this
mode:

> pytest apis/python/tests/ --hypothesis-profile=expensive

In this mode, tests will run significantly longer (very roughly, 100X longer than the default) and cover
many more test conditions. Because each invocation of the test starts with a unique random seed, you
can repeat this invocation until you are satisfied with your test coverage.

## Configuration

The `_ht_test_config.py` file is used to configure the tests. The most common use case is a config flag
which enables a defect workaround, while the issue is being resolved.

## For More Information

See the [Hypothesis documentation].


[Hypothesis documentation]: https://hypothesis.readthedocs.io/
