# Get the Default SOMA Context

Retrieve the current default [`SOMAContext`](SOMAContext.md) used by
TileDB-SOMA operations.

## Usage

``` r
get_default_context()
```

## Value

The context that will be used for TileDB-SOMA API when no context is
provided by the user.

## Details

This function returns the context that was either:

- Explicitly set via [`set_default_context`](set_default_context.md), or

- Automatically created when a SOMA object was first created

An error is raised if no default context is set.
