# Set the Default Global Context

Configure a default [`SOMAContext`](SOMAContext.md) to be used by all
TileDB-SOMA operations when no explicit context is provided.

## Usage

``` r
set_default_context(config = NULL, replace = FALSE)
```

## Arguments

- config:

  ...

- replace:

  Allow replacing the existing default context with new configuration
  parameters.

## Value

Invisibly, the default default context object.

## Details

This function should be called once at the beginning of your session
before opening any SOMA objects if you want to customize the TileDB
context parameters that will apply to all subsequent operations.
Otherwise, a default context will be created automatically with standard
parameters when you first open a SOMA object.

If the global context was already set, an error will be raised unless
`replace=True`. Setting a new default context will not change the
context for TileDB-SOMA objects that were already created.
