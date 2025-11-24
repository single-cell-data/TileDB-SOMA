# Create and cache a SOMA Context Object

Create and cache a SOMA Context Object

## Usage

``` r
soma_context(config)
```

## Arguments

- config:

  A named character vector with “key” and “value” pairs defining the
  configuration setting.

## Value

An external pointer object containing a shared pointer instance of
`SOMAContext`.

## Examples

``` r
head(cfgvec <- as.vector(tiledb::tiledb_config())) # TileDB config as a vector
#> config.env_var_prefix config.logging_format  config.logging_level 
#>             "TILEDB_"             "DEFAULT"                   "0" 
#> filestore.buffer_size           profile_dir          profile_name 
#>           "104857600"                    ""                    "" 
(sctx <- soma_context(cfgvec))
#> <pointer: 0x556b9fa0c8a0>
```
