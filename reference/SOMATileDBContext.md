# SOMA TileDB Context

Context map for TileDB-backed SOMA objects

## Super classes

[`tiledbsoma::MappingBase`](MappingBase.md) -\>
[`tiledbsoma::ScalarMap`](ScalarMap.md) -\>
[`tiledbsoma::SOMAContextBase`](SOMAContextBase.md) -\>
`SOMATileDBContext`

## Methods

### Public methods

- [`SOMATileDBContext$new()`](#method-SOMATileDBContext-new)

- [`SOMATileDBContext$keys()`](#method-SOMATileDBContext-keys)

- [`SOMATileDBContext$items()`](#method-SOMATileDBContext-items)

- [`SOMATileDBContext$length()`](#method-SOMATileDBContext-length)

- [`SOMATileDBContext$get()`](#method-SOMATileDBContext-get)

- [`SOMATileDBContext$set()`](#method-SOMATileDBContext-set)

- [`SOMATileDBContext$to_tiledb_context()`](#method-SOMATileDBContext-to_tiledb_context)

- [`SOMATileDBContext$context()`](#method-SOMATileDBContext-context)

- [`SOMATileDBContext$clone()`](#method-SOMATileDBContext-clone)

Inherited methods

- [`tiledbsoma::MappingBase$print()`](MappingBase.html#method-print)
- [`tiledbsoma::MappingBase$remove()`](MappingBase.html#method-remove)
- [`tiledbsoma::MappingBase$setv()`](MappingBase.html#method-setv)
- [`tiledbsoma::MappingBase$to_list()`](MappingBase.html#method-to_list)
- [`tiledbsoma::MappingBase$update()`](MappingBase.html#method-update)
- [`tiledbsoma::MappingBase$values()`](MappingBase.html#method-values)

------------------------------------------------------------------------

### Method `new()`

#### Usage

    SOMATileDBContext$new(config = NULL, cached = TRUE)

#### Arguments

- `config`:

  ...

- `cached`:

  Force new creation

#### Returns

An instantiated `SOMATileDBContext` object

------------------------------------------------------------------------

### Method `keys()`

#### Usage

    SOMATileDBContext$keys()

#### Returns

The keys of the map

------------------------------------------------------------------------

### Method `items()`

#### Usage

    SOMATileDBContext$items()

#### Returns

Return the items of the map as a list

------------------------------------------------------------------------

### Method [`length()`](https://rdrr.io/r/base/length.html)

#### Usage

    SOMATileDBContext$length()

#### Returns

The number of items in the map

------------------------------------------------------------------------

### Method [`get()`](https://rdrr.io/r/base/get.html)

#### Usage

    SOMATileDBContext$get(key, default = quote(expr = ))

#### Arguments

- `key`:

  Key to fetch

- `default`:

  Default value to fetch if `key` is not found; defaults to `NULL`

#### Returns

The value of `key` in the map, or `default` if `key` is not found

------------------------------------------------------------------------

### Method `set()`

#### Usage

    SOMATileDBContext$set(key, value)

#### Arguments

- `key`:

  Key to set

- `value`:

  Value to add for `key`, or `NULL` to remove the entry for `key`

#### Returns

\[chainable\] Invisibly returns `self` with `value` added as `key`

------------------------------------------------------------------------

### Method `to_tiledb_context()`

#### Usage

    SOMATileDBContext$to_tiledb_context()

#### Returns

A
[`tiledb_ctx`](https://tiledb-inc.github.io/TileDB-R/reference/tiledb_ctx.html)
object, dynamically constructed. Most useful for the constructor of this
class.

------------------------------------------------------------------------

### Method `context()`

#### Usage

    SOMATileDBContext$context()

#### Returns

A
[`tiledb_ctx`](https://tiledb-inc.github.io/TileDB-R/reference/tiledb_ctx.html)
object, which is a stored (and long-lived) result from
`to_tiledb_context`.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SOMATileDBContext$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
(ctx <- SOMATileDBContext$new())
#> Warning: SOMATileDBContext$new() was deprecated in tiledbsoma 2.3.0.
#> ℹ Use `SOMAContext` instead.
#> ℹ The deprecated feature was likely used in the R6 package.
#>   Please report the issue at <https://github.com/r-lib/R6/issues>.
#> <SOMATileDBContext>
#>   config.env_var_prefix: TILEDB_
#>   config.logging_format: DEFAULT
#>   config.logging_level: 0
#>   filestore.buffer_size: 104857600
#>   profile_dir: 
#>   profile_name: 
#>   rest.capnp_traversal_limit: 2147483648
#>   rest.curl.buffer_size: 524288
#>   rest.curl.retry_errors: true
#>   rest.curl.tcp_keepalive: true
#>   rest.curl.verbose: false
#>   rest.http_compressor: any
#>   rest.load_enumerations_on_array_open: false
#>   rest.load_enumerations_on_array_open_all_schemas: false
#>   rest.load_metadata_on_array_open: true
#>   rest.load_non_empty_domain_on_array_open: true
#>   rest.payer_namespace: 
#>   rest.retry_count: 25
#>   rest.retry_delay_factor: 1.25
#>   rest.retry_http_codes: 503
#>   rest.retry_initial_delay_ms: 500
#>   rest.server_address: https://api.tiledb.com
#>   rest.server_serialization_format: CAPNP
#>   rest.use_refactored_array_open: true
#>   rest.use_refactored_array_open_and_query_submit: true
#>   sm.allow_separate_attribute_writes: false
#>   sm.allow_updates_experimental: false
#>   sm.check_coord_dups: true
#>   sm.check_coord_oob: true
#>   sm.check_global_order: true
#>   sm.compute_concurrency_level: 4
#>   sm.consolidation.amplification: 1.0
#>   sm.consolidation.buffer_size: 50000000
#>   sm.consolidation.max_fragment_size: 18446744073709551615
#>   sm.consolidation.mode: fragments
#>   sm.consolidation.purge_deleted_cells: false
#>   sm.consolidation.step_max_frags: 4294967295
#>   sm.consolidation.step_min_frags: 4294967295
#>   sm.consolidation.step_size_ratio: 0.0
#>   sm.consolidation.steps: 4294967295
#>   sm.consolidation.timestamp_end: 18446744073709551615
#>   sm.consolidation.timestamp_start: 0
#>   sm.dedup_coords: false
#>   sm.enable_signal_handlers: true
#>   sm.encryption_key: 
#>   sm.encryption_type: NO_ENCRYPTION
#>   sm.enumerations_max_size: 10485760
#>   sm.enumerations_max_total_size: 52428800
#>   sm.fragment_info.preload_mbrs: false
#>   sm.group.timestamp_end: 18446744073709551615
#>   sm.group.timestamp_start: 0
#>   sm.io_concurrency_level: 4
#>   sm.max_tile_overlap_size: 314572800
#>   sm.mem.consolidation.buffers_weight: 1
#>   sm.mem.consolidation.reader_weight: 3
#>   sm.mem.consolidation.writer_weight: 2
#>   sm.mem.malloc_trim: true
#>   sm.mem.reader.sparse_global_order.ratio_array_data: 0.3
#>   sm.mem.reader.sparse_global_order.ratio_coords: 0.5
#>   sm.mem.reader.sparse_global_order.ratio_tile_ranges: 0.1
#>   sm.mem.reader.sparse_unordered_with_dups.ratio_array_data: 0.1
#>   sm.mem.reader.sparse_unordered_with_dups.ratio_coords: 0.5
#>   sm.mem.reader.sparse_unordered_with_dups.ratio_tile_ranges: 0.1
#>   sm.mem.tile_upper_memory_limit: 1073741824
#>   sm.mem.total_budget: 10737418240
#>   sm.memory_budget: 5368709120
#>   sm.memory_budget_var: 10737418240
#>   sm.merge_overlapping_ranges_experimental: true
#>   sm.partial_tile_offsets_loading: false
#>   sm.query.condition_evaluator: ast
#>   sm.query.dense.qc_coords_mode: false
#>   sm.query.dense.reader: refactored
#>   sm.query.sparse_global_order.preprocess_tile_merge: 32768
#>   sm.query.sparse_global_order.reader: refactored
#>   sm.query.sparse_unordered_with_dups.reader: refactored
#>   sm.read_range_oob: warn
#>   sm.skip_checksum_validation: false
#>   sm.skip_est_size_partitioning: false
#>   sm.skip_unary_partitioning_budget_check: false
#>   sm.vacuum.mode: fragments
#>   sm.var_offsets.bitsize: 64
#>   sm.var_offsets.extra_element: false
#>   sm.var_offsets.mode: bytes
#>   ssl.ca_file: 
#>   ssl.ca_path: 
#>   ssl.verify: true
#>   vfs.azure.blob_endpoint: 
#>   vfs.azure.block_list_block_size: 5242880
#>   vfs.azure.is_data_lake_endpoint: 
#>   vfs.azure.max_parallel_ops: 4
#>   vfs.azure.max_retries: 5
#>   vfs.azure.max_retry_delay_ms: 60000
#>   vfs.azure.retry_delay_ms: 800
#>   vfs.azure.storage_account_key: 
#>   vfs.azure.storage_account_name: 
#>   vfs.azure.storage_sas_token: 
#>   vfs.azure.use_block_list_upload: true
#>   vfs.file.posix_directory_permissions: 755
#>   vfs.file.posix_file_permissions: 644
#>   vfs.gcs.endpoint: 
#>   vfs.gcs.impersonate_service_account: 
#>   vfs.gcs.max_direct_upload_size: 10737418240
#>   vfs.gcs.max_parallel_ops: 4
#>   vfs.gcs.multi_part_size: 5242880
#>   vfs.gcs.project_id: 
#>   vfs.gcs.request_timeout_ms: 3000
#>   vfs.gcs.service_account_key: 
#>   vfs.gcs.use_multi_part_upload: true
#>   vfs.gcs.workload_identity_configuration: 
#>   vfs.log_operations: false
#>   vfs.max_batch_size: 104857600
#>   vfs.min_batch_gap: 512000
#>   vfs.min_batch_size: 20971520
#>   vfs.min_parallel_size: 10485760
#>   vfs.read_ahead_cache_size: 10485760
#>   vfs.read_ahead_size: 102400
#>   vfs.read_logging_mode: 
#>   vfs.s3.aws_access_key_id: 
#>   vfs.s3.aws_external_id: 
#>   vfs.s3.aws_load_frequency: 
#>   vfs.s3.aws_role_arn: 
#>   vfs.s3.aws_secret_access_key: 
#>   vfs.s3.aws_session_name: 
#>   vfs.s3.aws_session_token: 
#>   vfs.s3.bucket_canned_acl: NOT_SET
#>   vfs.s3.ca_file: 
#>   vfs.s3.ca_path: 
#>   vfs.s3.config_source: auto
#>   vfs.s3.connect_max_tries: 5
#>   vfs.s3.connect_scale_factor: 25
#>   vfs.s3.connect_timeout_ms: 10800
#>   vfs.s3.endpoint_override: 
#>   vfs.s3.install_sigpipe_handler: true
#>   vfs.s3.logging_level: Off
#>   vfs.s3.max_parallel_ops: 4
#>   vfs.s3.multipart_part_size: 5242880
#>   vfs.s3.no_sign_request: false
#>   vfs.s3.object_canned_acl: NOT_SET
#>   vfs.s3.proxy_host: 
#>   vfs.s3.proxy_password: 
#>   vfs.s3.proxy_port: 0
#>   vfs.s3.proxy_scheme: http
#>   vfs.s3.proxy_username: 
#>   vfs.s3.region: 
#>   vfs.s3.request_timeout_ms: 3000
#>   vfs.s3.requester_pays: false
#>   vfs.s3.scheme: https
#>   vfs.s3.skip_init: false
#>   vfs.s3.sse: 
#>   vfs.s3.sse_kms_key_id: 
#>   vfs.s3.storage_class: NOT_SET
#>   vfs.s3.use_multipart_upload: true
#>   vfs.s3.use_virtual_addressing: true
#>   vfs.s3.verify_ssl: true
ctx$get("sm.mem.reader.sparse_global_order.ratio_array_data")
#> sm.mem.reader.sparse_global_order.ratio_array_data 
#>                                              "0.3" 

ctx$to_tiledb_context()
#> An object of class "tiledb_ctx"
#> Slot "ptr":
#> <pointer: 0x55cf3ab35330>
#> 
```
