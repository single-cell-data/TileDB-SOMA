
## DataFrame state machine

* [DONE] read/write/count and data ledger
* read with indexing
* read with value_filter
* nullable columns/attributes (not index)
  * What happens if index column write gets an array with an null position/validity buffer?
  * What happens if array has null, but schema doesn't (for soma_data or a DataFrame column)?
  * What happens if we read it back - do we get the right validity buffer?
* read result_order
* time-travel
* [DONE] variable-length types (string, binary)
* [DONE] enums/dicts
* Add/drop columns (schema evolution)
* update_foo methods

## DenseArray state machine

* [DONE] read/write
* [DONE] base/open/close/metadata
* [DONE] read result_order
* [DONE] read with indexing (current is "all")
* Resize is turned OFF due to bugs - re-enable
* time-travel
* nullable - what happens if we send array to write that contains an invalid/null value? (i.e. w/ validity buffer). Expect: error or no write?

## SparseNDArray

* [DONE] NaNs stored in data
* read with indexing (currently is "all")
* [DONE] maxshape?
* [BLOCKED] how does shape interact with time travel? (currently all reshapes happen at 'now')
* read result_order
* time-travel
* nullable - what happens if we send array to write that contains an invalid/null value? (i.e. w/ validity buffer). Expect: error or no write?

## Array metadata

* [DONE] are nulls allowed in key or value?  NOT ALLOWED IN EITHER
* [DONE] non-ASCII code points in key or value?  ALLOWED IN BOTH
* time travel (use a ledger)
* [DONE] NaNs stored in metadata (currently disabled simply for ease of testing equality)

## Other

* [DONE] DataFrame state machine
* [DONE] DenseNDArray state machine
* [DONE] Generalize/abstract state machine (ArrayStateMachine) which is sub-classed by others
* [DONE] make a PR
* [DONE] CI
* Groups
* fuzz value_filter parser
* [BLOCKED] Time travel for schema, metadata and data
* soma.io to/from and fuzzing (perfect case for just round-tripping and testing equality)
* schema evolution (beyond shape)
* out of order time travel (write earlier than most recent write)

## Fastercsx

* to_numpy/dense
* from_soma

## soma.io

* [DONE] to/from round-trip
