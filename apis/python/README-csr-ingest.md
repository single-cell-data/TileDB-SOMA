# Overview

As of 2022-05-03 we offer a chunked ingestor for `X` data within larger data files. This enables
ingestion of larger data files within RAM limits of available hardware. It entails a
read-performance penalty under certain query modes, as this note will articulate.

# What it and isn't streamed/chunked

Input `.h5ad` files are read into memory using `anndata.read_h5ad` -- if the input file is 5GB, say,
about that much RAM will be needed. We don't have a way to stream the contents of the `.h5ad` file
itself -- it's read all at once.

Often `X` data is in CSR format; occasionally CSC, or dense (`numpy.ndarray`).  Suppose for the rest
of this note that we're looking at CSR `X` data -- a similar analysis will hold (_mutatis mutandis_)
for CSC data.

Given CSR `X` data, we find that an all-at-once `x.toarray()` can involve a huge explosion of memory
requirement if the input data is sparse -- for this reason, we don't do that; we use `x.tocoo()`. In
summary, `.toarray()` offers a variably huge explosion from on-disk `X` size to in-memory densified,
and we don't do this.

Given CSR `X` data, we find that an all-at-once `x.tocoo()` involves about a 2x or 2.5x expansion in
RSS as revealed by `htop` -- CSR data from disk (and in RAM) is a list of contiguous row
sub-sequences with row-subsequence values all spelled out per array cell, but only bounding column
dimensions written out; COO data is a list of `(i,j,v)` tuples with the `i,j` written out
individually -- which of course takes more memory. In summary, all-at-once `.tocoo()` has a memory
increase from on-disk size to in-memory COO-ified, but with a lower (and more predictable)
multiplication factor.

The alternative discussed here applies for `.h5ad` data files which are small enough to read into
RAM, but for which the 2.5x-or-so inflation from CSR to COO results in a COO matrix which is too big
for RAM.

# Sketch, and relevant features of TileDB storage

What we will do is take chunks of the CSR -- a few rows at a time -- and convert each CSR submatrix
to COO, writing each "chunk" as a TileDB fragment.  This way the 2.5x memory expansion is paid only
from CSR submatrix to COO submatrix, and we can lower the memory footprint needed for the ingestion
into TileDB.

Some facts about this:

* In the `.h5ad` we have `obs`/`var` names mapping from string to int, and integer-indexed sparse/dense `X` matrices.
* In TileDB, by contrast, we have the `obs`/`var` names being _themselves_ string indices into sparse `X` matrices.
* TileDB storage orders its dims. That means that if you have an input matrix as on the left, with `obs_id=A,B,C,D` and `var_id=S,T,U,V`, then it will be stored as on the right:

```

Input row labels from AnnData `obs`:
[ C A B D ]
  0 1 2 3
 
Input column labels from AnnData `var`:
[ T V S U ]
  0 1 2 3
  
   Input CSR                                   TileDB storage
from AnnData `X`
  ---------                                    --------------  all one fragment in this example
    0 1 2 3                                      S T U V       -- for larger fragmented X, see below
  0 1 2 . .                                    A 4 . . 3
  1 . 3 4 .                                    B: 5 . 6 .
  2 . . 5 6                                    C . 1 . 2
  3 7 . 8 .                                    D 8 7 . .
```

* TileDB storage is 3-level: _fragments_ (corresponding to different timestamped writes); _tiles_; and _cells_.
* Fragments and tiles both have MBRs. For this example (suppose for the moment that is it's written all at once in a single fragment) the fragment MBR is `A..D` in the `obs_id` dimension and `S..V` in the `var_id` dimension.
* Query modes: we expect queries by `obs_id,var_id` pairs, or by `obs_id`, or by `var_id`. Given the above representation, since tiles within the fragment are using ordered `obs_id` and `var_id`, then all three query modes will be efficient:
  * there's one fragment
  * Queries on `obs_id,var_id` will locate only one tile within the fragment
  * Queries on `obs_id` will locate one row of files within the fragment
  * Queries on `var_id` will locate one column of files within the fragment

```
  TileDB storage
  --------------  all one fragment
    S T : U V
  A 4 . : . 3
  B: 5 .:  6 .
 .......: ...... tile boundary
  C . 1 : . 2
  D 8 7 : . .
```

# Problem statement by example

## Cursor-sort of rows

We next look at what we need to be concerned about when we write multiple fragments using the chunked-CSR reader.

Suppose the input `X` array is in CSR format as above:

```
    T V S U
  C 1 2 . .
  A . 3 4 .
  B . . 5 6
  D 7 . 8 .
```

And suppose we want to write it in two chunks of two rows each.

We must cursor-sort row labels so (with zero copy) the matrix will effectively look like this

```
    T V S U
  A . 3 4 .
  B . . 5 6
  ---------- chunk / fragment boundary
  C 1 2 . .
  D 7 . 8 .
```

This is necessary, since otherwise every fragment would have the same MBRs in both dimensions and all queries -- whether by `obs_id,var_id`, or by `obs_id`, or by `var_id` -- would need to consult all fragments.

* Chunk 1 (written as fragment 1) gets these COOs:
  * `A,V,3`
  * `A,S,4`
  * `B,S,5`
  * `B,U,6`
* Chunk 2 (written as fragment 2) gets these COOs:
  * `C,T,1`
  * `C,V,2`
  * `D,T,7`
  * `D,S,8`
* Fragment 1 MBR is `[A..B, S..V]`
* Fragment 2 MBR is `[C..D, S..V]`
* TileDB guarantees sorting on both dims within the fragment

Here's the performance concern:

* Queries on `obs_id,var_id` will locate only one fragment, since a given `obs_id` can only be in one fragment
* Queries on `obs_id` will locate one fragment, since a given `obs_id` can only be in one fragment
* Queries on `var_id` will locate _all_ fragments. (Note, however, this is the same amount of data as when the TileDB array was all in one fragment.)

## Cursor-sort of columns

Suppose we were to column-sort the CSR too -- it would look like this:

```
    S T U V
  A 4 . . 3
  B: 5 . 6 .
  ---------- chunk boundary
  C . 1 . 2
  D 8 7 . .
```

* Chunk 1 (written as fragment 1) gets these COOs:
  * `A,S,4`
  * `A,V,3`
  * `B,S,5`
  * `B,U,6`
* Chunk 2 (written as fragment 2) gets these COOs:
  * `C,T,1`
  * `C,V,2`
  * `D,S,8`
  * `D,T,7`
* Fragment 1 MBR is `[A..B, S..V]` same as before
* Fragment 2 MBR is `[C..D, S..V]` same as before
* TileDB guarantees sorting on both dims within the fragment

But the performance concern is _identical_ to the situation without cursor-sort of columns: in fact,
cursor-sorting the columns provides no benefit since TileDB is already sorting by both dimensions
within fragments, and the `var_id` slot of the fragment MBRs are `S..V` in both cases.

## Checkerboarding

Another option is to cursor-sort by both dimensions and then checkerboard:

```
    S T | U V
  A 4 . | . 3
  B: 5 .|  6 .
  ------+----- chunk boundary
  C . 1 | . 2
  D 8 7 | . .
```

* Fragment 1 gets these COOs:
  * `A,S,4`
  * `B,S,5`
* Fragment 2 gets these COOs:
  * `A,V,3`
  * `B,U,6`
* Fragment 3 gets these COOs:
  * `C,T,1`
  * `D,S,8`
  * `D,T,7`
* Fragment 4)gets these COOs:
  * `C,V,2`

* Fragment 1 MBR is `[A..B, S..T]`
* Fragment 2 MBR is `[A..B, U..V]`
* Fragment 3 MBR is `[C..D, S..T]`
* Fragment 4 MBR is `[C..D, U..V]`

* A query for `obs_id==D` will have to look at fragments 3 and 4
* A query for `var_id==T` will have to look at fragments 1 and 3
* We still cannot achieve having only one fragment for a given `obs_id`, and only one fragment for a
  given `var_id` -- we'd need
  to have a 'block diagional matrix' _even when the row & column labels are sorted_ which is not
  reasonable to expect.

## Global-order writes

See also [Python API docs](https://tiledb-inc-tiledb.readthedocs-hosted.com/en/1.6.3/tutorials/writing-sparse.html#writing-in-global-layout).

Idea:

* Write in global order (sorted by `obs_id` then `var_id`)
* Given the above example, we'd write
  * Fragment 1 gets these COOs:
    * `A,S,4`
    * `A,V,3`
    * `B,S,5`
    * `B,U,6`
  * Fragment 2 gets these COOs:
    * `C,V,2`
    * `C,T,1`
    * `D,S,8`
    * `D,T,7`
* Easy to do in Python at the row-chunk level
* Then:
  * Fragment writes will be faster.
  * Fragments will be auto-concatenated so they won't need consolidation at all.
  * Feature exists and is well-supported in C++.
  * Not yet present in the Python API.

# Suggested approach

* Use row-based chunking (checkerboard is not implemented as of 2022-05-03).
* Given that queries on `obs_id,var_id` or on `obs_id` will be efficient, but that queries on `var_id` will require consulting multiple fragments, ingest larger arrays as row-chunked CSR but consolidate them afterward.
* As of TileDB core 2.8.2, we cannot consolidate arrays with col-major tile order: so we write `X` with row-major tile order.
* Read-performance impact should be measured explicitly.
* Global-order writes need to be looked into.
