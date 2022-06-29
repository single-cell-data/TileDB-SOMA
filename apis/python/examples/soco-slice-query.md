Here we show an example of doing a _slice query_ across a `SOMACollection` -- we extract a
relatively small subset out of the full collection for analysis.

A key point is that these data (shown here on local disk) can likewise be stored on object stores
like S3.  This means you can _query_ large datasets without having to first _download_ large
datasets.

:::{.callout-tip}
Best S3 read performance is obtained by querying from an in-region EC2 instance, or a TileDB Cloud
notebook -- this is preferred to laptop-to-S3 performance.
:::

Another key point is that the _out-of-core processing_ showing here allows you to slice data out of
a collection which is far larger than fits in RAM.

## Do the slice query

Using [soco-slice-query.py](soco-slice-query.py)

```
TWO-SIDED QUERY
Wrote mini-atlas-two-sided.h5ad (8524, 1)
Wrote mini-atlas-two-sided (8524, 1)

OBS-ONLY QUERY
Wrote mini-atlas-obs-sided.h5ad (8524, 21648)
Wrote mini-atlas-obs-sided (8524, 21648)

VAR-ONLY QUERY
Wrote mini-atlas-var-sided.h5ad (181544, 1)
Wrote mini-atlas-var-sided (181544, 1)

OBS-ONLY QUERY
Wrote cell-ontology-236.h5ad (8524, 21648)
Wrote cell-ontology-236 (8524, 21648)
```

## Examine the results

```
$ peek-soma mini-atlas-two-sided

johnkerl@Kerl-MBP[prod][python]$ peek-soma mini-atlas-obs-sided
>>> soma.obs.df()
                                 assay_ontology_term_id cell_type_ontology_term_id  ...      sex tissue
obs_id                                                                              ...
AAACCCACACCCAATA                            EFO:0009922                 CL:0000236  ...     male  blood
AAACCCAGTTCCACAA                            EFO:0009922                 CL:0000236  ...     male  blood
AAACCCATCCCTCATG                            EFO:0009922                 CL:0000236  ...     male  blood
AAACCCATCGAAGAAT                            EFO:0009922                 CL:0000236  ...     male  blood
AAACGAAAGAATTTGG                            EFO:0009922                 CL:0000236  ...     male  blood
...                                                 ...                        ...  ...      ...    ...
batch4_5p_rna|TTTGTCAAGACTGTAA-1            EFO:0011025                 CL:0000236  ...  unknown  blood
batch4_5p_rna|TTTGTCAAGGATGGTC-1            EFO:0011025                 CL:0000236  ...  unknown  blood
batch4_5p_rna|TTTGTCACATCGATTG-1            EFO:0011025                 CL:0000236  ...  unknown  blood
batch4_5p_rna|TTTGTCAGTATGAAAC-1            EFO:0011025                 CL:0000236  ...  unknown  blood
batch4_5p_rna|TTTGTCAGTCGCATAT-1            EFO:0011025                 CL:0000236  ...  unknown  blood

[8524 rows x 16 columns]

>>> soma.var.df()
Empty DataFrame
Columns: []
Index: [ENSG00000000003, ENSG00000000419, ENSG00000000457, ...]

[21648 rows x 0 columns]

```
