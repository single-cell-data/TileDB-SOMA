Next we show an example of doing a _batch query_ across a `SOMACollection`.  Here we reduce
statistics from across the entire SOMA collection: specifically, we batch-process the mean of
`X/data`, grouping by `obs['cell_type_ontology_term_id']`.

A key point is that the _out-of-core processing_ showing here allows you to do multi-pass queries on
data from a collection which is far larger than fits in RAM.

## Do the batch query

Using [soco-batch-query.py](soco-batch-query.py) (this takes a few minutes on a 2.2GB atlas and
is amenable to parallelization) looping over the distinct values of `cell_type_ontology_term_id` in the collection:

```
...   0.00% done CL:1000449
...   2.63% done CL:0000037
...   5.26% done CL:0000232
...
...  92.11% done CL:0000789
...  94.74% done CL:0000940
...  97.37% done CL:0000897
```

We have averaged over cells, so the result is a matrix with `var_id` along rows and
`cell_type_ontology_term_id` along columns.


```
----------------------------------------------------------------
                                               value
var_id          cell_type_ontology_term_id
ENSG00000000003 CL:0000784                  0.000000
ENSG00000000419 CL:0000784                  0.060000
ENSG00000000457 CL:0000784                  0.033846
ENSG00000000460 CL:0000784                  0.004615
ENSG00000000938 CL:0000784                  0.338462
...                                              ...
ENSG00000283096 CL:0000003                  0.000000
ENSG00000283103 CL:0000003                  0.034409
ENSG00000283117 CL:0000003                  0.000000
ENSG00000283118 CL:0000003                  0.000000
ENSG00000283125 CL:0000003                  0.000000

[926384 rows x 1 columns]

----------------------------------------------------------------

As dense with zero-fill:
cell_type_ontology_term_id  CL:0000003  CL:0000037  CL:0000084  ...  CL:0001054  CL:0002396  CL:1000449
var_id                                                          ...
ENSG00000000003               0.079570    0.014815    0.001179  ...    0.000211    0.000000    3.615132
ENSG00000000005               0.178495    0.000000    0.000000  ...    0.000000    0.000000    0.030428
ENSG00000000419               0.070968    0.414815    0.705362  ...    0.112399    0.180967    1.176809
ENSG00000000457               0.032258    0.122222    0.131998  ...    0.023342    0.045762    0.387336
ENSG00000000460               0.015054    0.051852    0.283441  ...    0.005117    0.017681    0.319079
...                                ...         ...         ...  ...         ...         ...         ...
ENSSASG00005000009            0.000000    0.000000    0.000000  ...    0.000000    0.000000    0.000000
ENSSASG00005000010            0.000000    0.000000    0.000000  ...    0.000000    0.000000    0.000000
ENSSASG00005000011            0.000000    0.000000    0.000000  ...    0.000000    0.000000    0.000000
ENSSASG00005000012            0.000000    0.000000    0.000000  ...    0.000000    0.000000    0.000000
ENSSASG00005000013            0.000000    0.000000    0.000000  ...    0.000000    0.000000    0.000000

[33462 rows x 38 columns]

----------------------------------------------------------------
Trimmed of all-zeroes rows:
cell_type_ontology_term_id  CL:0000003  CL:0000037  CL:0000084  ...  CL:0001054  CL:0002396  CL:1000449
var_id                                                          ...
ENSG00000000003               0.079570    0.014815    0.001179  ...    0.000211    0.000000    3.615132
ENSG00000000005               0.178495    0.000000    0.000000  ...    0.000000    0.000000    0.030428
ENSG00000000419               0.070968    0.414815    0.705362  ...    0.112399    0.180967    1.176809
ENSG00000000457               0.032258    0.122222    0.131998  ...    0.023342    0.045762    0.387336
ENSG00000000460               0.015054    0.051852    0.283441  ...    0.005117    0.017681    0.319079
...                                ...         ...         ...  ...         ...         ...         ...
ENSG00000285480               0.000000    0.000000    0.000000  ...    0.000000    0.000000    0.000000
ENSG00000285486               0.000000    0.000000    0.000000  ...    0.000000    0.000000    0.000000
ENSG00000285492               0.000000    0.000000    0.002357  ...    0.009176    0.026521    0.000000
ENSG00000285505               0.000000    0.000000    0.000000  ...    0.000000    0.000000    0.000000
ENSG00000285509               0.000000    0.000000    0.000000  ...    0.000127    0.000000    0.000000

[29296 rows x 38 columns]
```
