Here we show an example of doing a _slice query_ across a `SOMACollection`.

## Populate the collection

Here we use a few small sample files included in this repository.

```
import tiledbsc
import tiledbsc.io
import os
import shutil

soco_path = './soco-attribute-filter'
if os.path.exists(soco_path):
    shutil.rmtree(soco_path)

soco = tiledbsc.SOMACollection(soco_path)
if not soco.exists():
    soco._create()

for name, h5ad in [
    ('subset-soma-01', './anndata/subset-soma-01.h5ad'),
    ('subset-soma-02', './anndata/subset-soma-02.h5ad'),
    ('subset-soma-03', './anndata/subset-soma-03.h5ad'),
    ('subset-soma-04', './anndata/subset-soma-04.h5ad'),
]:
    soma_path = os.path.join(soco_path, name)
    soma = tiledbsc.SOMA(soma_path)
    tiledbsc.io.from_h5ad(soma, h5ad)
    soco.add(soma)
```

## Do the slice query

```
import tiledbsc
import os
import shutil
from typing import List, Dict

# ----------------------------------------------------------------
def soco_attribute_filter_prototype(
    soco: tiledbsc.SOMACollection,
    output_h5ad_path: str,
    output_soma_path: str,
    obs_attr_names: List[str] = [],
    obs_query_string: str = None,
    var_attr_names: List[str] = [],
    var_query_string: str = None,
) -> None:

    soma_slices = []
    for soma in soco:
        # E.g. querying for 'cell_type == "blood"' but this SOMA doesn'tiledbsc have a cell_type column in
        # its obs at all.
        if not soma.obs.has_attr_names(obs_attr_names):
            continue
        # E.g. querying for 'feature_name == "MT-CO3"' but this SOMA doesn'tiledbsc have a feature_name
        # column in its var at all.
        if not soma.var.has_attr_names(var_attr_names):
            continue

        soma_slice = soma.attribute_filter(obs_query_string, var_query_string)
        if soma_slice != None:
            print("Slice SOMA from", soma.name, soma.X.data.shape(), 'to', soma_slice.ann.X.shape)
            soma_slices.append(soma_slice)

    result_soma_slice = tiledbsc.SOMASlice.concat(soma_slices)
    if result_soma_slice is None:
        print("Empty slice")
    else:
        a = result_soma_slice.to_anndata()
        a.write_h5ad(output_h5ad_path)
        print("Wrote", output_h5ad_path, a.X.shape)

        if os.path.exists(output_soma_path):
            shutil.rmtree(output_soma_path)
        soma = tiledbsc.SOMA.from_soma_slice(result_soma_slice, output_soma_path, verbose=False)
        print("Wrote", output_soma_path, soma.X.data.shape())

# ----------------------------------------------------------------
soco_path = './soco-attribute-filter'
soco = tiledbsc.SOMACollection(soco_path)
soco_attribute_filter_prototype(
    soco,
    obs_attr_names=["tissue"],
    obs_query_string='tissue == "blood"',
    var_attr_names=["feature_name"],
    var_query_string='feature_name == "MT-CO3"',
)
```

## Examine the results

```
$ peek-soma slice-query-results
>>> soma.obs.df()
                 assay_ontology_term_id cell_type_ontology_term_id development_stage_ontology_term_id disease_ontology_term_id ethnicity_ontology_term_id  ... ethnicity      organism   sex tissue is_primary_data
obs_id                                                                                                                                                     ...
AAACCCAAGACGGTTG            EFO:0009922                 CL:0000814                     HsapDv:0000143            MONDO:0100096                    unknown  ...   unknown  Homo sapiens  male  blood               1
AAACCCACACAATGTC            EFO:0009922                 CL:0000814                     HsapDv:0000143            MONDO:0100096                    unknown  ...   unknown  Homo sapiens  male  blood               1
AAACCCACACCCAATA            EFO:0009922                 CL:0000236                     HsapDv:0000143            MONDO:0100096                    unknown  ...   unknown  Homo sapiens  male  blood               1
AAACCCACACGTACTA            EFO:0009922                 CL:0000763                     HsapDv:0000143            MONDO:0100096                    unknown  ...   unknown  Homo sapiens  male  blood               1
AAACCCACACTTCTCG            EFO:0009922                 CL:0000763                     HsapDv:0000143            MONDO:0100096                    unknown  ...   unknown  Homo sapiens  male  blood               1
...                                 ...                        ...                                ...                      ...                        ...  ...       ...           ...   ...    ...             ...
ACTATCTGTACTCGTA            EFO:0009922                 CL:0000233                     HsapDv:0000143            MONDO:0100096                    unknown  ...   unknown  Homo sapiens  male  blood               1
ACTATCTGTATTAAGG            EFO:0009922                 CL:0000814                     HsapDv:0000143            MONDO:0100096                    unknown  ...   unknown  Homo sapiens  male  blood               1
ACTATCTTCTAGACCA            EFO:0009922                 CL:0000763                     HsapDv:0000143            MONDO:0100096                    unknown  ...   unknown  Homo sapiens  male  blood               1
ACTATGGAGGATATGT            EFO:0009922                 CL:0000763                     HsapDv:0000143            MONDO:0100096                    unknown  ...   unknown  Homo sapiens  male  blood               1
ACTATGGGTGGAGAAA            EFO:0009922                 CL:0000763                     HsapDv:0000143            MONDO:0100096                    unknown  ...   unknown  Homo sapiens  male  blood               1

[400 rows x 17 columns]
>>> soma.var.df()
Empty DataFrame
Columns: []
Index: [ENSG00000198938]
>>> soma.X.data.df()
                                  value
obs_id           var_id
AAACCCAAGACGGTTG ENSG00000198938  190.0
AAACCCACACAATGTC ENSG00000198938  118.0
AAACCCACACCCAATA ENSG00000198938  151.0
AAACCCACACGTACTA ENSG00000198938   29.0
AAACCCACACTTCTCG ENSG00000198938  139.0
...                                 ...
ACTATCTGTACTCGTA ENSG00000198938   23.0
ACTATCTGTATTAAGG ENSG00000198938   48.0
ACTATCTTCTAGACCA ENSG00000198938  124.0
ACTATGGAGGATATGT ENSG00000198938   78.0
ACTATGGGTGGAGAAA ENSG00000198938   37.0

[399 rows x 1 columns]
```

## More queries

The `obs` or `var` sides can be omitted, to only slice on the other:

```
print("TWO-SIDED QUERY")
soco_attribute_filter_prototype(
    soco = tiledbsc.SOMACollection('mini-corpus/atlas'),
    output_h5ad_path = "mini-atlas-two-sided.h5ad",
    output_soma_path = "mini-atlas-two-sided",
    obs_attr_names=["cell_type"],
    obs_query_string='cell_type == "B cell"',
    var_attr_names=["feature_name"],
    var_query_string='feature_name == "MT-CO3"',
)

print()
print("OBS-ONLY QUERY")
soco_attribute_filter_prototype(
    soco = tiledbsc.SOMACollection('mini-corpus/atlas'),
    output_h5ad_path = "mini-atlas-obs-sided.h5ad",
    output_soma_path = "mini-atlas-obs-sided",
    obs_attr_names=["cell_type"],
    obs_query_string='cell_type == "B cell"',
)

print()
print("VAR-ONLY QUERY")
soco_attribute_filter_prototype(
    soco = tiledbsc.SOMACollection('mini-corpus/atlas'),
    output_h5ad_path = "mini-atlas-var-sided.h5ad",
    output_soma_path = "mini-atlas-var-sided",
    var_attr_names=["feature_name"],
    var_query_string='feature_name == "MT-CO3"',
)
```

```
TWO-SIDED QUERY
Slice SOMA from acute-covid19-cohort (59506, 24004) to (6131, 1)
Slice SOMA from 4056cbab-2a32-4c9e-a55f-c930bc793fb6 (4603, 33244) to (306, 1)
Slice SOMA from longitudinal-profiling-49 (11574, 33244) to (1453, 1)
Slice SOMA from 0cfab2d4-1b79-444e-8cbe-2ca9671ca85e (3726, 33244) to (124, 1)
Slice SOMA from autoimmunity-pbmcs (97499, 22224) to (510, 1)
Wrote mini-atlas-two-sided.h5ad (8524, 1)
Wrote mini-atlas-two-sided (8524, 1)

OBS-ONLY QUERY
Slice SOMA from acute-covid19-cohort (59506, 24004) to (6131, 24004)
Slice SOMA from 4056cbab-2a32-4c9e-a55f-c930bc793fb6 (4603, 33244) to (306, 33244)
Slice SOMA from longitudinal-profiling-49 (11574, 33244) to (1453, 33244)
Slice SOMA from 0cfab2d4-1b79-444e-8cbe-2ca9671ca85e (3726, 33244) to (124, 33244)
Slice SOMA from autoimmunity-pbmcs (97499, 22224) to (510, 22224)
Wrote mini-atlas-obs-sided.h5ad (8524, 21648)
Wrote mini-atlas-obs-sided (8524, 21648)

VAR-ONLY QUERY
Slice SOMA from acute-covid19-cohort (59506, 24004) to (59506, 1)
Slice SOMA from 4056cbab-2a32-4c9e-a55f-c930bc793fb6 (4603, 33244) to (4603, 1)
Slice SOMA from longitudinal-profiling-49 (11574, 33244) to (11574, 1)
Slice SOMA from 0cfab2d4-1b79-444e-8cbe-2ca9671ca85e (3726, 33244) to (3726, 1)
Slice SOMA from human-kidney-tumors-wilms (4636, 32922) to (4636, 1)
Slice SOMA from autoimmunity-pbmcs (97499, 22224) to (97499, 1)
Wrote mini-atlas-var-sided.h5ad (181544, 1)
Wrote mini-atlas-var-sided (181544, 1)
```
