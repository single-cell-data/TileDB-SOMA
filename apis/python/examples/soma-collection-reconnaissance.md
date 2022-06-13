Next, let's do some cross-cutting queries over schemas of all SOMAs in the collection. The goal is
-- in preparation for a collection-level query -- to find out which `obs` columns, and which values
in those columns, are most likely to be promising in terms of yielding results given our
mini-corpus.

## Cell-counts

The mini-corpus we prepared is 29 SOMAs, 26GB total:

```
$ du -hs /mini-corpus/tiledb-data
 26G  /mini-corpus/tiledb-data

$ ls /mini-corpus/tiledb-data | wc -l
      29
```

This collection includes data on about 2.4 million cells:

```
import tiledbsc

soco = tiledbsc.SOMACollection('/mini-corpus/soco')

print("TOTAL CELL COUNT:")
print(soco.cell_count())
```

```
TOTAL CELL COUNT:
2464363
```

```
print()
print([soma.cell_count() for soma in soco])
```

```
[264824, 4636, 6288, 2223, 59506, 100, 2638, 982538, 385, 67794, 2638, 104148, 44721, 3799, 11574, 1679, 3589, 700, 584884, 16245, 4603, 3726, 4636, 7348, 3589, 40268, 12971, 4232, 80, 82478, 97499, 38024]
```

```
tabula-sapiens-stromal                                       82478
Puck_200903_10                                               38024
autoimmunity-pbmcs                                           97499
pbmc-small                                                   80
vieira19_Alveoli_and_parenchyma_anonymised.processed         12971
af9d8c03-696c-4997-bde8-8ef00844881b                         4232
d4db74ad-a129-4b1a-b9da-1b30db86bbe4-issue-74                3589
single-cell-transcriptomes                                   40268
local2                                                       7348
human-kidney-tumors-wilms                                    4636
0cfab2d4-1b79-444e-8cbe-2ca9671ca85e                         3726
issue-74                                                     3589
10x_pbmc68k_reduced                                          700
integrated-human-lung-cell-atlas                             584884
4056cbab-2a32-4c9e-a55f-c930bc793fb6                         4603
adult-mouse-cortical-cell-taxonomy                           1679
tabula-sapiens-epithelial                                    104148
Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection 44721
longitudinal-profiling-49                                    11574
azimuth-meta-analysis                                        982538
developmental-single-cell-atlas-of-the-murine-lung           67794
local3                                                       385
pbmc3k-krilow                                                2638
pbmc3k_processed                                             2638
subset_100_100                                               100
tabula-sapiens-immune                                        264824
brown-adipose-tissue-mouse                                   2223
acute-covid19-cohort                                         59506
issue-69                                                     6288
```

## Datasets having all three of obs.cell_type, obs.tissue, and obs.feature_name

```
names = sorted([
  soma.name for soma in soco
    if 'cell_type' in soma.obs.keys() and 'tissue' in soma.obs.keys() and 'feature_name' in soma.var.keys()
])
for name in names: print(name)
```

```
0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
4056cbab-2a32-4c9e-a55f-c930bc793fb6
Puck_200903_10
acute-covid19-cohort
adult-mouse-cortical-cell-taxonomy
af9d8c03-696c-4997-bde8-8ef00844881b
autoimmunity-pbmcs
azimuth-meta-analysis
brown-adipose-tissue-mouse
developmental-single-cell-atlas-of-the-murine-lung
human-kidney-tumors-wilms
integrated-human-lung-cell-atlas
local2
local3
longitudinal-profiling-49
single-cell-transcriptomes
tabula-sapiens-epithelial
tabula-sapiens-immune
tabula-sapiens-stromal
```

## Show counts of obs_ids and var_ids across the collection

```
def show_obs_id_counts(soco: t.SOMACollection) -> None:
    counts = {}
    for soma in soco:
        for oid in soma.obs.ids():
            if oid in counts:
                counts[oid] += 1
            else:
                counts[oid] = 1
    df = pandas.DataFrame.from_dict(
        {"obs_id": counts.keys(), "counts": counts.values()}
    )
    # print(df.head())
    print(df)


def show_var_id_counts(soco: t.SOMACollection) -> None:
    counts = {}
    for soma in soco:
        for oid in soma.var.ids():
            if oid in counts:
                counts[oid] += 1
            else:
                counts[oid] = 1
    df = pandas.DataFrame.from_dict(
        {"var_id": counts.keys(), "counts": counts.values()}
    )
    # print(df.head())
    print(df)

print("OBS_ID COUNTS")
show_obs_id_counts(soco)

print("VAR_ID COUNTS")
show_var_id_counts(soco)
```

```
OBS_ID COUNTS
                                                    obs_id  counts
0        AAACCCAAGAAACTGT_TSP7_LymphNodes_Inguinal_10X_1_1       1
1                   AAACCCAAGAAGGTAG_TSP10_Skin_NA_10X_1_1       1
2                 AAACCCAAGAAGGTAG_TSP14_SI_Distal_10X_1_1       1
3               AAACCCAAGAAGTCAT_TSP2_LymphNode_NA_10X_2_1       1
4              AAACCCAAGAATAACC_TSP14_Skin_Abdomen_10X_1_1       1
...                                                    ...     ...
2374963                                          Z11041427       1
2374964                                          Z11041428       1
2374965                                          Z11041429       1
2374966                                          Z11041430       1
2374967                                          Z11041431       1

[2374968 rows x 2 columns]

VAR_ID COUNTS
                    var_id  counts
0          ENSG00000000003      17
1          ENSG00000000005      15
2          ENSG00000000419      17
3          ENSG00000000457      17
4          ENSG00000000460      17
...                    ...     ...
172538  ENSMUSG00000119805       1
172539  ENSMUSG00000119882       1
172540  ENSMUSG00000119887       1
172541  ENSMUSG00000119908       1
172542  ENSMUSG00000119931       1

[172543 rows x 2 columns]
```

## Show unique values for obs.cell_type and var.feature_name

```
def show_obs_column_unique_values(soco: t.SOMACollection, col_name: str) -> None:
    for soma in soco:
        print()
        print(soma.uri)
        if col_name in soma.obs.keys():
            print(soma.obs.df()[col_name].unique())

def show_var_column_unique_values(soco: t.SOMACollection, col_name: str) -> None:
    for soma in soco:
        print()
        print(soma.uri)
        if col_name in soma.var.keys():
            print(soma.var.df()[col_name].unique())

print("UNIQUE VALUES FOR OBS.CELL_TYPE")
show_obs_column_unique_values(soco, "cell_type")

print("UNIQUE VALUES FOR VAR.FEATURE_NAME")
show_var_column_unique_values(soco, "feature_name")
```

```
OBS UNIQUE VALUES FOR CELL_TYPE

file:///Users/testuser/mini-corpus/tiledb-data/tabula-sapiens-immune
[b'B cell' b'macrophage' b'CD4-positive, alpha-beta T cell'
 b'naive thymus-derived CD4-positive, alpha-beta T cell'
 b'mature NK T cell' b'erythrocyte' b'plasma cell' b'DN3 thymocyte'
 b'effector CD8-positive, alpha-beta T cell' b'granulocyte'
 b'effector CD4-positive, alpha-beta T cell' b'neutrophil' b'monocyte'
 b'T cell' b'DN1 thymic pro-T cell' b'intermediate monocyte'
 b'CD8-positive, alpha-beta T cell' b'naive B cell' b'leukocyte'
 b'CD8-positive, alpha-beta memory T cell' b'CD4-positive helper T cell'
 b'mesenchymal stem cell' b'classical monocyte' b'memory B cell'
 b'regulatory T cell' b'basophil' b'innate lymphoid cell'
 b'CD4-positive, alpha-beta memory T cell' b'dendritic cell'
 b'CD8-positive, alpha-beta cytokine secreting effector T cell'
 b'microglial cell' b'CD141-positive myeloid dendritic cell'
 b'type I NK T cell' b'naive regulatory T cell' b'hematopoietic stem cell'
 b'T follicular helper cell' b'thymocyte' b'erythroid progenitor cell'
 b'myeloid cell' b'non-classical monocyte'
 b'naive thymus-derived CD8-positive, alpha-beta T cell' b'mast cell'
 b'liver dendritic cell' b'CD8-positive, alpha-beta cytotoxic T cell'
 b'common myeloid progenitor' b'CD1c-positive myeloid dendritic cell'
 b'Langerhans cell' b'platelet' b'plasmacytoid dendritic cell'
 b'DN4 thymocyte' b'plasmablast' b'immature natural killer cell'
 b'mature conventional dendritic cell' b'erythroid lineage cell'
 b'double-positive, alpha-beta thymocyte' b'myeloid dendritic cell']

file:///Users/testuser/mini-corpus/tiledb-data/integrated-human-lung-cell-atlas
[b'capillary endothelial cell' b'plasmacytoid dendritic cell'
 b'alveolar macrophage' b'natural killer cell' b'type II pneumocyte'
 b'non-classical monocyte' b'elicited macrophage' b'plasma cell'
 b'CD8-positive, alpha-beta T cell' b'vein endothelial cell'
 b'lung macrophage' b'CD4-positive, alpha-beta T cell'
 b'CD1c-positive myeloid dendritic cell'
 b'pulmonary artery endothelial cell'
 b'tracheobronchial smooth muscle cell' b'bronchial goblet cell'
 b'classical monocyte' b'fibroblast of lung' b'B cell' b'pericyte cell'
 b'conventional dendritic cell' b'mast cell'
 b'ciliated columnar cell of tracheobronchial tree' b'dendritic cell'
 b'epithelial cell of lung' b'type I pneumocyte'
 b'nasal mucosa goblet cell' b'club cell' b'respiratory basal cell'
 b'endothelial cell of lymphatic vessel' b'mesothelial cell'
 b'bronchus fibroblast of lung' b'T cell'
 b'multi-ciliated epithelial cell' b'fibroblast' b'myofibroblast cell'
 b'stromal cell' b'acinar cell' b'mucus secreting cell'
 b'tracheobronchial serous cell' b'ionocyte' b'smooth muscle cell'
 b'brush cell of trachebronchial tree' b'tracheobronchial goblet cell'
 b'lung neuroendocrine cell' b'serous secreting cell']

...

OBS UNIQUE VALUES FOR FEATURE_NAME

file:///Users/testuser/mini-corpus/tiledb-data/tabula-sapiens-immune
[b'TSPAN6' b'TNMD' b'DPM1' ... b'XXyac-YX60D10.3' b'CTD-2201E18.6'
 b'RP11-444B5.1']

file:///Users/testuser/mini-corpus/tiledb-data/integrated-human-lung-cell-atlas
[b'TSPAN6' b'TNMD' b'DPM1' ... b'MGC4859' b'RP11-107E5.4' b'RP11-299P2.2']

...
```

## Show value counts for obs.cell_type and obs.tissue

```
def show_obs_value_counts(soco: t.SOMACollection, obs_labels: List[str]) -> None:

    for obs_label in obs_labels:
        counts = {}

        for soma in soco:
            print("...", soma.name)
            # print("\n".join(sorted(soma.obs.attr_names())))
            obs = soma.obs.df()

            if not obs_label in obs:
                continue

            obs_label_values = sorted(list(set(obs[obs_label])))
            for obs_label_value in obs_label_values:
                if obs_label_value in counts:
                    counts[obs_label_value] += 1
                else:
                    counts[obs_label_value] = 1

        print(
            "----------------------------------------------------------------",
            obs_label,
        )
        for k, v in dict(sorted(counts.items(), key=lambda item: item[1])).items():
            print(k, v)

print("OBS VALUE COUNTS FOR CELL_TYPE AND TISSUE")
show_obs_value_counts(soco, ["cell_type", "tissue"])
```

```
OBS VALUE COUNTS FOR CELL_TYPE AND TISSUE
---------------------------------------------------------------- cell_type
b'CD141-positive myeloid dendritic cell' 1
b'CD4-positive helper T cell' 1
b'CD8-positive, alpha-beta cytokine secreting effector T cell' 1
b'CD8-positive, alpha-beta cytotoxic T cell' 1
b'DN1 thymic pro-T cell' 1
b'DN3 thymocyte' 1
b'DN4 thymocyte' 1
b'Langerhans cell' 1
b'T follicular helper cell' 1
b'common myeloid progenitor' 1
b'double-positive, alpha-beta thymocyte' 1
b'effector CD4-positive, alpha-beta T cell' 1
...
b'astrocyte' 4
b'CD8-positive, alpha-beta T cell' 5
b'T cell' 5
b'endothelial cell' 5
b'CD4-positive, alpha-beta T cell' 6
b'mature NK T cell' 6
b'native cell' 6
b'dendritic cell' 7
b'plasmacytoid dendritic cell' 7
b'natural killer cell' 7
b'myeloid cell' 8
b'neutrophil' 8
b'plasma cell' 8
b'platelet' 8
b'B cell' 12

---------------------------------------------------------------- tissue
b'bone marrow' 1
b'inguinal lymph node' 1
b'spleen' 1
b'conjunctiva' 1
b'lung parenchyma' 1
b'nose' 1
b'respiratory airway' 1
b'fovea centralis' 1
b'macula lutea proper' 1
b'peripheral region of retina' 1
...
b'large intestine' 3
b'liver' 3
b'mammary gland' 3
b'myometrium' 3
b'parotid gland' 3
b'posterior part of tongue' 3
b'prostate gland' 3
b'retinal neural layer' 3
b'sclera' 3
b'skin of abdomen' 3
b'skin of body' 3
b'skin of chest' 3
b'small intestine' 3
b'submandibular gland' 3
b'thymus' 3
b'tongue' 3
b'trachea' 3
b'uterus' 3
b'kidney' 4
b'lung' 5
b'blood' 9
```

## Find datasets having specified values

```
def show_somas_having(
    soco: tiledbsc.SOMACollection,
    obs_labels_to_values: Dict[str, List],
    var_labels_to_values: Dict[str, List],
) -> None:

    for soma in soco:
        print("...", soma.uri)

        obs = soma.obs.df()
        for obs_label in obs_labels_to_values:
            if not obs_label in obs:
                continue
            soma_obs_label_values = sorted(list(set(obs[obs_label])))
            soma_obs_label_values = [e.decode() for e in soma_obs_label_values]
            for sought_obs_label_value in obs_labels_to_values[obs_label]:
                if sought_obs_label_value in soma_obs_label_values:
                    print("  found obs", sought_obs_label_value)

        var = soma.var.df()
        for var_label in var_labels_to_values:
            if not var_label in var:
                continue
            soma_var_label_values = sorted(list(set(var[var_label])))
            soma_var_label_values = [e.decode() for e in soma_var_label_values]
            for sought_var_label_value in var_labels_to_values[var_label]:
                if sought_var_label_value in soma_var_label_values:
                    print("  found var", sought_var_label_value)

show_somas_having(
    soco,
    {"cell_type": ["B cell", "T cell"], "tissue": ["blood", "lung"]},
    {"feature_name": ["MT-CO3"]},
)
```

```
... file:///mini-corpus/tiledb-data/tabula-sapiens-stromal
  found obs lung
  found var MT-CO3
... file:///mini-corpus/tiledb-data/Puck_200903_10
  found var MT-CO3
... file:///mini-corpus/tiledb-data/autoimmunity-pbmcs
  found obs B cell
  found obs T cell
  found obs blood
  found var MT-CO3
... file:///mini-corpus/tiledb-data/pbmc-small
... file:///mini-corpus/tiledb-data/vieira19_Alveoli_and_parenchyma_anonymised.processed
... file:///mini-corpus/tiledb-data/af9d8c03-696c-4997-bde8-8ef00844881b
  found obs B cell
  found obs blood
  found var MT-CO3
... file:///mini-corpus/tiledb-data/d4db74ad-a129-4b1a-b9da-1b30db86bbe4-issue-74
... file:///mini-corpus/tiledb-data/single-cell-transcriptomes
  found obs B cell
  found var MT-CO3
... file:///mini-corpus/tiledb-data/local2
  found var MT-CO3
... file:///mini-corpus/tiledb-data/human-kidney-tumors-wilms
  found var MT-CO3
... file:///mini-corpus/tiledb-data/0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
  found obs B cell
  found obs blood
  found var MT-CO3
... file:///mini-corpus/tiledb-data/issue-74
... file:///mini-corpus/tiledb-data/10x_pbmc68k_reduced
... file:///mini-corpus/tiledb-data/integrated-human-lung-cell-atlas
  found obs B cell
  found obs T cell
  found var MT-CO3
... file:///mini-corpus/tiledb-data/4056cbab-2a32-4c9e-a55f-c930bc793fb6
  found obs B cell
  found obs blood
  found var MT-CO3
... file:///mini-corpus/tiledb-data/adult-mouse-cortical-cell-taxonomy
... file:///mini-corpus/tiledb-data/tabula-sapiens-epithelial
  found obs lung
  found var MT-CO3
... file:///mini-corpus/tiledb-data/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection
  found obs B cell
  found obs blood
... file:///mini-corpus/tiledb-data/longitudinal-profiling-49
  found obs B cell
  found obs blood
  found var MT-CO3
... file:///mini-corpus/tiledb-data/azimuth-meta-analysis
  found obs B cell
  found obs T cell
  found obs blood
  found obs lung
... file:///mini-corpus/tiledb-data/developmental-single-cell-atlas-of-the-murine-lung
  found obs lung
... file:///mini-corpus/tiledb-data/local3
  found var MT-CO3
... file:///mini-corpus/tiledb-data/pbmc3k-krilow
... file:///mini-corpus/tiledb-data/pbmc3k_processed
... file:///mini-corpus/tiledb-data/subset_100_100
... file:///mini-corpus/tiledb-data/tabula-sapiens-immune
  found obs B cell
  found obs T cell
  found obs blood
  found obs lung
  found var MT-CO3
... file:///mini-corpus/tiledb-data/brown-adipose-tissue-mouse
  found obs B cell
  found obs T cell
... file:///mini-corpus/tiledb-data/acute-covid19-cohort
  found obs B cell
  found obs blood
  found var MT-CO3
... file:///mini-corpus/tiledb-data/issue-69
>>>
```

## Conclusion

From these we conclude that `obs.cell_type == "B cell"` and `obs.tissue == "blood"`, and
`var.feature_name == "MT-CO3"` (acquired similarly but not shown here) are likeliest to produce the
largest result set, given our local-disk mini-corpus.
