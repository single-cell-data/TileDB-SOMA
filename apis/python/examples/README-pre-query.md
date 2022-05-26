# Example corpus setup and pre-query reconnaissance

## Overview

Output here is presented/formatted from [pre-query.py](pre-query.py). The idea is:

* Obtain a collection of single-cell data
* Ingest it into TileDB SOMA format
* Before running a query (a subject for a separate note), use the SOMA API to look over the collection and see which datasets have which columns in their `obs` and `var`, and for selected columns, see which values occur.

## Obtaining AnnData files

Visit [https://cellxgene.cziscience.com](https://cellxgene.cziscience.com) and select from among various choices there and download.

Files used for this example:

<details>

```
$ ls -lh /mini-corpus/anndata
-rw-r--r--@ 1 testuser  staff    51M May 13 22:34 0cfab2d4-1b79-444e-8cbe-2ca9671ca85e.h5ad
-rw-r--r--  1 testuser  staff   4.3M May 10 18:12 10x_pbmc68k_reduced.h5ad
-rw-r--r--@ 1 testuser  staff    55M May 13 22:33 4056cbab-2a32-4c9e-a55f-c930bc793fb6.h5ad
-rw-r--r--@ 1 testuser  staff    32M May 13 22:17 Puck_200903_10.h5ad
-rw-r--r--@ 1 testuser  staff   230M Apr 29 17:13 Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection.h5ad
-rw-r--r--@ 1 testuser  staff   376M Apr 25 22:48 acute-covid19-cohort.h5ad
-rw-r--r--@ 1 testuser  staff   117M May 13 22:20 adult-mouse-cortical-cell-taxonomy.h5ad
-rw-r--r--@ 1 testuser  staff    27M May 13 22:34 af9d8c03-696c-4997-bde8-8ef00844881b.h5ad
-rw-r--r--@ 1 testuser  staff   686M Apr 25 22:50 autoimmunity-pbmcs.h5ad
-rw-r--r--@ 1 testuser  staff   6.6G May 13 22:40 azimuth-meta-analysis.h5ad
-rw-r--r--@ 1 testuser  staff    40M May 13 22:21 brown-adipose-tissue-mouse.h5ad
-rw-r--r--@ 1 testuser  staff    30M May 13 22:45 d4db74ad-a129-4b1a-b9da-1b30db86bbe4-issue-74.h5ad
-rw-r--r--@ 1 testuser  staff   712M May 13 22:22 developmental-single-cell-atlas-of-the-murine-lung.h5ad
-rw-r--r--@ 1 testuser  staff    70M Apr 25 22:45 human-kidney-tumors-wilms.h5ad
-rw-r--r--@ 1 testuser  staff   5.6G Apr 25 23:04 integrated-human-lung-cell-atlas.h5ad
-rw-r--r--@ 1 testuser  staff   357M May 10 20:08 issue-69.h5ad
-rw-r--r--@ 1 testuser  staff    99M May 11 09:57 issue-71.h5ad
-rw-r--r--@ 1 testuser  staff    30M May 13 23:25 issue-74.h5ad
-rw-r--r--@ 1 testuser  staff   329M Apr 25 11:13 local2.h5ad
-rw-r--r--@ 1 testuser  staff    36M Apr 25 11:13 local3.h5ad
-rw-r--r--@ 1 testuser  staff   119M Apr 25 22:47 longitudinal-profiling-49.h5ad
-rw-r--r--  1 testuser  staff   230K May 11 08:08 pbmc-small.h5ad
-rw-r--r--  1 testuser  staff    47M Apr 28 15:15 pbmc3k-krilow.h5ad
-rw-r--r--  1 testuser  staff    38M May 11 08:08 pbmc3k_processed.h5ad
-rw-r--r--@ 1 testuser  staff   221M Apr 25 22:46 single-cell-transcriptomes.h5ad
-rw-r--r--  1 testuser  staff    34K Apr 25 11:13 subset_100_100.h5ad
-rw-r--r--@ 1 testuser  staff   3.2G May 13 22:30 tabula-sapiens-epithelial.h5ad
-rw-r--r--@ 1 testuser  staff   5.7G May 13 22:38 tabula-sapiens-immune.h5ad
-rw-r--r--@ 1 testuser  staff   2.5G May 13 22:35 tabula-sapiens-stromal.h5ad
-rw-r--r--@ 1 testuser  staff   231M Apr 25 11:13 vieira19_Alveoli_and_parenchyma_anonymised.processed.h5ad
```

</details>

## Ingesting into SOMAs

Use the [ingestor script](../tools/ingestor) to ingest these into SOMAs:

```
tools/ingestor -o /mini-corpus/tiledb-data -n /mini-corpus/anndata/0cfab2d4-1b79-444e-8cbe-2ca9671ca85e.h5ad
tools/ingestor -o /mini-corpus/tiledb-data -n /mini-corpus/anndata/10x_pbmc68k_reduced.h5ad
...
```

Note this takes many hours. The benefit of using an optimized storage solution (with admittedly
non-negligible ingest time) is that all subsequent queries benefit from that optimized storage. In
particular, all data queries from here on down in this note completed in a minute total.

## Populate a SOMA collection

Use the [populator script](../tools/populate-soco) to mark these as members of a SOMA collection:

```
populate-soco -o /mini-corpus/soco -a /mini-corpus/tiledb-data/0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
populate-soco -o /mini-corpus/soco -a /mini-corpus/tiledb-data/10x_pbmc68k_reduced
...
```

or simply

```
populate-soco -o /mini-corpus/soco -a /mini-corpus/tiledb-data/*
```

Note this is quite quick.

## Introspect the data: overview

The rest of this note interleaves code and output from the following:

```
$ examples/pre-query.py /mini-corpus/soco
```

## Looking at single SOMA schemas

Before we do full-collection traversals, let's first take a look at a single SOMA.

```
import tiledbsc
soma = tiledbsc.SOMA('/mini-corpus/tiledb-data/tabula-sapiens-epithelial')
```

```
>>> soma.obs.shape()
(104148, 26)

>>> soma.var.shape()
(58559, 12)
```

```
>>> soma.obs.df()
                                                   tissue_in_publication  ...        development_stage
obs_id                                                                    ...
AAACCCAAGAACTCCT_TSP14_Lung_Distal_10X_1_1                          Lung  ...  59-year-old human stage
AAACCCAAGAGGGTAA_TSP8_Prostate_NA_10X_1_1                       Prostate  ...  56-year-old human stage
AAACCCAAGCCACTCG_TSP14_Prostate_NA_10X_1_2                      Prostate  ...  59-year-old human stage
AAACCCAAGCCGGAAT_TSP14_Liver_NA_10X_1_1                            Liver  ...  59-year-old human stage
AAACCCAAGCCTTGAT_TSP7_Tongue_Posterior_10X_1_1                    Tongue  ...  69-year-old human stage
...                                                                  ...  ...                      ...
TTTGTTGTCTACGGTA_TSP5_Eye_NA_10X_1_2                                 Eye  ...  40-year-old human stage
TTTGTTGTCTATCGGA_TSP2_Lung_proxmedialdistal_10X...                  Lung  ...  61-year-old human stage
TTTGTTGTCTCTCAAT_TSP2_Kidney_NA_10X_1_2                           Kidney  ...  61-year-old human stage
TTTGTTGTCTGCCTGT_TSP4_Mammary_NA_10X_1_2                         Mammary  ...  38-year-old human stage
TTTGTTGTCTGTAACG_TSP14_Prostate_NA_10X_1_1                      Prostate  ...  59-year-old human stage

[104148 rows x 26 columns]

>>> soma.var.df()
                    feature_type           ensemblid  ...        feature_name  feature_reference
var_id                                                ...
ENSG00000000003  Gene Expression  ENSG00000000003.14  ...           b'TSPAN6'     NCBITaxon:9606
ENSG00000000005  Gene Expression   ENSG00000000005.6  ...             b'TNMD'     NCBITaxon:9606
ENSG00000000419  Gene Expression  ENSG00000000419.12  ...             b'DPM1'     NCBITaxon:9606
ENSG00000000457  Gene Expression  ENSG00000000457.14  ...            b'SCYL3'     NCBITaxon:9606
ENSG00000000460  Gene Expression  ENSG00000000460.17  ...         b'C1orf112'     NCBITaxon:9606
...                          ...                 ...  ...                 ...                ...
ENSG00000286268  Gene Expression   ENSG00000286268.1  ...  b'LL0XNC01-30I4.1'     NCBITaxon:9606
ENSG00000286269  Gene Expression   ENSG00000286269.1  ...     b'RP11-510D4.1'     NCBITaxon:9606
ENSG00000286270  Gene Expression   ENSG00000286270.1  ...  b'XXyac-YX60D10.3'     NCBITaxon:9606
ENSG00000286271  Gene Expression   ENSG00000286271.1  ...    b'CTD-2201E18.6'     NCBITaxon:9606
ENSG00000286272  Gene Expression   ENSG00000286272.1  ...     b'RP11-444B5.1'     NCBITaxon:9606

[58559 rows x 12 columns]
```

```
>>> soma.obs.keys()
['tissue_in_publication', 'assay_ontology_term_id', 'donor', 'anatomical_information',
'n_counts_UMIs', 'n_genes', 'cell_ontology_class', 'free_annotation', 'manually_annotated',
'compartment', 'sex_ontology_term_id', 'is_primary_data', 'organism_ontology_term_id',
'disease_ontology_term_id', 'ethnicity_ontology_term_id', 'development_stage_ontology_term_id',
'cell_type_ontology_term_id', 'tissue_ontology_term_id', 'cell_type', 'assay', 'disease',
'organism', 'sex', 'tissue', 'ethnicity', 'development_stage']

>>> soma.var.attr_names_to_types()
{'feature_type': dtype('<U'), 'ensemblid': dtype('<U'), 'highly_variable': dtype('uint8'), 'means':
dtype('float64'), 'dispersions': dtype('float64'), 'dispersions_norm': dtype('float32'), 'mean':
dtype('float64'), 'std': dtype('float64'), 'feature_biotype': dtype('<U'), 'feature_is_filtered':
dtype('uint8'), 'feature_name': dtype('S'), 'feature_reference': dtype('<U')}
>>>
```

## Names and URIs

Next let's start taking a look across the collection.

```
def print_names_and_uris(soco: t.SOMACollection) -> None:
    for soma in soco:
        print("%-40s %s" % (soma.name, soma.uri))

print("NAMES AND URIS")
print_names_and_uris(soco)
```

```
NAMES AND URIS
tabula-sapiens-immune                    file:///mini-corpus/tiledb-data/tabula-sapiens-immune
tabula-sapiens-epithelial                file:///mini-corpus/tiledb-data/tabula-sapiens-epithelial
integrated-human-lung-cell-atlas         file:///mini-corpus/tiledb-data/integrated-human-lung-cell-atlas
af9d8c03-696c-4997-bde8-8ef00844881b     file:///mini-corpus/tiledb-data/af9d8c03-696c-4997-bde8-8ef00844881b
subset_100_100                           file:///mini-corpus/tiledb-data/subset_100_100
pbmc3k_processed                         file:///mini-corpus/tiledb-data/pbmc3k_processed
pbmc3k-krilow                            file:///mini-corpus/tiledb-data/pbmc3k-krilow
issue-74                                 file:///mini-corpus/tiledb-data/issue-74
local2                                   file:///mini-corpus/tiledb-data/local2
developmental-single-cell-atlas-of-the-murine-lung file:///mini-corpus/tiledb-data/developmental-single-cell-atlas-of-the-murine-lung
single-cell-transcriptomes               file:///mini-corpus/tiledb-data/single-cell-transcriptomes
tabula-sapiens-stromal                   file:///mini-corpus/tiledb-data/tabula-sapiens-stromal
azimuth-meta-analysis                    file:///mini-corpus/tiledb-data/azimuth-meta-analysis
vieira19_Alveoli_and_parenchyma_anonymised.processed file:///mini-corpus/tiledb-data/vieira19_Alveoli_and_parenchyma_anonymised.processed
local3                                   file:///mini-corpus/tiledb-data/local3
Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection file:///mini-corpus/tiledb-data/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection
4056cbab-2a32-4c9e-a55f-c930bc793fb6     file:///mini-corpus/tiledb-data/4056cbab-2a32-4c9e-a55f-c930bc793fb6
longitudinal-profiling-49                file:///mini-corpus/tiledb-data/longitudinal-profiling-49
human-kidney-tumors-wilms                file:///mini-corpus/tiledb-data/human-kidney-tumors-wilms
issue-69                                 file:///mini-corpus/tiledb-data/issue-69
autoimmunity-pbmcs                       file:///mini-corpus/tiledb-data/autoimmunity-pbmcs
brown-adipose-tissue-mouse               file:///mini-corpus/tiledb-data/brown-adipose-tissue-mouse
d4db74ad-a129-4b1a-b9da-1b30db86bbe4-issue-74 file:///mini-corpus/tiledb-data/d4db74ad-a129-4b1a-b9da-1b30db86bbe4-issue-74
pbmc-small                               file:///mini-corpus/tiledb-data/pbmc-small
10x_pbmc68k_reduced                      file:///mini-corpus/tiledb-data/10x_pbmc68k_reduced
Puck_200903_10                           file:///mini-corpus/tiledb-data/Puck_200903_10
0cfab2d4-1b79-444e-8cbe-2ca9671ca85e     file:///mini-corpus/tiledb-data/0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
acute-covid19-cohort                     file:///mini-corpus/tiledb-data/acute-covid19-cohort
adult-mouse-cortical-cell-taxonomy       file:///mini-corpus/tiledb-data/adult-mouse-cortical-cell-taxonomy
```

## Names of obs and var columns

```
def show_obs_names(soco: t.SOMACollection) -> None:
    for soma in soco:
        print(soma.uri)
        for attr_name in soma.obs.keys():
            print("  obs", attr_name)

def show_var_names(soco: t.SOMACollection) -> None:
    for soma in soco:
        print(soma.uri)
        for attr_name in soma.var.keys():
            print("  var", attr_name)

print("OBS NAMES")
show_obs_names(soco)

print("VAR NAMES")
show_var_names(soco)
```

<details>

```
OBS NAMES
file:///mini-corpus/tiledb-data/tabula-sapiens-immune
  obs tissue_in_publication
  obs assay_ontology_term_id
  obs donor
  obs anatomical_information
  obs n_counts_UMIs
  obs n_genes
  obs cell_ontology_class
  obs free_annotation
  obs manually_annotated
  obs compartment
  obs sex_ontology_term_id
  obs is_primary_data
  obs organism_ontology_term_id
  obs disease_ontology_term_id
  obs ethnicity_ontology_term_id
  obs development_stage_ontology_term_id
  obs cell_type_ontology_term_id
  obs tissue_ontology_term_id
  obs cell_type
  obs assay
  obs disease
  obs organism
  obs sex
  obs tissue
  obs ethnicity
  obs development_stage
file:///mini-corpus/tiledb-data/integrated-human-lung-cell-atlas
  obs is_primary_data
  obs assay_ontology_term_id
  obs cell_type_ontology_term_id
  obs development_stage_ontology_term_id
  obs disease_ontology_term_id
  obs ethnicity_ontology_term_id
  obs tissue_ontology_term_id
  obs organism_ontology_term_id
  obs sex_ontology_term_id
  obs sample
  obs study
  obs subject_ID
  obs smoking_status
  obs BMI
  obs condition
  obs subject_type
  obs sample_type
  obs 3'_or_5'
  obs sequencing_platform
  obs cell_ranger_version
  obs fresh_or_frozen
  obs dataset
  obs anatomical_region_level_2
  obs anatomical_region_level_3
  obs anatomical_region_highest_res
  obs age
  obs ann_highest_res
  obs n_genes
  obs size_factors
  obs log10_total_counts
  obs mito_frac
  obs ribo_frac
  obs original_ann_level_1
  obs original_ann_level_2
  obs original_ann_level_3
  obs original_ann_level_4
  obs original_ann_level_5
  obs original_ann_nonharmonized
  obs scanvi_label
  obs leiden_1
  obs leiden_2
  obs leiden_3
  obs anatomical_region_ccf_score
  obs entropy_study_leiden_3
  obs entropy_dataset_leiden_3
  obs entropy_subject_ID_leiden_3
  obs entropy_original_ann_level_1_leiden_3
  obs entropy_original_ann_level_2_clean_leiden_3
  obs entropy_original_ann_level_3_clean_leiden_3
  obs entropy_original_ann_level_4_clean_leiden_3
  obs entropy_original_ann_level_5_clean_leiden_3
  obs leiden_4
  obs reannotation_type
  obs leiden_5
  obs ann_finest_level
  obs ann_level_1
  obs ann_level_2
  obs ann_level_3
  obs ann_level_4
  obs ann_level_5
  obs ann_coarse_for_GWAS_and_modeling
  obs cell_type
  obs assay
  obs disease
  obs organism
  obs sex
  obs tissue
  obs ethnicity
  obs development_stage
...

VAR NAMES
file:///mini-corpus/tiledb-data/tabula-sapiens-immune
  var feature_type
  var ensemblid
  var highly_variable
  var means
  var dispersions
  var dispersions_norm
  var mean
  var std
  var feature_biotype
  var feature_is_filtered
  var feature_name
  var feature_reference
file:///mini-corpus/tiledb-data/integrated-human-lung-cell-atlas
  var n_cells
  var highly_variable
  var means
  var dispersions
  var feature_biotype
  var feature_is_filtered
  var feature_name
  var feature_reference
...
```

</details>

## Datasets having all three of obs.cell_type, obs.tissue, and obs.feature_name

```
def show_somas_with_all_three(soco: t.SOMACollection) -> None:
    for soma in soco:
        if "cell_type" in soma.obs.attr_names():
            if "tissue" in soma.obs.attr_names():
                if "feature_name" in soma.var.attr_names():
                    print(soma.uri)

print("SOMAS HAVING ALL THREE")
show_somas_with_all_three(soco)
```

```
SOMAS HAVING ALL THREE
file:///mini-corpus/tiledb-data/tabula-sapiens-immune
file:///mini-corpus/tiledb-data/tabula-sapiens-epithelial
file:///mini-corpus/tiledb-data/integrated-human-lung-cell-atlas
file:///mini-corpus/tiledb-data/af9d8c03-696c-4997-bde8-8ef00844881b
file:///mini-corpus/tiledb-data/local2
file:///mini-corpus/tiledb-data/developmental-single-cell-atlas-of-the-murine-lung
file:///mini-corpus/tiledb-data/single-cell-transcriptomes
file:///mini-corpus/tiledb-data/tabula-sapiens-stromal
file:///mini-corpus/tiledb-data/azimuth-meta-analysis
file:///mini-corpus/tiledb-data/local3
file:///mini-corpus/tiledb-data/4056cbab-2a32-4c9e-a55f-c930bc793fb6
file:///mini-corpus/tiledb-data/longitudinal-profiling-49
file:///mini-corpus/tiledb-data/human-kidney-tumors-wilms
file:///mini-corpus/tiledb-data/autoimmunity-pbmcs
file:///mini-corpus/tiledb-data/brown-adipose-tissue-mouse
file:///mini-corpus/tiledb-data/Puck_200903_10
file:///mini-corpus/tiledb-data/0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
file:///mini-corpus/tiledb-data/acute-covid19-cohort
file:///mini-corpus/tiledb-data/adult-mouse-cortical-cell-taxonomy
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
    soco: t.SOMACollection,
    obs_labels_to_values: Dict[str, List],
    var_labels_to_values: Dict[str, List],
) -> None:

    for soma in soco:
        print(soma.uri)

        obs = soma.obs.df()
        for obs_label in obs_labels_to_values:
            if not obs_label in obs:
                continue
            soma_obs_label_values = sorted(list(set(obs[obs_label])))
            for sought_obs_label_value in obs_labels_to_values[obs_label]:
                if sought_obs_label_value in soma_obs_label_values:
                    print("  found obs", sought_obs_label_value)

        var = soma.var.df()
        for var_label in var_labels_to_values:
            if not var_label in var:
                continue
            soma_var_label_values = sorted(list(set(var[var_label])))
            for sought_var_label_value in var_labels_to_values[var_label]:
                if sought_var_label_value in soma_var_label_values:
                    print("  found var", sought_var_label_value)

print("SHOW SOMAS HAVING")
show_somas_having(
    soco,
    {"cell_type": ["B cell", "T cell"], "tissue": ["blood", "lung"]},
    {"feature_name": ["MT-CO3"]},
)
```

XXX FIX ME

```
SHOW SOMAS HAVING
file:///Users/testuser/mini-corpus/tiledb-data/tabula-sapiens-stromal
file:///Users/testuser/mini-corpus/tiledb-data/Puck_200903_10
file:///Users/testuser/mini-corpus/tiledb-data/autoimmunity-pbmcs
file:///Users/testuser/mini-corpus/tiledb-data/pbmc-small
file:///Users/testuser/mini-corpus/tiledb-data/vieira19_Alveoli_and_parenchyma_anonymised.processed
file:///Users/testuser/mini-corpus/tiledb-data/af9d8c03-696c-4997-bde8-8ef00844881b
file:///Users/testuser/mini-corpus/tiledb-data/d4db74ad-a129-4b1a-b9da-1b30db86bbe4-issue-74
file:///Users/testuser/mini-corpus/tiledb-data/single-cell-transcriptomes
file:///Users/testuser/mini-corpus/tiledb-data/local2
file:///Users/testuser/mini-corpus/tiledb-data/human-kidney-tumors-wilms
file:///Users/testuser/mini-corpus/tiledb-data/0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
file:///Users/testuser/mini-corpus/tiledb-data/issue-74
file:///Users/testuser/mini-corpus/tiledb-data/10x_pbmc68k_reduced
file:///Users/testuser/mini-corpus/tiledb-data/integrated-human-lung-cell-atlas
file:///Users/testuser/mini-corpus/tiledb-data/4056cbab-2a32-4c9e-a55f-c930bc793fb6
file:///Users/testuser/mini-corpus/tiledb-data/adult-mouse-cortical-cell-taxonomy
file:///Users/testuser/mini-corpus/tiledb-data/tabula-sapiens-epithelial
file:///Users/testuser/mini-corpus/tiledb-data/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection
file:///Users/testuser/mini-corpus/tiledb-data/longitudinal-profiling-49
file:///Users/testuser/mini-corpus/tiledb-data/azimuth-meta-analysis
file:///Users/testuser/mini-corpus/tiledb-data/developmental-single-cell-atlas-of-the-murine-lung
file:///Users/testuser/mini-corpus/tiledb-data/local3
file:///Users/testuser/mini-corpus/tiledb-data/pbmc3k-krilow
file:///Users/testuser/mini-corpus/tiledb-data/pbmc3k_processed
file:///Users/testuser/mini-corpus/tiledb-data/subset_100_100
file:///Users/testuser/mini-corpus/tiledb-data/tabula-sapiens-immune
file:///Users/testuser/mini-corpus/tiledb-data/brown-adipose-tissue-mouse
file:///Users/testuser/mini-corpus/tiledb-data/acute-covid19-cohort
file:///Users/testuser/mini-corpus/tiledb-data/issue-69
```

## Conclusion

From these we conclude that `obs.cell_type == "B cell"` and `obs.tissue == "blood"`, and
`var.feature_name == "MT-CO3"` (acquired similarly but not shown here) are likeliest to produce the
largest result set, given our local-disk mini-corpus.
