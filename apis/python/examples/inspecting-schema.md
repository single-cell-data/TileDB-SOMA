## Inspecting single SOMA schemas

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

## Names of obs and var columns

```
print("OBS NAMES")
for soma in soco:
    print(soma.uri)
    for attr_name in soma.obs.keys():
        print("  obs", attr_name)

print("VAR NAMES")
for soma in soco:
    print(soma.uri)
    for attr_name in soma.var.keys():
        print("  var", attr_name)
```

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
