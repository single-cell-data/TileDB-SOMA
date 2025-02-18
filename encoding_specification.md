---
title: Encoding Specification
---

## SOMA Encoding Version 1.1.0

### SOMACollection

| Metadata key            | Datatype           | Value            | Description           | Required |
| :---------------------- | :----------------- | :--------------- | :-------------------- | :------- |
| "soma_datatype"         | TILEDB_STRING_UTF8 | "SOMACollection" | SOMA type name        | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0"          | SOMA encoding version | required |

### SOMADataFrame

| Metadata key            | Datatype           | Value           | Description           | Required |
| :---------------------- | :----------------- | :-------------- | :-------------------- | :------- |
| "soma_datatype"         | TILEDB_STRING_UTF8 | "SOMADataFrame" | SOMA type name        | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0"         | SOMA encoding version | required |

### SOMADenseNDArray

| Metadata key            | Datatype           | Value              | Description           | Required |
| :---------------------- | :----------------- | :----------------- | :-------------------- | :------- |
| "soma_datatype"         | TILEDB_STRING_UTF8 | "SOMADenseNDArray" | SOMA type name        | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0"            | SOMA encoding version | required |

### SOMAExperiment

| Metadata key            | Datatype           | Value            | Description           | Required |
| :---------------------- | :----------------- | :--------------- | :-------------------- | :------- |
| "soma_datatype"         | TILEDB_STRING_UTF8 | "SOMAExperiment" | SOMA type name        | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0"          | SOMA encoding version | required |

### SOMAMultiscaleImage

#### SOMAMultiscaleImage Spatial Encoding Version 0.2.0

| Metadata key                    | Datatype           | Value                                             | Description                                       | Required |
| :------------------------------ | :----------------- | :------------------------------------------------ | :------------------------------------------------ | :------- |
| "soma_datatype"                 | TILEDB_STRING_UTF8 | "SOMAMultiscaleImage"                             | SOMA type name                                    | required |
| "soma_encoding_version"         | TILEDB_STRING_UTF8 | "1.1.0"                                           | general SOMA encoding version                     | required |
| "soma_spatial_encoding_version" | TILEDB_STRING_UTF8 | "0.2.0"                                           | type-specific SOMA encoding version               | required |
| "soma_coordinate_space"         | TILEDB_STRING_UTF8 | JSON serialization of the coordinate space        | JSON serialization of the coordinate space        | required |
| "soma_multiscale_image_schema"  | TILEDB_STRING_UTF8 | JSON serialization fo the multiscale image schema | JSON serialization fo the multiscale image schema | required |

<!-- TODO: Add description of multiscale image schema metadata-->

#### SOMAMultiscaleImage Spatial Encoding Version 0.1.0

| Metadata key                    | Datatype           | Value                                             | Description                                       | Required |
| :------------------------------ | :----------------- | :------------------------------------------------ | :------------------------------------------------ | :------- |
| "soma_datatype"                 | TILEDB_STRING_UTF8 | "SOMAMultiscaleImage"                             | SOMA type name                                    | required |
| "soma_encoding_version"         | TILEDB_STRING_UTF8 | "1.1.0"                                           | general SOMA encoding version                     | required |
| "soma_spatial_encoding_version" | TILEDB_STRING_UTF8 | "0.1.0"                                           | type-specific SOMA encoding version               | required |
| "soma_coordinate_space"         | TILEDB_STRING_UTF8 | JSON serialization of the coordinate space        | JSON serialization of the coordinate space        | required |
| "soma_multiscale_image_schema"  | TILEDB_STRING_UTF8 | JSON serialization fo the multiscale image schema | JSON serialization fo the multiscale image schema | required |

<!-- TODO: Add description of multiscale image schema metadata -->

### SOMAMeasurement

| Metadata key            | Datatype           | Value             | Description           | Required |
| :---------------------- | :----------------- | :---------------- | :-------------------- | :------- |
| "soma_datatype"         | TILEDB_STRING_UTF8 | "SOMAMeasurement" | SOMA type name        | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0"           | SOMA encoding version | required |

### SOMAPointCloudDataFrame

#### SOMAPointCloudDataFram Spatial Encoding Version 0.2.0

| Metadata key                    | Datatype           | Value                                      | Description                                                   | Required |
| :------------------------------ | :----------------- | :----------------------------------------- | :------------------------------------------------------------ | :------- |
| "soma_datatype"                 | TILEDB_STRING_UTF8 | "SOMAPointCloud"                           | SOMA type name                                                | required |
| "soma_encoding_version"         | TILEDB_STRING_UTF8 | "1.1.0"                                    | general SOMA encoding version                                 | required |
| "soma_spatial_encoding_version" | TILEDB_STRING_UTF8 | "0.2.0"                                    | type-specific SOMA encoding version                           | required |
| "soma_coordinate_space"         | TILEDB_STRING_UTF8 | JSON serialization of the coordinate space | JSON serialization of the coordinate space                    | required |
| "soma_geometry_type"            | TILEDB_STRING_UTF8 | enum with possible values "radius"         | optional shape of the stored points                           | optional |
| "soma_geometry"                 | TILEDB_FLOAT64     | floating-point value                       | the radius of circles centered at the points in the dataframe | optional |

#### SOMAPointCloudDataFram Spatial Encoding Version 0.1.0

| Metadata key                    | Datatype           | Value                                      | Description                                                   | Required |
| :------------------------------ | :----------------- | :----------------------------------------- | :------------------------------------------------------------ | :------- |
| "soma_datatype"                 | TILEDB_STRING_UTF8 | "SOMAPointCloud"                           | SOMA type name                                                | required |
| "soma_encoding_version"         | TILEDB_STRING_UTF8 | "1.1.0"                                    | general SOMA encoding version                                 | required |
| "soma_spatial_encoding_version" | TILEDB_STRING_UTF8 | "0.1.0"                                    | type-specific SOMA encoding version                           | required |
| "soma_coordinate_space"         | TILEDB_STRING_UTF8 | JSON serialization of the coordinate space | JSON serialization of the coordinate space                    | required |
| "soma_geometry_type"            | TILEDB_STRING_UTF8 | enum with possible values "radius"         | optional shape of the stored points                           | optional |
| "soma_geometry"                 | TILEDB_FLOAT64     | floating-point value                       | the radius of circles centered at the points in the dataframe | optional |

### SOMAScene

#### SOMAScene Spatial Encoding Version 0.2.0

Scene Group Metadata

| Metadata key                    | Datatype           | Value                                      | Description                                | Required |
| :------------------------------ | :----------------- | :----------------------------------------- | :----------------------------------------- | :------- |
| "soma_datatype"                 | TILEDB_STRING_UTF8 | "SOMAMultiscaleImage"                      | SOMA type name                             | required |
| "soma_encoding_version"         | TILEDB_STRING_UTF8 | "1.1.0"                                    | general SOMA encoding version              | required |
| "soma_spatial_encoding_version" | TILEDB_STRING_UTF8 | "0.2.0"                                    | type-specific SOMA encoding version        | required |
| "soma_coordinate_space"         | TILEDB_STRING_UTF8 | JSON serialization of the coordinate space | JSON serialization of the coordinate space | optional |

Scene Transformation Metadata

<!-- TODO Define the format for the transformation data -->

#### SOMAScene Spatial Encoding Version 0.1.0

| Metadata key                    | Datatype           | Value                                      | Description                                | Required |
| :------------------------------ | :----------------- | :----------------------------------------- | :----------------------------------------- | :------- |
| "soma_datatype"                 | TILEDB_STRING_UTF8 | "SOMAMultiscaleImage"                      | SOMA type name                             | required |
| "soma_encoding_version"         | TILEDB_STRING_UTF8 | "1.1.0"                                    | general SOMA encoding version              | required |
| "soma_spatial_encoding_version" | TILEDB_STRING_UTF8 | "0.1.0"                                    | type-specific SOMA encoding version        | required |
| "soma_coordinate_space"         | TILEDB_STRING_UTF8 | JSON serialization of the coordinate space | JSON serialization of the coordinate space | optional |

### SOMASparseNDArray

| Metadata key            | Datatype           | Value               | Description           | Required |
| :---------------------- | :----------------- | :------------------ | :-------------------- | :------- |
| "soma_datatype"         | TILEDB_STRING_UTF8 | "SOMASparseNDArray" | SOMA type name        | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0"             | SOMA encoding version | required |

### JSON Serialized CoordinateSpace

<!-- TODO Add serialization of CoordinateSpace -->

## SOMA Encoding Version 1.0.0

### SOMACollection

| Metadata key            | Datatype           | Value            | Description           | Required |
| :---------------------- | :----------------- | :--------------- | :-------------------- | :------- |
| "soma_datatype"         | TILEDB_STRING_UTF8 | "SOMACollection" | SOMA type name        | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1"              | SOMA encoding version | required |

### SOMADataFrame

| Metadata key            | Datatype           | Value           | Description           | Required |
| :---------------------- | :----------------- | :-------------- | :-------------------- | :------- |
| "soma_datatype"         | TILEDB_STRING_UTF8 | "SOMADataFrame" | SOMA type name        | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1"             | SOMA encoding version | required |

### SOMADenseNDArray

| Metadata key            | Datatype           | Value              | Description           | Required |
| :---------------------- | :----------------- | :----------------- | :-------------------- | :------- |
| "soma_datatype"         | TILEDB_STRING_UTF8 | "SOMADenseNDArray" | SOMA type name        | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1"                | SOMA encoding version | required |

### SOMAExperiment

| Metadata key            | Datatype           | Value            | Description           | Required |
| :---------------------- | :----------------- | :--------------- | :-------------------- | :------- |
| "soma_datatype"         | TILEDB_STRING_UTF8 | "SOMAExperiment" | SOMA type name        | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1"              | SOMA encoding version | required |

### SOMAMeasurement

| Metadata key            | Datatype           | Value             | Description           | Required |
| :---------------------- | :----------------- | :---------------- | :-------------------- | :------- |
| "soma_datatype"         | TILEDB_STRING_UTF8 | "SOMAMeasurement" | SOMA type name        | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1"               | SOMA encoding version | required |

### SOMASparseNDArray

| Metadata key            | Datatype           | Value               | Description           | Required |
| :---------------------- | :----------------- | :------------------ | :-------------------- | :------- |
| "soma_datatype"         | TILEDB_STRING_UTF8 | "SOMASparseNDArray" | SOMA type name        | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1"                 | SOMA encoding version | required |
