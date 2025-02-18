---
title: Encoding Specification
---

## SOMA Encoding Version 2.0.0

SOMA encoding version 2.0.0 introduces per type versioning. T

### `SOMACollection`

#### `SOMACollection` Type Encoding Version 1.0

| Metadata key               | Datatype           | Value            | Description                         | Required |
| :------------------------- | :----------------- | :--------------- | :---------------------------------- | :------- |
| soma_datatype              | TILEDB_STRING_UTF8 | "SOMACollection" | SOMA type name                      | required |
| soma_encoding_version      | TILEDB_STRING_UTF8 | "2.0.0"          | general SOMA encoding version       | required |
| soma_type_encoding_version | TILEDB_STRING_UTF8 | "1.0"            | type-specific SOMA encoding version | required |

### `SOMAMultiscaleImage`

#### `SOMAMultiscaleImage` Type Encoding Version 1.0

| Metadata key                   | Datatype           | Value                                             | Description                                       | Required |
| :----------------------------- | :----------------- | :------------------------------------------------ | :------------------------------------------------ | :------- |
| "soma_datatype"                | TILEDB_STRING_UTF8 | "SOMAMultiscaleImage"                             | SOMA type name                                    | required |
| "soma_encoding_version"        | TILEDB_STRING_UTF8 | "2.0.0"                                           | general SOMA encoding version                     | required |
| "soma_type_encoding_version"   | TILEDB_STRING_UTF8 | "0.1"                                             | type-specific SOMA encoding version               | required |
| "soma_coordinate_space"        | TILEDB_STRING_UTF8 | JSON serialization of the coordinate space        | JSON serialization of the coordinate space        | required |
| "soma_multiscale_image_schema" | TILEDB_STRING_UTF8 | JSON serialization fo the multiscale image schema | JSON serialization fo the multiscale image schema | required |

Coordinate Space JSON

<!-- TODO: Add description of coordinate space metadata -->

Multiscale Image Schema JSON

<!-- TODO: Add description of multiscale image metadata -->

## SOMA Encoding Version 1.1.0

### SOMACollection

| Metadata key            | Datatype           | Value            | Description           | Required |
| :---------------------- | :----------------- | :--------------- | :-------------------- | :------- |
| "soma_datatype"         | TILEDB_STRING_UTF8 | "SOMACollection" | SOMA type name        | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1"              | SOMA encoding version | required |

### `SOMAMultiscaleImage`

#### `SOMAMultiscaleImage` Spatial Encoding Version 0.2.0

| Metadata key                    | Datatype           | Value                                             | Description                                       | Required |
| :------------------------------ | :----------------- | :------------------------------------------------ | :------------------------------------------------ | :------- |
| "soma_datatype"                 | TILEDB_STRING_UTF8 | "SOMAMultiscaleImage"                             | SOMA type name                                    | required |
| "soma_encoding_version"         | TILEDB_STRING_UTF8 | "1.1.0"                                           | general SOMA encoding version                     | required |
| "soma_spatial_encoding_version" | TILEDB_STRING_UTF8 | "0.2.0"                                           | type-specific SOMA encoding version               | required |
| "soma_coordinate_space"         | TILEDB_STRING_UTF8 | JSON serialization of the coordinate space        | JSON serialization of the coordinate space        | required |
| "soma_multiscale_image_schema"  | TILEDB_STRING_UTF8 | JSON serialization fo the multiscale image schema | JSON serialization fo the multiscale image schema | required |

Coordinate Space JSON

<!-- TODO: Add description of coordinate space metadata -->

Multiscale Image Schema JSON

<!-- TODO: Add description of multiscale image metadata -->

#### `SOMAMultiscaleImage` Spatial Encoding Version 0.1.0

| Metadata key                    | Datatype           | Value                                             | Description                                       | Required |
| :------------------------------ | :----------------- | :------------------------------------------------ | :------------------------------------------------ | :------- |
| "soma_datatype"                 | TILEDB_STRING_UTF8 | "SOMAMultiscaleImage"                             | SOMA type name                                    | required |
| "soma_encoding_version"         | TILEDB_STRING_UTF8 | "1.1.0"                                           | general SOMA encoding version                     | required |
| "soma_spatial_encoding_version" | TILEDB_STRING_UTF8 | "0.1.0"                                           | type-specific SOMA encoding version               | required |
| "soma_coordinate_space"         | TILEDB_STRING_UTF8 | JSON serialization of the coordinate space        | JSON serialization of the coordinate space        | required |
| "soma_multiscale_image_schema"  | TILEDB_STRING_UTF8 | JSON serialization fo the multiscale image schema | JSON serialization fo the multiscale image schema | required |

Coordinate Space JSON

<!-- TODO: Add description of coordinate space metadata -->

Multiscale Image Schema JSON

<!-- TODO: Add description of multiscale image metadata -->

## SOMA Encoding Version 1.0.0

### SOMACollection

| Metadata key            | Datatype           | Value            | Description           | Required |
| :---------------------- | :----------------- | :--------------- | :-------------------- | :------- |
| "soma_datatype"         | TILEDB_STRING_UTF8 | "SOMACollection" | SOMA type name        | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1"              | SOMA encoding version | required |
