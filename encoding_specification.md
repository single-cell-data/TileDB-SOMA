# Encoding specification

## TileDB-SOMA encoding version 1.1.0

### SOMA abstract specification support

- TileDB-SOMA encoding version 1.1.0 with spatial encoding version 0.1.0 satisfies the requirements for [SOMA abstract specification](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md) 0.3.0-dev.

- TileDB-SOMA encoding version 1.1.0 with spatial encoding version 0.2.0 satifies the requirements for [SOMA abstract specification](https://github.com/single-cell-data/SOMA/blob/main/abstract_specification.md) 0.3.1-dev.

### Collection encoding version 1.1.0

A SOMA Collection is stored as a TileDB group. The group has the following metadata.

| Metadata key | Datatype | Value | Description | Required |
| :---------------------- | :----------------- | :--------------- | :-------------------- | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMACollection" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0" | SOMA encoding version | required |

### DataFrame encoding version 1.1.0

A SOMA DataFrame is stored as a sparse TileDB array. Index columns are stored as TileDB dimensions. All other columns are stored as TileDB attributes. The array has the following metadata.

| Metadata key | Datatype | Value | Description | Required |
| :---------------------- | :----------------- | :-------------- | :-------------------- | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMADataFrame" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0" | SOMA encoding version | required |

### DenseNDArray encoding version 1.1.0

A SOMA DenseNDArray is stored as a dense TileDB array. The array has the following metadata.

| Metadata key | Datatype | Value | Description | Required |
| :---------------------- | :----------------- | :----------------- | :-------------------- | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMADenseNDArray" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0" | SOMA encoding version | required |

### Experiment encoding version 1.1.0

A SOMA Experiment is stored as a TileDB group. The group has the following metadata.

| Metadata key | Datatype | Value | Description | Required |
| :---------------------- | :----------------- | :--------------- | :-------------------- | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMAExperiment" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0" | SOMA encoding version | required |

### MultiscaleImage encoding version 1.1.0

A SOMA MultiscaleImage is stored as a TileDB group that contains a SOMA DenseNDArray for each resolution level. The group has the following metadata:

| Metadata key | Datatype | Value | Description | Required |
| :------------------------------ | :----------------- | :------------------------------------------------ | :------------------------------------------------ | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMAMultiscaleImage" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0" | general SOMA encoding version | required |
| "soma_spatial_encoding_version" | TILEDB_STRING_UTF8 | "0.2.0" | type-specific SOMA encoding version | required |
| "soma_coordinate_space" | TILEDB_STRING_UTF8 | JSON serialization of the coordinate space | JSON serialization of the coordinate space | required |
| "soma_multiscale_image_schema" | TILEDB_STRING_UTF8 | JSON serialization fo the multiscale image schema | JSON serialization of the multiscale image schema | required |

JSON serialization of the multiscale image schema is given by

```
{
    "type": "object",
    "properties": {
        "data_axis_permutation": {
            "type": "array",
            "items": {"type": "integer"}
        },
        "has_channel_axis": "boolean",
        "shape": {
            "type": "array",
            "items": {"type": "integer"}
        },
        "datatype": {"type": "string"},
    }
}
```

| Serialization Key | Description |
| :---------------------- | :------------------------------------------------------------------------------------------------------------------------ |
| "data_axis_permutation" | The permutation of the coordinate space axis with the channel stored last to the order of the channels as stored on disk. |
| "has_channel_axis" | Whether the image has a dimension specifically for channel information. |
| "shape" | The shape of the level 0 image. |
| "datatype" | The datatype of the image data stored using the PyArrow C-convention. |

### Measurement encoding version 1.1.0

A SOMA Experiment is stored as a TileDB group. The group has the following metadata.

| Metadata key | Datatype | Value | Description | Required |
| :---------------------- | :----------------- | :---------------- | :-------------------- | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMAMeasurement" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0" | SOMA encoding version | required |

### PointCloudDataFrame encoding version 1.1.0

The SOMA PointCloudDataFrame is stored as a sparse TileDB array. Index columns are stored as TileDB dimensions. All other columns are stored as TileDB attributes. The array has the following metadata.

| Metadata key | Datatype | Value | Description | Required |
| :------------------------------ | :----------------- | :----------------------------------------- | :------------------------------------------------------------ | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMAPointCloudDataFrame" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0" | general SOMA encoding version | required |
| "soma_spatial_encoding_version" | TILEDB_STRING_UTF8 | "0.2.0" | type-specific SOMA encoding version | required |
| "soma_coordinate_space" | TILEDB_STRING_UTF8 | JSON serialization of the coordinate space | JSON serialization of the coordinate space | required |
| "soma_geometry_type" | TILEDB_STRING_UTF8 | string with possible values "radius" | optional shape of the stored points | optional |
| "soma_geometry" | TILEDB_FLOAT64 | floating-point value | the radius of circles centered at the points in the dataframe | optional |

### Scene encoding version 1.1.0

The SOMA Scene is stored as a TileDB group. It has the following group metadata:

| Metadata key | Datatype | Value | Description | Required |
| :------------------------------ | :----------------- | :----------------------------------------- | :----------------------------------------- | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMAMultiscaleImage" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0" | general SOMA encoding version | required |
| "soma_spatial_encoding_version" | TILEDB_STRING_UTF8 | "0.2.0" | type-specific SOMA encoding version | required |
| "soma_coordinate_space" | TILEDB_STRING_UTF8 | JSON serialization of the coordinate space | JSON serialization of the coordinate space | optional |

Collections inside the SOMA Scene may store transformation data from the Scene to the elements inside the collections. In addition to the standard group metadata, these collections have the following metadata:

| Metadata key | Datatype | Value | Description | Required |
| :---------------------------------- | :----------------- | :--------------------------------------------- | :----------------------------------------------------------------------------------------------------- | :------- |
| "soma*scene_registry*{object_name}" | TILEDB_STRING_UTF8 | JSON serialization of the coordinate transform | The coordinate transform from a SOMA Scene to the element names `{object_name}` inside the collection. | optional |

The SOMA CoordinateTransform are serialized using the following JSON schemas:

- IdentityTransform serialization

  ```
  {
      "type": "object",
      "properties": {
          "transform_type": {"const": "IdentityTransform"},
          "input_axes": {
              "type": "array",
              "items": {"type": "string"},
          }
          "output_axes": {
              "type": "array",
              "items": {"type": "string"},
          }
      }
  }
  ```

- UniformScaleTransform serialization

  ```
  {
      "type": "object",
      "properties": {
          "transform_type": {"const": "UniformScaleTransform"},
          "input_axes": {
              "type": "array",
              "items": {"type": "string"},
          },
          "output_axes": {
              "type": "array",
              "items": {"type": "string"},
          },
          "scale": {"type": "number"},
      }
  }
  ```

- ScaleTransform serialization

  ```
  {
      "type": "object",
      "properties": {
          "transform_type": {"const": "ScaleTransform"},
          "input_axes": {
              "type": "array",
              "items": {"type": "string"},
          },
          "output_axes": {
              "type": "array",
              "items": {"type": "string"},
          },
          "scale_factors"" {
              "type": "array",
              "items": {"type": "number"},
          },
      }
  }
  ```

- AffineTransform serializaton

  ```
  {
      "type": "object",
      "properties": {
          "transform_type": {"const": "AffineTransform"},
          "input_axes": {
              "type": "array",
              "items": {"type": "string"},
          },
          "output_axes": {
              "type": "array",
              "items": {"type": "string"},
          }
          "matrix": {
              "type": "array",
          }
      }
  }
  ```

| Serialization Key | Description |
| :---------------- | :----------------------------------------------------------------------------------------------------------------------------------------------- |
| "transform_type" | The type of the coordinate space transform. Must be one of "IdentityTransform", "UniformScaleTransform", "ScaleTransform", or "AffineTransform". |
| "input_axes" | The names of the axes of the input coordinate space to the transform. |
| "output_axes" | The names of the axes of the output coordinate space to the transform. |
| "scale" | The scale factor for the "UniformScaleTransform". |
| "scale_factors" | An array of the scale factors for the "ScaleTransform" in order or the axis they are applied to." |
| "maxtrix" | The matrix representation of the "AffineTransform", flattened in row-major order. |

### SparseNDArray encoding version 1.1.0

A SOMA SparseNDArray is stored as a sparse TileDB array. The array has the following metadata.

| Metadata key | Datatype | Value | Description | Required |
| :---------------------- | :----------------- | :------------------ | :-------------------- | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMASparseNDArray" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1.1.0" | SOMA encoding version | required |

### CoordinateSpace serialization encoding version 1.1.0

The SOMA CoordinateSpace is serialized to JSON as an array of axes, stored in order. It satisfies the the following schema.

```
{
    "type": "array",
    "items": {
        "type": "object",
        "properities": {
            "name": {"type": "string"},
            "unit": {"type": ["string", "null"]}
        },
        "required": ["name", "unit"]
    }
}
```

| Serialization Key | Description |
| :---------------- | :--------------------------------------------------------------------- |
| "name" | Name of the axis. |
| "unit" | The name of the units for the axis, or null if no units are specified. |

## TileDB-SOMA encoding version 1.0.0

### SOMA abstract specification

TileDB-SOMA encoding version 1.0.0 satisfies the SOMA abstract specification 0.1.0-dev and 0.2.0-dev.

### Collection encoding version 1.0.0

A SOMA Collection is stored as a TileDB group. The group has the following metadata.

| Metadata key | Datatype | Value | Description | Required |
| :---------------------- | :----------------- | :--------------- | :-------------------- | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMACollection" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1" | SOMA encoding version | required |

### DataFrame encoding version 1.0.0

A SOMA DataFrame is stored as a sparse TileDB array. Index columns are stored as TileDB dimensions. All other columns are stored as TileDB attributes. The array has the following metadata.

| Metadata key | Datatype | Value | Description | Required |
| :---------------------- | :----------------- | :-------------- | :-------------------- | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMADataFrame" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1" | SOMA encoding version | required |

### DenseNDArray encoding version 1.0.0

A SOMA DenseNDArray is stored as a dense TileDB array. The array has the following metadata.

| Metadata key | Datatype | Value | Description | Required |
| :---------------------- | :----------------- | :----------------- | :-------------------- | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMADenseNDArray" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1" | SOMA encoding version | required |

### Experiment encoding version 1.0.0

A SOMA Experiment is stored as a TileDB group. The group has the following metadata.

| Metadata key | Datatype | Value | Description | Required |
| :---------------------- | :----------------- | :--------------- | :-------------------- | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMAExperiment" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1" | SOMA encoding version | required |

### Measurement encoding version 1.0.0

A SOMA Experiment is stored as a TileDB group. The group has the following metadata.

| Metadata key | Datatype | Value | Description | Required |
| :---------------------- | :----------------- | :---------------- | :-------------------- | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMAMeasurement" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1" | SOMA encoding version | required |

### SparseNDArray encoding version 1.0.0

A SOMA SparseNDArray is stored as a sparse TileDB array. The array has the following metadata.

| Metadata key | Datatype | Value | Description | Required |
| :---------------------- | :----------------- | :------------------ | :-------------------- | :------- |
| "soma_object_type" | TILEDB_STRING_UTF8 | "SOMASparseNDArray" | SOMA type name | required |
| "soma_encoding_version" | TILEDB_STRING_UTF8 | "1" | SOMA encoding version | required |
