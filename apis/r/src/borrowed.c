#include "tiledbsoma/nanoarrow/nanoarrow.h"

// We need three entry points from nanoarrow that are declared 'static' in the official
// (now vendored) C library so we bring them here
void ArrowSchemaReleaseInternal(struct ArrowSchema* schema) {
  if (schema->format != NULL) ArrowFree((void*)schema->format);
  if (schema->name != NULL) ArrowFree((void*)schema->name);
  if (schema->metadata != NULL) ArrowFree((void*)schema->metadata);

  // This object owns the memory for all the children, but those
  // children may have been generated elsewhere and might have
  // their own release() callback.
  if (schema->children != NULL) {
    for (int64_t i = 0; i < schema->n_children; i++) {
      if (schema->children[i] != NULL) {
        if (schema->children[i]->release != NULL) {
          ArrowSchemaRelease(schema->children[i]);
        }

        ArrowFree(schema->children[i]);
      }
    }

    ArrowFree(schema->children);
  }

  // This object owns the memory for the dictionary but it
  // may have been generated somewhere else and have its own
  // release() callback.
  if (schema->dictionary != NULL) {
    if (schema->dictionary->release != NULL) {
      ArrowSchemaRelease(schema->dictionary);
    }

    ArrowFree(schema->dictionary);
  }

  // private data not currently used
  if (schema->private_data != NULL) {
    ArrowFree(schema->private_data);
  }

  schema->release = NULL;
}

void ArrowArrayReleaseInternal(struct ArrowArray* array) {
  // Release buffers held by this array
  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;
  if (private_data != NULL) {
    ArrowBitmapReset(&private_data->bitmap);
    ArrowBufferReset(&private_data->buffers[0]);
    ArrowBufferReset(&private_data->buffers[1]);
    ArrowFree(private_data);
  }

  // This object owns the memory for all the children, but those
  // children may have been generated elsewhere and might have
  // their own release() callback.
  if (array->children != NULL) {
    for (int64_t i = 0; i < array->n_children; i++) {
      if (array->children[i] != NULL) {
        if (array->children[i]->release != NULL) {
          ArrowArrayRelease(array->children[i]);
        }

        ArrowFree(array->children[i]);
      }
    }

    ArrowFree(array->children);
  }

  // This object owns the memory for the dictionary but it
  // may have been generated somewhere else and have its own
  // release() callback.
  if (array->dictionary != NULL) {
    if (array->dictionary->release != NULL) {
      ArrowArrayRelease(array->dictionary);
    }

    ArrowFree(array->dictionary);
  }

  // Mark released
  array->release = NULL;
}

ArrowErrorCode ArrowArraySetStorageType(struct ArrowArray* array,
                                               enum ArrowType storage_type) {
  switch (storage_type) {
    case NANOARROW_TYPE_UNINITIALIZED:
    case NANOARROW_TYPE_NA:
      array->n_buffers = 0;
      break;

    case NANOARROW_TYPE_FIXED_SIZE_LIST:
    case NANOARROW_TYPE_STRUCT:
    case NANOARROW_TYPE_SPARSE_UNION:
      array->n_buffers = 1;
      break;

    case NANOARROW_TYPE_LIST:
    case NANOARROW_TYPE_LARGE_LIST:
    case NANOARROW_TYPE_MAP:
    case NANOARROW_TYPE_BOOL:
    case NANOARROW_TYPE_UINT8:
    case NANOARROW_TYPE_INT8:
    case NANOARROW_TYPE_UINT16:
    case NANOARROW_TYPE_INT16:
    case NANOARROW_TYPE_UINT32:
    case NANOARROW_TYPE_INT32:
    case NANOARROW_TYPE_UINT64:
    case NANOARROW_TYPE_INT64:
    case NANOARROW_TYPE_HALF_FLOAT:
    case NANOARROW_TYPE_FLOAT:
    case NANOARROW_TYPE_DOUBLE:
    case NANOARROW_TYPE_DECIMAL128:
    case NANOARROW_TYPE_DECIMAL256:
    case NANOARROW_TYPE_INTERVAL_MONTHS:
    case NANOARROW_TYPE_INTERVAL_DAY_TIME:
    case NANOARROW_TYPE_INTERVAL_MONTH_DAY_NANO:
    case NANOARROW_TYPE_FIXED_SIZE_BINARY:
    case NANOARROW_TYPE_DENSE_UNION:
      array->n_buffers = 2;
      break;

    case NANOARROW_TYPE_STRING:
    case NANOARROW_TYPE_LARGE_STRING:
    case NANOARROW_TYPE_BINARY:
    case NANOARROW_TYPE_LARGE_BINARY:
      array->n_buffers = 3;
      break;

    default:
      return EINVAL;

      return NANOARROW_OK;
  }

  struct ArrowArrayPrivateData* private_data =
      (struct ArrowArrayPrivateData*)array->private_data;
  private_data->storage_type = storage_type;
  return NANOARROW_OK;
}
