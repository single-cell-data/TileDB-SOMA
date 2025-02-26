/**
 * @file   enums.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines enumerated values.
 */

#ifndef SOMA_ENUMS
#define SOMA_ENUMS

/** Defines whether the SOMAObject should be opened in read or write mode */
enum class OpenMode { read = 0, write };

/** Defines whether the result should be opened in row-major or column-major
 * order */
enum class ResultOrder { automatic = 0, rowmajor, colmajor, unordered, global };

/** Defines whether the SOMAGroup URI is absolute or relative */
enum class URIType { automatic = 0, absolute, relative };

typedef enum {
    SOMA_COLUMN_DIMENSION = 0,
    SOMA_COLUMN_ATTRIBUTE = 1,
    SOMA_COLUMN_GEOMETRY = 2
} soma_column_datatype_t;

// This enables some code deduplication between core domain, core current
// domain, and core non-empty domain.
enum class Domainish {
    kind_core_domain = 0,
    kind_core_current_domain = 1,
    kind_non_empty_domain = 2
};

#endif  // SOMA_ENUMS
