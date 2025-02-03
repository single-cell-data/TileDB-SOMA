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
enum class ResultOrder { automatic = 0, rowmajor, colmajor };

/** Defines whether the SOMAGroup URI is absolute or relative */
enum class URIType { automatic = 0, absolute, relative };

#endif  // SOMA_ENUMS
