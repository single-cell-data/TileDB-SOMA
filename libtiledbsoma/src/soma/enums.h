/**
 * @file   enums.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2023 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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