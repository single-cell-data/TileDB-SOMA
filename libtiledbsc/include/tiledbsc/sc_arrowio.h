#ifndef TILEDBSC_ARROW_H
#define TILEDBSC_ARROW_H


/**      -*-C++-*-
 * vim: set ft=cpp:
 * @file   arrowio
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2020-2021 TileDB, Inc.
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
 * This file defines the experimental TileDB interoperation with Apache Arrow.
 */

#include <memory>

#include <tiledbsc/query_result.h>
#include <tiledbsc/carrow.h>


namespace tiledbsc {

// forward declare
class QueryResult;
class BufferSet;

namespace arrow {

class TILEDBSC_EXPORT ArrowExporter;

struct TILEDBSC_EXPORT ArrowPair {
    ArrowSchema* schema;
    ArrowArray* array;
    bool own_;

    ArrowPair() : own_(true) {
      schema = (ArrowSchema*)malloc(sizeof(ArrowSchema));
      array = (ArrowArray*)malloc(sizeof(ArrowArray));

      if (!schema || !array)
        throw std::runtime_error("Failed to allocate ArrowSchema and ArrowArray structs");
    }

    ~ArrowPair() {
      if (own_) {
        free(schema);
        free(array);
      }
    }

    void disown() {
      own_ = false;
    }
};

/**
 * Adapter to export TileDB (read) Query results to Apache Arrow buffers
 * and import Arrow buffers into a TileDB (write) Query.
 *
 * This adapter exports buffers conforming to the Arrow C Data Interface
 * as documented at:
 *
 *   https://arrow.apache.org/docs/format/CDataInterface.html
 *
 */
class TILEDBSC_EXPORT ArrowAdapter {
public:
  /* Constructs an ArrowAdapter wrapping the given QueryResult */
  ArrowAdapter(tiledbsc::QueryResult& qr);
  ~ArrowAdapter();

  /**
   * Exports named Query buffer to ArrowArray/ArrowSchema struct pair,
   * as defined in the Arrow C Data Interface.
   *
   * @param name The name of the buffer to export.
   * @param arrow_array Pointer to pre-allocated ArrowArray struct
   * @param arrow_schema Pointer to pre-allocated ArrowSchema struct
   * @throws tiledb::TileDBError with error-specific message.
   */
  void export_array(const char* name, void* arrow_array, void* arrow_schema);

   /**
   * Exports given BufferSet to ArrowArray/ArrowSchema struct pair,
   * as defined in the Arrow C Data Interface.
   *
   * @param BufferSet The BufferSet to export.
   * @param arrow_array Pointer to pre-allocated ArrowArray struct
   * @param arrow_schema Pointer to pre-allocated ArrowSchema struct
   * @throws tiledb::TileDBError with error-specific message.
   */
  static void export_buffer(BufferSet& buffer_set, ArrowArray* arrow_array, ArrowSchema* arrow_schema);

  /**
   * Exports *all* results in QueryResults to ArrowArray/ArrowSchema
   *
   * @param arrow_array Pointer to pre-allocated ArrowArray struct
   * @param arrow_schema Pointer to pre-allocated ArrowSchema struct
   * @throws tiledb::TileDBError with error-specific message.
   */
  void export_table(void* arrow_array, void* arrow_schema);

  /**
   * Set named Query buffer from ArrowArray/ArrowSchema struct pair
   * representing external data buffers conforming to the
   * Arrow C Data Interface.
   *
   * @param name The name of the buffer to export.
   * @param arrow_array Pointer to pre-allocated ArrowArray struct
   * @param arrow_schema Pointer to pre-allocated ArrowSchema struct
   * @throws tiledb::TileDBError with error-specific message.
   */
  //void import_buffer(const char* name, void* arrow_array, void* arrow_schema);

private:
  //ArrowImporter* importer_;
  ArrowExporter* exporter_;
};

}  // end namespace arrow
}  // end namespace tiledb

//#include "arrow_io_impl.h"

#endif // TILEDBSC_ARROW_H