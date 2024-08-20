/**
 * @file   soma_object.h
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
 * This file defines the SOMAObject class. SOMAObject is an abstract base
 * class for all public SOMA classes: SOMAObject, SOMAExperiment,
 * SOMAMeasurement, SOMA{Sparse,Dense}NdArray, and SOMADataFrame.
 */

#ifndef SOMA_OBJECT
#define SOMA_OBJECT

#include <filesystem>
#include <map>
#include <string>
#include <tiledb/tiledb>

#include "soma_context.h"

namespace tiledbsoma {

using namespace tiledb;
class SOMAObject {
   public:
    //===================================================================
    //= public non-static
    //===================================================================
    virtual ~SOMAObject() = default;

    static std::unique_ptr<SOMAObject> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt,
        std::optional<std::string> soma_type = std::nullopt);

    /**
     * @brief Return a constant string describing the type of the object.
     */
    const std::optional<std::string> type();

    /**
     * @brief Get URI of the SOMAObject.
     *
     * @return std::string URI
     */
    virtual const std::string uri() const = 0;

    /**
     * @brief Get the context associated with the SOMAObject.
     *
     * @return SOMAContext
     */
    virtual std::shared_ptr<SOMAContext> ctx() = 0;

    /**
     * Get whether the SOMAObject was open in read or write mode.
     *
     * @return OpenMode
     */
    virtual OpenMode mode() const = 0;

    /**
     * @brief Close the SOMAObject.
     */
    virtual void close() = 0;

    /**
     * @brief Check if the SOMAObject is open.
     */
    virtual bool is_open() const = 0;

    /**
     * Set metadata key-value items to a SOMAObject. The SOMAObject must
     * opened in WRITE mode, otherwise the function will error out.
     *
     * @param key The key of the metadata item to be added. UTF-8 encodings
     *     are acceptable.
     * @param value_type The datatype of the value.
     * @param value_num The value may consist of more than one items of the
     *     same datatype. This argument indicates the number of items in the
     *     value component of the metadata.
     * @param value The metadata value in binary form.
     * @param force A boolean toggle to suppress internal checks, defaults to
     *     false.
     *
     * @note The writes will take effect only upon closing the array.
     */
    virtual void set_metadata(
        const std::string& key,
        tiledb_datatype_t value_type,
        uint32_t value_num,
        const void* value,
        bool force = false) = 0;

    /**
     * Delete a metadata key-value item from an open SOMAObject. The
     * SOMAObject must be opened in WRITE mode, otherwise the function will
     * error out.
     *
     * @param key The key of the metadata item to be deleted.
     *
     * @note The writes will take effect only upon closing the group.
     *
     * @note If the key does not exist, this will take no effect
     *     (i.e., the function will not error out).
     */
    virtual void delete_metadata(const std::string& key) = 0;

    /**
     * @brief Given a key, get the associated value datatype, number of
     * values, and value in binary form.
     *
     * The value may consist of more than one items of the same datatype.
     Keys
     * that do not exist in the metadata will return std::nullopt for the value.
     *
     * **Example:**
     * @code{.cpp}
     * // Open the group for reading
     * tiledbsoma::SOMAGroup soma_group = SOMAGroup::open(TILEDB_READ,
     "s3://bucket-name/group-name");
     * tiledbsoma::MetadataValue meta_val = soma_group->get_metadata("key");
     * std::string key = std::get<MetadataInfo::key>(meta_val);
     * tiledb_datatype_t dtype = std::get<MetadataInfo::dtype>(meta_val);
     * uint32_t num = std::get<MetadataInfo::num>(meta_val);
     * const void* value = *((const
     int32_t*)std::get<MetadataInfo::value>(meta_val));
     * @endcode
     *
     * @param key The key of the metadata item to be retrieved. UTF-8
     encodings
     *     are acceptable.
     * @return std::optional<MetadataValue>
     */
    virtual std::optional<MetadataValue> get_metadata(
        const std::string& key) = 0;

    /**
     * Get a mapping of all metadata keys with its associated value datatype,
     * number of values, and value in binary form.
     *
     * @return std::map<std::string, MetadataValue>
     */
    virtual std::map<std::string, MetadataValue> get_metadata() = 0;

    /**
     * Check if the key exists in metadata from an open SOMAObject.
     *
     * @param key The key of the metadata item to be checked. UTF-8 encodings
     *     are acceptable.
     * @return true if the key exists, else false.
     */
    virtual bool has_metadata(const std::string& key) = 0;

    /**
     * Return then number of metadata items in an open SOMAObject. The group
     * must be opened in READ mode, otherwise the function will error out.
     */
    virtual uint64_t metadata_num() const = 0;
};
}  // namespace tiledbsoma

#endif  // SOMA_OBJECT
