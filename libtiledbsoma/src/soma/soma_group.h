/**
 * @file   soma_group.h
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
 *   This declares the SOMAGroup class.
 */

#ifndef SOMA_GROUP
#define SOMA_GROUP

#include <future>
#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'
#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

namespace tiledbsoma {
using namespace tiledb;

class SOMAGroup {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMAGroup
     *
     * * **Example:**
     * @code{.cpp}
     * tiledb::Group::create(ctx, "s3://bucket-name/group-name");
     * @endcode
     *
     * @param ctx tiledb context
     * @param uri URI where group will be created.
     */
    static void create(std::shared_ptr<Context> ctx, const std::string& uri);

    /**
     * @brief Open a group at the specified URI and return SOMAGroup
     * object.
     *
     * @param mode TILEDB_READ or TILEDB_WRITE
     * @param uri URI of the array
     * @param platform_config Config parameter dictionary
     * @return std::unique_ptr<SOMAGroup> SOMAGroup
     */
    static std::unique_ptr<SOMAGroup> open(
        tiledb_query_type_t mode,
        std::string_view uri,
        std::map<std::string, std::string> platform_config = {});

    /**
     * @brief Open a group at the specified URI and return SOMAGroup
     * object.
     *

     * @return std::unique_ptr<SOMAGroup> SOMAGroup
     */
    static std::unique_ptr<SOMAGroup> open(
        tiledb_query_type_t mode,
        std::shared_ptr<Context> ctx,
        std::string_view uri);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMACollection object.
     *
     * @param mode TILEDB_READ or TILEDB_WRITE
     * @param uri URI of the array
     * @param ctx TileDB context
     */
    SOMAGroup(
        tiledb_query_type_t mode,
        std::string_view uri,
        std::shared_ptr<Context> ctx);

    SOMAGroup(const SOMAGroup&) = default;
    SOMAGroup(SOMAGroup&&) = default;
    SOMAGroup& operator=(const SOMAGroup&) = default;
    SOMAGroup& operator=(SOMAGroup&&) = default;

    /**
     * Opens the SOMAGroup object.
     */
    void open(tiledb_query_type_t);

    /**
     * Closes the SOMAGroup object.
     */
    void close();

    /**
     * Get the SOMAGroup URI.
     */
    std::string uri() const;

    /**
     * Get the Context associated with the SOMAGroup.
     *
     * @return std::shared_ptr<Context>
     */
    std::shared_ptr<Context> ctx();

    /**
     * Get a member from the SOMAGroup given the index.
     *
     * @param index of member
     */
    tiledb::Object get_member(uint64_t) const;

    /**
     * Get a member from the SOMAGroup given the name.
     *
     * @param name of member
     */
    tiledb::Object get_member(const std::string&) const;

    /**
     * Check if the SOMAGroup contains the given name.
     *
     * @param name of member
     */
    bool has_member(const std::string&);

    /**
     * Add a named member to a SOMAGroup.
     *
     * @param uri of member to add
     * @param relative is the URI relative to the SOMAGroup location
     * @param name of member
     */
    void add_member(const std::string&, bool, const std::string&);

    /**
     * Get the number of members in the SOMAGroup.
     */
    uint64_t get_length() const;

    /**
     * Remove a named member from the SOMAGroup.
     *
     * @param name of member
     */
    void remove_member(const std::string&);

    /**
     * Returns a SOMAGroup member to URI mapping.
     *
     * @return std::map<std::string, std::string
     */
    std::map<std::string, std::string> member_to_uri_mapping() const;

    /**
     * Set a metadata key-value item to an open group. The array must
     * be opened in WRITE mode, otherwise the function will error out.
     *
     * @param key The key of the metadata item to be added. UTF-8 encodings
     *     are acceptable.
     * @param value_type The datatype of the value.
     * @param value_num The value may consist of more than one items of the
     *     same datatype. This argument indicates the number of items in the
     *     value component of the metadata.
     * @param value The metadata value in binary form.
     *
     * @note The writes will take effect only upon closing the array.
     */
    void set_metadata(
        const std::string&, tiledb_datatype_t, uint32_t, const void*);

    /**
     * It deletes a metadata key-value item from an open group. The array must
     * be opened in WRITE mode, otherwise the function will error out.
     *
     * @param key The key of the metadata item to be deleted.
     *
     * @note The writes will take effect only upon closing the array.
     *
     * @note If the key does not exist, this will take no effect
     *     (i.e., the function will not error out).
     */
    void delete_metadata(const std::string&);

    /**
     * @brief Get the value of a metadata key-value item given the key.
     *
     * @param key The key of the metadata item to be retrieved. UTF-8 encodings
     *     are acceptable.
     * @param value_type The datatype of the value.
     * @param value_num The value may consist of more than one items of the
     *     same datatype. This argument indicates the number of items in the
     *     value component of the metadata. Keys with empty values are indicated
     *     by value_num == 1 and value == NULL.
     * @param value The metadata value in binary form.
     *
     * @note If the key does not exist, then `value` will be NULL.
     */
    void get_metadata(
        const std::string&, tiledb_datatype_t*, uint32_t*, const void**);

    /**
     * It gets a metadata item from an open group using an index.
     * The array must be opened in READ mode, otherwise the function will
     * error out.
     *
     * @param index The index used to get the metadata.
     * @param key The metadata key.
     * @param value_type The datatype of the value.
     * @param value_num The value may consist of more than one items of the
     *     same datatype. This argument indicates the number of items in the
     *     value component of the metadata. Keys with empty values are indicated
     *     by value_num == 1 and value == NULL.
     * @param value The metadata value in binary form.
     */
    void get_metadata_from_index(
        uint64_t, std::string*, tiledb_datatype_t*, uint32_t*, const void**);

    /**
     * Checks if key exists in metadata from an open group. The array must
     * be opened in READ mode, otherwise the function will error out.
     *
     * @param key The key of the metadata item to be retrieved. UTF-8 encodings
     *     are acceptable.
     * @param value_type The datatype of the value associated with the key (if
     * any).
     * @return true if the key exists, else false.
     * @note If the key does not exist, then `value_type` will not be modified.
     */
    bool has_metadata(const std::string&, tiledb_datatype_t*);

    /**
     * Returns then number of metadata items in an open group. The array must
     * be opened in READ mode, otherwise the function will error out.
     */
    uint64_t metadata_num() const;

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // TileDB context
    std::shared_ptr<Context> ctx_;

    // SOMAGroup
    std::unique_ptr<Group> group_;
};

}  // namespace tiledbsoma

#endif  // SOMA_GROUP
