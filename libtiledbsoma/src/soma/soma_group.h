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
    //= public non-static
    //===================================================================

    SOMAGroup(
        std::shared_ptr<Context>, const std::string&, tiledb_query_type_t);

    SOMAGroup(const SOMAGroup&) = default;
    SOMAGroup(SOMAGroup&&) = default;
    SOMAGroup& operator=(const SOMAGroup&) = default;
    SOMAGroup& operator=(SOMAGroup&&) = default;

    /**
     * Closes the SOMAGroup object.
     */
    void close();

    /**
     * Get the SOMAGroup URI.
     */
    std::string uri() const;

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
    void add_member(const std::string&, const bool&, const std::string&);

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
