#ifndef SOMA_COLLECTION_H
#define SOMA_COLLECTION_H

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

namespace tiledbsc {
using namespace tiledb;

class SOMACollection {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Open a SOMACollection at the specified URI and return a
     * SOMACollection object.
     *
     * @param uri URI of SOMACollection
     * @return SOMACollection object
     */
    static SOMACollection open(std::string_view uri);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMACollection object
     *
     * @param uri URI of SOMACollection
     */
    SOMACollection(std::string_view uri);

    /**
     * @brief Return a map of hierarchical SOMA names to SOMA URIs for all
     * SOMAs in the SOMACollection. This includes nested SOMACollections.
     *
     * @return std::unordered_map<std::string, std::string> Map of SOMA name to
     * SOMA URI
     */
    std::unordered_map<std::string, std::string> list_somas();

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // TileDB context
    Context ctx_;

    // SOMACollection URI
    std::string uri_;

    // Map of hierarchical SOMA name to SOMA URI
    std::unordered_map<std::string, std::string> soma_uri_map_;

    /**
     * @brief Walk the TileDB group tree to discover SOMAs and populate the
     * SOMA URI map. This function is called recursively in order to discover
     * all SOMAs in all subgroups.
     *
     * @param group TileDB group to inspect for SOMAs and subgroups.
     * @param parent Hierarchical group name of the group's parent.
     */
    void build_uri_map(Group& group, std::string_view parent = "");

    /**
     * @brief Check if the provided URI is a TileDB Cloud URI.
     *
     * @param uri URI to check
     * @return true URI is a TileBD Cloud URI
     * @return false URI is not a TileBD Cloud URI
     */
    // TODO: move this to utils
    bool is_tiledb_uri(std::string_view uri) {
        return uri.find("tiledb://") == 0;
    }
};

};  // namespace tiledbsc

#endif
