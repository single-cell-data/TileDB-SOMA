#ifndef SOMA_H
#define SOMA_H

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

namespace tiledbsc {
using namespace tiledb;

class SOMA {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Open a SOMA at the specified URI and return a SOMA object.
     *
     * @param uri URI of SOMA
     * @return SOMA object
     */
    static SOMA open(std::string_view uri);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMA object
     *
     * @param uri URI of SOMA
     */
    SOMA(std::string_view uri);

    /**
     * @brief Return a map of hierarchical array names to array URIs for all
     * arrays in the SOMA.
     *
     * NOTE: If the SOMA URI is *not* a TileDB Cloud URI and the array URIs are
     * TileDB Cloud URIs, the array URIs are converted to relative URIs.
     *
     * @return std::unordered_map<std::string, std::string> Map of array name to
     * array URI
     */
    std::unordered_map<std::string, std::string> list_arrays();

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // TileDB context
    Context ctx_;

    // SOMA URI
    std::string uri_;

    // Map of array name to array URI
    std::unordered_map<std::string, std::string> array_uri_map_;

    // Flag that is true if TileDB Cloud URIs were converted to relative URIs
    bool group_uri_override_ = false;

    /**
     * @brief Walk the TileDB group tree to discover arrays and populate the
     * array URI map. This function is called recursively in order to discover
     * all arrays in all subgroups.
     *
     * @param group TileDB group to inspect for arrays and subgroups.
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
