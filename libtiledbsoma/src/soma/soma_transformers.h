#ifndef SOMA_TRANSFORMERS_H
#define SOMA_TRANSFORMERS_H

#include "../utils/transformer.h"
#include "soma_coordinates.h"

namespace tiledbsoma {

class OutlineTransformer : public Transformer {
   public:
    OutlineTransformer(
        SOMACoordinateSpace coordinate_space,
        const std::vector<std::pair<std::string, bool>>& columns);

    virtual ~OutlineTransformer();

    ArrowTable apply(
        std::unique_ptr<ArrowArray>, std::unique_ptr<ArrowSchema>) override;

   private:
    /**
     * @brief Cast an array containing the outer rings of polygons to an Arrow
     * array holding the WKB-encoded polygons and generate the additional
     * index-column arrays based on the spatial axes.
     */
    std::pair<
        std::vector<std::unique_ptr<ArrowArray>>,
        std::vector<std::unique_ptr<ArrowSchema>>>
    _cast_polygon_vertex_list_to_wkb(
        ArrowArray* array,
        ArrowSchema* schema,
        const SOMACoordinateSpace& coordinate_space,
        bool generate_dimensions);

    tiledbsoma::SOMACoordinateSpace coordinate_space;
    std::vector<std::pair<std::string, bool>> columns;
};
}  // namespace tiledbsoma
#endif