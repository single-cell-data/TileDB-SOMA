#ifndef SOMA_TRANSFORMERS_H
#define SOMA_TRANSFORMERS_H

#include "../utils/transformer.h"
#include "common/arrow/utils.h"
#include "soma_coordinates.h"

namespace tiledbsoma {

class OutlineTransformer : public Transformer {
   public:
    OutlineTransformer(SOMACoordinateSpace coordinate_space);

    virtual ~OutlineTransformer();

    common::arrow::ArrowTable apply(
        common::arrow::managed_unique_ptr<ArrowArray>, common::arrow::managed_unique_ptr<ArrowSchema>) override;

   private:
    /**
     * @brief Cast an array containing the outer rings of polygons to an Arrow
     * array holding the WKB-encoded polygons and generate the additional
     * index-column arrays based on the spatial axes.
     */
    std::pair<
        std::vector<common::arrow::managed_unique_ptr<ArrowArray>>,
        std::vector<common::arrow::managed_unique_ptr<ArrowSchema>>>
    _cast_polygon_vertex_list_to_wkb(ArrowArray* array, const SOMACoordinateSpace& coordinate_space);

    tiledbsoma::SOMACoordinateSpace coordinate_space;
};
}  // namespace tiledbsoma
#endif