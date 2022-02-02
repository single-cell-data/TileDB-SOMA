// external includes
#include <string>

// local includes
#include "tiledbsc_export.h"

namespace tiledbsc {

    class TILEDBSC_EXPORT SCDataset {
        public:
            SCDataset();

        public:
            // variables
            std::string name_;
    };

    class TILEDBSC_EXPORT SCDatasetLoader {};
};