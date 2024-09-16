#include <iostream>

#include "common.h"

#include <tiledbsoma/tiledbsoma>

int main()
{
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-dev";
    auto [schema, index_columns] = helper::create_arrow_schema(1000, false);

    tiledbsoma::SOMAGeometryDataFrame::create(uri, std::move(schema), ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)), std::vector<std::string>(), ctx);

    std::cout << "Running v1234" << std::endl;
    return 0;
}
