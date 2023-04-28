
#define SPDLOG_FMT_EXTERNAL_HO 1

#include <tiledb/tiledb>
#include <tiledbsoma/tiledbsoma>

int main() {
    auto ver = tiledb::version();
    std::cout << "Seeing "
              << std::get<0>(ver) << "." << std::get<1>(ver) << "." << std::get<2>(ver)
              << " and "
              << tiledbsoma::version::as_string()
              << std::endl;
}
