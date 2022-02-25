#include <tiledbsc/tiledbsc.h>

using namespace tiledbsc;

int main(int argc, char** argv) {
    (void)argc;

    std::string uri(argv[1]);
    auto sc = SCGroup(uri);

    std::cout << "&sc: " << &sc << std::endl;

    return 0;
};