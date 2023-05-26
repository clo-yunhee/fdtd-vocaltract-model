#include <cstdlib>

#include "types.hh"
// #include "routines/find_cell_types.hh"

void talkingTube();

int main(int argc, char** argv) {
    af::setBackend(AF_BACKEND_CUDA);
    af::info();

    talkingTube();

    return EXIT_SUCCESS;
}