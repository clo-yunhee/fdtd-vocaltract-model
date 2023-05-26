include(ExternalProject)

set(EXTERNAL_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(EXTERNAL_BUILD_ROOT "${CMAKE_CURRENT_BINARY_DIR}")

### Colormap

set(COLORMAP_INSTALL_DIR "${EXTERNAL_BUILD_ROOT}/colormap-install/")

FetchContent_Declare(
    colormap
    GIT_REPOSITORY https://github.com/jgreitemann/colormap
    GIT_TAG        40f3a7d8cb648ecb881e0666450f647846627200
    SOURCE_DIR "${EXTERNAL_SOURCE_ROOT}/colormap-src"
    BINARY_DIR "${EXTERNAL_BUILD_ROOT}/colormap-build"
    INSTALL_DIR "${COLORMAP_INSTALL_DIR}"
    CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    -DCMAKE_BUILD_TYPE=Release
)
FetchContent_MakeAvailable(colormap)

add_library(colormap::colormap ALIAS colormap)