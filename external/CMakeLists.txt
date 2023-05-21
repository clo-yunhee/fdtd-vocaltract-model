include(ExternalProject)

set(EXTERNAL_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(EXTERNAL_BUILD_ROOT "${CMAKE_CURRENT_BINARY_DIR}")

### Eigen

set(EIGEN_INSTALL_DIR "${EXTERNAL_BUILD_ROOT}/eigen-install/")

FetchContent_Declare(
    eigen
    URL      https://gitlab.com/libeigen/eigen/-/archive/3026fc0d3cf7becebd0761d8213536071eb1aa50/eigen-3026fc0d3cf7becebd0761d8213536071eb1aa50.tar.gz
    URL_HASH SHA256=4CF5F2D6ABB45F9E043B99BAE7416DBA136FAD7DAF6A54C0A2B8BD74EEB628B1
    SOURCE_DIR "${EXTERNAL_SOURCE_ROOT}/eigen-src"
    BINARY_DIR "${EXTERNAL_BUILD_ROOT}/eigen-build"
    INSTALL_DIR "${EIGEN_INSTALL_DIR}"
    CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    -DCMAKE_BUILD_TYPE=Release
)
FetchContent_MakeAvailable(eigen)

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