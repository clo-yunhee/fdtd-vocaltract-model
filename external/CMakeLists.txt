include(ExternalProject)

### Eigen

set(EIGEN_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen-install/")

ExternalProject_Add(
    eigen
    URL      https://gitlab.com/libeigen/eigen/-/archive/3026fc0d3cf7becebd0761d8213536071eb1aa50/eigen-3026fc0d3cf7becebd0761d8213536071eb1aa50.tar.gz
    URL_HASH SHA256=4CF5F2D6ABB45F9E043B99BAE7416DBA136FAD7DAF6A54C0A2B8BD74EEB628B1
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen-src"
    BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen-build"
    INSTALL_DIR "${EIGEN_INSTALL_DIR}"
    CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    -DCMAKE_BUILD_TYPE=Release
)

file(MAKE_DIRECTORY ${EIGEN_INSTALL_DIR}/include/eigen3)  # avoid race condition

add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
add_dependencies(Eigen3::Eigen eigen)
target_compile_features(Eigen3::Eigen INTERFACE cxx_std_14)

set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${EIGEN_INSTALL_DIR}/include/eigen3
)

### Colormap

set(COLORMAP_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/colormap-install/")

ExternalProject_Add(
    colormap
    GIT_REPOSITORY https://github.com/jgreitemann/colormap
    GIT_TAG        40f3a7d8cb648ecb881e0666450f647846627200
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/colormap-src"
    BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/colormap-build"
    INSTALL_DIR "${COLORMAP_INSTALL_DIR}"
    CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    -DCMAKE_BUILD_TYPE=Release
)

file(MAKE_DIRECTORY ${COLORMAP_INSTALL_DIR}/include)  # avoid race condition

add_library(colormap::colormap INTERFACE IMPORTED GLOBAL)
add_dependencies(colormap::colormap colormap)
target_compile_features(colormap::colormap INTERFACE cxx_std_14)

set_target_properties(colormap::colormap PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${COLORMAP_INSTALL_DIR}/include
)