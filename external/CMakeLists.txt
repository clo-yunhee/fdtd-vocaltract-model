include(ExternalProject)

set(EIGEN_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen-install/")

ExternalProject_Add(
    eigen
    URL  https://gitlab.com/libeigen/eigen/-/archive/3026fc0d3cf7becebd0761d8213536071eb1aa50/eigen-3026fc0d3cf7becebd0761d8213536071eb1aa50.tar.gz
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