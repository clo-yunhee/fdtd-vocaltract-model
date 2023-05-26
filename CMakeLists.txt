cmake_minimum_required(VERSION 3.23...3.26)

# Project name and a few settings.
project(FDTD-VocalTract-Model
    VERSION 0.1
    DESCRIPTION "Finite difference time domain 2.5D vocal tract model"
    LANGUAGES CXX
)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules" ${CMAKE_MODULE_PATH})

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

set(default_build_type "RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    # Let's nicely support folders in IDEs
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # Set /MTd on MSVC
    if(MSVC)
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()

endif()

include(FetchContent)
include(PackageWithFallback)

# Formatting library
FetchContent_Declare(
    fmtlib
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 10.0.0
)
FetchContent_MakeAvailable(fmtlib)

# Sound file I/O library
find_package(SndFile REQUIRED)

# PNG library
find_package(PNG REQUIRED)

# Boost library
find_package(Boost REQUIRED)

find_package(CUDAToolkit 12.0 EXACT REQUIRED)
set(CMAKE_CUDA_ARCHITECTURES all-major)

# ArrayFire library
find_package(ArrayFire REQUIRED)

# External libraries
add_subdirectory(external)

# Project code
add_subdirectory(src)
