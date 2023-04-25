cmake_minimum_required(VERSION 3.13)

project(FTDT-VocalTract-Model LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "CMAKE_BUILD_TYPE was not set: defaults to RelWithDebInfo")
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    add_definitions(-D_WIN32_WINNT=0x0601)
endif()

add_subdirectory(external)
add_subdirectory(src)
        