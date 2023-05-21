# Find a dependency using find_package, fallback to using pkg-config otherwise.
function (FIND_PACKAGE_WITH_FALLBACK _FIND_CMAKE_NAME)

    cmake_parse_arguments(
        _FIND
        "QUIET;REQUIRED"
        "PC_NAME"
        ""
        ${ARGN}
    )

    if (NOT DEFINED _FIND_CMAKE_NAME)
        message(WARNING
            "First argument wasn't provided to find_package_with_fallback function.")
    endif ()

    # Try using CMake find_package first. 
    find_package(${_FIND_CMAKE_NAME} QUIET NO_CMAKE_PATH)

    # Use pkg_config if find_package didn't work or wasn't requested.
    if (DEFINED _FIND_PC_NAME AND NOT ${_FIND_CMAKE_NAME}_FOUND) 
        find_package(PkgConfig)
        if (PkgConfig_FOUND)
            set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH ON)
            pkg_search_module(${_FIND_CMAKE_NAME} ${_FIND_PC_NAME})

            # Mark as found by pkg-config
            if (${_FIND_CMAKE_NAME}_FOUND)
                set(${_FIND_CMAKE_NAME}_PC_FOUND TRUE CACHE BOOL "${_FIND_CMAKE_NAME} found with pkg-config.")
            endif()
        endif ()
    endif ()
    
    if (NOT ${_FIND_CMAKE_NAME}_FOUND AND NOT _FIND_QUIET)
        if (_FIND_REQUIRED)
            message(FATAL_ERROR "Could not find ${_FIND_CMAKE_NAME}.")
        else ()
            message(WARNING "Could not find ${_FIND_CMAKE_NAME}. Ignoring because not marked as REQUIRED.")
        endif ()
    endif ()

endfunction (FIND_PACKAGE_WITH_FALLBACK)

# Create an alias target from pkg-config variables.
function (CREATE_ALIAS_TARGET _ALIAS_NAME _TARGET_NAME)

    add_library(lib_${_TARGET_NAME} INTERFACE)

    target_include_directories(lib_${_TARGET_NAME} INTERFACE
        ${${_TARGET_NAME}_INCLUDE_DIRECTORIES})
    target_link_libraries(lib_${_TARGET_NAME} INTERFACE
        ${${_TARGET_NAME}_LIBRARIES})
    target_link_options(lib_${_TARGET_NAME} INTERFACE
        ${${_TARGET_NAME}_LDFLAGS_OTHER})
    target_compile_options(lib_${_TARGET_NAME} INTERFACE
        ${${_TARGET_NAME}_CFLAGS})
    
    add_library(${_ALIAS_NAME} ALIAS lib_${_TARGET_NAME})

endfunction (CREATE_ALIAS_TARGET)