find_package_with_fallback(PNG REQUIRED PC_NAME libpng)

if (PNG_PC_FOUND)
    create_alias_target(PNG::PNG PNG)
endif ()