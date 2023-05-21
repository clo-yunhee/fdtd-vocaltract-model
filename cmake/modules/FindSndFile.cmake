find_package_with_fallback(SndFile REQUIRED PC_NAME sndfile)

if (SndFile_PC_FOUND)
    create_alias_target(SndFile::sndfile SndFile)
endif ()