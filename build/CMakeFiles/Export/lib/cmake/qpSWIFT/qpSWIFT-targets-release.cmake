#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "qpSWIFT::qpSWIFT-static" for configuration "Release"
set_property(TARGET qpSWIFT::qpSWIFT-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(qpSWIFT::qpSWIFT-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libqpSWIFT.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS qpSWIFT::qpSWIFT-static )
list(APPEND _IMPORT_CHECK_FILES_FOR_qpSWIFT::qpSWIFT-static "${_IMPORT_PREFIX}/lib/libqpSWIFT.a" )

# Import target "qpSWIFT::qpSWIFT-shared" for configuration "Release"
set_property(TARGET qpSWIFT::qpSWIFT-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(qpSWIFT::qpSWIFT-shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libqpSWIFT.so"
  IMPORTED_SONAME_RELEASE "libqpSWIFT.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS qpSWIFT::qpSWIFT-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_qpSWIFT::qpSWIFT-shared "${_IMPORT_PREFIX}/lib/libqpSWIFT.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
