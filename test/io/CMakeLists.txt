
  include_directories(
    ${MOAB_SOURCE_DIR}/src
    ${MOAB_BINARY_DIR}/src
    ${MOAB_SOURCE_DIR}/src/parallel
    ${MOAB_BINARY_DIR}/src/parallel
    ${MOAB_SOURCE_DIR}/src/moab
    ${MOAB_BINARY_DIR}/src/moab
    ${MOAB_SOURCE_DIR}/src/parallel/moab
    ${MOAB_BINARY_DIR}/src/parallel/moab
  )

  add_executable( tqdcfr ${MOAB_SOURCE_DIR}/src/io/Tqdcfr.cpp )
  set_source_files_properties( ${MOAB_SOURCE_DIR}/src/io/Tqdcfr.cpp
    COMPILE_FLAGS "-DTEST_TQDCFR -DIS_BUILDING_MB ${MOAB_DEFINES}" )
  target_link_libraries( tqdcfr MOAB )
  add_test( TestTQDCFR ${EXECUTABLE_OUTPUT_PATH}/tqdcfr )
