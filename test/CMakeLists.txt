  add_subdirectory(io)

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

  add_executable( moab_test MBTest.cpp )
  target_link_libraries( moab_test MOAB )
  add_test( TestMOAB ${EXECUTABLE_OUTPUT_PATH}/moab_test -nostress )
  # FIXME: moab_test depends on: test/mb_big_test.g test/cell1.gen test/cell2.gen

  add_executable( seq_man_test TestTypeSequenceManager.cpp )
  target_link_libraries( seq_man_test MOAB )
  set_source_files_properties( TestTypeSequenceManager.cpp
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}" )
   add_test( TestTypeSequenceManager ${EXECUTABLE_OUTPUT_PATH}/seq_man_test )

  add_executable( homxform_test ${MOAB_SOURCE_DIR}/src/HomXform.cpp )
  set_source_files_properties( ${MOAB_SOURCE_DIR}/src/HomXform.cpp
    COMPILE_FLAGS "-DTEST ${MOAB_DEFINES}" )
  target_link_libraries( homxform_test MOAB )
  add_test( TestHomXform ${EXECUTABLE_OUTPUT_PATH}/homxform_test )

  add_executable( scdseq_test scdseq_test.cpp )
  set_source_files_properties( scdseq_test.cpp
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}" )
  target_link_libraries( scdseq_test MOAB )
  add_test( TestSCDSeq ${EXECUTABLE_OUTPUT_PATH}/scdseq_test )

  add_executable( tqdcfr ${MOAB_SOURCE_DIR}/src/io/Tqdcfr.cpp )
  set_source_files_properties( ${MOAB_SOURCE_DIR}/src/io/Tqdcfr.cpp
    COMPILE_FLAGS "-DTEST_TQDCFR -DIS_BUILDING_MB ${MOAB_DEFINES}" )
  target_link_libraries( tqdcfr MOAB )
  add_test( TestTQDCFR ${EXECUTABLE_OUTPUT_PATH}/tqdcfr )

  add_executable( bsp_tree_test bsp_tree_test.cpp )
  set_source_files_properties( bsp_tree_test.cpp
    COMPILE_FLAGS "-DTEST ${MOAB_DEFINES}" )
  target_link_libraries( bsp_tree_test MOAB )
  add_test( TestBSPTree ${EXECUTABLE_OUTPUT_PATH}/bsp_tree_test )


