  add_subdirectory(io)

  include_directories(
    ${MOAB_BINARY_DIR}
    ${MOAB_BINARY_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${MOAB_SOURCE_DIR}/src
    ${MOAB_BINARY_DIR}/src/parallel
    ${MOAB_SOURCE_DIR}/src/parallel
  )

  add_custom_command(OUTPUT mb_big_test.g
                     COMMAND /bin/zcat ${MOAB_SOURCE_DIR}/MeshFiles/unittest/mb_big_test.g.gz > mb_big_test.g
                     DEPENDS ${MOAB_SOURCE_DIR}/MeshFiles/unittest/mb_big_test.g.gz
                     )
  add_custom_command(OUTPUT cell1.gen
                     COMMAND /bin/zcat ${MOAB_SOURCE_DIR}/MeshFiles/unittest/cell1.gen.gz > cell1.gen
                     DEPENDS ${MOAB_SOURCE_DIR}/MeshFiles/unittest/cell1.gen.gz
                     )
  add_custom_command(OUTPUT cell2.gen
                     COMMAND /bin/zcat ${MOAB_SOURCE_DIR}/MeshFiles/unittest/cell2.gen.gz > cell2.gen
                     DEPENDS ${MOAB_SOURCE_DIR}/MeshFiles/unittest/cell2.gen.gz
                     )
  
  add_custom_target(moab_test_files DEPENDS mb_big_test.g cell1.gen cell2.gen)
  set_source_files_properties( mBTest.cpp
    COMPILE_FLAGS "-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR} -DMESHDIR=${MOAB_SOURCE_DIR}/MeshFiles/unittest" )
  add_executable( moab_test MBTest.cpp )
  add_dependencies(moab_test moab_test_files)
  target_link_libraries( moab_test MOAB )
  add_test( TestMOAB ${EXECUTABLE_OUTPUT_PATH}/moab_test )
  # FIXME: moab_test depends on: test/mb_big_test.g test/cell1.gen test/cell2.gen

  add_executable( tag_test TagTest.cpp )
  target_link_libraries( tag_test MOAB )
  add_test( TagTest ${EXECUTABLE_OUTPUT_PATH}/tag_test )

  add_executable( seq_man_test TestTypeSequenceManager.cpp )
  target_link_libraries( seq_man_test MOAB )
  set_source_files_properties( TestTypeSequenceManager.cpp
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES} -DMESHDIR=${MOAB_SOURCE_DIR}/MeshFiles/unittest" )
  add_test( TestTypeSequenceManager ${EXECUTABLE_OUTPUT_PATH}/seq_man_test )

  add_executable( homxform_test ${MOAB_SOURCE_DIR}/src/HomXform.cpp )
  set_source_files_properties( ${MOAB_SOURCE_DIR}/src/HomXform.cpp
    COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} -DMESHDIR=${MOAB_SOURCE_DIR}/MeshFiles/unittest" )
  target_link_libraries( homxform_test MOAB )
  add_test( TestHomXform ${EXECUTABLE_OUTPUT_PATH}/homxform_test )

  add_executable( scdseq_test scdseq_test.cpp )
  set_source_files_properties( scdseq_test.cpp
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES} -DMESHDIR=${MOAB_SOURCE_DIR}/MeshFiles/unittest" )
  target_link_libraries( scdseq_test MOAB )
  add_test( TestSCDSeq ${EXECUTABLE_OUTPUT_PATH}/scdseq_test )

  add_executable( bsp_tree_test bsp_tree_test.cpp )
  set_source_files_properties( bsp_tree_test.cpp
    COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} -DMESHDIR=${MOAB_SOURCE_DIR}/MeshFiles/unittest" )
  target_link_libraries( bsp_tree_test MOAB )
  add_test( TestBSPTree ${EXECUTABLE_OUTPUT_PATH}/bsp_tree_test )

  add_executable( geomutiltests GeomUtilTests.cpp )
  set_source_files_properties( GeomUtilTests.cpp
    COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} -DMESHDIR=${MOAB_SOURCE_DIR}/MeshFiles/unittest" )
  target_link_libraries( geomutiltests MOAB )
  add_test( TestGeomUtil ${EXECUTABLE_OUTPUT_PATH}/geomutiltests )

  add_executable( range_test range_test.cpp )
  set_source_files_properties( range_test.cpp
    COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} -DMESHDIR=${MOAB_SOURCE_DIR}/MeshFiles/unittest" )
  target_link_libraries( range_test MOAB )
  add_test( TestRange ${EXECUTABLE_OUTPUT_PATH}/range_test )

if ( MOAB_USE_MPI AND MPI_FOUND )

  add_executable ( mbparallelcomm_test parallel/mbparallelcomm_test.cpp )
  target_link_libraries( mbparallelcomm_test MOAB )
  set_source_files_properties( parallel/mbparallelcomm_test.cpp
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES} -DMESHDIR=${MOAB_SOURCE_DIR}/MeshFiles/unittest" )
  add_test( TestParallelComm-BcastDelete
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} 0 ${MOAB_SOURCE_DIR}/parallel/ptest.cub )
  add_test( TestParallelComm-ReadDelete
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} -1 ${MOAB_SOURCE_DIR}/parallel/ptest.cub )
  add_test( TestParallelComm-ReadParallel
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} -2 ${MOAB_SOURCE_DIR}/parallel/ptest.cub )
  add_test( TestParallelComm-Broadcast
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} -3 ${MOAB_SOURCE_DIR}/parallel/ptest.cub )

  if ( MOAB_USE_HDF )
    include_directories(
      ${HDF5_INCLUDE_DIR}
      ${MOAB_SOURCE_DIR}/src/io/mhdf/include
    )
    set_source_files_properties( parallel/mhdf_parallel.c
      COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} -DMESHDIR=${MOAB_SOURCE_DIR}/MeshFiles/unittest" )
    add_executable( mhdf_parallel parallel/mhdf_parallel.c )
    target_link_libraries( mhdf_parallel mhdf MOAB MOABpar )
    add_test( TestMHDFParallel ${EXECUTABLE_OUTPUT_PATH}/mhdf_parallel )
  endif ( MOAB_USE_HDF )

  set_source_files_properties( parallel/parallel_unit_tests.cpp
    COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} -DMESHDIR=${MOAB_SOURCE_DIR}/MeshFiles/unittest" )
  add_executable ( parallel_unit_tests parallel/parallel_unit_tests.cpp )
  target_link_libraries( parallel_unit_tests MOAB )
  add_test( TestParallel
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/parallel_unit_tests ${MPIEXEC_POSTFLAGS} -p 2 ${MOAB_SOURCE_DIR}/MeshFiles/unittest/ptest.cub )
  
  set_source_files_properties( parallel/pcomm_serial.cpp
    COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} -DMESHDIR=${MOAB_SOURCE_DIR}/MeshFiles/unittest" )
  add_executable ( pcomm_serial parallel/pcomm_serial.cpp )
  target_link_libraries( pcomm_serial MOAB )
  add_test( TestPCommSerial
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/pcomm_serial ${MPIEXEC_POSTFLAGS} 1 ${MOAB_SOURCE_DIR}/MeshFiles/unittest/ptest.cub )

endif ( MOAB_USE_MPI AND MPI_FOUND )
