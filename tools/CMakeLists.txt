
  # Documentation
  # Optional tools
  option ( MOAB_BUILD_MBCONVERT    "Build the MOAB mesh converter tool?"  ON )
  option ( MOAB_BUILD_HEXMODOPS    "Build the hexahedral modops tool?"    ON )
  option ( MOAB_BUILD_MBCHACO      "Build MOAB Chaco interface?"          OFF )
  option ( MOAB_BUILD_MBPERF       "Build the MOAB performance tool?"     ON )
  option ( MOAB_BUILD_QVDUAL       "Build the qvdual application?"        OFF )
  option ( MOAB_BUILD_MBSIZE       "Build the mbsize tool?"               ON )
  option ( MOAB_BUILD_MBMEM        "Build the mbmem tool?"                ON )
  option ( MOAB_BUILD_MBSKIN       "Build the mbskin tool?"               ON )
  option ( MOAB_BUILD_MBDEPTH      "Build the mbdepth tool?"              ON )
  option ( MOAB_BUILD_MBTAGPROP    "Build the mbtagprop tool?"            ON )
  option ( MOAB_BUILD_MBGSETS      "Build the mbgsets tool?"              ON )
  option ( MOAB_BUILD_MBCONVERT    "Build the mbconvert tool?"            ON )
  option ( MOAB_BUILD_SPHEREDECOMP "Build the sphere decomposition tool?" ON )
  option ( MOAB_BUILD_MBSURFPLOT   "Build the mbsurfplot application?"    ON )
  option ( MOAB_BUILD_MBZOLTAN     "Build the mbzoltan partitioner?"      OFF )
  option ( MOAB_BUILD_DAGMC        "Build the dagmc tool?"                ON )
  option ( MOAB_BUILD_MBCOUPLER    "Build the mesh coupler tool?"         ON )

  if ( MOAB_USE_MPI AND MPI_FOUND )
    add_subdirectory( refiner )
  endif ( MOAB_USE_MPI AND MPI_FOUND )

  include_directories(
    ${MOAB_BINARY_DIR}
    ${MOAB_BINARY_DIR}/src
    ${MOAB_BINARY_DIR}/src/parallel
    ${MOAB_SOURCE_DIR}/src
    ${MOAB_SOURCE_DIR}/src/parallel
    ${CMAKE_CURRENT_SOURCE_DIR}
  )

  # MBSIZE
  if ( MOAB_BUILD_MBSIZE )
    add_executable( mbsize size.cpp measure.cpp)
    set_source_files_properties( size.cpp measure.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
    target_link_libraries( mbsize MOAB )
  endif ( MOAB_BUILD_MBSIZE )

  # MBMEM
  if ( MOAB_BUILD_MBMEM )
    add_executable( mbmem mbmem.cpp)
    set_source_files_properties( mbmem.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
    target_link_libraries( mbmem MOAB )
  endif ( MOAB_BUILD_MBMEM )

  # MBCONVERT
  if ( MOAB_BUILD_MBCONVERT )
    add_executable( mbconvert convert.cpp)
    set_source_files_properties( convert.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
    target_link_libraries( mbconvert MOAB )
  endif ( MOAB_BUILD_MBCONVERT )

  # MBMEM
  if ( MOAB_BUILD_MBDEPTH )
    add_executable( mbdepth depth.cpp)
    set_source_files_properties( depth.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
    target_link_libraries( mbdepth MOAB )
  endif ( MOAB_BUILD_MBDEPTH )

  # MBSKIN
  if ( MOAB_BUILD_MBSKIN )
    add_executable( mbskin skin.cpp)
    set_source_files_properties( skin.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
    target_link_libraries( mbskin MOAB )
  endif ( MOAB_BUILD_MBSKIN )

  # MBSURFPLOT
  if ( MOAB_BUILD_MBSURFPLOT )
    add_executable( mbsurfplot surfplot.cpp )
    set_source_files_properties( surfplot.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
    target_link_libraries( mbsurfplot MOAB )
  endif ( MOAB_BUILD_MBSURFPLOT )

  # MBTAGPROP
  if ( MOAB_BUILD_MBTAGPROP )
    add_executable( mbtagprop parse.cpp propagate_tags.cpp  )
    set_source_files_properties( parse.cpp propagate_tags.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
    target_link_libraries( mbtagprop MOAB )
  endif ( MOAB_BUILD_MBTAGPROP )

  # MBGSETS
  if ( MOAB_BUILD_MBGSETS )
    add_executable( mbgsets gsets.cc )
    set_source_files_properties( gsets.cc COMPILE_FLAGS "${MOAB_DEFINES}" )
    target_link_libraries( mbgsets MOAB )
  endif ( MOAB_BUILD_MBGSETS )

  # CGM
  if ( MOAB_BUILD_DAGMC )
    add_subdirectory( dagmc )
  endif ( MOAB_BUILD_DAGMC )

  # Zoltan
  if ( MOAB_BUILD_MBZOLTAN )
    add_subdirectory( mbzoltan )
  endif ( MOAB_BUILD_MBZOLTAN )

  # Chaco

  # QVDual
  if ( MOAB_BUILD_QVDUAL )
    add_subdirectory( qvdual )
  endif ( MOAB_BUILD_QVDUAL )

  # MBCoupler
  if ( MOAB_USE_MPI AND MPI_FOUND AND MOAB_BUILD_MBCOUPLER )
#  if ( MOAB_BUILD_MBCOUPLER )
    add_subdirectory( mbcoupler )
  endif ( MOAB_USE_MPI AND MPI_FOUND AND MOAB_BUILD_MBCOUPLER )
