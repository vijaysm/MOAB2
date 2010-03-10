  
  # Documentation
  # Optional tools
  option ( MOAB_BUILD_MBCONVERT    "Build the MOAB mesh converter tool?"  ON )
  option ( MOAB_BUILD_HEXMODOPS    "Build the hexahedral modops tool?"    ON )
  option ( MOAB_BUILD_MBCHACO      "Build MOAB Chaco interface?"          OFF )
  option ( MOAB_BUILD_MBPERF       "Build the MOAB performance tool?"     ON )
  option ( MOAB_BUILD_QVDUAL       "Build the qvdual application?"        OFF )
  option ( MOAB_BUILD_MBSIZE       "Build the mbsize tool?"               ON )
  option ( MOAB_BUILD_MBSKIN       "Build the mbskin tool?"               ON )
  option ( MOAB_BUILD_MBTAGPROP    "Build the mbtagprop tool?"            ON )
  option ( MOAB_BUILD_SPHEREDECOMP "Build the sphere decomposition tool?" ON )
  option ( MOAB_BUILD_MBSURFPLOT   "Build the mbsurfplot application?"    ON )
  option ( MOAB_BUILD_MBZOLTAN     "Build the mbzoltan partitioner?"      OFF )
  option ( MOAB_BUILD_DAGMC        "Build the dagmc tool?"                ON )
  option ( MOAB_BUILD_MBCOUPLER    "Build the mesh coupler tool?"         ON )

  if ( MOAB_USE_MPI AND MPI_FOUND )
    add_subdirectory( refiner )
    include_directories( ${MOAB_SOURCE_DIR}/refiner )
    target_link_libraries( MOAB MOABrefiner )
  endif ( MOAB_USE_MPI AND MPI_FOUND )


  # Zoltan
  # Chaco

  # QVDual
  if ( MOAB_BUILD_QVDUAL )
    add_subdirectory( qvdual )
  endif ( MOAB_BUILD_QVDUAL )

  # CGM
  if ( MOAB_BUILD_DAGMC )
    add_subdirectory( dagmc )
  endif ( MOAB_BUILD_DAGMC )

  # MBCoupler
  if ( MOAB_USE_MPI AND MPI_FOUND AND MOAB_BUILD_MBCOUPLER )
    add_subdirectory( mbcoupler )
  endif ( MOAB_USE_MPI AND MPI_FOUND AND MOAB_BUILD_MBCOUPLER )
