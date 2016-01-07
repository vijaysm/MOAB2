# Documentation
# Optional tools
option ( MOAB_BUILD_MBCONVERT    "Build the MOAB mesh converter tool?"           ON )
option ( MOAB_BUILD_HEXMODOPS    "Build the hexahedral modops tool?"             ON )
option ( MOAB_BUILD_MBCHACO      "Build MOAB Chaco interface?"                  OFF )
option ( MOAB_BUILD_MBSIZE       "Build the mbsize tool?"                        ON )
option ( MOAB_BUILD_MBMEM        "Build the mbmem tool?"                         ON )
option ( MOAB_BUILD_MBSKIN       "Build the mbskin tool?"                        ON )
option ( MOAB_BUILD_MBDEPTH      "Build the mbdepth tool?"                       ON )
option ( MOAB_BUILD_MBTAGPROP    "Build the mbtagprop tool?"                     ON )
option ( MOAB_BUILD_MBGSETS      "Build the mbgsets tool?"                       ON )
option ( MOAB_BUILD_SPHEREDECOMP "Build the sphere decomposition tool?"          ON )
option ( MOAB_BUILD_MBSURFPLOT   "Build the mbsurfplot application?"             ON )
option ( MOAB_BUILD_MBPART       "Build the mbpart (Zoltan-based) partitioner?" OFF )
option ( MOAB_BUILD_DAGMC        "Build the dagmc tool?"                         ON )
option ( MOAB_BUILD_MBCOUPLER    "Build the mesh coupler tool?"                  ON )

if ( MOAB_HAVE_MPI AND MPI_FOUND )
  add_subdirectory( refiner )
endif ( )

if (ENABLE_MESQUITE AND MOAB_HAVE_CGM)
  add_subdirectory(mesquite)
endif(ENABLE_MESQUITE AND MOAB_HAVE_CGM)

include_directories(
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src
  ${CMAKE_BINARY_DIR}/src/parallel
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/parallel
  ${CMAKE_CURRENT_SOURCE_DIR}
)

set( TOOLS_LINK_LIBS MOAB )
if ( MOAB_HAVE_CGM )
  set( TOOLS_LINK_LIBS ${TOOLS_LINK_LIBS} ${CGM_LIBRARIES} )
endif ( MOAB_HAVE_CGM )

# MBSIZE
if ( MOAB_BUILD_MBSIZE )
  add_executable( mbsize size.cpp measure.cpp)
  set_source_files_properties( size.cpp measure.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbsize ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBSIZE 1 )
  install(
    TARGETS mbsize
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( )

# MBMEM
if ( MOAB_BUILD_MBMEM )
  add_executable( mbmem mbmem.cpp)
  set_source_files_properties( mbmem.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbmem ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBMEM 1 )
  install(
    TARGETS mbmem
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( )

# MBCONVERT
if ( MOAB_BUILD_MBCONVERT )
  add_executable( mbconvert convert.cpp)
  set_source_files_properties( convert.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbconvert ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBCONVERT 1 )
  install(
    TARGETS mbconvert
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( )

# MBDEPTH
if ( MOAB_BUILD_MBDEPTH )
  add_executable( mbdepth depth.cpp)
  set_source_files_properties( depth.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbdepth ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBDEPTH 1 )
  install(
    TARGETS mbdepth
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( )

# MBSKIN
if ( MOAB_BUILD_MBSKIN )
  add_executable( mbskin skin.cpp)
  set_source_files_properties( skin.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbskin ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBSKIN 1 )
  install(
    TARGETS mbskin
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( )

# MBSURFPLOT
if ( MOAB_BUILD_MBSURFPLOT )
  add_executable( mbsurfplot surfplot.cpp )
  set_source_files_properties( surfplot.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbsurfplot ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBSURFPLOT 1 )
  install(
    TARGETS mbsurfplot
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( )

# MBTAGPROP
if ( MOAB_BUILD_MBTAGPROP )
  add_executable( mbtagprop parse.cpp propagate_tags.cpp  )
  set_source_files_properties( parse.cpp propagate_tags.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbtagprop ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBTAGPROP 1 )
  install(
    TARGETS mbtagprop
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( )

# MBGSETS
if ( MOAB_BUILD_MBGSETS )
  add_executable( mbgsets gsets.cc )
  set_source_files_properties( gsets.cc COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbgsets ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBGSETS 1 )
  install(
    TARGETS mbgsets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif( )

# CGM + DagMC
if ( MOAB_BUILD_DAGMC )
  add_subdirectory( dagmc )
endif( )

# Zoltan
if ( MOAB_BUILD_MBPART AND ( MOAB_HAVE_ZOLTAN OR MOAB_HAVE_METIS ) )
  # Build the mbpart executable
  add_executable( mbpart mbpart.cpp )
  include_directories( ${ZOLTAN_INCLUDES} ${PARMETIS_INCLUDES} ${METIS_INCLUDES} )
  if (MOAB_HAVE_CGM)
    include_directories( ${CGM_INCLUDE_DIRS} )
  endif (MOAB_HAVE_CGM)
  set_source_files_properties( mbpart.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbpart ${TOOLS_LINK_LIBS} ${ZOLTAN_LIBRARIES} ${PARMETIS_LIBRARIES} ${METIS_LIBRARIES} )
  set( MOAB_HAVE_MBPART 1 )
  install(
    TARGETS mbpart
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( MOAB_BUILD_MBPART AND ( MOAB_HAVE_ZOLTAN OR MOAB_HAVE_METIS ) )

# MBCoupler
if ( ENABLE_IMESH AND MOAB_HAVE_MPI AND MPI_FOUND AND MOAB_BUILD_MBCOUPLER )
  add_subdirectory( mbcoupler )
endif ( )

# hexmodops
if( MOAB_BUILD_HEXMODOPS )
  add_executable( hexmodops makeops.cpp )
  target_link_libraries( hexmodops ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_HEXMODOPS 1 )
  install(
    TARGETS hexmodops
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif()
