ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/FBiGeom_protos.h
  COMMAND ${CMAKE_COMMAND}
          -DDEF:STRING=MOAB
          -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/FBiGeom.h
          -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/FBiGeom_protos.h
          -Dprefix:STRING=FBiGeom
          -P ${CMAKE_SOURCE_DIR}/config/ITAPSFortranMangling.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../imesh/iMesh.h )

set( FBiGEOM_SOURCE FBiGeom_MOAB.hpp FBiGeom_MOAB.cpp )
set(
  FBiGEOM_INCLUDE_HEADERS
  FBiGeom.h
  ${CMAKE_CURRENT_BINARY_DIR}/FBiGeom_protos.h
  FBiGeom_f.h
  )

add_library(
  FBiGeomMOAB
  ${FBiGEOM_SOURCE}
  ${FBiGEOM_INCLUDE_HEADERS}
  )

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/itaps
  ${CMAKE_BINARY_DIR}/itaps
  ${CMAKE_SOURCE_DIR}/itaps/imesh
  ${CMAKE_BINARY_DIR}/itaps/imesh
 )

if ( MOAB_HAVE_MPI )
  include_directories( ${CMAKE_SOURCE_DIR}/src/parallel
                       ${CMAKE_BINARY_DIR}/src/parallel)
endif ()

target_link_libraries( FBiGeomMOAB MOAB )

set_target_properties(
  FBiGeomMOAB
  PROPERTIES
  VERSION ${MOAB_VERSION}
  SOVERSION ${MOAB_VERSION_MAJOR}
  )
