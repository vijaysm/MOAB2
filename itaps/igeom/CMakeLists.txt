ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/FBiGeom_protos.h
  COMMAND ${CMAKE_COMMAND}
          -DDEF:STRING=MOAB
          -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/FBiGeom.h
          -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/FBiGeom_protos.h
          -Dprefix:STRING=FBiGeom
          -P ${CMAKE_SOURCE_DIR}/config/ITAPSFortranMangling.cmake
  )

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
  ${CMAKE_CURRENT_BINARY_DIR}  # for FBiGeom_protos.h
  ${CMAKE_BINARY_DIR}/src  # for MOAB_FCDefs.h
  ${CMAKE_SOURCE_DIR}/src  # for moab/Interface.hpp
  ${CMAKE_SOURCE_DIR}/itaps/imesh  # for iMesh.h
  ${CMAKE_BINARY_DIR}/itaps/imesh  # for iMesh_protos.h
  )

target_link_libraries(
  FBiGeomMOAB
  MOAB iMesh
  )

set_target_properties(
  FBiGeomMOAB
  PROPERTIES
  VERSION ${MOAB_VERSION}
  SOVERSION ${MOAB_VERSION_MAJOR}
  )
