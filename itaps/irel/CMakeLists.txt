project(Lasso CXX)

if (MOAB_HAS_MPI)
  find_package(MPI REQUIRED)
  add_definitions(-DUSE_MPI)
endif ()

if (CGM_FOUND)
  add_definitions(-DENABLE_IGEOM)
endif ()

if (ENABLE_IGEOM)
  add_definitions(-DENABLE_FBIGEOM)
endif ()
if (ENABLE_IMESH)
  add_definitions(-DENABLE_IMESH)
endif ()

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/iBase_FCDefs.h.in" fcdefs)
string(REPLACE "FC_FUNC" "IREL_FC_FUNC" fcdefs "${fcdefs}")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/iRel_FCDefs.h" "${fcdefs}")

add_custom_command(
  OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/iRel_protos.h"
  COMMAND "${CMAKE_COMMAND}"
          "-Doutput=${CMAKE_CURRENT_BINARY_DIR}/iRel_protos.h"
          "-Dinput=${CMAKE_CURRENT_SOURCE_DIR}/iRel.h"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/generate_irel_protos.cmake"
  DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/iRel.h"
          "${CMAKE_CURRENT_SOURCE_DIR}/generate_irel_protos.cmake"
  COMMENT "Generating ${CMAKE_CURRENT_BINARY_DIR}/itaps/irel/iRel_protos.h")

set(lasso_srcs
  AssocPair.cpp
  iRel_Lasso.cpp
  Lasso.cpp)

set(lasso_includes
  ArrayManager.hpp
  AssocPair.hpp
  AssocPairSide.hpp
  iRel_Lasso.hpp
  Lasso.hpp

  iRel.h
  iRel_f.h
  "${CMAKE_CURRENT_BINARY_DIR}/iRel_protos.h"
  "${CMAKE_CURRENT_BINARY_DIR}/iRel_FCDefs.h")

set(link_libs)

if (CGM_FOUND)
  list(APPEND lasso_srcs
    GeomAssocPairSide.cpp)
  list(APPEND lasso_includes
    GeomAssocPairSide.hpp)
  list(APPEND link_libs
    iGeom)
endif ()

if (ENABLE_FBIGEOM)
  list(APPEND lasso_srcs
    FBGeomAssocPairSide.cpp)
  list(APPEND lasso_includes
    FBGeomAssocPairSide.hpp)
  list(APPEND link_libs
    FBiGeomMOAB)
endif ()

if (ENABLE_IMESH)
  list(APPEND lasso_srcs
    MeshAssocPairSide.cpp)
  list(APPEND lasso_includes
    MeshAssocPairSide.hpp)
  list(APPEND link_libs
    iMesh)
endif ()

include_directories(
  ${CGM_INCLUDE_DIRS}
  ${MOAB_INCLUDE_DIRS}
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_BINARY_DIR}")

add_library(iRel
  ${lasso_srcs}
  ${lasso_includes})
target_link_libraries(iRel
  ${link_libs})
target_include_directories(iRel
  PUBLIC
    "$<INSTALL_INTERFACE:include>")
include(GenerateExportHeader)
generate_export_header(iRel)
list(APPEND lasso_includes
  "${CMAKE_CURRENT_BINARY_DIR}/irel_export.h")

set_target_properties(iRel
  PROPERTIES
    VERSION   ${MOAB_VERSION}
    SOVERSION ${MOAB_VERSION_MAJOR}
  )

install(
  TARGETS   iRel
  EXPORT    MOABTargets
  RUNTIME   DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY   DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE   DESTINATION ${CMAKE_INSTALL_LIBDIR}
  COMPONENT Runtime)
install(
  FILES       ${lasso_includes}
  DESTINATION include
  COMPONENT   Development)
