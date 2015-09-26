
# Generate Fortran name mangling headers
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iRel_protos.h
  COMMAND ${CMAKE_COMMAND}
          -DDEF:STRING=MOAB
          -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/iRel.h
          -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/iRel_protos.h
          -Dprefix:STRING=iRel
          -P ${CMAKE_SOURCE_DIR}/config/ITAPSFortranMangling.cmake
          DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/iRel.h )

set ( MOAB_IREL_SRCS
          ArrayManager.hpp
          AssocPair.cpp
          AssocPair.hpp
          AssocPairSide.hpp
          iRel_Lasso.hpp
          iRel_Lasso.cpp
          Lasso.hpp
          Lasso.cpp
    )
set ( MOAB_IREL_LIB_INCLUDES
          iRel.h
          iRel_f.h
          ${CMAKE_CURRENT_BINARY_DIR}/iRel_protos.h
    )

set(IREL_LINK_LIBS )
if (MOAB_HAVE_CGM)
  include_directories( ${CGM_INCLUDE_DIRS} )
  set(IREL_LINK_LIBS ${IGEOM_LIBS} ${CGM_LIBRARIES})

  set(MOAB_IREL_SRCS ${MOAB_IREL_SRCS}
                    GeomAssocPairSide.cpp
                    GeomAssocPairSide.hpp
      )
endif()

set(IREL_LINK_LIBS MOAB ${IREL_LINK_LIBS})
if (MOAB_HAVE_FBIGEOM)
  include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src/itaps/fbigeom ${CMAKE_CURRENT_BINARY_DIR}/src/itaps/fbigeom )
  set(IREL_LINK_LIBS FBiGeomMOAB ${IREL_LINK_LIBS})
  set(MOAB_IREL_SRCS ${MOAB_IREL_SRCS}
                     FBGeomAssocPairSide.cpp
                     FBGeomAssocPairSide.hpp
      )
endif()

if (MOAB_HAVE_IMESH)
  include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src/itaps/imesh ${CMAKE_CURRENT_BINARY_DIR}/src/itaps/imesh )
  set(IREL_LINK_LIBS iMesh ${IREL_LINK_LIBS})
  set(MOAB_IREL_SRCS ${MOAB_IREL_SRCS}
                     MeshAssocPairSide.cpp
                     MeshAssocPairSide.hpp
    )
endif()

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/itaps
  ${CMAKE_BINARY_DIR}/itaps
  )

if ( MOAB_HAVE_HDF5 AND HDF5_FOUND )
  include_directories( ${HDF5_INCLUDE_DIRECTORIES}
                       ${CMAKE_SOURCE_DIR}/src/io/mhdf/include )
endif( )

set_source_files_properties(
  ${MOAB_IREL_SRCS} COMPILE_FLAGS "${MOAB_DEFINES}"
  )

add_library(
  iRel
  ${MOAB_IREL_SRCS}
  ${MOAB_IREL_LIB_INCLUDES}
  )

target_link_libraries( iRel ${IREL_LINK_LIBS} )

set_target_properties(
  iRel
  PROPERTIES
  VERSION ${MOAB_VERSION}
  SOVERSION ${MOAB_VERSION_MAJOR}
  )

install(
  TARGETS iRel
  EXPORT MOABTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
set( MOAB_INSTALL_TARGETS ${MOAB_INSTALL_TARGETS} iRel)

install(
  FILES ${MOAB_IREL_LIB_INCLUDES}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

