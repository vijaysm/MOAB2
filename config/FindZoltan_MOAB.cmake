#
# Find Zoltan include directories and libraries
#
# ZOLTAN_INCLUDES            - list of include paths to find netcdf.h
# ZOLTAN_LIBRARIES           - list of libraries to link against when using Zoltan
# ZOLTAN_FOUND               - Do not attempt to use Zoltan if "no", "0", or undefined.

set (ZOLTAN_DIR "" CACHE PATH "Path to search for Zoltan header and library files" )
set (ZOLTAN_FOUND NO CACHE INTERNAL "Found Zoltan components successfully." )

# First try to find netCDF with default CMake logic.
# Search twice, <http://cmake.org/cmake/help/v3.0/command/find_package.html>.
find_package(
  Trilinos COMPONENTS Zoltan
  PATHS ${ZOLTAN_DIR}
  NO_DEFAULT_PATH
  )
find_package(Trilinos COMPONENTS Zoltan)

if (${Trilinos_Zoltan_FOUND})
  set(ZOLTAN_FOUND ${Trilinos_Zoltan_FOUND})
  set(ZOLTAN_LIBRARIES ${Zoltan_LIBRARIES})
  set(ZOLTAN_INCLUDES ${Zoltan_INCLUDE_DIRS})
else()
  # try to find zoltan manually
  find_path( ZOLTAN_INCLUDE_DIR zoltan.h
    ${ZOLTAN_DIR}
    ${ZOLTAN_DIR}/include
    /usr/local/include
    /usr/include
    )

  find_library( ZOLTAN_LIBRARY
    NAMES zoltan
    HINTS ${ZOLTAN_DIR}
    ${ZOLTAN_DIR}/lib64
    ${ZOLTAN_DIR}/lib
    /usr/local/lib64
    /usr/lib64
    /usr/lib64/zoltan
    /usr/local/lib
    /usr/lib
    /usr/lib/zoltan
    )


  macro (ZOLTAN_GET_VARIABLE makefile name var)
    set (${var} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
    execute_process (COMMAND ${CMAKE_BUILD_TOOL} -f ${${makefile}} show VARIABLE=${name}
      OUTPUT_VARIABLE ${var}
      RESULT_VARIABLE zoltan_return)
  endmacro (ZOLTAN_GET_VARIABLE)

  macro (ZOLTAN_GET_ALL_VARIABLES)
    if (NOT zoltan_config_current)
      # A temporary makefile to probe this Zoltan components's configuration
      # The current inspection is based on Zoltan-3.6 installation
      set (zoltan_config_makefile "${CMAKE_CURRENT_BINARY_DIR}/Makefile.zoltan")
      file (WRITE ${zoltan_config_makefile}
        "## This file was autogenerated by FindZoltan.cmake
        include ${ZOLTAN_INCLUDE_DIR}/Makefile.export.zoltan
        include ${ZOLTAN_INCLUDE_DIR}/Makefile.export.zoltan.macros
        show :
        -@echo -n \${\${VARIABLE}}"
        )
      ZOLTAN_GET_VARIABLE (zoltan_config_makefile ZOLTAN_CPPFLAGS    zoltan_extra_cppflags)
      ZOLTAN_GET_VARIABLE (zoltan_config_makefile ZOLTAN_EXTRA_LIBS  zoltan_extra_libs)
      ZOLTAN_GET_VARIABLE (zoltan_config_makefile ZOLTAN_LDFLAGS     zoltan_ldflags)

      file (REMOVE ${zoltan_config_makefile})
      SET(tmp_incs "-I${ZOLTAN_INCLUDE_DIR} ${zoltan_extra_cppflags}")
      resolve_includes(ZOLTAN_INCLUDES ${tmp_incs})
      SET(tmp_libs "${ZOLTAN_LIBRARY} ${zoltan_ldflags} ${zoltan_extra_libs}")
      resolve_libraries (ZOLTAN_LIBRARIES "${tmp_libs}")
    endif ()
  endmacro (ZOLTAN_GET_ALL_VARIABLES)

  IF (NOT ZOLTAN_FOUND)
    if ( ZOLTAN_INCLUDE_DIR AND ZOLTAN_LIBRARY )
      set( ZOLTAN_FOUND YES )
      if(EXISTS ${ZOLTAN_INCLUDE_DIR}/Makefile.export.zoltan)
        include (ResolveCompilerPaths)
        ZOLTAN_GET_ALL_VARIABLES()
      else(EXISTS ${ZOLTAN_INCLUDE_DIR}/Makefile.export.zoltan)
        SET(ZOLTAN_INCLUDES ${ZOLTAN_INCLUDE_DIR})
        SET(ZOLTAN_LIBRARIES ${ZOLTAN_LIBRARY})
      endif(EXISTS ${ZOLTAN_INCLUDE_DIR}/Makefile.export.zoltan)
    else ( ZOLTAN_INCLUDE_DIR AND ZOLTAN_LIBRARY )
      set( ZOLTAN_FOUND NO )
      message("finding Zoltan failed, please try to set the var ZOLTAN_DIR")
    endif ( ZOLTAN_INCLUDE_DIR AND ZOLTAN_LIBRARY )
  ENDIF (NOT ZOLTAN_FOUND)

endif()

mark_as_advanced(
  ZOLTAN_DIR
  ZOLTAN_INCLUDES
  ZOLTAN_LIBRARIES
  )

message (STATUS "---   Zoltan Configuration ::")
message (STATUS "        INCLUDES  : ${ZOLTAN_INCLUDES}")
message (STATUS "        LIBRARIES : ${ZOLTAN_LIBRARIES}")

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (
  Zoltan "Zoltan not found, check environment variables ZOLTAN_DIR"
  ZOLTAN_INCLUDES
  ZOLTAN_LIBRARIES
  )
