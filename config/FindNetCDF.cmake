#
# Find NetCDF include directories and libraries
#
# NetCDF_INCLUDES            - list of include paths to find netcdf.h
# NetCDF_LIBRARIES           - list of libraries to link against when using NetCDF
# NetCDF_FOUND               - Do not attempt to use NetCDF if "no", "0", or undefined.

set (NetCDF_DIR "" CACHE PATH "Path to search for NetCDF header and library files" )
set (NetCDF_FOUND NO CACHE INTERNAL "Found NetCDF components successfully." )

if (EXISTS "${NetCDF_DIR}/lib/pkgconfig/netcdf.pc")
  execute_process(COMMAND pkg-config --cflags-only-I "${NetCDF_DIR}/lib/pkgconfig/netcdf.pc" OUTPUT_VARIABLE NETCDF_INCLUDES)
  execute_process(COMMAND pkg-config --libs-only-L "${NetCDF_DIR}/lib/pkgconfig/netcdf.pc" OUTPUT_VARIABLE NETCDF_LDLAGS)
  execute_process(COMMAND pkg-config --libs-only-l "${NetCDF_DIR}/lib/pkgconfig/netcdf.pc" OUTPUT_VARIABLE NETCDF_LIBRARY)
  execute_process(COMMAND pkg-config --libs "${NetCDF_DIR}/lib/pkgconfig/netcdf.pc" OUTPUT_VARIABLE NETCDF_LIBS)

  IF (NOT NetCDF_FOUND)
    if ( NETCDF_INCLUDES AND NETCDF_LIBS )
      string(STRIP "${NETCDF_INCLUDES}" NETCDF_INCLUDES)
      string(STRIP "${NETCDF_LIBS}" NETCDF_LIBS)
      include (ResolveCompilerPaths)
      RESOLVE_INCLUDES(NetCDF_INCLUDES "${NETCDF_INCLUDES}")
      #RESOLVE_LIBRARIES(NetCDF_LIBRARIES "${NETCDF_LIBS}")
      SET(NetCDF_LIBRARIES ${NETCDF_LIBS})
      SET(NetCDF_FOUND YES )
    else ( NETCDF_INCLUDES AND NETCDF_LIBS )
      set( NetCDF_FOUND NO )
    endif ( NETCDF_INCLUDES AND NETCDF_LIBS )
  ENDIF (NOT NetCDF_FOUND)

else (EXISTS "${NetCDF_DIR}/lib/pkgconfig/netcdf.pc")

  find_path( NetCDF_INCLUDE_DIR netcdf.h
    HINTS ${NetCDF_DIR}
    ${NetCDF_DIR}/include
    ENV CPLUS_INCLUDE_PATH
    NO_DEFAULT_PATH
  )

  find_library( NetCDF_C_LIBRARY
    NAMES libnetcdf.a netcdf
    HINTS ${NetCDF_DIR}
    ${NetCDF_DIR}/lib64
    ${NetCDF_DIR}/lib
    NO_DEFAULT_PATH
  )

  find_library( NetCDF_CXX_LIBRARY
    NAMES netcdf_c++
    HINTS ${NetCDF_DIR}
    ${NetCDF_DIR}/lib64
    ${NetCDF_DIR}/lib
    NO_DEFAULT_PATH
  )

  find_library( NetCDF_FORTRAN_LIBRARY
    NAMES netcdf_g77 netcdff netcdf_ifc netcdf_x86_64
    HINTS ${NetCDF_DIR}
    ${NetCDF_DIR}/lib64
    ${NetCDF_DIR}/lib
    NO_DEFAULT_PATH
  )

  IF (NOT NetCDF_FOUND)
    if ( NetCDF_INCLUDE_DIR AND NetCDF_C_LIBRARY )
      set(NetCDF_FOUND YES )
      set(NetCDF_INCLUDES "${NetCDF_INCLUDE_DIR}")
      set(NetCDF_LIBRARIES ${NetCDF_C_LIBRARY})
      if ( NetCDF_CXX_LIBRARY )
        set(NetCDF_LIBRARIES ${NetCDF_LIBRARIES} ${NetCDF_CXX_LIBRARY})
      endif ( NetCDF_CXX_LIBRARY )
      if ( NetCDF_FORTRAN_LIBRARY )
        set(NetCDF_LIBRARIES ${NetCDF_LIBRARIES} ${NetCDF_FORTRAN_LIBRARY})
      endif ( NetCDF_FORTRAN_LIBRARY )
    else ( NetCDF_INCLUDE_DIR AND NetCDF_C_LIBRARY )
      set( NetCDF_FOUND NO )
      message("finding NetCDF failed, please try to set the var NetCDF_DIR")
    endif ( NetCDF_INCLUDE_DIR AND NetCDF_C_LIBRARY )
  ENDIF (NOT NetCDF_FOUND)

endif (EXISTS "${NetCDF_DIR}/lib/pkgconfig/netcdf.pc")

message (STATUS "---   NetCDF Configuration ::")
message (STATUS "        INCLUDES  : ${NetCDF_INCLUDES}")
message (STATUS "        LIBRARIES : ${NetCDF_LIBRARIES}")

mark_as_advanced(
  NetCDF_DIR
  NetCDF_INCLUDES
  NetCDF_LIBRARIES
)

IF (MOAB_USE_MPI)
  set (PNetCDF_DIR "" CACHE PATH "Path to search for parallel NetCDF header and library files" )
  set (PNetCDF_FOUND NO CACHE INTERNAL "Found parallel NetCDF components successfully." )

  find_path( PNetCDF_INCLUDES pnetcdf.h
    ${PNetCDF_DIR}
    ${PNetCDF_DIR}/include
    ENV CPLUS_INCLUDE_PATH
    NO_DEFAULT_PATH
  )

  find_library( PNetCDF_LIBRARIES
    NAMES pnetcdf
    HINTS ${PNetCDF_DIR}
    ${PNetCDF_DIR}/lib64
    ${PNetCDF_DIR}/lib
    NO_DEFAULT_PATH
  )

  IF (NOT PNetCDF_FOUND)
    if ( PNetCDF_INCLUDES AND PNetCDF_LIBRARIES )
      set( PNetCDF_FOUND YES )
      message (STATUS "---   PNetCDF Configuration ::")
      message (STATUS "        INCLUDES  : ${PNetCDF_INCLUDES}")
      message (STATUS "        LIBRARIES : ${PNetCDF_LIBRARIES}")
    else ( PNetCDF_INCLUDES AND PNetCDF_LIBRARIES )
      set( NetCDF_FOUND NO )
      message("finding PNetCDF failed, please try to set the var PNetCDF_DIR")
    endif ( PNetCDF_INCLUDES AND PNetCDF_LIBRARIES )
  ENDIF (NOT PNetCDF_FOUND)

  mark_as_advanced(
    PNetCDF_DIR
    PNetCDF_INCLUDES
    PNetCDF_LIBRARIES
  )
ELSE (MOAB_USE_MPI)
  message (STATUS "Not configuring with PNetCDF since MPI installation not specified")
ENDIF()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NetCDF "NetCDF not found, check environment variables NetCDF_DIR"
  NetCDF_DIR NetCDF_INCLUDES NetCDF_LIBRARIES)
find_package_handle_standard_args (PNetCDF "PNetCDF not found, check environment variables PNetCDF_DIR"
  PNetCDF_DIR PNetCDF_INCLUDES PNetCDF_LIBRARIES)
