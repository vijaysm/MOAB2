#
# Find the native HDF5 includes and library
#
# HDF5_INCLUDES    - where to find hdf5.h, H5public.h, etc.
# HDF5_LIBRARIES   - List of fully qualified libraries to link against when using hdf5.
# HDF5_FOUND       - Do not attempt to use hdf5 if "no" or undefined.

set( HDF5_DIR "" CACHE PATH "Path to search for HDF5 header and library files" )
set (HDF5_FOUND NO CACHE INTERNAL "Found HDF5 components successfully." )

if(EXISTS "${HDF5_DIR}/share/cmake/hdf5/hdf5-config.cmake")
  include(${HDF5_DIR}/share/cmake/hdf5/hdf5-config.cmake)
else()

FIND_PATH(HDF5_INCLUDE_DIR
  NAMES hdf5.h H5public.h
  PATHS ${HDF5_DIR}/include
  /usr/local/include
  /usr/include
  /opt/local/include
)

foreach (VARIANT dl m z )
  FIND_LIBRARY(hdf5_deplibs_${VARIANT} ${VARIANT}
    PATHS /lib /usr/local/lib /usr/lib /opt/local/lib
  )
  list(APPEND HDF5_DEP_LIBRARIES ${hdf5_deplibs_${VARIANT}})
endforeach()

FIND_LIBRARY(HDF5_BASE_LIBRARY hdf5 hdf5d)

FIND_LIBRARY(HDF5_BASE_LIBRARY NAMES hdf5 hdf5d
  PATHS ${HDF5_DIR}/lib /usr/local/lib /usr/lib /opt/local/lib
)
FIND_LIBRARY(HDF5_HLBASE_LIBRARY hdf5_hl hdf5_hld
  PATHS ${HDF5_DIR}/lib /usr/local/lib /usr/lib /opt/local/lib
)

IF (NOT HDF5_FOUND)
  IF (HDF5_INCLUDE_DIR AND HDF5_BASE_LIBRARY)
    FIND_LIBRARY(HDF5_CXX_LIBRARY hdf5_cxx
      PATHS ${HDF5_DIR}/lib /usr/local/lib /usr/lib /opt/local/lib
    )
    FIND_LIBRARY(HDF5_HLCXX_LIBRARY hdf5_hl_cxx
      PATHS ${HDF5_DIR}/lib /usr/local/lib /usr/lib /opt/local/lib
    )
    FIND_LIBRARY(HDF5_FORT_LIBRARY hdf5_fortran
      PATHS ${HDF5_DIR}/lib /usr/local/lib /usr/lib /opt/local/lib
    )
    FIND_LIBRARY(HDF5_HLFORT_LIBRARY
      NAMES hdf5hl_fortran hdf5_hl_fortran
      PATHS ${HDF5_DIR}/lib /usr/local/lib /usr/lib /opt/local/lib
    )
    SET( HDF5_INCLUDES "${HDF5_INCLUDE_DIR}" )
    if (HDF5_FORT_LIBRARY)
      FIND_PATH(HDF5_FORT_INCLUDE_DIR
        NAMES hdf5.mod
        PATHS ${HDF5_DIR}/include
        ${HDF5_DIR}/include/fortran
        /usr/local/include
        /usr/include
        /opt/local/include
      )
      if (HDF5_FORT_INCLUDE_DIR AND NOT ${HDF5_FORT_INCLUDE_DIR} STREQUAL ${HDF5_INCLUDE_DIR})
        SET( HDF5_INCLUDES "${HDF5_INCLUDES} ${HDF5_FORT_INCLUDE_DIR}" )
      endif (HDF5_FORT_INCLUDE_DIR AND NOT ${HDF5_FORT_INCLUDE_DIR} STREQUAL ${HDF5_INCLUDE_DIR})
      unset(HDF5_FORT_INCLUDE_DIR CACHE)
    endif (HDF5_FORT_LIBRARY)
    # Add the libraries based on availability
    foreach (VARIANT CXX FORT BASE )
      if (HDF5_HL${VARIANT}_LIBRARY)
        list(APPEND HDF5_LIBRARIES ${HDF5_HL${VARIANT}_LIBRARY})
      endif (HDF5_HL${VARIANT}_LIBRARY)
      if (HDF5_${VARIANT}_LIBRARY)
        list(APPEND HDF5_LIBRARIES ${HDF5_${VARIANT}_LIBRARY})
      endif (HDF5_${VARIANT}_LIBRARY)
      unset(HDF5_HL${VARIANT}_LIBRARY CACHE)
      unset(HDF5_${VARIANT}_LIBRARY CACHE)
    endforeach()
    list(APPEND HDF5_LIBRARIES ${HDF5_DEP_LIBRARIES})
    SET( HDF5_FOUND YES )
    message (STATUS "---   HDF5 Configuration ::")
    message (STATUS "        INCLUDES  : ${HDF5_INCLUDES}")
    message (STATUS "        LIBRARIES : ${HDF5_LIBRARIES}")
  ELSE (HDF5_INCLUDE_DIR AND HDF5_BASE_LIBRARY)
    set( HDF5_FOUND NO )
    message("finding HDF5 failed, please try to set the var HDF5_DIR")
  ENDIF(HDF5_INCLUDE_DIR AND HDF5_BASE_LIBRARY)
ENDIF (NOT HDF5_FOUND)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (HDF5 "HDF5 not found, check environment variables HDF5_DIR"
  HDF5_DIR HDF5_INCLUDES HDF5_LIBRARIES)

#now we create fake targets to be used
if(EXISTS ${HDF5_DIR}/share/cmake/hdf5/hdf5-targets.cmake)
  include(${HDF5_DIR}/share/cmake/hdf5/hdf5-targets.cmake)
endif()

endif(EXISTS "${HDF5_DIR}/share/cmake/hdf5/hdf5-config.cmake")
