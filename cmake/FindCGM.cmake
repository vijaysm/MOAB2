# Find the common geometry module libraries

if ( NOT CGM_SOURCE_DIR )
  find_path( CGM_SOURCE_DIR
    NAMES geom/CubitGeomConfigure.h.in
    PATHS
      ${CGM_INCLUDE_DIRECTORIES}
  )
endif ( NOT CGM_SOURCE_DIR )

if ( NOT CGM_BINARY_DIR )
  if ( CGM_LIBRARIES )
    set( CGM_BUILD_SEARCH "" )
    foreach ( LIB ${CGM_LIBRARIES} )
      get_filename_component( PATH DIR ${LIB} )
      set( CGM_BUILD_SEARCH "${CGM_BUILD_SEARCH};${DIR}/.." )
    endforeach ( DIR )
    set( DIR )
  endif ( CGM_LIBRARIES )
  find_path( CGM_BINARY_DIR
    NAMES geom/CubitGeomConfigure.h
    PATHS
      ${CGM_BUILD_SEARCH}
  )
endif ( NOT CGM_BINARY_DIR )

if ( NOT CGM_LIBRARIES )
  set( CGM_LIBSEARCH
    ${CGM_BINARY_DIR}/init
    ${CGM_BINARY_DIR}/geom/virtual
    ${CGM_BINARY_DIR}/geom/facetbool
    ${CGM_BINARY_DIR}/geom/facet
    ${CGM_BINARY_DIR}/geom/Cholla
    ${CGM_BINARY_DIR}/geom
    ${CGM_BINARY_DIR}/util
    /usr/local/lib64
    /usr/local/bin64
    /usr/lib64
    /usr/bin64
    /usr/local/lib
    /usr/local/bin
    /usr/lib
    /usr/bin
  )
  find_library( CGM_UTIL_LIBRARY
    NAMES cubit_util
    PATHS ${CGM_LIBSEARCH}
  )

  find_library( CGM_GEOM_LIBRARY
    NAMES cubit_geom
    PATHS ${CGM_LIBSEARCH}
  )

  find_library( CGM_CHOLLA_LIBRARY
    NAMES cholla
    PATHS ${CGM_LIBSEARCH}
  )

  find_library( CGM_FACET_LIBRARY
    NAMES cubit_facet
    PATHS ${CGM_LIBSEARCH}
  )

  find_library( CGM_FACETBOOL_LIBRARY
    NAMES cubit_facetboolstub
    PATHS ${CGM_LIBSEARCH}
  )

  find_library( CGM_VIRTUAL_LIBRARY
    NAMES cubit_virtual
    PATHS ${CGM_LIBSEARCH}
  )

  find_library( CGM_INIT_LIBRARY
    NAMES cgma_init
    PATHS ${CGM_LIBSEARCH}
  )
endif ( NOT CGM_LIBRARIES )

if ( CGM_SOURCE_DIR AND CGM_BINARY_DIR AND NOT CGM_INCLUDE_DIRECTORIES )
  set( cubit_geom_SOURCE_DIR "${CGM_SOURCE_DIR}/geom" )
  set( cubit_geom_BINARY_DIR "${CGM_BINARY_DIR}/geom" )
  set( cubit_util_SOURCE_DIR "${CGM_SOURCE_DIR}/util" )
  set( cubit_util_BINARY_DIR "${CGM_BINARY_DIR}/util" )
  set( cgma_init_SOURCE_DIR "${CGM_SOURCE_DIR}/init" )
  set( cgma_init_BINARY_DIR "${CGM_BINARY_DIR}/init" )
  set( CGM_INCLUDE_DIRECTORIES
    "${cubit_util_SOURCE_DIR}"
    "${cubit_util_BINARY_DIR}"
    "${cubit_geom_SOURCE_DIR}"
    "${cubit_geom_SOURCE_DIR}/ACIS"
    "${cubit_geom_SOURCE_DIR}/Cholla"
    "${cubit_geom_SOURCE_DIR}/facet"
    "${cubit_geom_SOURCE_DIR}/facetbool"
    "${cubit_geom_SOURCE_DIR}/OCC"
    "${cubit_geom_SOURCE_DIR}/parallel"
    "${cubit_geom_SOURCE_DIR}/SolidWorks"
    "${cubit_geom_SOURCE_DIR}/virtual"
    "${cubit_geom_BINARY_DIR}"
    "${cgma_init_SOURCE_DIR}"
    "${cgma_init_BINARY_DIR}"
    CACHE PATH "Paths to util/CubitUtilConfigure.h.in, util/CubitUtilConfigure.h, geom/CubitGeomConfigure.h.in, AND geom/CubitGeomConfigure.h files." FORCE
  )
endif ( CGM_SOURCE_DIR AND CGM_BINARY_DIR AND NOT CGM_INCLUDE_DIRECTORIES )

if ( CGM_UTIL_LIBRARY AND CGM_GEOM_LIBRARY AND CGM_INIT_LIBRARY AND NOT CGM_LIBRARIES )
  # Link to libdl in case shared libraries are used.
  set( CGM_LIBRARIES
    "${CGM_INIT_LIBRARY}"
    "${CGM_CHOLLA_LIBRARY}"
    "${CGM_FACETBOOL_LIBRARY}"
    "${CGM_FACET_LIBRARY}"
    "${CGM_VIRTUAL_LIBRARY}"
    "${CGM_GEOM_LIBRARY}"
    "${CGM_UTIL_LIBRARY}"
    dl
    CACHE FILEPATH "Full paths to cubit_util AND cubit_geom libraries." FORCE
  )
endif ( CGM_UTIL_LIBRARY AND CGM_GEOM_LIBRARY AND CGM_INIT_LIBRARY AND NOT CGM_LIBRARIES )

mark_as_advanced(
  CGM_SOURCE_DIR
  CGM_BINARY_DIR
  CGM_INCLUDE_DIRECTORIES
  CGM_LIBRARIES
  CGM_UTIL_LIBRARY
  CGM_GEOM_LIBRARY
  CGM_INIT_LIBRARY
  CGM_CHOLLA_LIBRARY
  CGM_FACETBOOL_LIBRARY
  CGM_FACET_LIBRARY
  CGM_VIRTUAL_LIBRARY
)

if ( NOT CGM_INCLUDE_DIRECTORIES OR NOT CGM_LIBRARIES )
  set( CGM_INCLUDE_DIRECTORIES "" CACHE PATH "Paths to util/CubitUtilConfigure.h.in, util/CubitUtilConfigure.h, geom/CubitGeomConfigure.h.in, AND geom/CubitGeomConfigure.h files." )
  set( CGM_LIBRARIES           "" CACHE PATH "Paths to cubit_util AND cubit_geom libraries." )
  set( CGM_FOUND 0 )
  if ( CGM_FIND_REQUIRED )
    message( FATAL_ERROR "CGM is required. Please set CGM_INCLUDE_DIRECTORIES and CGM_LIBRARIES or set CGM_SOURCE_DIR and CGM_BINARY_DIR" )
  endif ( CGM_FIND_REQUIRED )
else ( NOT CGM_INCLUDE_DIRECTORIES OR NOT CGM_LIBRARIES )
  set( CGM_FOUND 1 )
endif ( NOT CGM_INCLUDE_DIRECTORIES OR NOT CGM_LIBRARIES )

