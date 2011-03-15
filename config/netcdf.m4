#######################################################################################
# Check values of NC_MAX_DIMS and NC_MAX_VARS for ExodusII compatability.
# 
# Arguments are: 
#  1) required value for NC_MAX_DIMS, may be emtpy
#  2) required value for NC_MAX_VARS, may be emtpy
#  3) name of header in which to check for values
#  4) variable in which to store result (set to yes if
#        limits are sufficient, and on otherwise.)
#######################################################################################
AC_DEFUN([FATHOM_CHECK_NETCDF_LIMITS],[
  MIN_NC_MAX_DIMS="$1"
  MIN_NC_MAX_VARS="$2"
  $4=yes

  if test "x" != "x$MIN_NC_MAX_DIMS"; then
    AC_MSG_CHECKING([if NC_MAX_DIMS is at least ${MIN_NC_MAX_DIMS}])
    AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM([#include <$3>],
                       [[int arr[1 + (int)(NC_MAX_DIMS) - (int)(${MIN_NC_MAX_DIMS})];]])],
      [AC_MSG_RESULT([yes])],
      [AC_MSG_RESULT([no]); $4=no])
  fi
  if test "x" != "x$MIN_NC_MAX_VARS"; then
    AC_MSG_CHECKING([if NC_MAX_VARS is at least ${MIN_NC_MAX_VARS}])
    AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM([#include <netcdf.h>],
                       [[int arr[1 + (int)(NC_MAX_VARS) - (int)(${MIN_NC_MAX_VARS})];]])],
      [AC_MSG_RESULT([yes])],
      [AC_MSG_RESULT([no]); $4=no])
  fi
])

#######################################################################################
# Check for NetCDF library
# Sets HAVE_NETCDF to 'yes' or 'no'
# If HAVE_NETCDF == yes, then exports:
#   NETCDF_CPPFLAGS
#   NETCDF_LDFLAGS
#   NETCDF_LIBS
#   NETCDF_SUFFICIENT_DIMS_VARS
#
# This macro has two optional arguments:  a minimum value for
# NC_MAX_DIMS and a minimum value for NC_MAX_VARS.  If either or
# both of these are specified, NETCDF_SUFFICIENT_DIMS_VARS will
# be set to yes if the NetCDF library is built with limits that
# are at least the passed miminums.  It will be set to no if 
# either limit is less than the passed minimum.
#######################################################################################
AC_DEFUN([FATHOM_CHECK_NETCDF],[

AC_MSG_CHECKING([if NetCDF support is enabled])
AC_ARG_WITH(netcdf, 
[AC_HELP_STRING([--with-netcdf@<:@=DIR@:>@], [Specify NetCDF library to use for ExodusII file format])
AC_HELP_STRING([--without-netcdf], [Disable support for ExodusII file format])],
[NETCDF_ARG=$withval
DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-netcdf=\"${withval}\""
]
, [NETCDF_ARG=])
if test "xno" != "x$NETCDF_ARG"; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

 # if NetCDF support is not disabled
HAVE_NETCDF=no
if test "xno" != "x$NETCDF_ARG"; then
  HAVE_NETCDF=yes
  
    # if a path is specified, update LIBS and INCLUDES accordingly
  if test "xyes" != "x$NETCDF_ARG" && test "x" != "x$NETCDF_ARG"; then
    if test -d "${NETCDF_ARG}/lib"; then
      NETCDF_LDFLAGS="-L${NETCDF_ARG}/lib"
    elif test -d "${NETCDF_ARG}"; then
      NETCDF_LDFLAGS="-L${NETCDF_ARG}"
    else
      AC_MSG_ERROR("$NETCDF_ARG is not a directory.")
    fi
    if test -d "${NETCDF_ARG}/include"; then
      NETCDF_CPPFLAGS="-I${NETCDF_ARG}/include"
    elif test -d "${NETCDF_ARG}/inc"; then
      NETCDF_CPPFLAGS="-I${NETCDF_ARG}/inc"
    else
      NETCDF_CPPFLAGS="-I${NETCDF_ARG}"
    fi
  fi
  
  old_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$NETCDF_CPPFLAGS $CPPFLAGS"
  old_LDFLAGS="$LDFLAGS"
  LDFLAGS="$NETCDF_LDFLAGS $HDF5_LDFLAGS $LDFLAGS"
  
   # Check for C library
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS( [netcdf.h], 
                    [FATHOM_CHECK_NETCDF_LIMITS([$1],[$2],[netcdf.h],[NETCDF_SUFFICIENT_DIM_VARS])], 
                    [AC_MSG_WARN([[NetCDF header not found.]]); HAVE_NETCDF=no] )
  
      # Check if netcdf is usable by itself
  AC_CHECK_LIB( [netcdf], [nc_create], [NETCDF_LIBS="-lnetcdf"], [
      # Check if netcdf is usable with HDF5
    unset ac_cv_lib_netcdf
    unset ac_cv_lib_netcdf_nc_create
      # If we haven't already looked for HDF5 libraries, again now incase
      # they're in the NetCDF lib directory.
    FATHOM_DETECT_HDF5_LIBS
    LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
    AC_CHECK_LIB( [netcdf], [nc_create], [NETCDF_LIBS="-lnetcdf -lhdf5_hl"], [
      # Try one more time with HDF5 and libcurl
      unset ac_cv_lib_netcdf
      unset ac_cv_lib_netcdf_nc_create
      AC_CHECK_LIB( [netcdf], [nc_create], [NETCDF_LIBS="-lnetcdf -lhdf5_hl -lcurl"], 
        [HAVE_NETCDF=no], [-lhdf5_hl $HDF5_LIBS -lcurl] )],
      [-lhdf5_hl $HDF5_LIBS] )],
    )
  
  CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$old_LDFLAGS"
  AC_LANG_POP([C])

  if test "x$HAVE_NETCDF" = "xno"; then
    if test "x$NETCDF_ARG" != "x"; then 
      AC_MSG_ERROR("NetCDF not found or not working")
    else
      AC_MSG_WARN("NetCDF support disabled")
    fi
    NETCDF_CPPFLAGS=
    NETCDF_LDFLAGS=
  fi
fi

]) # FATHOM_HAVE_NETCDF

AC_DEFUN([FATHOM_CHECK_PNETCDF],[

AC_ARG_WITH(pnetcdf, 
[AC_HELP_STRING([--with-pnetcdf@<:@=DIR@:>@], [Specify Pnetcdf library to use])
AC_HELP_STRING([--without-pnetcdf], [Disable support for Pnetcdf-based file formats])],
[PNETCDF_ARG=$withval
DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-pnetcdf=\"${withval}\""
]
, [PNETCDF_ARG=])

  # PNETCDF requires MPI too
if test "xyes" != "x$WITH_MPI"; then
  if test "x" == "x$PNETCDF_ARG"; then
    PNETCDF_ARG=no
  elif test "xno" != "xPNETCDF_ARG"; then
    AC_MSG_ERROR([Pnetcdf requires --with-mpi])
  fi
  HAVE_PNETCDF=no
fi

AC_MSG_CHECKING([if Pnetcdf support is enabled])
if test "xno" != "x$PNETCDF_ARG"; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

 # if Pnetcdf support is not disabled
HAVE_PNETCDF=no
if test "xno" != "x$PNETCDF_ARG"; then
  HAVE_PNETCDF=yes
  
    # if a path is specified, update LIBS and INCLUDES accordingly
  if test "xyes" != "x$PNETCDF_ARG" && test "x" != "x$PNETCDF_ARG"; then
    if test -d "${PNETCDF_ARG}/lib"; then
      PNETCDF_LDFLAGS="-L${PNETCDF_ARG}/lib"
    elif test -d "${PNETCDF_ARG}"; then
      PNETCDF_LDFLAGS="-L${PNETCDF_ARG}"
    else
      AC_MSG_ERROR("$PNETCDF_ARG is not a directory.")
    fi
    if test -d "${PNETCDF_ARG}/include"; then
      PNETCDF_CPPFLAGS="-I${PNETCDF_ARG}/include"
    elif test -d "${PNETCDF_ARG}/inc"; then
      PNETCDF_CPPFLAGS="-I${PNETCDF_ARG}/inc"
    else
      PNETCDF_CPPFLAGS="-I${PNETCDF_ARG}"
    fi
  fi
  
  old_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$PNETCDF_CPPFLAGS $CPPFLAGS"
  old_LDFLAGS="$LDFLAGS"
  LDFLAGS="$PNETCDF_LDFLAGS $LDFLAGS"
  
   # Check for C library
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS( [pnetcdf.h], 
                    [FATHOM_CHECK_NETCDF_LIMITS([$1],[$2],[pnetcdf.h],[PNETCDF_SUFFICIENT_DIM_VARS])], 
                    [HAVE_PNETCDF=no] )

  AC_CHECK_LIB( [pnetcdf], [ncmpi_create], [PNETCDF_LIBS="-lpnetcdf"], [HAVE_PNETCDF=no] )

  
  AC_LANG_POP([C])
  CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$old_LDFLAGS"

  if test "x$HAVE_PNETCDF" = "xno"; then
    if test "x$PNETCDF_ARG" != "x"; then 
      AC_MSG_ERROR("Pnetcdf not found or not working")
    fi
    PNETCDF_CPPFLAGS=
    PNETCDF_LDFLAGS=
  fi
fi

]) # FATHOM_HAVE_PNETCDF
