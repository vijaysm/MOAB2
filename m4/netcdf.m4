#######################################################################################
# Check for NetCDF library ((C++)
# Sets HAVE_NETCDF to 'yes' or 'no'
# If HAVE_NETCDF == yes, then exports:
#   NETCDF_CPPFLAGS
#   NETCDF_LDFLAGS
#   NETCDF_LIBS
#######################################################################################
AC_DEFUN([FATHOM_CHECK_NETCDF],[

AC_MSG_CHECKING([if NetCDF support is enabled])
AC_ARG_WITH(netcdf, 
[AC_HELP_STRING([--with-netcdf=DIR], [Specify NetCDF library to use for ExodusII file format])
AC_HELP_STRING([--without-netcdf], [Disable support for ExodusII file format])],
[NETCDF_ARG=$withval
DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-netcdf=\"${withval}\""
]
, [NETCDF_ARG=])
if test "xno" = "x$NETCDF_ARG"; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([yes])
fi

 # if NetCDF support is not disabled
HAVE_NETCDF=no
if test "xno" != "x$NETCDF_ARG"; then
  HAVE_NETCDF=yes

    # Check for stream headers and set STRSTREAM_H_SPEC accordingly
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_CHECK_HEADER( [strstream.h], [NETCDF_DEF="<strstream.h>"], [
    AC_CHECK_HEADER( [sstream.h], [NETCDF_DEF="<sstream.h>"], [
      AC_CHECK_HEADER( [strstream], [NETCDF_DEF="<strstream>"], [
        AC_CHECK_HEADER( [sstream], [NETCDF_DEF="<sstream>"] )
  ] ) ] ) ] )
  AC_LANG_RESTORE
  
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
  LDFLAGS="$NETCDF_LDFLAGS $LDFLAGS"
  
   # Check for C library
  AC_CHECK_HEADERS( [netcdf.h], [], [AC_MSG_WARN([[NetCDF header not found.]]); HAVE_NETCDF=no] )
  AC_MSG_CHECKING([for netcdf.hh])
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  HAVE_NETCDF_HH=no
  AC_TRY_COMPILE( 
[#include "netcdf.hh"], [], [HAVE_NETCDF_HH=yes; NETCDF_DEF=], [
    AC_TRY_COMPILE( 
[#define STRSTREAM_H_SPEC $NETCDF_DEF
 #include "netcdf.hh"], [], [HAVE_NETCDF_HH=yes], [NAVE_NETCDF_HH=no])])
  AC_MSG_RESULT([$HAVE_NETCDF_HH])
  if test $HAVE_NETCDF_HH != yes; then
    AC_MSG_WARN([NetCDF C++ header not found])
    HAVE_NETCDF=no
  fi
  if test "x$NETCDF_DEF" != "x"; then
    NETCDF_CPPFLAGS="$NETCDF_CPPFLAGS -DSTRSTREAM_H_SPEC=$NETCDF_DEF"
    CPPFLAGS="$CPPFLAGS -DSTRSTREAM_H_SPEC=$NETCDF_DEF"
  fi
  AC_MSG_CHECKING([[for netcdf_c++ library]])
  old_LIBS="$LIBS"
  LIBS="$LIBS -lnetcdf_c++ -lnetcdf"
  AC_TRY_LINK(
    [#include <netcdf.hh>], [NcFile ncf("foo",NcFile::ReadOnly);],
    [AC_MSG_RESULT([yes]); NETCDF_LIBS="-lnetcdf_c++ -lnetcdf"], 
    [AC_MSG_RESULT([no]);
     AC_MSG_CHECKING([for netcdf_c++ library requiring HDF5-high-level])
     LIBS="$LIBS -lhdf5_hl $HDF5_LIBS"
     AC_TRY_LINK(
           [#include <netcdf.hh>], [NcFile ncf("foo",NcFile::ReadOnly);],
           [AC_MSG_RESULT([yes]); NETCDF_LIBS="-lnetcdf_c++ -lnetcdf -lhdf5_hl"], 
           [AC_MSG_RESULT([no]); HAVE_NETCDF=no] )    
     ])
  LIBS="$old_LIBS"
  AC_LANG_RESTORE
  CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$old_LDFLAGS"

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
