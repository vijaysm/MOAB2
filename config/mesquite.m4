#######################################################################################
# Check for MESQUITE library ((C++)
# Sets HAVE_MESQUITE to 'yes' or 'no'
# If HAVE_MESQUITE == yes, then exports:
#   MESQUITE_CPPFLAGS
#   MESQUITE_LDFLAGS
#   MESQUITE_LIBS
#######################################################################################
AC_DEFUN([FATHOM_CHECK_MESQUITE],[

AC_MSG_CHECKING([if MESQUITE support is enabled])
AC_ARG_WITH(mesquite,
[AS_HELP_STRING([--with-mesquite=DIR], [Specify MESQUITE location])
AS_HELP_STRING([--without-mesquite], [Disable support for MESQUITE file format])],
[MESQUITE_ARG=$withval
DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-mesquite=\"${withval}\""
]
, [MESQUITE_ARG=])
if test "xno" != "x$MESQUITE_ARG"; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

 # if MESQUITE support is not disabled
AC_MSG_CHECKING([if MESQUITE support available])
AC_MSG_RESULT([])
HAVE_MESQUITE=no
if test "xno" != "x$MESQUITE_ARG"; then
  HAVE_MESQUITE=yes

    # if a path is specified, update LIBS and INCLUDES accordingly
  if test "xyes" != "x$MESQUITE_ARG" && test "x" != "x$MESQUITE_ARG"; then
    if test -d "${MESQUITE_ARG}/lib"; then
      MESQUITE_LDFLAGS="-L${MESQUITE_ARG}/lib"
    elif test -d "${MESQUITE_ARG}"; then
      MESQUITE_LDFLAGS="-L${MESQUITE_ARG}/lib"
    elif test -d "${MESQUITE_ARG}/lib"; then
      MESQUITE_LDFLAGS="-L${MESQUITE_ARG}"
    else
      AC_MSG_ERROR("$MESQUITE_ARG is not a directory.")
    fi
    if test -d "${MESQUITE_ARG}/include"; then
      MESQUITE_CPPFLAGS="-I${MESQUITE_ARG}/include"
    else
      MESQUITE_CPPFLAGS="-I${MESQUITE_ARG}"
    fi
  fi
  
  old_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$MESQUITE_CPPFLAGS $CPPFLAGS"
  old_LDFLAGS="$LDFLAGS"
  LDFLAGS="$MESQUITE_LDFLAGS $LDFLAGS"
  
   # Check for C library
  AC_LANG_PUSH([C++])
  AC_CHECK_HEADERS( [MeshInterface.hpp], [], [HAVE_MESQUITE=no] )
  AC_CHECK_HEADERS( [MsqVertex.hpp], [], [HAVE_MESQUITE=no] )
  AC_CHECK_HEADERS( [Mesquite.hpp], [], [HAVE_MESQUITE=no] )
  CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$old_LDFLAGS"

  if test "x$HAVE_MESQUITE" = "xno"; then
    if test "x$MESQUITE_ARG" != "x"; then
      AC_MSG_ERROR("MESQUITE not found or not working")
    else
      AC_MSG_CHECKING([unsuccessful, MESQUITE support disabled])
      AC_MSG_RESULT([])
    fi
    MESQUITE_CPPFLAGS=
    MESQUITE_LDFLAGS=
  else
    MESQUITE_LIBS="-lmesquite"
  fi
fi

]) # FATHOM_HAVE_MESQUITE
