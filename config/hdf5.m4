AC_DEFUN([FATHOM_HDF5_LIBS_HELPER],[
if test $HAVE_LIB_HDF5 = no; then
   unset "ac_cv_lib_${HDF5_LIBNAME}_H5Fopen"
   unset "ac_cv_lib_${HDF5_LIBNAME}___H5Fopen"
   AC_CHECK_LIB( [${HDF5_LIBNAME}], [H5Fopen], [HAVE_LIB_HDF5=yes; HDF5_LIBS="$HDF5_LIBS $1"], [], [$1] )
fi
])

#######################################################################################
# Helper function for FATHOM_CHECK_HDF5 and FATHOM_CHECK_NETCDF
# If HAVE_LIB_HDF5 == yes, then does nothing.
# Otherwise sets HAVE_LIB_HDF5 to yes or no depending on whether or not
# HDF5 libraries are detected and sets HDF5_LIBS
# Respects caller's LDFLAGS, CPPFLAGS, and LIBS.  Caller should set appropriately.
# Respections optional HDF5_LIBNAME uses to specify alternate library name.  If
# not specified, will be set to -lhdf5
#######################################################################################
AC_DEFUN([FATHOM_DETECT_HDF5_LIBS],[

 # if we've already done this check, then don't do it again
if test "xyes" != "x$HAVE_LIB_HDF5"; then
    # Check for IBM parallel IO library
  if test "x$WITH_MPI" != "xno"; then
    AC_CHECK_LIB([gpfs],[gpfs_stat],[HDF5_LIBS="-lgpfs $HDF5_LIBS"])
  fi
  
  test "x" != "x$HDF5_LIBNAME" || HDF5_LIBNAME=hdf5
  
  HAVE_LIB_HDF5=no
  FATHOM_HDF5_LIBS_HELPER
  if test $HAVE_ZLIB = yes; then
    FATHOM_HDF5_LIBS_HELPER([-lz])
    FATHOM_HDF5_LIBS_HELPER([-lz -lpthread])
  fi
  if test $HAVE_SZIP = yes; then
    FATHOM_HDF5_LIBS_HELPER([-lsz])
    FATHOM_HDF5_LIBS_HELPER([-lsz -lpthread])
    if test $HAVE_ZLIB = yes; then
      FATHOM_HDF5_LIBS_HELPER([-lsz -lz])
      FATHOM_HDF5_LIBS_HELPER([-lsz -lz -lpthread])
    fi
  fi
fi
])

#######################################################################################
# Check for HDF5 library and related stuff
# Sets HAVE_HDF5 to 'yes' or 'no'
# If HAVE_HDF5 == yes, then sets:
#  HDF5_CPPFLAGS
#  HDF5_LDFLAGS
#  HDF5_LIBS
#######################################################################################
AC_DEFUN([FATHOM_CHECK_HDF5],[

  # CLI option for linking zlib
AC_ARG_WITH(zlib,
  [AC_HELP_STRING([--with-zlib=DIR],[HDF5 requires zlib, and zlib can be found at...])],
  [WITH_ZLIB=$withval
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-zlib=\"${withval}\""
  ],[WITH_ZLIB=])
case "x$WITH_ZLIB" in
  xyes|xno|x)
    ;;
  *)
    if ! test -d  ${WITH_ZLIB}/lib; then
      AC_MSG_ERROR([Not a directory: ${WITH_ZLIB}/lib])
    fi
    HDF5_LDFLAGS="$HDF5_LDFLAGS -L${WITH_ZLIB}/lib"
    ;;
esac
HAVE_ZLIB=no
if test "x$WITH_ZLIB" != "xno"; then
  old_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
  AC_CHECK_LIB([z],[deflate],[HAVE_ZLIB=yes],
    [if test "x$WITH_ZLIB" != "x"; then AC_MSG_ERROR([Could not find zlib]); fi])
  LDFLAGS="$old_LDFLAGS"
fi

  # CLI option for linking szip
AC_ARG_WITH(szip,
  [AC_HELP_STRING([--with-szip=DIR],[HDF5 requires szip, and szip an be found at...])],
  [WITH_SZIP=$withval
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-szip=\"${withval}\""
  ],[WITH_SZIP=])
case "x$WITH_SZIP" in
  xyes|xno|x)
    ;;
  *)
    if ! test -d  ${WITH_SZIP}/lib; then
      AC_MSG_ERROR([Not a directory: ${WITH_SZIP}/lib])
    fi
    HDF5_LDFLAGS="$HDF5_LDFLAGS -L${WITH_SZIP}/lib"
    ;;
esac
HAVE_SZIP=no
if test "x$WITH_SZIP" != "xno"; then
  old_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
  AC_CHECK_LIB([sz],[SZ_Decompress],[HAVE_SZIP=yes],
    [if test "x$WITH_SZIP" != "x"; then AC_MSG_ERROR([Could not find libsz]); fi])
  LDFLAGS="$old_LDFLAGS"
fi

  # CLI option for extra HDF5 link flags
AC_ARG_WITH([hdf5-ldflags],[AC_HELP_STRING([--with-hdf5-ldflags=...],
 [Extra LDFLAGS required for HDF5 library (e.g. parallel IO lib)])],
 [HDF5_LDFLAGS_WITHVAL="$withval"; 
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-hdf5-ldflags=\"${withval}\""
],[HDF5_LDFLAGS_WITHVAL=])
case "x$HDF5_LDFLAGS_WITHVAL" in
  xno)
    AC_MSG_ERROR("Invalid argument: --without-hdf5-ldflags")
    ;;
  xyes)
    AC_MSG_ERROR("Nonsensical argument:  --with-hdf5-ldflags without any specified flags")
    ;;
  *)
    HDF5_LDFLAGS="$HDF5_LDFLAGS $HDF5_LDFLAGS_WITHVAL"
    ;;
esac


  # CLI option for HDF5
AC_MSG_CHECKING([if HDF5 support is enabled])
AC_ARG_WITH(hdf5, 
[AC_HELP_STRING([--with-hdf5@<:@=DIR@:>@], [Specify HDF5 library to use for native file format])
AC_HELP_STRING([--without-hdf5], [Disable support for native HDF5 file format])],
[HDF5_ARG=$withval
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-hdf5=\"${withval}\""
], [HDF5_ARG=])
if test "xno" = "x$HDF5_ARG"; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([yes])
fi

HAVE_HDF5=no
if test "xno" != "x$HDF5_ARG"; then
  HAVE_HDF5=yes

    # if a path is specified, update LIBS and INCLUDES accordingly
  if test "xyes" != "x$HDF5_ARG" && test "x" != "x$HDF5_ARG"; then
    if test -d "${HDF5_ARG}/dll"; then
      HDF5_LDFLAGS="$HDF5_LDFLAGS -L${HDF5_ARG}/dll"
      HDF5_LIBNAME=hdf5dll
    elif test -d "${HDF5_ARG}/lib"; then
      HDF5_LDFLAGS="$HDF5_LDFLAGS -L${HDF5_ARG}/lib"
    elif test -d "${HDF5_ARG}"; then
      HDF5_LDFLAGS="$HDF5_LDFLAGS -L${HDF5_ARG}"
    else
      AC_MSG_ERROR("$HDF5_ARG is not a directory.")
    fi
    if test -d "${HDF5_ARG}/include"; then
      HDF5_CPPFLAGS="$HDF5_CPPFLAGS -I${HDF5_ARG}/include"
      if test "x$GXX" = "xyes" && test "x$GCC" = "xyes"; then
        HDF5_CPPFLAGS="$HDF5_CPPFLAGS -isystem ${HDF5_ARG}/include"
      fi
    else
      HDF5_CPPFLAGS="$HDF5_CPPFLAGS -I${HDF5_ARG}"
    fi
  fi
  
    # Add flag to defines

  old_CPPFLAGS="$CPPFLAGS"
  old_LDFLAGS="$LDFLAGS"
  old_LIBS="$LIBS"
  CPPFLAGS="$HDF5_CPPFLAGS $CPPFLAGS"
  LDFLAGS="$HDF5_LDFLAGS $LDFLAGS"
  LIBS="$HDF5_LIBS $LIBS"
  
    # check for libraries and headers
  AC_CHECK_HEADERS( [hdf5.h], [], [HAVE_HDF5=no] )

  HAVE_LIB_HDF5=no
  FATHOM_DETECT_HDF5_LIBS

  if test "x$HAVE_LIB_HDF5" = "xno"; then
    HAVE_HDF5=no
  fi
  
  if test "x$HAVE_HDF5" = "xno"; then
    if test "x" = "x$HDF5_ARG"; then
      AC_MSG_WARN([HDF5 library not found or not usable.])
    else
      AC_MSG_ERROR([HDF5 library not found or not usable.])
    fi
    HDF5_CPPFLAGS=
    HDF5_LDFLAGS=
    HDF5_LIBS=
  else
    HDF5_LIBS="-l$HDF5_LIBNAME $HDF5_LIBS"
  fi
  
  CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$old_LDFLAGS"
  LIBS="$old_LIBS"
fi

]) # FATHOM_CHECK_HDF5
