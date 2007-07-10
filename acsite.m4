#######################################################################################
# Get libtool configuration variable
# Arguments:
#  libtool config tag (e.g CXX)
#  libtool variable name
#  variable in which to store result
#######################################################################################
AC_DEFUN( [ITAPS_LIBTOOL_VAR], [
  $3=`./libtool --tag=$1 --config | sed -e 's/^$2=//p' -e 'd' | tr -d '"\n'`
])

#######################################################################################
# Implement checks for C and C++ compilers, with all corresponding options
#
# Sets the following variables:
#  CPP      - The C preprocessor
#  CC       - The C compiler
#  CXX      - The C++ compiler
#  CFLAGS   - C compiler flags
#  CXXFLAGS - C++ compiler flags
#  WITH_MPI - 'yes' if parallel support, 'no' otherwise
#
#######################################################################################
AC_DEFUN([SNL_CHECK_COMPILERS], [

  # Save these before calling AC_PROG_CC or AC_PROG_CXX
  # because those macros will modify them, and we want
  # the original user values, not the autoconf defaults.
USER_CXXFLAGS="$CXXFLAGS"
USER_CFLAGS="$CFLAGS"

  # Check for Parallel
  # Need to check this early so we can look for the correct compiler
AC_ARG_WITH( [mpi], AC_HELP_STRING([[--with-mpi(=DIR)]], [Enable parallel support]),
             [WITH_MPI=$withval],[WITH_MPI=no] )
case "x$WITH_MPI" in
  xno)
    CC_LIST="cc gcc cl egcs pgcc"
    CXX_LIST="CC aCC cxx xlC_r xlC c++ g++ pgCC gpp cc++ cl FCC KCC RCC"
    ;;
  xyes)
    CC_LIST="mpicc mpcc"
    CXX_LIST="mpiCC mpCC mpicxx"
    ;;
  x*)
    if test -z "$CC";then
      for prog in mpicc mpcc; do
        if test -x ${WITH_MPI}/bin/$prog; then
          CC="${WITH_MPI}/bin/$prog"
          CC_LIST="$prog"
        fi
      done
    else
      CC_LIST="$CC"
    fi
    if test -z "$CXX";then
      for prog in mpicxx mpiCC mpCC mpicxx; do
        if test -x ${WITH_MPI}/bin/$prog; then
          CXX="${WITH_MPI}/bin/$prog"
          CXX_LIST="$prog"
        fi
      done
    else
      CXX_LIST="$CXX"
    fi
    WITH_MPI=yes
    ;;
esac
AC_PROG_CC( [$CC_LIST] )
AC_PROG_CPP
AC_PROG_CXX( [$CXX_LIST] )
AC_PROG_CXXCPP

# Try to determine compiler-specific flags.  This must be done
# before setting up libtool so that it can override libtool settings.
SNL_CC_FLAGS
SNL_CXX_FLAGS
CFLAGS="$USER_CFLAGS $SNL_CC_SPECIAL"
CXXFLAGS="$USER_CXXFLAGS $SNL_CXX_SPECIAL"

# On IBM/AIX, the check for OBJEXT fails for the mpcc compiler.
# (Comment out this hack, it should be fixed correctly now)
#if test "x$OBJEXT" = "x"; then
#  OBJEXT=o
#fi

  # Check for debug flags
AC_ARG_ENABLE( debug, AC_HELP_STRING([--enable-debug],[Debug symbols (-g)]),
               [enable_debug=$enableval], [enable_debug=] )  
AC_ARG_ENABLE( optimize, AC_HELP_STRING([--enable-optimize],[Compile optimized (-O2)]),
               [enable_cxx_optimize=$enableval
                enable_cc_optimize=$enableval], 
               [enable_cxx_optimize=
                enable_cc_optimize=] )

# Do enable_optimize by default, unless user has specified
# custom CXXFLAGS or CFLAGS
if test "x$enable_debug" = "x"; then
  if test "x$enable_cxx_optimize" = "x"; then
    if test "x$USER_CXXFLAGS" = "x"; then
      enable_cxx_optimize=yes
    fi
  fi
  if test "x$enable_cc_optimize" = "x"; then
    if test "x$USER_CFLAGS" = "x"; then
      enable_cc_optimize=yes
    fi
  fi
fi

# Choose compiler flags from CLI args
if test "xyes" = "x$enable_debug"; then
  CXXFLAGS="$CXXFLAGS -g"
  CFLAGS="$CLFAGS -g"
fi
if test "xyes" = "x$enable_cxx_optimize"; then
  CXXFLAGS="$CXXFLAGS -O2 -DNDEBUG"
fi
if test "xyes" = "x$enable_cc_optimize"; then
  CFLAGS="$CFLAGS -O2 -DNDEBUG"
fi

  # Check for 32/64 bit.
  # This requires SNL_CXX_FLAGS and SNL_CC_FLAGS to have been called first
AC_ARG_ENABLE( 32bit, AC_HELP_STRING([--enable-32bit],[Force 32-bit objects]),
[
  if test "xyes" != "x$enableval"; then
    AC_MSG_ERROR([Unknown argument --enable-32bit=$enableval])
  elif test "x" = "x$SNL_CXX_32BIT"; then
    AC_MSG_ERROR([Don't know how to force 32-bit C++ on this platform.  Try setting CXXFLAGS manually])
  elif test "x" = "x$SNL_CC_32BIT"; then
    AC_MSG_ERROR([Don't know how to force 32-bit C on this platform.  Try setting CFLAGS manually])
  fi
  CXXFLAGS="$CXXFLAGS $SNL_CXX_32BIT"
  CFLAGS="$CFLAGS $SNL_CC_32BIT"
  enable_32bit=yes
])
# This requires SNL_CXX_FLAGS and SNL_CC_FLAGS to have been called first
AC_ARG_ENABLE( 64bit, AC_HELP_STRING([--enable-64bit],[Force 64-bit objects]),
[
  if test "xyes" != "x$enableval"; then
    AC_MSG_ERROR([Unknown argument --enable-64bit=$enableval])
  elif test "x" = "x$SNL_CXX_64BIT"; then
    AC_MSG_ERROR([Don't know how to force 64-bit C++ on this platform.  Try setting CXXFLAGS manually])
  elif test "x" = "x$SNL_CC_64BIT"; then
    AC_MSG_ERROR([Don't know how to force 64-bit C on this platform.  Try setting CFLAGS manually])
  elif test "xyes" = "x$enable_32bit"; then
    AC_MSG_ERROR([Cannot do both --enable-32bit and --enable-64bit])
  fi
  CXXFLAGS="$CXXFLAGS $SNL_CXX_64BIT"
  CFLAGS="$CFLAGS $SNL_CC_64BIT"
])

]) # SNL_CHECK_COMPILERS



#######################################################################################
# Extract the value of a variable from a makefile fragment.
# Arguments:
#   - The path of the makefile fragment
#   - The name of the makefile variable
#   - Action on success (the shell variable "make_val" will contain the value
#         of the variable)
#   - Action on failure
#######################################################################################
AC_DEFUN([SNL_MAKE_INC_VAR], [
make_val=
snl_makefile="snl_check.mak"
rm -f $snl_makefile

if test ! -f $1 ; then
  AC_MSG_WARN([File not found: $1])
  $4
else
cat >$snl_makefile <<SNL_END_OF_MAKEFILE
default:
	@echo "\$($2)"

include $1
SNL_END_OF_MAKEFILE
if make -f $snl_makefile > /dev/null 2>&1; then
  make_val=`make -f $snl_makefile`
  rm -f $snl_makefile
  $3
else
  rm -f $snl_makefile
  $4
fi
fi
])


#######################################################################################
# Check for existance of new no-file-extension C++ headers.
# Conditinionally sets the following preprocessor macros:
#   CANT_USE_STD
#    - The new no-extension versions of the c++ headers do not
#      exist or do not define symbols to be in the "std" namespace.
#      For example: if this is defined #include <map.h> rather than <map>.
#   CANT_USE_STD_IO
#    - The new no-extension versions of the c++ IO headers do not
#      exist or do not define symbols to be in the "std" namespace.
#      For example: if this is defined, #include <iostrea.h> rather than <iostream>
#######################################################################################
AC_DEFUN([SNL_CANT_USE_STD], [
AC_LANG_SAVE
AC_LANG_CPLUSPLUS

AC_MSG_CHECKING([for C++-standard header files])
AC_TRY_COMPILE([#include <vector>
                #include <string>
                #include <map>
                #include <algorithm>
               ],
               [std::vector<std::string> v;
                std::map<std::string,std::string> m;
               ],
               [AC_MSG_RESULT(yes)],
               [AC_MSG_RESULT(no); CANT_USE_STD=-DCANT_USE_STD])

AC_MSG_CHECKING([for C++-standard I/O header files])
AC_TRY_COMPILE([#include <iosfwd>
                #include <iostream>
                #include <ostream>
                #include <sstream>
               ],
               [std::cout << std::endl;],
               [AC_MSG_RESULT(yes)],
               [AC_MSG_RESULT(no); CANT_USE_STD_IO=-DCANT_USE_STD_IO])

AC_LANG_RESTORE
]) # SNL_CANT_USE_STD




#######################################################################################
# Check if template definitions (.cpp files) must be
# included (in the .hpp files).
# Sets TEMPLATE_DEFS_INCLUDED=1
#######################################################################################
AC_DEFUN([SNL_TEMPLATE_DEFS_INCLUDED], [
AC_LANG_SAVE
AC_LANG_CPLUSPLUS

AC_MSG_CHECKING([if template definitions should be included])

src=conftest.cc
templ=conftemp.cc
exe=conftest

echo "template <typename T> class MyTempl { public: T get() const; };" >$templ
echo "template <typename T> T MyTempl<T>::get() const { return 0; }"  >>$templ
echo "template <typename T> class MyTempl { public: T get() const; };" >$src
echo "int main( ) { MyTempl<int> c; return c.get(); }"                >>$src
if $CXX $CXXFLAGS $LDFLAGS $src $templ -o $exe >/dev/null 2>/dev/null; then
  TEMPLATE_DEFS_INCLUDED=
  AC_MSG_RESULT(no)
else
  TEMPLATE_DEFS_INCLUDED=-DTEMPLATE_DEFS_INCLUDED
  AC_MSG_RESULT(yes)
fi

rm -f $src $templ $exe
AC_LANG_RESTORE
]) # SNL_TEMPLATE_DEFS_INCLUDED




#######################################################################################
# Check for HDF5 library and related stuff
# Sets HAVE_HDF5 to 'yes' or 'no'
# If HAVE_HDF5 == yes, then updates INCLUDES and LIBS accordingly.
#######################################################################################
AC_DEFUN([SNL_CHECK_HDF5],[

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
    LIBS="$LIBS -L${WITH_ZLIB}/lib"
    ;;
esac
HAVE_ZLIB=no
if test "x$WITH_ZLIB" != "xno"; then
  AC_CHECK_LIB([z],[deflate],[HAVE_ZLIB=yes],
    [if test "x$WITH_ZLIB" != "x"; then AC_MSG_ERROR([Could not find zlib]); fi])
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
    LIBS="$LIBS -L${WITH_SZIP}/lib"
    ;;
esac
HAVE_SZIP=no
if test "x$WITH_SZIP" != "xno"; then
  AC_CHECK_LIB([sz],[SZ_Decompress],[HAVE_SZIP=yes],
    [if test "x$WITH_SZIP" != "x"; then AC_MSG_ERROR([Could not find libsz]); fi])
fi

  # CLI option for extra HDF5 link flags
AC_ARG_WITH([--with-hdf5-ldflags],[AC_HELP_STRING([--with-hdf5-ldflags=...],
 [Extra LDFLAGS required for HDF5 library (e.g. parallel IO lib)])],
 [HDF5_LDFLAGS="$withval"
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-hdf5-ldflags=\"${withval}\""
],[HDF5_LDFLAGS=])
case "x$HDF5_LDFLAGS" in
  xno)
    AC_MSG_ERROR("Invalid argument: --without-hdf5-ldflags")
    ;;
  xyes)
    AC_MSG_ERROR("Nonsensical argument:  --with-hdf5-ldflags without any specified flags")
    ;;
esac


  # CLI option for HDF5
AC_MSG_CHECKING([if HDF5 support is enabled])
AC_ARG_WITH(hdf5, 
[AC_HELP_STRING([--with-hdf5=DIR], [Specify HDF5 library to use for native file format])
AC_HELP_STRING([--without-hdf5], [Disable support for native HDF5 file format])],
[HDF5_ARG=$withval
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-hdf5=\"${withval}\""
], [HDF5_ARG=yes])
if test "xno" = "x$HDF5_ARG"; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([yes])
fi

 # if HDF5 support is not disabled, check to make sure we have HDF5
HAVE_HDF5=no
if test "xno" != "x$HDF5_ARG"; then
  HAVE_HDF5=yes
    # Check for IBM parallel IO library
  if test "x$WITH_MPI" != "xno"; then
    AC_CHECK_LIB([gpfs],[gpfs_stat],[LIBS="-lgpfs $LIBS"])
  fi

    # if a path is specified, update LIBS and INCLUDES accordingly
  if test "xyes" != "x$HDF5_ARG"; then
    if test -d "${HDF5_ARG}/lib"; then
      HDF5_LDFLAGS="$HDF5_LDFLAGS -L${HDF5_ARG}/lib"
    elif test -d "${HDF5_ARG}"; then
      HDF5_LDFLAGS="$HDF5_LDFLAGS -L${HDF5_ARG}"
    else
      AC_MSG_ERROR("$HDF5_ARG is not a directory.")
    fi
    if test -d "${HDF5_ARG}/include"; then
      INCLUDES="$INCLUDES -I${HDF5_ARG}/include"
    else
      INCLUDES="$INCLUDES -I${HDF5_ARG}"
    fi
  fi
  
    # Add flag to defines
  old_LIBS="$LIBS"
  LIBS="$LIBS $HDF5_LDFLAGS"
  
    # check for libraries and headers
  old_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $INCLUDES"
  AC_CHECK_HEADERS( [hdf5.h], [], [HAVE_HDF5=no] )
  CPPFLAGS="$old_CPPFLAGSS"
  
  HAVE_LIB_HDF5=no
  AC_CHECK_LIB( [hdf5], [H5Fopen], [HAVE_LIB_HDF5=yes] )
  if test $HAVE_LIB_HDF5 = no; then
    if test $HAVE_ZLIB = yes; then
      unset ac_cv_lib_hdf5_H5Fopen
      AC_CHECK_LIB( [hdf5], [H5Fopen], [HAVE_LIB_HDF5=yes; old_LIBS="-lz $old_LIBS"], [], [-lz] )
    fi
  fi
  if test $HAVE_LIB_HDF5 = no; then
    if test $HAVE_SZIP = yes; then
      unset ac_cv_lib_hdf5_H5Fopen
      AC_CHECK_LIB( [hdf5], [H5Fopen], [HAVE_LIB_HDF5=yes; old_LIBS="-lsz $old_LIBS"], [], [-lsz] )
    fi
  fi
  if test $HAVE_LIB_HDF5 = no; then
    if test $HAVE_SZIP = yes; then
      if test $HAVE_ZLIB = yes; then
        unset ac_cv_lib_hdf5_H5Fopen
        AC_CHECK_LIB( [hdf5], [H5Fopen], [HAVE_LIB_HDF5=yes; old_LIBS="-lsz -lz $old_LIBS"], [], [-lz -lsz] )
      fi
    fi
  fi
  if test "x$HAVE_LIB_HDF5" = "xno"; then
    AC_MSG_WARN([Could not find HDF5 library (-lhdf5)])
    HAVE_HDF5=no
  fi
  
  if test "x$HAVE_HDF5" = "xyes"; then
    LIBS="$HDF5_LDFLAGS -lhdf5 $old_LIBS"
  else
    LIBS="$old_LIBS"
  fi
fi


]) # SNL_CHECK_HDF5




#######################################################################################
# Check for NetCDF library ((C++)
# Sets HAVE_NETCDF to 'yes' or 'no'
# If HAVE_NETCDF == yes, then updates INCLUDES and LIBS accordingly.
#######################################################################################
AC_DEFUN([SNL_CHECK_NETCDF],[

AC_MSG_CHECKING([if NetCDF support is enabled])
AC_ARG_WITH(netcdf, 
[AC_HELP_STRING([--with-netcdf=DIR], [Specify NetCDF library to use for ExodusII file format])
AC_HELP_STRING([--without-netcdf], [Disable support for ExodusII file format])],
[NETCDF_ARG=$withval
DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-netcdf=\"${withval}\""
]
, [NETCDF_ARG=yes])
if test "xno" = "x$NETCDF_ARG"; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([yes])
fi

 # if NetCDF support is not disabled
HAVE_NETCDF=no
if test "xno" != "x$NETCDF_ARG"; then
  HAVE_NETCDF=yes
    # Add flag to defines
  DEFINES="$DEFINES -DNETCDF_FILE"

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
  if test "xyes" != "x$NETCDF_ARG"; then
    if test -d "${NETCDF_ARG}/lib"; then
      LIBS="$LIBS -L${NETCDF_ARG}/lib"
    elif test -d "${NETCDF_ARG}"; then
      LIBS="$LIBS -L${NETCDF_ARG}"
    else
      AC_MSG_ERROR("$NETCDF_ARG is not a directory.")
    fi
    if test -d "${NETCDF_ARG}/include"; then
      INCLUDES="$INCLUDES -I${NETCDF_ARG}/include"
    elif test -d "${NETCDF_ARG}/inc"; then
      INCLUDES="$INCLUDES -I${NETCDF_ARG}/inc"
    else
      INCLUDES="$INCLUDES -I${NETCDF_ARG}"
    fi
  fi
  
  old_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $INCLUDES"
  old_CXXFLAGS="$CXXFLAGS"
  CXXFLAGS="$CXXFLAGS $INCLUDES"
  old_LIBS="$LIBS"
  
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
  CPPFLAGS="$CPPFLAGS -DSTRSTREAM_H_SPEC=$NETCDF_DEF"
  CXXFLAGS="$CXXFLAGS -DSTRSTREAM_H_SPEC=$NETCDF_DEF"
  fi
  AC_MSG_CHECKING([[for netcdf_c++ library]])
  LIBS="$LIBS -lnetcdf_c++ -lnetcdf"
  AC_TRY_LINK(
    [#include <netcdf.hh>], [NcFile ncf("foo",NcFile::ReadOnly);],
    [AC_MSG_RESULT([yes])], 
    [AC_MSG_RESULT([no]); 
     AC_MSG_ERROR([NetCDF C++ API not found])
     HAVE_NETCDF=no] )
  AC_LANG_RESTORE
  CPPFLAGS="$old_CPPFLAGS"
  CXXFLAGS="$old_CXXFLAGS"
  if test "x$NETCDF_DEF" != "x"; then
    DEFINES="$DEFINES '-DSTRSTREAM_H_SPEC=$NETCDF_DEF'"
  fi
  if test "x$HAVE_NETCDF" = "xno"; then
    AC_MSG_WARN("NetCDF support disabled")
    LIBS="$old_LIBS"
  fi
fi

]) # SNL_HAVE_NETCDF


#######################################################################################
# Check for Boost libraryies
# Sets HAVE_BOOST to 'yes' or 'no' and adds any user-specified path to INCLUDES
# Unless user specified --without-boost, will check for any passed headers.
# For any found headers, -DHAVE_[header] will be added to DEFINES
#######################################################################################
AC_DEFUN([SNL_CHECK_BOOST],[
AC_MSG_CHECKING([if boost library is enabled])
SNL_BOOST_OPT_HEADER_LIST="$1"
AC_ARG_WITH(boost, 
[AC_HELP_STRING([--with-boost=DIR], [Specify directory where boost is installed])
AC_HELP_STRING([--without-boost], [Disable support for boost libraries])],
[BOOST_ARG=$withval
DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-boost=\"${withval}\""
]
, [BOOST_ARG=yes])
if test "xno" = "x$BOOST_ARG"; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([yes])
  # If user-specified directory, add to includes
  if test "xyes" != "x$BOOST_ARG"; then
    INCLUDES="$INCLUDES -I$BOOST_ARG"
  fi
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_CHECK_HEADERS( [$SNL_BOOST_OPT_HEADER_LIST],[def=`echo "$ac_header" | $as_tr_cpp`; DEFINES="$DEFINES -DHAVE_$def"] )
  AC_LANG_RESTORE
fi
]) # SNL_CHECK_BOOST

#######################################################################################
# Macro to set the following variables:
# ACIS_VERSION
# ACIS_LIB_DIR
# ACIS_DEBUG_LIB_DIR
# ACIS_PLATFORM
# This macro expects ACIS_DIR to be set.
# If ACIS_SYSTEM is set, ACIS_LIB_DIR will
# be constructed using that value rather than
# via a trial-and-error search.
# If ACIS_VERSION is set, the corresponding test
# will be skipped.
#######################################################################################
AC_DEFUN([SNL_ACIS_ENV], [

AC_CHECK_FILE([$ACIS_DIR/include/acis.hxx], [], [AC_MSG_ERROR("Invalid ACIS path")])
 
if test "x" = "x$ACIS_SYSTEM"; then
  dir_list="$ACIS_DIR/bin/*"
else
  dir_list="$ACIS_DIR/bin/$ACIS_SYSTEM"
fi

for dir in $dir_list; do
    # first try non-debug libraries
  case "$dir" in
    *_debug)
      ;;
    *)
      if ! test "$ACIS_VERSION" -gt "600"; then
        SNL_ACIS_VERSION( [$dir] )
        ACIS_LIB_DIR="$dir"
      fi
      ;;
  esac
    # next try debug libraries
  case "$dir" in
    *_debug)
      if ! test "$ACIS_VERSION" -gt "600"; then
        SNL_ACIS_VERSION( [$dir] )
        ACIS_DEBUG_LIB_DIR="$dir"
      fi
      ;;
  esac
done

if ! test "$ACIS_VERSION" -gt "600"; then
  AC_MSG_ERROR([ACIS configuration failed.  Verify ACIS directory.  Try specifying ACIS system and version])
fi

# If either ACIS_LIB_DIR or ACIS_DEBUG_LIB_DIR is not defined,
# make both the same
if test "x" = "x$ACIS_LIB_DIR"; then
  ACIS_LIB_DIR="$ACIS_DEBUG_LIB_DIR"
elif test "x" = "x$ACIS_DEBUG_LIB_DIR"; then
  ACIS_DEBUG_LIB_DIR="$ACIS_LIB_DIR"
fi

# Determine the ACIS platform name from the lib dir
AC_MSG_CHECKING([ACIS platform name])
case "$ACIS_LIB_DIR" in
  $ACIS_DIR/bin/aix*)  
    ACIS_PLATFORM=aix
    ;;
  $ACIS_DIR/bin/hp700*)
    ACIS_PLATFORM=hp700
    ;;
  $ACIS_DIR/bin/linux*)
    ACIS_PLATFORM=linux
    ;;
  $ACIS_DIR/bin/mac*)
    ACIS_PLATFORM=mac
    ;;
  $ACIS_DIR/bin/osf1*)
    ACIS_PLATFORM=osf1
    ;;
  $ACIS_DIR/bin/sgi*)
    ACIS_PLATFORM=sgi
    ;;
  $ACIS_DIR/bin/solaris*)
    ACIS_PLATFORM=solaris
    ;;
  *)
    AC_MSG_ERROR([Cannot determine ACIS platform name.])
    ;;
esac
AC_MSG_RESULT([$ACIS_PLATFORM])
])



#######################################################################################
# Macro to check for ACIS translators.
#######################################################################################
AC_DEFUN([SNL_ACIS_TRANSLATOR], [
AC_REQUIRE([SNL_ACIS_ENV])
case "$ACIS_VERSION" in
  11??)
    ACIS_XLATE_LIBS='-lxacis2k -lxcore2k -lxutil'
    ;;
  1[23]??)
    ACIS_XLATE_LIBS='-lSPAXAssemblyRep -lSPAXInterop -lSPAXBase -lxacis2k -lxcore2k -lxutil'
    ;;
  14??)
    ACIS_XLATE_LIBS='-lSPAXAssemblyRep -lSPAXInterop -lSPAXAcisBase -lSPAXDefaultGeometryRep -lSPAXGeometryRepresentation -lSPAXPropertiesBRepImporter -lSPAXPropertiesBase -lSPAXBase -lxacis2k -lxcore2k'
    ;;
  *)
    ACIS_XLATE_LIBS='-lSPAXEBOMBase -lSPAXEBOMNewAssemblyExporter -lSPAXEBOMNewAssemblyImporter -lSPAXXMLTk -lpthread -lSPAXPropertiesNewAssemblyImporter -lSPAXAssemblyRep -lSPAXInterop -lSPAXAcisBase -lSPAXDefaultGeometryRep -lSPAXGeometryRepresentation -lSPAXPropertiesBRepImporter -lSPAXPropertiesBase -lSPAXBase -lxacis2k -lxcore2k'
    ;;
esac
old_LIBS="$LIBS"
LIBS="$LIBS -L$ACIS_LIB_DIR $ACIS_XLATE_LIBS $ACIS_BASE_LIBS"
AC_CHECK_LIB([xstep],[main],[ACIS_STEP_TRANSLATOR=-DACIS_STEP_TRANSLATOR])
AC_CHECK_LIB([xiges],[main],[ACIS_IGES_TRANSLATOR=-DACIS_IGES_TRANSLATOR])
LIBS="$old_LIBS"
])


#######################################################################################
# *******************************************************************************
# **************************** INTERNAL STUFF ***********************************
# *******************************************************************************
#######################################################################################

#################################################################################
# Check if the compiler defines a specific preprocessor macro
# Arguments:
#  - preprocessor define to check for
#  - action upon success
#  - action upon failure
#################################################################################
AC_DEFUN([SNL_TRY_COMPILER_DEFINE], [
AC_COMPILE_IFELSE([
AC_LANG_PROGRAM( [[#ifndef $1
  choke me
#endif]], []) ],
[$2],[$3])
])


#######################################################################################
# Check for compiler-specific flags.
# Sets the following environmental variables:
#   SNL_CXX_SPECIAL : Any compiler-specific flags which must be specified
#   SNL_CXX_32BIT   : Flag to force compilation of 32-bit code
#   SNL_CXX_64BIT   : Flag to force compilation of 64-bit code
#######################################################################################
AC_DEFUN([SNL_CXX_FLAGS], [
AC_REQUIRE([AC_PROG_CXX])

# Detect compiler 
AC_MSG_CHECKING([for known c++ compilers])
# Autoconf does G++ for us
if test x$GXX = xyes; then
  cxx_compiler=GNU
# Search for other compiler types
# For efficiency, limit checks to relevant OSs
else
  cxx_compiler=unknown
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  case "$target_os" in
    aix*)
      SNL_TRY_COMPILER_DEFINE([__IBMCPP__],[cxx_compiler=VisualAge])
      ;;
    solaris*|sunos*)
      SNL_TRY_COMPILER_DEFINE([__SUNPRO_CC],[cxx_compiler=SunWorkshop])
      ;;
    irix*)
      SNL_TRY_COMPILER_DEFINE([__sgi],[cxx_compiler=MIPSpro])
      ;;
    linux*)
      SNL_TRY_COMPILER_DEFINE([__INTEL_COMPILER],[cxx_compiler=Intel])
      SNL_TRY_COMPILER_DEFINE([__IBMCPP__],[cxx_compiler=VisualAge])
      SNL_TRY_COMPILER_DEFINE([__DECCXX_VER],[cxx_compiler=Compaq])
      SNL_TRY_COMPILER_DEFINE([__SUNPRO_CC],[cxx_compiler=SunWorkshop])
      SNL_TRY_COMPILER_DEFINE([__PGI],[cxx_cmopiler=PortlandGroup])
      ;;
    hpux*)
      SNL_TRY_COMPILER_DEFINE([__HP_aCC],[cxx_compiler=HP])
      ;;
    windows*)
      SNL_TRY_COMPILER_DEFINE([__MSC_VER],[cxx_compiler=VisualStudio])
      SNL_TRY_COMPILER_DEFINE([__INTEL_COMPILER],[cxx_compiler=Intel])
      SNL_TRY_COMPILER_DEFINE([__DECCXX_VER],[cxx_compiler=Compaq])
      SNL_TRY_COMPILER_DEFINE([__BORLANDC__],[cxx_compiler=Borland])
      SNL_TRY_COMPILER_DEFINE([__CYGWIN__],[cxx_compiler=Cygwin])
      SNL_TRY_COMPILER_DEFINE([__MINGW32__],[cxx_compiler=MinGW])
      ;;
    *)
      SNL_TRY_COMPILER_DEFINE([__PGI],[cc_cmopiler=PortlandGroup])
      ;;
  esac
  AC_LANG_RESTORE
fi
AC_MSG_RESULT([$cxx_compiler])
if test "x$cxx_compiler" = "xunknown"; then
  AC_MSG_WARN([Unrecognized C++ compiler: $CXX])
fi

# Now decide special compiler flags using compiler/cpu combination
AC_MSG_CHECKING([for known compiler/OS combinations])
case "$cxx_compiler:$host_cpu" in
  GNU:sparc*)
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    SNL_CXX_SPECIAL="-Wall -pipe"
    ;;
  GNU:powerpc*)
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    SNL_CXX_SPECIAL="-Wall -pipe"
    ;;
  GNU:i?86|GNU:x86_64)
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    SNL_CXX_SPECIAL="-Wall -pipe"
    ;;
  GNU:mips*)
    SNL_CXX_32BIT="-mips32 -mabi=32"
    SNL_CXX_64BIT="-mips64 -mabi=64"
    SNL_CXX_SPECIAL="-Wall -pipe"
    ;;
  VisualAge:*)
    SNL_CXX_32BIT=-q32
    SNL_CXX_64BIT=-q64
    # Do V5.0 namemangling for compatibility with ACIS, and enable RTTI
    SNL_CXX_SPECIAL="-qrtti=all -qnamemangling=v5"
    AR="ar -X 32_64"
    NM="nm -B -X 32_64"
    ;;
  MIPSpro:mips)
    SNL_CXX_32BIT=-n32
    SNL_CXX_64BIT=-64
    SNL_CXX_SPECIAL=-LANG:std
    ;;
  MIPSpro:*)
    SNL_CXX_SPECIAL=-LANG:std
    ;;
  SunWorkshop:sparc*)
    SNL_CXX_32BIT=-xarch=generic
    SNL_CXX_64BIT=-xarch=generic64
    ;;
  SunWorkshop:i?86|SunWorkshop:x86_64)
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    ;;
  *)
    ;;
esac
AC_MSG_RESULT([$cxx_compiler:$host_cpu])
]) # end SNL_CXX_FLAGS

#######################################################################################
# Check for compiler-specific flags.
# Sets the following environmental variables:
#   SNL_CC_SPECIAL : Any compiler-specific flags which must be specified
#   SNL_CC_32BIT   : Flag to force compilation of 32-bit code
#   SNL_CC_64BIT   : Flag to force compilation of 64-bit code
#######################################################################################
AC_DEFUN([SNL_CC_FLAGS], [
AC_REQUIRE([AC_PROG_CC])

# Detect compiler 

AC_MSG_CHECKING([for known C compilers])
# Autoconf does gcc for us
if test x$GCC = xyes; then
  cc_compiler=GNU
# Search for other compiler types
# For efficiency, limit checks to relevant OSs
else
  cc_compiler=unknown
  case "$target_os" in
    aix*)
      SNL_TRY_COMPILER_DEFINE([__IBMC__],[cc_compiler=VisualAge])
      ;;
    solaris*|sunos*)
      SNL_TRY_COMPILER_DEFINE([__SUNPRO_C],[cc_compiler=SunWorkshop])
      ;;
    irix*)
      SNL_TRY_COMPILER_DEFINE([__sgi],[cc_compiler=MIPSpro])
      ;;
    linux*)
      SNL_TRY_COMPILER_DEFINE([__INTEL_COMPILER],[cc_compiler=Intel])
      SNL_TRY_COMPILER_DEFINE([__IBMC__],[cc_compiler=VisualAge])
      SNL_TRY_COMPILER_DEFINE([__DECC_VER],[cc_compiler=Compaq])
      SNL_TRY_COMPILER_DEFINE([__SUNPRO_C],[cc_compiler=SunWorkshop])
      SNL_TRY_COMPILER_DEFINE([__PGI],[cc_cmopiler=PortlandGroup])
      ;;
    hpux*)
      SNL_TRY_COMPILER_DEFINE([__HP_cc],[cc_compiler=HP])
      ;;
    windows*)
      SNL_TRY_COMPILER_DEFINE([__MSC_VER],[cc_compiler=VisualStudio])
      SNL_TRY_COMPILER_DEFINE([__INTEL_COMPILER],[cc_compiler=Intel])
      SNL_TRY_COMPILER_DEFINE([__DECC_VER],[cc_compiler=Compaq])
      SNL_TRY_COMPILER_DEFINE([__BORLANDC__],[cc_compiler=Borland])
      SNL_TRY_COMPILER_DEFINE([__TURBOC__],[cc_compiler=TurboC])
      SNL_TRY_COMPILER_DEFINE([__CYGWIN__],[cc_compiler=Cygwin])
      SNL_TRY_COMPILER_DEFINE([__MINGW32__],[cc_compiler=MinGW])
      ;;
    *)
      SNL_TRY_COMPILER_DEFINE([__PGI],[cc_cmopiler=PortlandGroup])
      ;;
  esac
fi
AC_MSG_RESULT([$cc_compiler])
if test "x$cc_compiler" = "xunknown"; then
  AC_MSG_WARN([Unrecognized C compiler: $CXX])
fi

# Now decide special compiler flags using compiler/cpu combination
AC_MSG_CHECKING([for known compiler/OS combinations])
case "$cc_compiler:$host_cpu" in
  GNU:sparc*)
    SNL_CC_32BIT=-m32
    SNL_CC_64BIT=-m64
    SNL_CC_SPECIAL="-Wall -pipe"
    ;;
  GNU:powerpc*)
    SNL_CC_32BIT=-m32
    SNL_CC_64BIT=-m64
    SNL_CC_SPECIAL="-Wall -pipe"
    ;;
  GNU:i?86|GNU:x86_64)
    SNL_CC_32BIT=-m32
    SNL_CC_64BIT=-m64
    SNL_CC_SPECIAL="-Wall -pipe"
    ;;
  GNU:mips*)
    SNL_CC_32BIT="-mips32 -mabi=32"
    SNL_CC_64BIT="-mips64 -mabi=64"
    SNL_CC_SPECIAL="-Wall -pipe"
    ;;
  VisualAge:*)
    SNL_CC_32BIT=-q32
    SNL_CC_64BIT=-q64
    AR="ar -X 32_64"
    NM="nm -B -X 32_64"
    ;;
  MIPSpro:mips)
    SNL_CC_32BIT=-n32
    SNL_CC_64BIT=-64
    SNL_CC_SPECIAL=-LANG:std
    ;;
  MIPSpro:*)
    SNL_CC_SPECIAL=-LANG:std
    ;;
  SunWorkshop:sparc*)
    SNL_CC_32BIT=-xarch=generic
    SNL_CC_64BIT=-xarch=generic64
    ;;
  SunWorkshop:i?86|SunWorkshop:x86_64)
    SNL_CC_32BIT=-m32
    SNL_CC_64BIT=-m64
    ;;
  *)
    ;;
esac
AC_MSG_RESULT([$cc_compiler:$host_cpu])
]) # end SNL_CC_FLAGS


#######################################################################################
# Macro to get ACIS_VERSION
# Expected arguments: ACIS library path
#######################################################################################
AC_DEFUN([SNL_ACIS_VERSION], [
AC_REQUIRE([AC_PROG_LIBTOOL])
AC_LANG_SAVE
AC_LANG_CPLUSPLUS
AC_MSG_CHECKING([ACIS version of $1])
src=conftest.cc
exe=conftest

cat <<END_SNL_ACIS_VERSION_SRC  >$src
#include <stdio.h>
int get_major_version();
int get_minor_version();
int get_point_version();
int main(int,char**) { printf("%d\n", 
100*get_major_version() +
 10*get_minor_version() +
    get_point_version()); 
return 0;}

END_SNL_ACIS_VERSION_SRC

FLAGS="$CXXFLAGS $LDFLAGS -L$1 -R$1"
if ./libtool --mode=link $CXX $FLAGS $src -lSpaBase -o $exe >/dev/null 2>/dev/null; then
  ACIS_VERSION=`./$exe || echo 0`
  AC_MSG_RESULT($ACIS_VERSION)
else
  ACIS_VERSION=0
  AC_MSG_RESULT(failed)
fi

#rm -f $src $exe
AC_LANG_RESTORE
]) # SNL_ACIS_VERSION
