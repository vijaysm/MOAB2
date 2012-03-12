################################################################################
#                           Standard Stuff
################################################################################
AC_INIT(MOAB, 4.5.0)
AC_CONFIG_SRCDIR([moab.make.in])
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_AUX_DIR([config])
AC_CANONICAL_SYSTEM
m4_ifdef([AM_SILENT_RULES],
 [AM_INIT_AUTOMAKE([color-tests parallel-tests])],
 [AM_INIT_AUTOMAKE])

AC_CHECK_PROG( [SED], [sed], [sed], [true] )


################################################################################
#                           Compilers
################################################################################

AC_ARG_ENABLE([fortran],[AC_HELP_STRING([--disable-fortran],
  [No Fortran name mangling in ITAPS/MBCN C headers])],
  [ENABLE_FORTRAN=$enableval],[ENABLE_FORTRAN=yes])

EXTRA_GNU_FLAGS='-Wall -pipe -pedantic -Wno-long-long'
EXTRA_INTEL_FLAGS='-Wall'
FATHOM_CHECK_COMPILERS([yes],[yes],[$ENABLE_FORTRAN])
AM_CONDITIONAL(PARALLEL,[test "x$WITH_MPI" != "xno"])
AM_CONDITIONAL(ENABLE_FORTRAN,[test "x$ENABLE_FORTRAN" != "xno"])
LIBS="$LIBS -lm"

AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_DISABLE_SHARED
AC_PROG_LIBTOOL
FATHOM_COMPILER_FLAGS
m4_ifdef([AM_SILENT_RULES],[
if test "x$DEBUG" = "xyes"; then
  AM_SILENT_RULES(yes)
else
  AM_SILENT_RULES(no)
fi
])

# Check if platform is BlueGene
AC_MSG_CHECKING([if platform is IBM BlueGene])
FATHOM_TRY_COMPILER_DEFINE([__bg__],
  [AM_CPPFLAGS="$AM_CPPFLAGS -DBLUEGENE"
   AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])])

AC_CHECK_PROG([ZCAT],[gunzip],[gunzip -c],[])
AC_CHECK_PROG([ZCAT],[gzip],[gzip -cd],[])
AC_CHECK_PROG([ZCAT],[zcat],[zcat],[])
if test "x" = "x$ZCAT"; then
  AC_MSG_WARN([[Cannot run tests: no way to uncompress input files.]])
fi
AC_SUBST([ZCAT])
AM_CONDITIONAL([HAVE_ZCAT],[test "x" = "x$ZCAT"])

FATHOM_TEMPLATE_SPECIALIZATION
FATHOM_TEMPLATE_FUNC_SPECIALIZATION
AM_CPPFLAGS="$AM_CPPFLAGS $TEMPLATE_SPECIALIZATION $TEMPLATE_FUNC_SPECIALIZATION"

FATHOM_CHECK_CXX_WORKS([], [AC_MSG_ERROR([Cannot build without C++ compiler])])
if test "xyes" = "x$ENABLE_FORTRAN" && test "x" != "x$FC"; then
  AC_FC_WRAPPERS
fi

################################################################################
#                           Check for need for extra flags to support cray pointers
################################################################################
if test "xyes" = "x$ENABLE_FORTRAN" && test "x" != "x$F77"; then
  PAC_PROG_F77_CRAY_POINTER([
    FFLAGS="$FFLAGS $CRAYPTR_FFLAGS"
  ])
fi

if test "xyes" = "x$ENABLE_FORTRAN" && test "x" != "x$FC"; then
  PAC_PROG_FC_CRAY_POINTER([
    FCFLAGS="$FCFLAGS $CRAYPTR_FCFLAGS"
  ])
fi

################################################################################
#                           Basic Portability Stuff
################################################################################
AC_CHECK_FUNC([vsnprintf],[AM_CPPFLAGS="$AM_CPPFLAGS -DHAVE_VSNPRINTF"])

# Use 64-bit off_t value for file offsets.  W/out this, the stat/fseek calls
# in serial_open_file and other places will fail on 32-bit systems if the file 
# is larger that 2GB.
AM_CPPFLAGS="$AM_CPPFLAGS -D_FILE_OFFSET_BITS=64"

FATHOM_VECTOR_TEMPLATE_INSERT
FATHOM_OLD_STD_COUNT
AM_CPPFLAGS="$AM_CPPFLAGS $NO_VECTOR_TEMPLATE_INSERT $OLD_STD_COUNT"


################################################################################
#                           No-Namepsace (3.0) Headers
################################################################################
AC_ARG_ENABLE([old-headers],
 [AC_HELP_STRING([--disable-old-headers],
     [Do not install no-namespace headers for compatability with MOAB 3.0 API])],
 [OLD_HEADERS=$enableval],[OLD_HEADERS=yes])
AM_CONDITIONAL([OLD_HEADERS],[test "xyes" = "x$OLD_HEADERS"])


################################################################################
#                              HANDLE SIZE
################################################################################
AC_CHECK_HEADER([inttypes.h], [AC_DEFINE(MOAB_HAVE_INTTYPES_H,  [], [MOAB qualified HAVE_INTTYPES_H]) ])
AC_CHECK_HEADER([stdint.h],   [AC_DEFINE(MOAB_HAVE_STDINT_H,    [], [MOAB qualified HAVE_STDINT_H])   ])
AC_CHECK_HEADER([sys/types.h],[AC_DEFINE(MOAB_HAVE_SYS_TYPES_H, [], [MOAB qualified HAVE_SYS_TYPES_H])])
AC_CHECK_HEADER([stddef.h],   [AC_DEFINE(MOAB_HAVE_STDDEF_H,    [], [MOAB qualified HAVE_STDDEF_H])   ])
AC_CHECK_HEADER([stdlib.h],   [AC_DEFINE(MOAB_HAVE_STDLIB_H,    [], [MOAB qualified HAVE_STDLIB_H])   ])
SIZEOF_MBENTITYHANDLE=0
AC_CHECK_SIZEOF([void *])
SIZEOF_VOID_P=$ac_cv_sizeof_void_p
AC_SUBST(SIZEOF_VOID_P)
AC_ARG_ENABLE([64bit-handles],
  [AC_HELP_STRING([--enable-64bit-handles],[Force EntityHandle to be 64 bits])],
  [AC_CHECK_TYPE([uint64_t],
    [AC_DEFINE(MOAB_FORCE_64_BIT_HANDLES,[],[Use int64_t for handles])
     SIZEOF_MBENTITYHANDLE=8
     DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --enable-64bit-handles"],
    [AC_MSG_ERROR([Platform does not have uint64_t. Cannot force 64-bit handles])])
])
AC_ARG_ENABLE([32bit-handles],
  [AC_HELP_STRING([--enable-32bit-handles],[Force EntityHandle to be 32 bits])],
  [AC_CHECK_TYPE([uint32_t],
    [AC_DEFINE(MOAB_FORCE_32_BIT_HANDLES,[],[Use int32_t for handles])
     SIZEOF_MBENTITYHANDLE=4
     DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --enable-32bit-handles"],
    [AC_MSG_ERROR([Platform does not have uint32_t. Cannot force 32-bit handles])])
])
if test "x$SIZEOF_MBENTITYHANDLE" = "x0"; then
  SIZEOF_MBENTITYHANDLE=$SIZEOF_VOID_P
  AC_CHECK_TYPE([size_t], 
   [AC_DEFINE(MOAB_HAVE_SIZE_T,[],[System provides size_t typedef])],
   [AC_CHECK_SIZEOF([unsigned long])
    SIZEOF_UNSIGNGED_LONG=$ac_cv_sizeof_unsigned_long
    if test "$SIZEOF_UNSIGNED_LONG" -ne "$SIZEOF_VOID_P"; then
      AC_MSG_ERROR([Cannot define handle with sizeof void*: size_t not found and long is of incorrect size.])
    fi]
  )
  AC_CHECK_TYPE([ptrdiff_t], 
   [AC_DEFINE(MOAB_HAVE_PTRDIFF_T,[],[System provides ptrdiff_t typedef])],
   [AC_CHECK_SIZEOF([long])
    AC_CHECK_SIZEOF([void *])
    SIZEOF_VOID_P=$ac_cv_sizeof_void_p
    SIZEOF_LONG=$ac_cv_sizeof_long
    if test "$SIZEOF_LONG" -ne "$SIZEOF_VOID_P"; then
      AC_MSG_ERROR([Cannot define handle with sizeof void*: ptrdiff_t not found and long is of incorrect size.])
    fi]
  )
fi


################################################################################
#                           iMesh Part 1 of 2
################################################################################

FORCE_IMESH=no
AC_ARG_ENABLE( [imesh],
[AC_HELP_STRING([--disable-imesh],[Do not build support for iMesh interface.])],
[FORCE_IMESH=yes],
[enableval=yes] )

case $enableval in
  yes)
    ENABLE_imesh=yes
    ;;
  nofortran|NoFortran|noFORTRAN|NoFORTRAN|Nofortran)
    AC_MSG_ERROR([Depricated value for --enable-imesh: \"$enableval\".  Try \"--disable-fortran\"])
    ;;
  no)
    ENABLE_imesh=no
    ;;
  *)
    AC_MSG_ERROR([Invalid argument to --enable-imesh : $enableval])
    ;;
esac

################################################################################
#                           iGeom Part 1 of 2
################################################################################

FORCE_IGEOM=no
AC_ARG_ENABLE( [igeom],
[AC_HELP_STRING([--enable-igeom],[Build mesh-based iGeom implementation])],
[FORCE_IGEOM=yes],[enableval=no] )

case $enableval in
  yes)
    ENABLE_igeom=yes
    ;;
  nofortran|NoFortran|noFORTRAN|NoFORTRAN|Nofortran)
    AC_MSG_ERROR([Depricated value for --enable-igeom: \"$enableval\".  Try \"FC= F77=\"])
    ;;
  no)
    ENABLE_igeom=no
    ;;
  *)
    AC_MSG_ERROR([Invalid argument to --enable-igeom : $enableval])
    ;;
esac

################################################################################
#                           iMesh/iGeom Part 2 of 2
################################################################################

# The iMesh option must be checked before configuring the compilers
# so that we know if we need Fortran.  This stuff has to be done
# after configuring the compilers so that we know what size of 
# various types are.

AC_MSG_CHECKING([size of EntityHandle])
AC_MSG_RESULT([$SIZEOF_MBENTITYHANDLE])
ITAPS_VALID_HANDLE_SIZE=yes
if test $SIZEOF_VOID_P -ne $SIZEOF_MBENTITYHANDLE; then
  ITAPS_VALID_HANDLE_SIZE=no
  AC_MSG_WARN([Cannot build iMesh or iGeom API because sizeof(EntityHandle) != sizeof(void*)])
  if test "x$FORCE_IMESH" = "xno"; then
    ENABLE_imesh=no
  fi
  if test "x$FORCE_IGEOM" = "xno"; then
    ENABLE_igeom=no
  fi
fi

AC_MSG_CHECKING([if iMesh support is to be built])
AC_MSG_RESULT([$ENABLE_imesh])
AM_CONDITIONAL([ENABLE_imesh],[test "xyes" = "x$ENABLE_imesh"])
if test "x$ENABLE_imesh" = "xyes"; then
  if test "x$ITAPS_VALID_HANDLE_SIZE" != "xyes"; then
    AC_MSG_ERROR([Cannot build iMesh API with sizeof(EntityHandle) != sizeof(void*)])
  fi
fi

AC_MSG_CHECKING([if iGeom support is to be built])
AC_MSG_RESULT([$ENABLE_igeom])
AM_CONDITIONAL([ENABLE_igeom],[test "xyes" = "x$ENABLE_igeom"])
if test "x$ENABLE_igeom" = "xyes"; then
  if test "x$ITAPS_VALID_HANDLE_SIZE" != "xyes"; then
    AC_MSG_ERROR([Cannot build iGeom API with sizeof(EntityHandle) != sizeof(void*)])
  fi
fi

################################################################################
#                              Version Stuff
################################################################################

AC_DEFINE(MB_VERSION,["AC_PACKAGE_VERSION"],[MOAB Version])
VERSION_MAJOR=`expr AC_PACKAGE_VERSION : '\([[0-9]]*\)'`
VERSION_MINOR=`expr AC_PACKAGE_VERSION : '[[0-9]]*\.\([[0-9]]*\)'`
VERSION_PATCH=`expr AC_PACKAGE_VERSION : '[[0-9]]*\.[[0-9]]*\.\(.*\)'`
test "x" != "x$VERSION_MAJOR" || AC_MSG_ERROR("Invalid version string: AC_PACKAGE_VERSION")
test "x" != "x$VERSION_MINOR" || AC_MSG_ERROR("Invalid version string: AC_PACKAGE_VERSION")
AC_DEFINE_UNQUOTED(MB_VERSION_MAJOR,$VERSION_MAJOR,[MOAB Major Version])
AC_DEFINE_UNQUOTED(MB_VERSION_MINOR,$VERSION_MINOR,[MOAB Minor Version])
if test "x" != "x$VERSION_PATCH"; then
  AC_DEFINE_UNQUOTED(MB_VERSION_PATCH,$VERSION_PATCH,[MOAB Patch Level])
  VERSION_STRING="AC_PACKAGE_STRING"
elif test $VERSION_MINOR -eq 99; then
  VERSION_STRING="AC_PACKAGE_STRING (Alpha)"
else
  VERSION_STRING="AC_PACKAGE_STRING (Beta)"
fi
AC_DEFINE_UNQUOTED(MB_VERSION_STRING,"${VERSION_STRING}",[MOAB Version String])

################################################################################
#                Create libtool script
################################################################################
# We need the libtool script both for the calls to the ITAPS_LIBTOOL_VAR
# macro and for us in testing for ACIS libraries.  If we're using a newer
# version of libtool, the script normally isn't generated until AC_OUTPUT
# at the end of the configure script.  Ask that it be created now instead.
m4_ifdef([LT_OUTPUT],[LT_OUTPUT])


################################################################################
#                              Extract libtool config
################################################################################
FATHOM_LIBTOOL_VAR( [CXX], [compiler_lib_search_path], [MOAB_CXX_LINKFLAGS])
FATHOM_LIBTOOL_VAR( [CXX], [postdeps], [MOAB_CXX_LIBS_TMP])
MOAB_CXX_LIBS=
for lib in $MOAB_CXX_LIBS_TMP; do
  case $lib in
    -lgcc_s)
      ;;
    -lm)
      ;;
    -lc)
      ;;
    *)
      MOAB_CXX_LIBS="$MOAB_CXX_LIBS $lib"
      ;;
  esac
done
AC_SUBST(MOAB_CXX_LINKFLAGS)
AC_SUBST(MOAB_CXX_LIBS)


################################################################################
#                             System Headers
################################################################################
AC_CHECK_HEADER([ieeefp.h],[AM_CPPFLAGS="$AM_CPPFLAGS -DHAVE_IEEEFP_H"] )
if test "xyes" = "x$enable_debug"; then
  AC_CHECK_HEADER( [valgrind/memcheck.h], [CPPFLAGS="$CPPFLAGS -DVALGRIND"] )
fi

################################################################################
#                             BOOST OPTIONS
################################################################################

dnl FATHOM_CHECK_BOOST( [boost/pool/object_pool.hpp] )
MK_CHECK_UNORDERED_MAP([CPPFLAGS="$CPPFLAGS -DUNORDERED_MAP_NS=$result -DHAVE_UNORDERED_MAP=${incdir}unordered_map -DHAVE_UNORDERED_SET=${incdir}unordered_set"])

################################################################################
#                              MPI OPTIONS
################################################################################

FATHOM_CHECK_MPI
FATHOM_CONFIG_MPI_EXEC
test "xyes" != "x$WITH_MPI" || AM_CPPFLAGS="$AM_CPPFLAGS -DUSE_MPI"
AM_CONDITIONAL(USE_MPI, [test "xno" != "x$WITH_MPI"])
# Used to generate moab_mpi.h
if test "xyes" = "x$MPI_CXX_HELP_NEEDED"; then
  AC_DEFINE( [MB_MPI_CXX_CONFLICT], [1], 
    [MPICH_IGNORE_CXX_SEEK is not sufficient to avoid conflicts] )
  AC_DEFINE_UNQUOTED([MB_SEEK_SET],[$SEEK_SET],["Value of C SEEK_SET"])
  AC_DEFINE_UNQUOTED([MB_SEEK_CUR],[$SEEK_CUR],["Value of C SEEK_CUR"])
  AC_DEFINE_UNQUOTED([MB_SEEK_END],[$SEEK_END],["Value of C SEEK_END"])
fi

if test "xyes" = "x$WITH_MPI"; then
  FATHOM_CHECK_MPE
  AC_SUBST([WITH_MPE])
  AC_SUBST([MPE_LIBS])
  test "xyes" != "x$WITH_MPE" || AM_CPPFLAGS="$AM_CPPFLAGS -DUSE_MPE"
fi

################################################################################
#                              HDF5 OPTIONS
################################################################################

FATHOM_CHECK_HDF5
if test "xno" = "x$HAVE_HDF5"; then
  AC_MSG_WARN([Support for native HDF5 file format disabled])
else
  AM_CPPFLAGS="$AM_CPPFLAGS -DHDF5_FILE"
fi
AM_CONDITIONAL(HDF5_FILE, [test "xno" != "x$HAVE_HDF5"])
AM_CPPFLAGS="$HDF5_CPPFLAGS $AM_CPPFLAGS"
EXPORT_LDFLAGS="$EXPORT_LDFLAGS $HDF5_LDFLAGS"
AC_SUBST(HDF5_LIBS)

WARN_PARALLEL_HDF5=no
WARN_PARALLEL_HDF5_NO_COMPLEX=no
HAVE_HDF5_PARALLEL=no
if test "xno" != "x$HAVE_HDF5"; then
  if test "xno" != "x$WITH_MPI"; then
    old_LDFLAGS="$LDFLAGS"
    LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
    AC_CHECK_LIB( [hdf5], [H5Pset_fapl_mpio], [HAVE_HDF5_PARALLEL=yes],
      [WARN_PARALLEL_HDF5=yes; AC_MSG_WARN("libhdf5 library does not include parallel support.  Parallel HDF5 I/O disabled")],
      [$HDF5_LIBS]
      )
    LDFLAGS="$old_LDFLAGS"
  fi
fi
AM_CONDITIONAL(PARALLEL_HDF5, [test "xno" != "x$HAVE_HDF5_PARALLEL"])
if test "xno" != "x$HAVE_HDF5_PARALLEL"; then
  AM_CPPFLAGS="$AM_CPPFLAGS -DHDF5_PARALLEL"
  
  AC_MSG_CHECKING([for H5_MPI_COMPLEX_DERIVED_DATATYPE_WORKS])
  old_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS"
  AC_PREPROC_IFELSE([AC_LANG_PROGRAM([#include <H5pubconf.h>],[
#ifndef H5_MPI_COMPLEX_DERIVED_DATATYPE_WORKS
  choke me
#endif])],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no); WARN_PARALLEL_HDF5_NO_COMPLEX=yes])
  CPPFLAGS="$old_CPPFLAGS"
fi

################################################################################
#                             CCMIO OPTIONS
################################################################################

old_LDFLAGS="$LDFLAGS"
LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
FATHOM_CHECK_CCMIO
LDFLAGS="$old_LDFLAGS"
if test "xno" != "x$HAVE_CCMIO"; then
  AM_CPPFLAGS="$AM_CPPFLAGS -DCCMIO_FILE"
fi
AM_CONDITIONAL(CCMIO_FILE, [test "xno" != "x$HAVE_CCMIO"])
AM_CPPFLAGS="$CCMIO_CPPFLAGS $AM_CPPFLAGS"
EXPORT_LDFLAGS="$CCMIO_LDFLAGS $EXPORT_LDFLAGS"
AC_SUBST(CCMIO_LIBS)


################################################################################
#                             DAMSEL OPTIONS
################################################################################

old_LDFLAGS="$LDFLAGS"
LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
FATHOM_CHECK_DAMSEL
LDFLAGS="$old_LDFLAGS"
if test "xno" != "x$HAVE_DAMSEL"; then
  AM_CPPFLAGS="$AM_CPPFLAGS -DDAMSEL_FILE"
fi
AM_CONDITIONAL(DAMSEL_FILE, [test "xno" != "x$HAVE_DAMSEL"])
AM_CPPFLAGS="$DAMSEL_CPPFLAGS $AM_CPPFLAGS"
EXPORT_LDFLAGS="$DAMSEL_LDFLAGS $EXPORT_LDFLAGS"
AC_SUBST(DAMSEL_LIBS)


################################################################################
#                             NetCDF OPTIONS
################################################################################

old_LDFLAGS="$LDFLAGS"
LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
FATHOM_CHECK_NETCDF([65536],[524288])
LDFLAGS="$old_LDFLAGS"
if test "xyes" = "x$HAVE_NETCDF"; then
  AM_CPPFLAGS="$AM_CPPFLAGS -DNETCDF_FILE"
fi
AM_CONDITIONAL(NETCDF_FILE, [test "xno" != "x$HAVE_NETCDF"])
AM_CPPFLAGS="$NETCDF_CPPFLAGS $AM_CPPFLAGS"
EXPORT_LDFLAGS="$NETCDF_LDFLAGS $EXPORT_LDFLAGS"
AC_SUBST(NETCDF_LIBS)

  # now pnetcdf
FATHOM_CHECK_PNETCDF([65536],[524288])
if test "xyes" = "x$HAVE_PNETCDF"; then
  AM_CPPFLAGS="$AM_CPPFLAGS -DPNETCDF_FILE"
fi
AM_CONDITIONAL(PNETCDF_FILE, [test "xyes" = "x$HAVE_PNETCDF"])
AM_CPPFLAGS="$PNETCDF_CPPFLAGS $AM_CPPFLAGS"
EXPORT_LDFLAGS="$PNETCDF_LDFLAGS $EXPORT_LDFLAGS"
AC_SUBST(PNETCDF_LIBS)


#################################################################################
#                             Documentation
#################################################################################
AC_ARG_ENABLE([doxygen],
[AC_HELP_STRING([[--enable-doxygen@<:@=DIR@:>@]],[Specify directory where Doxygen program is installed])
AC_HELP_STRING([--disable-doxygen],[Do not generate API documentation (default)])],
                        [ENABLE_DOXYGEN="$enableval"],[ENABLE_DOXYGEN=no] )
if test "x$ENABLE_DOXYGEN" = "xyes"; then
  AC_PATH_PROG( [DOXYGEN], [doxygen], [no] )
elif test "x$ENABLE_DOXYGEN" != "xno"; then
  AC_PATH_PROG( [DOXYGEN], [doxygen], [no], [$ENABLE_DOXYGEN] )
fi
if test "x$ENABLE_DOXYGEN" != "xno"; then
  if test "x$DOXYGEN" = "xno"; then
    AC_MSG_ERROR("Doxygen executable not found.")
  fi
fi
AC_SUBST([DOXYGEN])
AM_CONDITIONAL([ENABLE_DOXYGEN],[test "x$ENABLE_DOXYGEN" != "xno"])


################################################################################
#                           Zoltan - Part 1 of 2
################################################################################
WITH_ZOLTAN=no
ZOLTAN_DIR=
ZOLTAN_ARCH=
WITH_PARAMETIS=no
PARMETIS_DIR=

ZOLTAN_LIBS=
ZOLTAN_LIB_FLAGS=
ZOLTAN_INC_FLAGS=

# If user specified zoltan-arch, but not zoltan, assume
# they want zoltan support, otherwise assume no zoltan
ZOLTAN_DEFAULT=no

# Some versions of Zoltan put libraries in directories such as Obj_Linux,
# where the part after the Obj_ prefix is the zoltan architecture.  Here
# we allow the user to specify this architecture.  We try to avoid requiring
# this option, but we allow the user to specify it in case something is 
# going wrong with the logic to choose one automatically.  The 'choose one
# automatically' logic is basically to add everything that matches Obj_*
# to the search path, with the assumption that the linker will ignore libraries
# for other architectures.
AC_ARG_WITH( [zoltan-arch],
             [AC_HELP_STRING([--with-zoltan-arch=OS],[Specify zoltan architecture])],
             [case "x$withval" in
               xyes|xno)
                 AC_MSG_ERROR([Expeted value for --with-zoltan-arch])
                 ;;
               x*)
                 ZOLTAN_ARCH="$withval"
                 ZOLTAN_DEFAULT=yes
                 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-zoltan-arch=\"${withval}\""
                 ;;
               esac])

# Enable support for Zoltan, and optionally specify directory in which
# to search for Zoltan libs and headers.  Specifying this also enables
# mbzoltan unless it is explicitly disabled.
AC_ARG_WITH( [zoltan], 
             [AC_HELP_STRING([--with-zoltan=DIR],[Enable Zoltan support, and optionally specify Zoltan location])],
             [case "x$withval" in
               xyes)
                 WITH_ZOLTAN=yes
                 ;;
               xno)
                 WITH_ZOLTAN=no
                 ;;
               x*)
                 WITH_ZOLTAN=yes
                 ZOLTAN_DIR="$withval"
                 test -d "$withval" || AC_MSG_ERROR([Not a directory: with-zoltan=$withval])
                 if test -d "$withval/include"; then
                   ZOLTAN_INC_FLAGS="-I$withval/include"
                 else
                   ZOLTAN_INC_FLAGS="-I$withval"
                 fi
                 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-zoltan=\"${withval}\""
                 ;;
               esac],[WITH_ZOLTAN=$ZOLTAN_DEFAULT] )

# If user specified with-parmetis, then assume they are trying
# to add an additional library search directory in which to check
# for the parmetis library.  We always check if zoltan requires
# parmetis, so the only thing the user would really need to tell
# us is if we need to look in some special place for it. 
AC_ARG_WITH([parmetis], 
            [AC_HELP_STRING([--with-parmetis=DIR],[Directory containing Parmetis library for Zoltan])],
             [case "x$withval" in
               xyes|xno)
                 # do nothing
                 ;;
               x*)
                 PARMETIS_DIR="$withval"
                 test -d "$withval" || AC_MSG_ERROR([Not a directory: with-parmetis=$withval])
                 if test -d "$withval/lib"; then
                   ZOLTAN_LIB_FLAGS="-L$withval/lib"
                 else
                   ZOLTAN_LIB_FLAGS="-L$withval"
                 fi
                 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-parmetis=\"${withval}\""
		 ;;
               esac] )



################################################################################
#                            Optional Libraries
################################################################################

AC_ARG_WITH([vtk],
  [AC_HELP_STRING([--with-vtk@<:@=DIR@:>@], [Enable VTK support, and specify VTK location])],
  [WITH_VTK=$withval
   DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-vtk=\"${withval}\""
  ],[WITH_VTK=])

VTKMOAB_DEFAULT=no
if test "x$WITH_VTK" != "xno" -a "x$WITH_VTK" != "x"; then
  VTKMOAB_DEFAULT=yes
fi
  

###############################################################################
#                           Optional Tools
###############################################################################
  # Define a macro to avoid typing this for each individual tool
  # Usage: MB_OPTIONAL_TOOL( name, default )
  #  name      - name of option
  #  default   - yes/no, $ENABLE_TOOLS overrides this if set
  # Actions:
  #  sets ENABLE_${tool} to 'yes' or 'no'
  #  creates ENABLE_${tool} automake conditional
AC_DEFUN([MB_OPTIONAL_TOOL],[
  mb_default_$1=$2
  if test "x" != "x$ENABLE_TOOLS"; then
    mb_default_$1="$ENABLE_TOOLS"
  fi
  AC_ARG_ENABLE( [$1],
[AC_HELP_STRING([--enable-$1],[Build tool: $1])
AC_HELP_STRING([--disable-$1],[Don't build $1])],
                 [ENABLE_$1=$enableval
DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --enable-$1=\"${enableval}\""
],[ENABLE_$1=${mb_default_$1}] )
  AM_CONDITIONAL([ENABLE_$1],[test "x${ENABLE_$1}" != "xno"])
  AC_MSG_CHECKING([if $1 is to be built])
  AC_MSG_RESULT([${ENABLE_$1}])
])
  # The default for all tools
AC_ARG_ENABLE( [tools], [
AC_HELP_STRING([--enable-tools],[Build all tools by default])
AC_HELP_STRING([--disable-tools],[Disable all tools by default])],
                         [ENABLE_TOOLS=$enableval],[ENABLE_TOOLS=] )
  # Individual tools
MB_OPTIONAL_TOOL([mbconvert],    [yes])
MB_OPTIONAL_TOOL([hexmodops],    [yes])
MB_OPTIONAL_TOOL([qvdual],       [no] )
MB_OPTIONAL_TOOL([vtkMOABReader],[${VTKMOAB_DEFAULT}] )
MB_OPTIONAL_TOOL([mbsize],       [yes])
MB_OPTIONAL_TOOL([mbskin],       [yes])
MB_OPTIONAL_TOOL([mbtagprop],    [yes])
MB_OPTIONAL_TOOL([mbmem],        [yes])
MB_OPTIONAL_TOOL([mcnpmit],      [no])
MB_OPTIONAL_TOOL([mbcoupler],    [$ENABLE_mcnpmit] )
MB_OPTIONAL_TOOL([spheredecomp], [yes])
MB_OPTIONAL_TOOL([mbsurfplot],   [yes])
MB_OPTIONAL_TOOL([mbzoltan],     [${WITH_ZOLTAN}] )
MB_OPTIONAL_TOOL([dagmc],        [yes])
MB_OPTIONAL_TOOL([gsets],        [yes])
MB_OPTIONAL_TOOL([mbdepth],      [yes])
MB_OPTIONAL_TOOL([refiner],      [no])
MB_OPTIONAL_TOOL([h5mtools],     [yes])

if test "xyes" = "x$ENABLE_refiner"; then
  if test "xyes" != "x$WITH_MPI"; then
    AC_MSG_ERROR([refiner requires MPI/parallel configuration])
  fi
fi

if test "xyes" = "x$ENABLE_vtkMOABReader"; then
  if test "xno" == "x$WITH_VTK"; then
    AC_MSG_ERROR([vtkMOABReader requires VTK])
  elif test "x" == "x$WITH_VTK"; then 
    WITH_VTK=yes
  fi
fi

if test "xyes" = "x$ENABLE_qvdual"; then
  if test "xno" == "x$WITH_VTK"; then
    AC_MSG_ERROR([qvdual requires VTK])
  elif test "x" == "x$WITH_VTK"; then 
    WITH_VTK=yes
  fi
fi

################################################################################
#                           Zoltan - Part 2 of 2
################################################################################

if test "xyes" = "x$ENABLE_mbzoltan"; then
  WITH_ZOLTAN=yes
else
  WITH_ZOLTAN=no
fi

# If user requested Zoltan support, verify that we can
# find the library and headers, and set the necessary
# build variables.
if test "xyes" = "x$WITH_ZOLTAN"; then

  if test "xyes" != "x$WITH_MPI"; then
    AC_MSG_ERROR([Cannot include Zoltan support without MPI.  See --with-mpi])
  fi

    # Test for zoltan headers
  AC_LANG_PUSH([C++])
  old_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $ZOLTAN_INC_FLAGS -DMPICH_IGNORE_CXX_SEEK"
  AC_CHECK_HEADER([zoltan_cpp.h],,[AC_MSG_ERROR([zoltan_cpp.h: header not found in $ZOLTAN_INC_FLAGS .])],
                  [#include <stdlib.h>])
  CPPFLAGS="$old_CPPFLAGS"
  AC_LANG_POP([C++])
  
    # Build LDFLAGS for Zoltan
  found=no
  if test "x" != "x$ZOLTAN_DIR"; then
    if test "x" != "x$ZOLTAN_ARCH" && test -d "$ZOLTAN_DIR/Obj_$ZOLTAN_ARCH"; then
      ZOLTAN_LIB_FLAGS="$ZOLTAN_LIB_FLAGS -L$ZOLTAN_DIR/Obj_$ZOLTAN_ARCH"
      found=yes
    elif test -d "$ZOLTAN_DIR/lib"; then
      ZOLTAN_LIB_FLAGS="$ZOLTAN_LIB_FLAGS -L$ZOLTAN_DIR/lib"
      found=yes
    else
      for d in $ZOLTAN_DIR/Obj_*; do
        if test -d "$d"; then
          found=yes
          ZOLTAN_LIB_FLAGS="$ZOLTAN_LIB_FLAGS -L$d"
        fi
      done
    fi
    if test "no" = "$found"; then
      ZOLTAN_LIB_FLAGS="$ZOLTAN_LIB_FLAGS -L$ZOLTAN_DIR"
    fi
  fi
  
    # Test for Zoltan library, possibly requiring parmetis
  old_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $ZOLTAN_LIB_FLAGS"
  AC_CHECK_LIB([zoltan],[Zoltan_LB_Partition],[ZOLTAN_LIBS='-lzoltan'],
    [unset ac_cv_lib_zoltan_Zoltan_LB_Partition
     AC_CHECK_LIB([zoltan],[Zoltan_LB_Partition],
                  [ZOLTAN_LIBS="-lzoltan -lparmetis -lmetis"],
                  [AC_MSG_ERROR([zoltan library not found or not usable.  Perhaps you need to specify --with-parmetis=@<:@DIR@:>@?])],
                  [-lparmetis -lmetis])])
  LDFLAGS="$old_LDFLAGS"
  
fi

AC_SUBST(ZOLTAN_LIBS)
AC_SUBST(ZOLTAN_LIB_FLAGS)
AC_SUBST(ZOLTAN_INC_FLAGS)

################################################################################
#                    VTK for qvdual and/or vtkMOABReader
################################################################################
AC_DEFUN([FATHOM_MERGE_LISTS],[
  $1="$2"
  # Reverse list two
  for item3 in $3; do 
    found=no
    for item2 in $2; do
      if test "x$item2" = "x$item3"; then found=yes; fi
    done
    if test "$found" = "no"; then $1="${$1} $item3"; fi
  done
  $1="${$1} $2"
])

if test "xno" != "x$WITH_VTK" -a "x" != "x$WITH_VTK"; then

  AC_PATH_X
  AC_PATH_XTRA
  if test "xyes" = "x$no_x"; then
    AC_MSG_ERROR([X11 not found])
  fi

  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS

  VTK_INCLUDES=  
  case "x$WITH_VTK" in
    xyes|x)
      VTK_LIB_DIR=
      VTK_INC_DIR=
      for dir in  /usr/local/include/vtk /usr/local/include/vtk-* /usr/include/vtk  /usr/include/vtk-*; do
        if test "x" = "x$VTK_INC_DIR"; then
          AC_CHECK_FILE([${dir}/vtkRenderingInstantiator.h],[VTK_INC_DIR="-I$dir"])
        fi
      done
      if test x = x"$VTK_INC_DIR"; then
        AC_CHECK_HEADER([vtkRenderingInstantiator.h],[],
                        [AC_MSG_ERROR([VTK includes not found])])
      else
        VTK_INCLUDES="$VTK_INC_DIR"
      fi
      old_LIBS="$LIBS"
      for dir in /usr/lib/vtk /usr/lib /usr/bin /usr; do
        if test "x" = "x$VTK_LIB_DIR"; then
          LIBS="$old_LIBS -L$dir"
          unset ac_cv_lib_vtkCommon
	  unset ac_cv_lib_vtkCommon_main
          AC_HAVE_LIBRARY([vtkCommon],[VTK_LIB_DIR=$dir])
        fi
      done
      if test "x" = "x$VTK_LIB_DIR"; then
        AC_MSG_ERROR([VTK libraries not found])
      fi
      LIBS="$old_LIBS"
      ;;
    x?*)
      VTK_LIB_DIR=
      VTK_INC_DIR=
      for dir in ${WITH_VTK}/include/vtk ${VTK_DIR}/include/vtk/Rendering ${WITH_VTK}/include ${WITH_VTK}/include/paraview-* ${WITH_VTK}/VTK ${WITH_VTK}; do
        if test "x" = "x$VTK_INC_DIR"; then
          AC_CHECK_FILE([${dir}/vtkRenderingInstantiator.h],[VTK_INC_DIR="-I$dir"])
        fi
      done
      if test "x" = "x$VTK_INC_DIR"; then
        # Check for non-installed VTK source
        AC_CHECK_FILE([${VTK_DIR}/Rendering/vtkRenderingInstantiator.h],[
          for subdir in . Common Filtering GenericFiltering Graphics Hybrid Imaging Parallel Patented Rendering Utilities Utilities/vtkexpat Utilities/expat Utilities/DICOMParser Utilties/vtkpng Utilities/vtkzlib Utilities/freetype/include/freetype Wrapping; do
  	  AC_MSG_CHECKING([for directory ${WITH_VTK}/$subdir])
  	  if test -d ${WITH_VTK}/$subdir; then
  	    VTK_INC_DIR="$VTK_INC_DIR -I${WITH_VTK}/$subdir"
  	    AC_MSG_RESULT([yes])
  	  else
  	    AC_MSG_RESULT([no])
  	  fi
  	done
  	],[AC_MSG_ERROR([VTK includes not found])]
        )
      fi
      VTK_INCLUDES="$VTK_INC_DIR"
      old_LIBS="$LIBS"
      for dir in ${WITH_VTK}/lib/vtk ${WITH_VTK}/lib ${VTK_DIR}/bin ${WITH_VTK}/bin ${WITH_VTK}; do
        if test "x" = "x$VTK_LIB_DIR"; then
          LIBS="$old_LIBS -L$dir"
          unset ac_cv_lib_vtkCommon
  	unset ac_cv_lib_vtkCommon_main
          AC_HAVE_LIBRARY([vtkCommon],[VTK_LIB_DIR=$dir])
        fi
      done
      if test "x" = "x$VTK_LIB_DIR"; then
        AC_MSG_ERROR([VTK libraries not found])
      fi
      LIBS="$old_LIBS"
      ;;
  esac
  
dnl Define a function that checks for the existence of a library
dnl and if it exists, prepends the corresponding link flag to a 
dnl sepecified variable
dnl Arguments:
dnl  - The library name
dnl  - The name of the variable to which to prepend the library link flag
dnl  - List of extra libraries required for link test
dnl  - Optional library path linker argument (e.g. -L/opt/foo)
  AC_DEFUN([FATHOM_CHECK_LIB_LIST],[
    old_LIBS="$LIBS"
    LIBS="$4 ${$2} $old_LIBS"
    AC_HAVE_LIBRARY([$1],[$2="-l$1 ${$2}"],[],[$3])
    LIBS="$old_LIBS"
  ])
  
   # These are ordered such that a given library depends only on those libraries
   # that preceed it in the list.
  FATHOM_CHECK_LIB_LIST([pthread],[vtk_sys_libs])
  FATHOM_CHECK_LIB_LIST([dl],     [vtk_sys_libs])
  FATHOM_CHECK_LIB_LIST([m],      [vtk_sys_libs])
  FATHOM_CHECK_LIB_LIST([gcc],    [vtk_sys_libs])
  x_libs_arg="$X_EXTRA_LIBS -lX11 $X_PRE_LIBS"
  FATHOM_CHECK_LIB_LIST([Xext],[vtk_x_libs],[$x_libs_arg],[-L$x_libraries]) 
  FATHOM_CHECK_LIB_LIST([ICE], [vtk_x_libs],[$x_libs_arg],[-L$x_libraries]) 
  FATHOM_CHECK_LIB_LIST([SM],  [vtk_x_libs],[$x_libs_arg],[-L$x_libraries]) 
  FATHOM_CHECK_LIB_LIST([Xt],  [vtk_x_libs],[$x_libs_arg],[-L$x_libraries]) 
  FATHOM_CHECK_LIB_LIST([GL],  [vtk_x_libs],[$x_libs_arg],[-L$x_libraries]) 
  FATHOM_CHECK_LIB_LIST([GLU], [vtk_x_libs],[$x_libs_arg],[-L$x_libraries]) 
  vtk_libs_arg="-L$x_libraries $vtk_x_libs $X_EXTRA_LIBS -lX11 $X_PRE_LIBS $vtk_sys_libs"
  FATHOM_CHECK_LIB_LIST([vtksys],        [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkCommon],     [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkDICOMParser],[vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkfreetype],   [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkexpat],      [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkzlib],       [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkpng],        [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkjpeg],       [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtktiff],       [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkftgl],       [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkFiltering],  [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkImaging],    [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([verdict],       [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkGraphics],   [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkNetCDF],     [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkexoIIc],     [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkmetaio],     [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkIO],         [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkRendering],  [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkParallel],   [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkHybrid],     [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  
   # Make sure we found at least vtkCommon and vtkRendering
  if ! echo "x $vtk_libs" | grep lvtkCommon >/dev/null; then
    AC_MSG_ERROR("VTK library 'vtkCommon' not found.")
  fi
  if ! echo "x $vtk_libs" | grep lvtkRendering >/dev/null; then
    AC_MSG_ERROR("VTK library 'vtkRendering' not found.")
  fi
  
  AC_LANG_RESTORE

  FATHOM_MERGE_LISTS([XLIBS_tmp],[$X_EXTRA_LIBS -lX11 $X_PRE_LIBS],[$vtk_x_libs])
  FATHOM_MERGE_LISTS([XLIBS],[$vtk_sys_libs],[$XLIBS_tmp])
  VTK_LIBS="$LIBS -L$VTK_LIB_DIR $vtk_libs $XLIBS"
  VTK_INCLUDES="$VTK_INCLUDES"
  AC_SUBST(VTK_LIBS)
  AC_SUBST(VTK_INCLUDES)
  AC_SUBST(VTK_LIB_DIR)
  AC_SUBST(VTK_INC_DIR)
fi

################################################################################
#                               qvdual
################################################################################

if test $ENABLE_qvdual = yes; then


################################################################################
#                               GraphViz
################################################################################
  AC_ARG_WITH([graphviz],
  [AC_HELP_STRING([--with-graphviz=DIR], [Specify directory containing GraphViz])],
  [GRAPH_DIR=$withval
   DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-graphviz=\"${withval}\""
  ],[GRAPH_DIR=])

  case "x$GRAPH_DIR" in
    xno)
      GRAPHVIZ_MISSING=yes
      ;;
    xyes|x)
      GRAPH_INC_DIR=
      GRAPH_DIR=
      for dir in /usr /usr/local; do
        if test "x" = "x$GRAPH_DIR"; then
          AC_CHECK_FILE([${dir}/include/graphviz/gvc.h],[GRAPH_DIR="$dir"])
        fi
      done
      AC_MSG_CHECKING([for GraphViz include dir])
      if test "x" != x"$GRAPH_DIR"; then
        GRAPH_INC_DIR="-I${GRAPH_DIR}/include/graphviz"
      fi
      AC_MSG_RESULT([$GRAPH_INC_DIR])
      AC_MSG_CHECKING([for GraphViz library dir])
      GRAPH_LIB_DIR=
      for dir in ${GRAPH_DIR}/lib64/graphviz ${GRAPH_DIR}/lib/graphviz /usr/lib64/graphviz /usr/lib/graphviz /usr/local/lib64/graphviz /usr/local/lib/graphviz /lib64/graphviz /lib/graphviz; do
        if test "x" = "x$GRAPH_LIB_DIR"; then
          if test -d $dir; then
            GRAPH_LIB_DIR=-L$dir
          fi
        fi
      done
      AC_MSG_RESULT([$GRAPH_LIB_DIR])
      ;;
    x?*)
      GRAPH_INC_DIR=
      for dir in ${GRAPH_DIR}/include/graphviz ${GRAPH_DIR}/include ${GRAPH_DIR}; do
        if test "x" = "x$GRAPH_INC_DIR"; then
          AC_CHECK_FILE([${dir}/gvc.h],[GRAPH_INC_DIR="-I$dir"])
        fi
      done
      AC_MSG_CHECKING([for GraphViz include dir])
      if test "x" = x"$GRAPH_INC_DIR"; then
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([Graphviz includes for qvdual not found])
      else
        AC_MSG_RESULT([$GRAPH_INC_DIR])
      fi
      AC_MSG_CHECKING([for GraphViz library dir])
      GRAPH_LIB_DIR=
      for dir in ${GRAPH_DIR}/lib/graphviz ${GRAPH_DIR}/lib ${GRAPH_DIR}; do
        if test "x" = "x$GRAPH_LIB_DIR"; then
          if test -d $dir; then
            GRAPH_LIB_DIR=-L$dir
          fi
        fi
      done
      if test "x" = x"$GRAPH_LIB_DIR"; then
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([Graphviz lib for qvdual not found])
      else
        AC_MSG_RESULT([$GRAPH_LIB_DIR])
      fi
      ;;
  esac

  GRAPH_LIBS="$GRAPH_LIB_DIR -lgvc -lgraph"
  old_CXXFLAGS="$CXXFLAGS"
  old_CPPFLAGS="$CPPFLAGS"
  CXXFLAGS="$GRAPH_INC_DIR"
  CPPFLAGS="$GRAPH_INC_DIR"
  old_LIBS="$LIBS"
  LIBS=$GRAPH_LIBS
  AC_CHECK_HEADERS([gvc.h],[],[AC_MSG_WARN([Missing GraphViz header]); GRAPHVIZ_MISSING=yes])
  AC_CHECK_LIB([gvc],[gvContext],[],[AC_MSG_WARN([Cannot find GraphViz library: -lgraph]); GRAPHVIZ_MISSING=yes])
  LIBS="$old_LIBS"
  CXXFLAGS="$old_CXXFLAGS"
  CPPFLAGS="$old_CPPFLAGS"
  
  
################################################################################
#                                   QT
################################################################################
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_ARG_VAR([QTDIR],[QT installation directory])
  if test x"$QTDIR" = x; then
    search_path="$PATH"
    qlibdir=
    incdir=
    for dir in /usr/include/qt3 /usr/local/include/qt3 /include/qt3; do
      if test x = x${incdir}; then
        if test -d $dir; then
          incdir=-I$dir
        fi
      fi
    done
  else
    search_path="${QTDIR}/bin"
    qlibdir="-L${QTDIR}/lib"
    incdir="-I${QTDIR}/include"
  fi
  
  AC_PATH_PROG([MOC],[moc],[],[$search_path])
  if test x = x"$MOC";then
    AC_MSG_WARN([Cannot find moc (Qt meta object compiler)])
    QT_MISSING=yes
  fi
  AC_SUBST(MOC)
  AC_PATH_PROG([UIC],[uic],[],[$search_path])
  if test x = x"$UIC";then
    AC_MSG_WARN([Cannot find uic (Qt user interface compiler)])
    QT_MISSING=yes
  fi
  AC_SUBST(UIC)
  
  old_CXXFLAGS="$CXXFLAGS"
  old_CPPFLAGS="$CPPFLAGS"
  CXXFLAGS="$incdir"
  CPPFLAGS="$incdir"
  AC_CHECK_HEADERS([qlineedit.h qevent.h qapplication.h \
                    qwidget.h qtimer.h qpixmap.h qwidgetplugin.h \
                    qobject.h qmetaobject.h],[],
                   [QT_MISSING=yes])
  CXXFLAGS="$old_CXXFLAGS"
  CPPFLAGS="$old_CPPFLAGS"
  QT_INCLUDES="$incdir"
  
  old_LIBS="$LIBS"
  LIBS="$LIBS $qlibdir"
  # Cannot use AC_HAVE_LIBRARY for qt-mt because of
  # autoconf bug with libraries that have a '-' in their name.
  #AC_HAVE_LIBRARY([qt-mt],[QT_LIBS="$qlibdir -lqt-mt"],[
  AC_MSG_CHECKING([for main in -lqt-mt])
  QT_LIBS="$qlibdir -lqt-mt"
  AC_TRY_LINK_FUNC(main,[AC_MSG_RESULT([yes])],[
    AC_MSG_RESULT([no])
    LIBS="$old_LIBS $qlibdir"
    AC_HAVE_LIBRARY([qt],[QT_LIBS="$qlibdir -lqt"],[
      QT_MISSING=yes
    ])
  ])
  LIBS="$VTK_LIBS $QT_LIBS"
  AC_HAVE_LIBRARY([QVTK],[QT_LIBS="-lQVTL $QT_LIBS"],
                  [AC_MSG_ERROR([QVTK library for qvdual not found])])
  LIBS="$old_LIBS"
  
################################################################################
#                                 X/VTK/QT
################################################################################

  QVDUAL_LIBS="$LIBS $GRAPH_LIBS $QT_LIBS $VTK_LIBS"
  QVDUAL_INCLUDES="$GRAPH_INC_DIR $VTK_INCLUDES $QT_INCLUDES"
  AC_SUBST(QVDUAL_LIBS)
  AC_SUBST(QVDUAL_INCLUDES)
fi

################################################################################
#                           CGM
################################################################################
AC_ARG_WITH(cgm, 
[AC_HELP_STRING([--with-cgm=DIR], [Specify directory containing CGM])],
[CGM_DIR=$withval
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-cgm=\"${withval}\""
 ], [CGM_DIR=])

case "x$CGM_DIR" in
  xno|x)
    CGM_MISSING=yes
    ;;
  xyes)
    AC_MSG_ERROR([--with-cgm requires an argument])
    ;;
  *)
    if ! test -d "$CGM_DIR"; then
      AC_MSG_ERROR([$CGM_DIR : not a directory.])
    fi
    CGM_MISSING=no
    CGM_MAKE=
    AC_CHECK_FILE([${CGM_DIR}/lib/cgm.make],
                  [CGM_MAKE="${CGM_DIR}/lib/cgm.make";
                   AM_CPPFLAGS="$AM_CPPFLAGS -DCGM"],
                  [AC_CHECK_FILE([${CGM_DIR}/cgm.make],
                                 [CGM_MAKE="${CGM_DIR}/cgm.make";
                   AM_CPPFLAGS="$AM_CPPFLAGS -DCGM"],
                                 [AC_MSG_ERROR([$CGM_DIR : not a configured CGM]); 
                                  CGM_MISSING=yes])
                  ])
    ;;
esac

if test "x$CGM_MISSING" = "xno"; then
  CGM_CONFIG_OPTIONS="include $CGM_MAKE"

  AC_MSG_CHECKING([Checking for GeometryQueryTool::ray_fire(RefFace*,...)])
  FATHOM_MAKE_INC_VAR([$CGM_MAKE],[CGM_DEFINES], [CGM_CPPFLAGS="$make_val"])   
  FATHOM_MAKE_INC_VAR([$CGM_MAKE],[CGM_INCLUDES],[CGM_CPPFLAGS="$CGM_CPPFLAGS $make_val"])   
  FATHOM_MAKE_INC_VAR([$CGM_MAKE],[CGM_LTFLAGS], [CGM_LTFLAGS="$make_val"])   
  FATHOM_MAKE_INC_VAR([$CGM_MAKE],[CGM_LDFLAGS], [CGM_LDFLAGS="$make_val"])   
  FATHOM_MAKE_INC_VAR([$CGM_MAKE],[CGM_LIBS],    [CGM_LIBS="$make_val"])   
  AC_LANG_PUSH([C++])
  old_CPPFLAGS="$CPPFLAGS"
  old_LDFLAGS="$LDFLAGS"
  old_LIBS="$LIBS"
  CPPFLAGS="$CGM_CPPFLAGS $CPPFLAGS"
  LDFLAGS="$CGM_LTFLAGS $CGM_LDFLAGS $CPPFLAGS"
  LIBS="$CGM_LIBS $LIBS"
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include "GeometryQueryTool.hpp"
                      #include "CubitVector.hpp"],
                    [CubitVector a,b;
                     DLIList<double> c;
                     GeometryQueryTool::instance()->fire_ray((RefFace*)0, a, b, c );
                    ])],
    [MOAB_CGM_DEFINES="-DHAVE_CGM_FIRE_RAY"; AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)])
  CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$old_LDFLAGS"
  LIBS="$old_LIBS"
#  LDFLAGS="$old_LDFLAGS"
#  LIBS="$old_LIBS"
  AC_LANG_POP([C++])
  EXPORT_LDFLAGS="$EXPORT_LDFLAGS $CGM_LDFLAGS"
  EXPORT_LTFLAGS="$EXPORT_LTFLAGS $CGM_LTFLAGS"
fi

AC_SUBST(MOAB_CGM_DEFINES)
AC_SUBST(CGM_CPPFLAGS)
AC_SUBST(CGM_LDFLAGS)
AC_SUBST(CGM_LTFLAGS)
AC_SUBST(CGM_LIBS)

AM_CONDITIONAL( HAVE_CGM, [test "x$CGM_MISSING" = "xno"] )

################################################################################
#                           Output Files
################################################################################

# Some old versions of autoconf don't define docdir.  
# Define it if autoconf did not.
if test "x" = "x$docdir"; then
  docdir='${datadir}/doc/moab'
  AC_SUBST(docdir)
fi

AC_SUBST([AM_CPPFLAGS])
AM_LDFLAGS="$AM_LDFLAGS $EXPORT_LTFLAGS $EXPORT_LDFLAGS"
AC_SUBST([AM_LDFLAGS])
AC_SUBST([EXPORT_LTFLAGS])
AC_SUBST([EXPORT_LDFLAGS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([DISTCHECK_CONFIGURE_FLAGS])

AC_ARG_VAR([FC], [FORTRAN compiler command])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_HEADERS([src/moab/Version.h
                   src/moab/EntityHandle.hpp
                   src/parallel/moab_mpi_config.h
                   src/FCDefs.h])
AC_CONFIG_FILES([Makefile 
                 moab.make                  
                 src/Makefile
                 src/io/Makefile
                 src/io/mhdf/Makefile
                 src/parallel/Makefile
                 src/oldinc/Makefile
                 test/Makefile
                 test/h5file/Makefile
                 test/dual/Makefile
                 test/obb/Makefile
                 test/perf/Makefile
                 test/io/Makefile
                 test/parallel/Makefile
                 test/oldinc/Makefile
		 itaps/iBase_f.h
                 itaps/Makefile
                 itaps/imesh/Makefile
                 itaps/imesh/iMesh-Defs.inc
		 itaps/igeom/Makefile
                 itaps/igeom/FBiGeom-Defs.inc
                 tools/Makefile
                 tools/refiner/Makefile
                 tools/mbcoupler/Makefile
                 tools/mbzoltan/Config.moab
                 tools/mbzoltan/Makefile
                 tools/mcnpmit/Makefile
                 tools/qvdual/Makefile
                 tools/vtkMOABReader/Makefile
                 tools/dagmc/Makefile
		 tools/vtkMOABReader/CMakeLists.txt
                 doc/Makefile
                 doc/user.dox
                 doc/config.tex
                 examples/Makefile
                 examples/examples.make
                 examples/simple/makefile
                 examples/itaps/Makefile
                 MeshFiles/Makefile
                 MeshFiles/unittest/Makefile
                 MeshFiles/unittest/io/Makefile
                 ])
AC_CONFIG_COMMANDS([src/MOAB_FCDefs.h],
  [sed -e "s/FC_FUNC/MOAB_FC_FUNC/" src/FCDefs.h >src/MOAB_FCDefs.h])

AC_OUTPUT

AC_MSG_RESULT([C:   $CC $CFLAGS $CPPFLAGS])
AC_MSG_RESULT([C++: $CXX $CXXFLAGS $CPPFLAGS])
if test "x" != "x$FC"; then
  AC_MSG_RESULT([Fortran: $FC $FCFLAGS])
fi
if test "x" != "x$F77"; then
  AC_MSG_RESULT([Fortran: $F77 $FFLAGS])
fi

if test "x$HAVE_HDF5" = "xno"; then
  AC_MSG_WARN([
*************************************************************************
*            MOAB has been configured w/out the HDF5 library.
*              Support for native file format is disabled!
*************************************************************************])
fi

if test "x$HAVE_NETCDF" = "xno"; then
  AC_MSG_WARN([
*************************************************************************
*           MOAB has been configured w/out the NetCDF library.
*     Support for ExodusII/Genesis and other file formats is disabled!
*************************************************************************])
elif test "xno" = "x$NETCDF_SUFFICIENT_DIMS_VARS"; then
  AC_MSG_WARN([
*************************************************************************
*  The configured NetCDF library does not have sufficient constant limits
*  to support reading of ExodusII with large numbers of blocks or sidesets
*  or nodesets.  Please consider re-compiling your NetCDF library after
*  increasing the constants NC_MAX_DIMS and NC_MAX_VARS in netcdf.h to
*  65536 and 524288, respectively.  See: 
*  http://cubit.sandia.gov/help-version11.1/finite_element_model/export/exporting_exodus2_file.htm
*************************************************************************])
fi

if test "x$WARN_PARALLEL_HDF5" = "xyes"; then
  AC_MSG_WARN([
*************************************************************************
*        MOAB has been configured with parallel and HDF5 support
*     but the configured HDF5 library does not support parallel IO.
*            Some parallel IO capabilities will be disabled.
*************************************************************************])
fi

if test "x$WARN_PARALLEL_HDF5_NO_COMPLEX" = "xyes"; then
  AC_MSG_WARN([
*************************************************************************
*     Your parallel HDF5 library is configured without 
*     H5_MPI_COMPLEX_DERIVED_DATATYPE_WORKS .  For the types of IO
*     patterns MOAB typically does this will result in degrading
*     collective IO calls to independent IO, which may have a very
*     significant impact on IO performance.
*************************************************************************])
fi

