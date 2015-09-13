dnl Process this file with autoconf 2.60 to produce a configure script
AC_DEFUN([CONFIGURE_MESQUITE],[

mesquite_ns="MOAB_Mesquite299"
AC_DEFINE(MESQUITE_NS,[mesquite_ns],"Mesquite namespace")
#AC_SUBST(MESQUITE_NS)
# AC_DEFINE_UNQUOTED(MESQUITE_NS,[$namespace],"Mesquite namespace")
AC_DEFINE(MESQUITE_NS_ALIAS,[1],"Mesquite namespace alias")

AC_DEFINE_UNQUOTED(MSQ_VERSION_MAJOR,[2],"Mesquite Major Version")
AC_DEFINE_UNQUOTED(MSQ_VERSION_MINOR,[99],"Mesquite Minor Version")
AC_DEFINE(MSQ_VERSION_STRING,["MOAB :: Mesquite 2.99"],"Mesquite Version String")

MSQ_DO_RELEASE="$enable_optimize"
MSQ_DO_DEBUG="$enable_debug"
MSQ_DO_OPTIMIZE="$enable_optimize"
MSQ_DO_FAST=""
MSQ_DEBUG_SYMBOLS="$enable_debug"
MSQ_DEBUG_ASSERTS=""
MSQ_DEBUG_OUT=""
MSQ_DO_TIMERS="no"
MSQ_TRAP_FPE=""
MSQ_QUITE_MAKE="yes"
MSQ_DO_32BIT="$enable_32bit"
MSQ_DO_64BIT="$enable_64bit"

#------------------------------------------------------------------------------
# Compile flag options -- need to do this *before* detecting the compiler
# otherwise cannot tell if user has CXXFLAGS defined already or if config
# set to a default value.
#------------------------------------------------------------------------------

# m4_ifdef([AM_SILENT_RULES],[
#   if test "xyes" = "x$MSQ_DO_DEBUG"; then
#     AM_SILENT_RULES(yes)
#   else
#     AM_SILENT_RULES(no)
#   fi
# ])



# If neather debug or release is specified, enable release
if test -z "$MSQ_DO_RELEASE" -a -z "$MSQ_DO_DEBUG"; then
  if test "x$have_user_CXXFLAGS" = "xno"; then
    MSQ_DO_RELEASE=yes
  fi
fi

# if release, then enable appropriate sub-options
if test "$MSQ_DO_RELEASE" = yes; then
  if test -z "$MSQ_DO_OPTIMIZE"; then
    MSQ_DO_OPTIMIZE=yes
  fi
  if test -z "$MSQ_DEBUG_ASSERTS"; then
    MSQ_DEBUG_ASSERTS=no
  fi
fi

# if debug, then enable appropriate sub-options
if test "$MSQ_DO_DEBUG" = "yes"; then
  if test "x$enable_dependency_tracking" = "x"; then
    enable_dependency_tracking=yes
  fi
  if test -z "$MSQ_DEBUG_SYMBOLS"; then
    MSQ_DEBUG_SYMBOLS=yes
  fi
  if test -z "$MSQ_DEBUG_OUT"; then
    MSQ_DEBUG_OUT=yes
  fi
  if test -z "$MSQ_TRAP_FPE"; then
    MSQ_TRAP_FPE=yes
  fi
  if test -z "$MSQ_QUITE_MAKE"; then
    MSQ_QUITE_MAKE=yes
  fi
fi


#------------------------------------------------------------------------------
# Construct compiler flags from options
#------------------------------------------------------------------------------

# Construct compiler options from above configure options
#SNL_DETECT_CXX
#SNL_REQUIRED_CXX_FLAGS
if test "$MSQ_DO_FAST" = "yes"; then
  SNL_CXX_COMPILE_FAST
elif test "$MSQ_DO_OPTIMIZE" = "yes"; then
  SNL_CXX_COMPILE_OPTIMIZED
fi
if test "$MSQ_DEBUG_SYMBOLS" = "yes"; then
  SNL_CXX_DEBUG_SYMBOLS
fi
if test "$MSQ_DEBUG_ASSERTS" = "no"; then
  MSQ_AM_CPPFLAGS="$MSQ_AM_CPPFLAGS -DNDEBUG"
fi
if test -n "$MSQ_DEBUG_OUT"; then
  if test "$MSQ_DEBUG_OUT" = "yes"; then 
    MSQ_AM_CPPFLAGS="$MSQ_AM_CPPFLAGS -DMSQ_ENABLE_DEBUG"
  elif test "$MSQ_DEBUG_OUT" != "no"; then
    MSQ_AM_CPPFLAGS="$MSQ_AM_CPPFLAGS -DMSQ_ENABLE_DEBUG=$MSQ_DEBUG_OUT"
  fi
fi
if test "$MSQ_TRAP_FPE" = "yes"; then
  MSQ_AM_CPPFLAGS="$MSQ_AM_CPPFLAGS -DMSQ_TRAP_FPE"
fi
if test "$MSQ_DO_TIMERS" = "yes"; then
  MSQ_AM_CPPFLAGS="MSQ_$AM_CPPFLAGS -DMSQ_USE_FUNCTION_TIMERS"
fi
if test "x$MSQ_QUITE_MAKE" = "xyes"; then
  MSQ_LIBTOOL_PREFIX='@echo "building $@...";'
  MSQ_LIBTOOL_FLAGS='--silent'
else
  MSQ_LIBTOOL_PREFIX=
  MSQ_LIBTOOL_FLAGS=
fi

# if test "x$MSQ_DO_32BIT" = "xyes"; then
#   SNL_CXX_COMPILE_32BIT
# fi
# if test "x$MSQ_DO_64BIT" = "xyes"; then
#   if test "x$MSQ_DO_32BIT" = "xyes"; then
#     AC_MSG_ERROR([--enable-32bit and --enable-64bit are mutually exclusive.])
#   fi
#   SNL_CXX_COMPILE_64BIT
# fi
# AC_LANG_CPLUSPLUS

#------------------------------------------------------------------------------
# LIBTOOL
#------------------------------------------------------------------------------
# Don't set up libtool until we're done messing with the compiler flags
# so we're sure libtool is set up correctly for the flags we want.
# Mostly need to make sure ar and nm are correct on IBM AIX, and that's
# handled inside our compiler detection code.

# AC_LIBTOOL_WIN32_DLL
# AM_PROG_LIBTOOL

#-----------------------------------------------------------------------------
# Check for required headers
#-----------------------------------------------------------------------------
#AC_HEADER_STDC


# Use C++ compiler because C allows undefined functions, so these
# checks don't achive much if compiled as C.  Also, if C and C++
# compilers are mis-matched, what works for one may not work for the
# other and the C++ one is what is actually used for this in Mesquite.
AC_LANG_PUSH(C++)
AC_MSG_CHECKING( for fpsetmask );
AC_TRY_COMPILE( [#include <ieeefp.h>],
                [fpsetmask(FP_X_INV|FP_X_OFL|FP_X_DZ);],
                [AM_CPPFLAGS="${AM_CPPFLAGS} -DHAVE_FPSETMASK"
                 AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)] )
AC_MSG_CHECKING( for feenableexcept );
AC_TRY_COMPILE( [#define _GNU_SOURCE
                 #include <fenv.h>  ],
                [feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);],
                [AM_CPPFLAGS="${AM_CPPFLAGS} -DHAVE_FEENABLEEXCEPT"
                 AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)] )
AC_LANG_POP(C++)

#-----------------------------------------------------------------------------
# Check C++ environment
#-----------------------------------------------------------------------------
# MSQ_CPLUSPLUS_FUNC

#-----------------------------------------------------------------------------
# Check for C++ features
#-----------------------------------------------------------------------------
AC_MSG_CHECKING( [for fd() in std::basic_file] )
AC_TRY_COMPILE([#include <fstream>],
               [using namespace std; ofstream f; int fd = f.rdbuf()->fd();],
               [CPPFLAGS="$CPPFLAGS -DFSTREAM_HAS_FD"; AC_MSG_RESULT(yes)],
               [AC_MSG_RESULT(no)])

#------------------------------------------------------------------------------
# Other build options
#------------------------------------------------------------------------------

# AC_ARG_WITH(exodus,
#   [AC_HELP_STRING([--with-exodus(=DIR)],
#                    [Enable exodusII support and specifiy directory])
# AC_HELP_STRING([--without-exodus],[Disable exodusII support (default)])],
#   [EXODUS_ARG=${withval}],[EXODUS_ARG=no])

# AC_ARG_WITH(netcdf,
#   [AC_HELP_STRING([--with-netcdf(=DIR)],
#                   [ExodusII requires NetCDF - defaults to values for --with-exodus])
# AC_HELP_STRING([--without-netcdf],[Skip NetCDF check])],
#   [NETCDF_ARG=$withval], [NETCDF_ARG=])

# AC_ARG_WITH(cppunit,
#   [AC_HELP_STRING([--with-cppunit(=DIR)],[Specify directory where CppUnit is installed.])
# AC_HELP_STRING([--without-cppunit],   [Disable CppUnit tests])],
#   [CPPUNIT_ARG=${withval}], [CPPUNIT_ARG=])


#-------------------------------------------------------------------------------
# Configure different options
#-------------------------------------------------------------------------------

# Configure Exodus and NetCDF
# AM_CONDITIONAL([WITH_EXODUS],[test "x$EXODUS_ARG" != "xno"])
# if test "x$EXODUS_ARG" != "xno"; then
#   old_CPPFLAGS="$CPPFLAGS"
#   old_LDFLAGS="$LDFLAGS"
#   old_LIBS="$LIBS"
#   AM_CPPFLAGS="$AM_CPPFLAGS -DMSQ_USING_EXODUS"
  
#     # If user specified path for Exodus, add to appropriate vars
#   if test "x$EXODUS_ARG" != "xyes"; then
#     EXODUS_INC="-I${EXODUS_ARG}/include"
#     EXODUS_LNK="-L${EXODUS_ARG}/lib"
#     CPPFLAGS="$CPPFLAGS $EXODUS_INC"
#     LDFLAGS="$LDFLAGS $EXODUS_LNK"
#     AM_CPPFLAGS="$AM_CPPFLAGS $EXODUS_INC"
#     AM_LDFLAGS="$AM_LDFLAGS $EXODUS_LNK"
#   fi
#   AM_LDFLAGS="$AM_LDFLAGS -lexoIIv2c"
    
#     # Check for netcdf unless user explicitly said --without-netcdf
#   if test "x$NETCDF_ARG" != "xno"; then
#       # If user specified path for NetCDF, add to appropriate vars
#     if test "x$NETCDF_ARG" != "xyes" -a "x$NETCDF_ARG" != "x" ; then
#       NETCDF_INC="-I${NETCDF_ARG}/inc -I${NETCDF_ARG}/include"
#       NETCDF_LNK="-L${NETCDF_ARG}/lib"
#       CPPFLAGS="$CPPFLAGS $NETCDF_INC"
#       LDFLAGS="$LDFLAGS $NETCDF_LNK"
#       AM_CPPFLAGS="$AM_CPPFLAGS $NETCDF_INC"
#       AM_LDFLAGS="$AM_LDFLAGS $NETCDF_LNK"
#     fi
#       # Check for NetCDF, but don't fail configure unless user
#       # explicitly specified that netcdf should be included
#     AC_CHECK_LIB( [netcdf], [nc_open], 
#       [LIBS="$LIBS -lnetcdf"; AM_LDFLAGS="$AM_LDFLAGS -lnetcdf"],
#       [if test "x$NETCDF_ARG" != "x"; then AC_MSG_ERROR([NetCDF library not found]); fi])
#     AC_CHECK_HEADER( [netcdf.h], [], 
#       [if test "x$NETCDF_ARG" != "x"; then AC_MSG_ERROR([netcdf.h not found]); fi ])
#   fi
  
#     # Check for ExodusII
#   AC_CHECK_LIB( [exoIIv2c], [ex_close], [], [AC_MSG_ERROR([ExodusII library not found])] )
#   AC_CHECK_HEADER( [exodusII.h], [], [AC_MSG_ERROR([exodusII.h not found])] )
#   CPPFLAGS="$old_CPPFLAGS"
#   LDFLAGS="$old_LDFLAGS"
#   LIBS="$old_LIBS"
# fi



# CPPUnit
HAVE_CPPUNIT="no"
# CPPUNIT_CPPFLAGS=
# CPPUNIT_LDFLAGS=
# if test "x$CPPUNIT_ARG" = "xno"; then
#   HAVE_CPPUNIT=no
# else
#     # Incorporate user-specified search paths
#   old_CPPFLAGS="$CPPFLAGS"
#   old_LDFLAGS="$LDFLAGS"
#   old_LIBS="$LIBS"
#   if test "x$CPPUNIT_ARG" != "x"; then
#     if test "x$CPPUNIT_ARG" != "xyes"; then
#       CPPUNIT_CPPFLAGS="-I$CPPUNIT_ARG/include"
#       CPPUNIT_LDFLAGS="-L$CPPUNIT_ARG/lib"
#       CPPFLAGS="$CPPFLAGS $CPPUNIT_CPPFLAGS"
#       LDFLAGS="$CPPFLAGS $CPPUNIT_LDFLAGS"
#     fi
#   fi
#     # Set some variables 
#   HAVE_CPPUNIT=yes
#   CPPUNIT_LDFLAGS="$CPPUNIT_LDFLAGS -lcppunit"
#     # Check that CppUnit exists
#   AC_CHECK_LIB( [dl], [dlopen], 
#     [LDFLAGS="$LDFLAGS -ldl"; CPPUNIT_LDFLAGS="$CPPUNIT_LDFLAGS -ldl"] )
#   AC_LANG_PUSH(C++)
#   AC_CHECK_HEADER( [cppunit/Test.h], [], [HAVE_CPPUNIT=no] )
#   AC_CHECK_LIB( [cppunit], [main], [], [HAVE_CPPUNIT=no] )
#   AC_LANG_POP(C++)
#     # If user explicitly specified --with-cppunit, fail if it was not found
#   if test "x$CPPUNIT_ARG" != "x"; then
#     if test "x$HAVE_CPPUNIT" == "xno"; then
#       AC_MSG_ERROR([CppUnit not found])
#     fi
#   fi
#     # restore some modified state
#   CPPFLAGS="$old_CPPFLAGS"
#   LDFLAGS="$old_LDFLAGS"
#   LIBS="$old_LIBS"
# fi

#------------------------------------------------------------------------------
# ITAPS
#------------------------------------------------------------------------------

# ITAPS_API([imesh], [IMESH], [iMesh])
# ITAPS_API([igeom], [IGEOM], [iGeom])
# ITAPS_API([irel],  [IREL],  [iRel],  [IMESH], [IGEOM])
# ITAPS_API([imeshp],[IMESHP],[iMeshP],[IMESH])

MSQ_AM_CPPFLAGS="${MSQ_AM_CPPFLAGS}"
# AM_CONDITIONAL([ENABLE_ITAPS],[test "x$ENABLE_ITAPS" = "xyes"])
# AC_SUBST(IBASE_INCL)

#------------------------------------------------------------------------------
# The End
#------------------------------------------------------------------------------
if test "x" = "x$docdir"; then
  docdir='${datadir}/doc/mesquite'
  AC_SUBST(docdir)
fi

MESQUITE_LIBS="$MESQUITE_LIBS -lm"
AC_SUBST(MSQ_AM_CXXFLAGS)
AC_SUBST(MSQ_AM_CPPFLAGS)
AC_SUBST(MSQ_AM_LDFLAGS)
AC_SUBST(MESQUITE_LIBS)
AC_SUBST(MSQ_LIBTOOL_PREFIX)
AC_SUBST(MSQ_LIBTOOL_FLAGS)
# AC_SUBST(CPPUNIT_CPPFLAGS)
# AC_SUBST(CPPUNIT_LDFLAGS)
AC_SUBST(MESQUITE_IMESH_MAKE_INCLUDE)

# We don't want to automatically generate mesquite_config.h.in
# automatically using autoheader.  It puts #defines in the file
# that are not appropriate for a library (will conflict with app's
# #defines).  Try to disable it
#AUTOHEADER="touch \$@"
# AUTOHEADER=":"
# AC_SUBST(AUTOHEADER)

TEST_MAKE_INCLUDE='include $(abs_top_builddir)/src/mesquite/msqcppflags.make'
UTIL_MAKE_INCLUDE='include $(abs_top_builddir)/src/mesquite/msqcppflags.make'
AC_SUBST(TEST_MAKE_INCLUDE)
AC_SUBST(UTIL_MAKE_INCLUDE)

# The first header is the header that autoheader creates.
# It should not be used because it contains stuff that is
# not restricted to the Mesquite namespace and therefore
# may conflict with the config.h of the application using
# mesquite.
AC_CONFIG_HEADERS([src/mesquite/include/mesquite_config.h
                   src/mesquite/include/mesquite_version.h])

AC_CONFIG_FILES([src/mesquite/Makefile
                 itaps/mesquite/Makefile
                 test/mesquite/Makefile
                 test/mesquite/meshfiles.h
                 test/mesquite/ActiveSetTest/Makefile
                 test/mesquite/algorithm_test/Makefile
                 test/mesquite/analytical_grad_3D/Makefile
                 test/mesquite/convert/Makefile
                 test/mesquite/benchmark_tests/Makefile
                 test/mesquite/paraboloid_domain_test/Makefile
                 test/mesquite/headers/Makefile
                 test/mesquite/idft_time/Makefile
                 test/mesquite/igeom/Makefile
                 test/mesquite/imesh/Makefile
                 test/mesquite/jacobi/Makefile
                 test/mesquite/laplacian_test/Makefile
                 test/mesquite/laplacian_polygon_test/Makefile
                 test/mesquite/nongradient_test/Makefile
                 test/mesquite/pyramid/Makefile
                 test/mesquite/simple_hybrid_test/Makefile
                 test/mesquite/test_1/Makefile
                 test/mesquite/transform/Makefile
                 test/mesquite/tutorial/Makefile
                 test/mesquite/tutorial/tutorial.make
                 test/mesquite/unit/Makefile
                 test/mesquite/untangle_test/Makefile
                 test/mesquite/parallel_untangle_shape/Makefile 
                 test/mesquite/parallel_smooth_laplace/Makefile
                 test/mesquite/wedge/Makefile
                 test/mesquite/wrapper_tests/Makefile
                 test/mesquite/2d_target/Makefile
                 test/mesquite/2d_metrics/Makefile
                 test/mesquite/2d_formulation/Makefile
                 test/mesquite/synchronous/Makefile
                 test/mesquite/high_aspect_ratio/Makefile
                 test/mesquite/higher_order/Makefile
                 test/mesquite/slaved/Makefile
                 tools/mesquite/Makefile
                 MeshFiles/mesquite/Makefile
                 ])
                 
# Remove any source depenency data so that it is regenerated.
# This isn't necessary for a correctly working configure script,
# but from a support POV, it is a easier to explain and remember
# that the configure script should be re-run without having to
# also remember to do "rm .deps/*" when a source file is removed
# from the build.
# AC_CONFIG_COMMANDS_POST([rm -f .deps/*])
                 
# echo "MSQ_CXXFLAGS = $MSQ_CXXFLAGS"

])
