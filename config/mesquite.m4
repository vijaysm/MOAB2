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

# AC_ARG_WITH(cppunit,
#   [AC_HELP_STRING([--with-cppunit(=DIR)],[Specify directory where CppUnit is installed.])
# AC_HELP_STRING([--without-cppunit],   [Disable CppUnit tests])],
#   [CPPUNIT_ARG=${withval}], [CPPUNIT_ARG=])

#-------------------------------------------------------------------------------
# Configure different options
#-------------------------------------------------------------------------------

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

MSQ_AM_CPPFLAGS="${MSQ_AM_CPPFLAGS}"

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

TEST_MAKE_INCLUDE='include $(abs_top_builddir)/src/mesquite/msqcppflags.make'
UTIL_MAKE_INCLUDE='include $(abs_top_builddir)/src/mesquite/msqcppflags.make'
AC_SUBST(TEST_MAKE_INCLUDE)
AC_SUBST(UTIL_MAKE_INCLUDE)

])
