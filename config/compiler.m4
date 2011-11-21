#######################################################################################
# Get libtool configuration variable
# Arguments:
#  libtool config tag (e.g CXX)
#  libtool variable name
#  variable in which to store result
#######################################################################################
AC_DEFUN([FATHOM_LIBTOOL_VAR], [
  $3=`./libtool --tag=$1 --config | sed -e 's/^$2=//p' -e 'd' | tr -d '"\n'`
])


#######################################################################################
# Check if the C++ compiler works.
# Arguments: action on success and action on failure
#######################################################################################
AC_DEFUN([FATHOM_CHECK_CXX_WORKS], [
  AC_LANG_PUSH([C++])
  AC_MSG_CHECKING([if $CXX works])
  AC_COMPILE_IFELSE(
   [AC_LANG_PROGRAM( [class Cl { int i; public: Cl(int x) : i(x) {} };] [Cl x(5);])],
   [AC_COMPILE_IFELSE( 
      [AC_LANG_PROGRAM( [],
[#ifdef __cplusplus
   choke me
 #endif])],
      [AC_MSG_RESULT([no (accepted invalid input)]); $2], [AC_MSG_RESULT([yes]); $1])],
    [AC_MSG_RESULT([no (failed on trival valid C++ source)]); $2])

  AC_LANG_POP([C++])
])


########## Helper function for FATHOM_CHECK_COMPILERS #############
# args: compiler variable, compiler list, path
AC_DEFUN([FATHOM_SET_MPI_COMPILER], [
  if test "x" = "x${$1}"; then
    if test "x" = "x$3"; then
      AC_CHECK_PROGS([$1], [$2], [false])
    else
      $1=false
      for prog in $2 ; do
        if test -x "$3/$prog"; then
          $1="$3/$prog"
          AC_SUBST($1)
          export $1
          break
        fi
      done
    fi
    
    if test "x${$1}" = "xfalse"; then
      AC_MSG_ERROR([Cannot find MPI compiler.  Try specifying \$$1])
    fi
  fi
])

#######################################################################################
# Implement checks for C and C++ compilers, with all corresponding options
#
# Sets the following variables:
#  CPP      - The C preprocessor
#  CC       - The C compiler
#  CXX      - The C++ compiler
#  FC       - The Fortran compiler
#  CFLAGS   - C compiler flags
#  CXXFLAGS - C++ compiler flags
#  WITH_MPI - 'yes' if parallel support, 'no' otherwise
#
# Arguments:  three strings that msut be either "yes" or "no".
#             - test for C compiler
#             - test for C++ compiler
#             - test for Fortran compiler
#######################################################################################
AC_DEFUN([FATHOM_CHECK_COMPILERS], [

CHECK_CC="$1"
CHECK_CXX="$2"
CHECK_FC="$3"

# If not specified or invalid value, change to yes.
test "xno" = "x$CHECK_CC" || CHECK_CC=yes 
test "xno" = "x$CHECK_CXX" || CHECK_CXX=yes 
test "xno" = "x$CHECK_FC" || CHECK_FC=yes 

  # Save these before calling AC_PROG_CC or AC_PROG_CXX
  # because those macros will modify them, and we want
  # the original user values, not the autoconf defaults.
USER_CXXFLAGS="$CXXFLAGS"
USER_CFLAGS="$CFLAGS"

  # Check for Parallel
  # Need to check this early so we can look for the correct compiler
AC_ARG_WITH( [mpi], AC_HELP_STRING([[--with-mpi@<:@=DIR@:>@]], [Enable parallel support]),
             [WITH_MPI=$withval],[WITH_MPI=no] )
if test "xno" != "x$WITH_MPI"; then

  CC_LIST="mpixlc mpicc mpcc"
  CXX_LIST="mpixlcxx mpiCC mpCC mpicxx"
  FC_LIST="mpixlf95 mpixlf90 mpif90"
  F77_LIST="mpixlf77 mpif77"
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-mpi=\"${withval}\""
  
  if test "xyes" == "x$WITH_MPI"; then
    FATHOM_SET_MPI_COMPILER([CC],  [$CC_LIST])
    FATHOM_SET_MPI_COMPILER([CXX],[$CXX_LIST])
    FATHOM_SET_MPI_COMPILER([FC],  [$FC_LIST])
    FATHOM_SET_MPI_COMPILER([F77],[$F77_LIST])
  else
    FATHOM_SET_MPI_COMPILER([CC],  [$CC_LIST],[${WITH_MPI}/bin])
    FATHOM_SET_MPI_COMPILER([CXX],[$CXX_LIST],[${WITH_MPI}/bin])
    FATHOM_SET_MPI_COMPILER([FC],  [$FC_LIST],[${WITH_MPI}/bin])
    FATHOM_SET_MPI_COMPILER([F77],[$F77_LIST],[${WITH_MPI}/bin])
    WITH_MPI=yes
  fi
fi

if test "xno" != "x$CHECK_CC"; then
  AC_PROG_CC
fi
AC_PROG_CPP
if test "xno" != "x$CHECK_CXX"; then
  AC_PROG_CXX
  AC_PROG_CXXCPP
fi
if test "xno" != "x$CHECK_FC"; then
  AC_PROG_FC
  AC_PROG_F77
fi

]) # FATHOM_CHECK_COMPILERS



#######################################################################################
# Implement checks for C and C++ compiler options
#
#  CFLAGS   - C compiler flags
#  CXXFLAGS - C++ compiler flags
#  DEBUG - yes if specified, no otherwise
#
#######################################################################################
AC_DEFUN([FATHOM_COMPILER_FLAGS], [


if test "xno" != "x$CHECK_CC"; then
  FATHOM_CC_FLAGS
fi
if test "xno" != "x$CHECK_CXX"; then
  FATHOM_CXX_FLAGS
fi


# Try to determine compiler-specific flags.  This must be done
# before setting up libtool so that it can override libtool settings.
CFLAGS="$USER_CFLAGS $FATHOM_CC_SPECIAL"
CXXFLAGS="$USER_CXXFLAGS $FATHOM_CXX_SPECIAL"
FFLAGS="$USER_FFLAGS $FATHOM_F77_SPECIAL"
FCFLAGS="$USER_FCFLAGS $FATHOM_FC_SPECIAL"

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
                enable_cc_optimize=$enableval
		enable_fc_optimize=$enableval
		], 
               [enable_cxx_optimize=
                enable_cc_optimize=
		enable_fc_optimize=
		] )

# Do enable_optimize by default, unless user has specified
# custom CXXFLAGS or CFLAGS
DEBUG=no
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
  if test "x$enable_fc_optimize" = "x"; then
    if test "x$USER_FCFLAGS" = "x"; then
      enable_fc_optimize=yes
    fi
  fi
  if test "x$enable_f77_optimize" = "x"; then
    if test "x$USER_FFLAGS" = "x"; then
      enable_f77_optimize=yes
    fi
  fi
fi

# Choose compiler flags from CLI args
if test "xyes" = "x$enable_debug"; then
  DEBUG=yes
  CXXFLAGS="$CXXFLAGS -g"
  CFLAGS="$CFLAGS -g"
  FCFLAGS="$FCFLAGS -g"
  FFLAGS="$FFLAGS -g"
fi
if test "xyes" = "x$enable_cxx_optimize"; then
  CXXFLAGS="$CXXFLAGS -O2 -DNDEBUG"
fi
if test "xyes" = "x$enable_cc_optimize"; then
  CFLAGS="$CFLAGS -O2 -DNDEBUG"
fi
if test "xyes" = "x$enable_fc_optimize"; then
  FCFLAGS="$FCFLAGS -O2"
fi
if test "xyes" = "x$enable_f77_optimize"; then
  FFLAGS="$FFLAGS -O2"
fi

  # Check for 32/64 bit.
  # This requires FATHOM_CXX_FLAGS and FATHOM_CC_FLAGS to have been called first
AC_ARG_ENABLE( 32bit, AC_HELP_STRING([--enable-32bit],[Force 32-bit objects]),
[
  if test "xyes" != "x$enableval"; then
    AC_MSG_ERROR([Unknown argument --enable-32bit=$enableval])
  elif test "x" = "x$FATHOM_CXX_32BIT"; then
    AC_MSG_ERROR([Don't know how to force 32-bit C++ on this platform.  Try setting CXXFLAGS manually])
  elif test "x" = "x$FATHOM_CC_32BIT"; then
    AC_MSG_ERROR([Don't know how to force 32-bit C on this platform.  Try setting CFLAGS manually])
  fi
  CXXFLAGS="$CXXFLAGS $FATHOM_CXX_32BIT"
  CFLAGS="$CFLAGS $FATHOM_CC_32BIT"
  enable_32bit=yes
])
# This requires FATHOM_CXX_FLAGS and FATHOM_CC_FLAGS to have been called first
AC_ARG_ENABLE( 64bit, AC_HELP_STRING([--enable-64bit],[Force 64-bit objects]),
[
  if test "xyes" != "x$enableval"; then
    AC_MSG_ERROR([Unknown argument --enable-64bit=$enableval])
  elif test "x" = "x$FATHOM_CXX_64BIT"; then
    AC_MSG_ERROR([Don't know how to force 64-bit C++ on this platform.  Try setting CXXFLAGS manually])
  elif test "x" = "x$FATHOM_CC_64BIT"; then
    AC_MSG_ERROR([Don't know how to force 64-bit C on this platform.  Try setting CFLAGS manually])
  elif test "xyes" = "x$enable_32bit"; then
    AC_MSG_ERROR([Cannot do both --enable-32bit and --enable-64bit])
  fi
  CXXFLAGS="$CXXFLAGS $FATHOM_CXX_64BIT"
  CFLAGS="$CFLAGS $FATHOM_CC_64BIT"
])

]) # FATHOM_COMPILER_FLAGS

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
AC_DEFUN([FATHOM_TRY_COMPILER_DEFINE], [
 AC_COMPILE_IFELSE([
 AC_LANG_PROGRAM( [[#ifndef $1
   choke me
 #endif]], []) ],
 [$2],[$3])
])

#################################################################################
# Check if the compiler defines a specific preprocessor macro with an integer
# value greater than or equal to the passed value
# Arguments:
#  - preprocessor define to check for
#  - numeric value to test
#  - action upon success
#  - action upon failure
#################################################################################
AC_DEFUN([FATHOM_TRY_COMPILER_DEFINE_GE], [
 AC_COMPILE_IFELSE([
 AC_LANG_PROGRAM( [[#if !defined($1) || $1 < $2
   choke me
 #endif]], []) ],
 [$3],[$4])
])


#######################################################################################
# Check for compiler-specific flags.
# Sets the following environmental variables:
#   FATHOM_CXX_SPECIAL : Any compiler-specific flags which must be specified
#   FATHOM_CXX_32BIT   : Flag to force compilation of 32-bit code
#   FATHOM_CXX_64BIT   : Flag to force compilation of 64-bit code
#######################################################################################
AC_DEFUN([FATHOM_CXX_FLAGS], [
AC_LANG_PUSH([C++])

# Detect compiler 
AC_MSG_CHECKING([for known c++ compilers])
# Autoconf does G++ for us
if test x$GXX = xyes; then
  cxx_compiler=GNU
  # Intel claims to be GCC, check for it here
  FATHOM_TRY_COMPILER_DEFINE([__INTEL_COMPILER],[cxx_compiler=Intel])
# Search for other compiler types
# For efficiency, limit checks to relevant OSs
else
  cxx_compiler=unknown
  case "$target_os" in
    aix*)
      FATHOM_TRY_COMPILER_DEFINE_GE([__IBMCPP__],[800],[cxx_compiler=VisualAge8],
        [FATHOM_TRY_COMPILER_DEFINE([__IBMCPP__],[cxx_compiler=VisualAge])])
      ;;
    solaris*|sunos*)
      FATHOM_TRY_COMPILER_DEFINE([__SUNPRO_CC],[cxx_compiler=SunWorkshop])
      ;;
    irix*)
      FATHOM_TRY_COMPILER_DEFINE([__sgi],[cxx_compiler=MIPSpro])
      ;;
    linux*)
      FATHOM_TRY_COMPILER_DEFINE([__INTEL_COMPILER],[cxx_compiler=Intel])
      FATHOM_TRY_COMPILER_DEFINE_GE([__IBMCPP__],[800],[cxx_compiler=VisualAge8],
        [FATHOM_TRY_COMPILER_DEFINE([__IBMCPP__],[cxx_compiler=VisualAge])])
      FATHOM_TRY_COMPILER_DEFINE([__DECCXX_VER],[cxx_compiler=Compaq])
      FATHOM_TRY_COMPILER_DEFINE([__SUNPRO_CC],[cxx_compiler=SunWorkshop])
      FATHOM_TRY_COMPILER_DEFINE([__PGI],[cxx_compiler=PortlandGroup])
      ;;
    hpux*)
      FATHOM_TRY_COMPILER_DEFINE([__HP_aCC],[cxx_compiler=HP])
      ;;
    windows*)
      FATHOM_TRY_COMPILER_DEFINE([__MSC_VER],[cxx_compiler=VisualStudio])
      FATHOM_TRY_COMPILER_DEFINE([__INTEL_COMPILER],[cxx_compiler=Intel])
      FATHOM_TRY_COMPILER_DEFINE([__DECCXX_VER],[cxx_compiler=Compaq])
      FATHOM_TRY_COMPILER_DEFINE([__BORLANDC__],[cxx_compiler=Borland])
      FATHOM_TRY_COMPILER_DEFINE([__CYGWIN__],[cxx_compiler=Cygwin])
      FATHOM_TRY_COMPILER_DEFINE([__MINGW32__],[cxx_compiler=MinGW])
      ;;
    *)
      FATHOM_TRY_COMPILER_DEFINE([__PGI],[cxx_compiler=PortlandGroup])
      ;;
  esac
fi

AC_LANG_POP([C++])
AC_MSG_RESULT([$cxx_compiler])
if test "x$cxx_compiler" = "xunknown"; then
  AC_MSG_WARN([Unrecognized C++ compiler: $CXX])
fi

# Now decide special compiler flags using compiler/cpu combination
AC_MSG_CHECKING([for known compiler/OS combinations])
case "$cxx_compiler:$host_cpu" in
  GNU:sparc*)
    FATHOM_CXX_32BIT=-m32
    FATHOM_CXX_64BIT=-m64
    FATHOM_CXX_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  GNU:powerpc*)
    FATHOM_CXX_32BIT=-m32
    FATHOM_CXX_64BIT=-m64
    FATHOM_CXX_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  GNU:i?86|GNU:x86_64)
    FATHOM_CXX_32BIT=-m32
    FATHOM_CXX_64BIT=-m64
    FATHOM_CXX_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  GNU:mips*)
    FATHOM_CXX_32BIT="-mips32 -mabi=32"
    FATHOM_CXX_64BIT="-mips64 -mabi=64"
    FATHOM_CXX_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  GNU:*)
    FATHOM_CXX_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  Intel:*)
    FATHOM_CXX_32BIT=-m32
    FATHOM_CXX_64BIT=-m64
    FATHOM_CXX_SPECIAL="$EXTRA_INTEL_FLAGS -wd981 -wd383 -wd1572 -wd2259"
    ;;
  VisualAge:*)
    FATHOM_CXX_32BIT=-q32
    FATHOM_CXX_64BIT=-q64
    FATHOM_CXX_SPECIAL="-qrtti=all"
    AR="ar -X 32_64"
    NM="nm -B -X 32_64"
    ;;
  VisualAge8:*)
    FATHOM_CXX_32BIT=-q32
    FATHOM_CXX_64BIT=-q64
    NM="nm -B -X 32_64"
    ;;
  MIPSpro:mips)
    FATHOM_CXX_32BIT=-n32
    FATHOM_CXX_64BIT=-64
    FATHOM_CXX_SPECIAL=-LANG:std
    ;;
  MIPSpro:*)
    FATHOM_CXX_SPECIAL=-LANG:std
    ;;
  SunWorkshop:sparc*)
    FATHOM_CXX_32BIT=-xarch=generic
    FATHOM_CXX_64BIT=-xarch=generic64
    ;;
  SunWorkshop:i?86|SunWorkshop:x86_64)
    FATHOM_CXX_32BIT=-m32
    FATHOM_CXX_64BIT=-m64
    ;;
  *)
    ;;
esac
AC_MSG_RESULT([$cxx_compiler:$host_cpu])
]) # end FATHOM_CXX_FLAGS

#######################################################################################
# Check for compiler-specific flags.
# Sets the following environmental variables:
#   FATHOM_CC_SPECIAL : Any compiler-specific flags which must be specified
#   FATHOM_CC_32BIT   : Flag to force compilation of 32-bit code
#   FATHOM_CC_64BIT   : Flag to force compilation of 64-bit code
#######################################################################################
AC_DEFUN([FATHOM_CC_FLAGS], [
AC_LANG_PUSH([C])

# Detect compiler 

AC_MSG_CHECKING([for known C compilers])
# Autoconf does gcc for us
if test x$GCC = xyes; then
  cc_compiler=GNU
  # Intel claims to be GCC, check for it here
  FATHOM_TRY_COMPILER_DEFINE([__INTEL_COMPILER],[cc_compiler=Intel])
# Search for other compiler types
# For efficiency, limit checks to relevant OSs
else
  cc_compiler=unknown
  case "$target_os" in
    aix*)
      FATHOM_TRY_COMPILER_DEFINE([__IBMC__],[cc_compiler=VisualAge])
      ;;
    solaris*|sunos*)
      FATHOM_TRY_COMPILER_DEFINE([__SUNPRO_C],[cc_compiler=SunWorkshop])
      ;;
    irix*)
      FATHOM_TRY_COMPILER_DEFINE([__sgi],[cc_compiler=MIPSpro])
      ;;
    linux*)
      FATHOM_TRY_COMPILER_DEFINE([__INTEL_COMPILER],[cc_compiler=Intel])
      FATHOM_TRY_COMPILER_DEFINE([__IBMC__],[cc_compiler=VisualAge])
      FATHOM_TRY_COMPILER_DEFINE([__DECC_VER],[cc_compiler=Compaq])
      FATHOM_TRY_COMPILER_DEFINE([__SUNPRO_C],[cc_compiler=SunWorkshop])
      FATHOM_TRY_COMPILER_DEFINE([__PGI],[cc_compiler=PortlandGroup])
      ;;
    hpux*)
      FATHOM_TRY_COMPILER_DEFINE([__HP_cc],[cc_compiler=HP])
      ;;
    windows*)
      FATHOM_TRY_COMPILER_DEFINE([__MSC_VER],[cc_compiler=VisualStudio])
      FATHOM_TRY_COMPILER_DEFINE([__INTEL_COMPILER],[cc_compiler=Intel])
      FATHOM_TRY_COMPILER_DEFINE([__DECC_VER],[cc_compiler=Compaq])
      FATHOM_TRY_COMPILER_DEFINE([__BORLANDC__],[cc_compiler=Borland])
      FATHOM_TRY_COMPILER_DEFINE([__TURBOC__],[cc_compiler=TurboC])
      FATHOM_TRY_COMPILER_DEFINE([__CYGWIN__],[cc_compiler=Cygwin])
      FATHOM_TRY_COMPILER_DEFINE([__MINGW32__],[cc_compiler=MinGW])
      ;;
    *)
      FATHOM_TRY_COMPILER_DEFINE([__PGI],[cc_compiler=PortlandGroup])
      ;;
  esac
fi
AC_LANG_POP([C])
AC_MSG_RESULT([$cc_compiler])
if test "x$cc_compiler" = "xunknown"; then
  AC_MSG_WARN([Unrecognized C compiler: $CXX])
fi

# Now decide special compiler flags using compiler/cpu combination
AC_MSG_CHECKING([for known compiler/OS combinations])
case "$cc_compiler:$host_cpu" in
  GNU:sparc*)
    FATHOM_CC_32BIT=-m32
    FATHOM_CC_64BIT=-m64
    FATHOM_CC_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  GNU:powerpc*)
    FATHOM_CC_32BIT=-m32
    FATHOM_CC_64BIT=-m64
    FATHOM_CC_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  GNU:i?86|GNU:x86_64)
    FATHOM_CC_32BIT=-m32
    FATHOM_CC_64BIT=-m64
    FATHOM_CC_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  Intel:*)
    FATHOM_CC_32BIT=-m32
    FATHOM_CC_64BIT=-m64
    FATHOM_CC_SPECIAL="$EXTRA_INTEL_FLAGS -wd981 -wd383 -wd1572"
    ;;
  GNU:mips*)
    FATHOM_CC_32BIT="-mips32 -mabi=32"
    FATHOM_CC_64BIT="-mips64 -mabi=64"
    FATHOM_CC_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  GNU:*)
    FATHOM_CC_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  VisualAge:*)
    case "$target_vendor" in
      bgp)
        FATHOM_CC_32BIT=-q32
        FATHOM_CC_64BIT=-q64
        AR="ar"
        NM="nm -B"
        ;;
      *)
        FATHOM_CC_32BIT=-q32
        FATHOM_CC_64BIT=-q64
        AR="ar -X 32_64"
        NM="nm -B -X 32_64"
        ;;
    esac
    ;;
  MIPSpro:mips)
    FATHOM_CC_32BIT=-n32
    FATHOM_CC_64BIT=-64
    FATHOM_CC_SPECIAL=-LANG:std
    ;;
  MIPSpro:*)
    FATHOM_CC_SPECIAL=-LANG:std
    ;;
  SunWorkshop:sparc*)
    FATHOM_CC_32BIT=-xarch=generic
    FATHOM_CC_64BIT=-xarch=generic64
    ;;
  SunWorkshop:i?86|SunWorkshop:x86_64)
    FATHOM_CC_32BIT=-m32
    FATHOM_CC_64BIT=-m64
    ;;
  *)
    ;;
esac
AC_MSG_RESULT([$cc_compiler:$host_cpu])
]) # end FATHOM_CC_FLAGS


#######################################################################################
# Check for Fortran compiler-specific flags.
#######################################################################################

dnl PAC_PROG_F77_CRAY_POINTER - Check if Fortran 77 supports Cray-style pointer.
dnl                             If so, set pac_cv_prog_f77_has_pointer to yes
dnl                             and find out if any extra compiler flag is
dnl                             needed and set it as CRAYPTR_FFLAGS.
dnl                             i.e. CRAYPTR_FFLAGS is meaningful only if
dnl                             pac_cv_prog_f77_has_pointer = yes.
dnl
dnl Synopsis:
dnl   PAC_PROG_F77_CRAY_POINTER([action-if-true],[action-if-false])
dnl D*/
AC_DEFUN([PAC_PROG_F77_CRAY_POINTER],[
AC_CACHE_CHECK([whether Fortran 77 supports Cray-style pointer],
pac_cv_prog_f77_has_pointer,[
AC_LANG_PUSH([Fortran 77])
AC_LANG_CONFTEST([
    AC_LANG_PROGRAM([],[
        integer M
        pointer (MPTR,M)
        data MPTR/0/
    ])
])
saved_FFLAGS="$FFLAGS"
pac_cv_prog_f77_has_pointer=no
CRAYPTR_FFLAGS=""
for ptrflag in '' '-fcray-pointer' ; do
    FFLAGS="$saved_FFLAGS $ptrflag"
    AC_COMPILE_IFELSE([], [
        pac_cv_prog_f77_has_pointer=yes
        CRAYPTR_FFLAGS="$ptrflag"
        break
    ])
done
dnl Restore FFLAGS first, since user may not want to modify FFLAGS
FFLAGS="$saved_FFLAGS"
dnl remove conftest after ac_lang_conftest
rm -f conftest.$ac_ext
AC_LANG_POP([Fortran 77])
])
if test "$pac_cv_prog_f77_has_pointer" = "yes" ; then
    AC_MSG_CHECKING([for Fortran 77 compiler flag for Cray-style pointer])
    if test "X$CRAYPTR_FFLAGS" != "X" ; then
        AC_MSG_RESULT([$CRAYPTR_FFLAGS])
    else
        AC_MSG_RESULT([none])
    fi
    ifelse([$1],[],[:],[$1])
else
    ifelse([$2],[],[:],[$2])
fi
])

dnl PAC_PROG_FC_CRAY_POINTER - Check if Fortran supports Cray-style pointer.
dnl                            If so, set pac_cv_prog_fc_has_pointer to yes
dnl                            and find out if any extra compiler flag is
dnl                            needed and set it as CRAYPTR_FCFLAGS.
dnl                            i.e. CRAYPTR_FCFLAGS is meaningful only if
dnl                            pac_cv_prog_fc_has_pointer = yes.
dnl
dnl Synopsis:
dnl   PAC_PROG_FC_CRAY_POINTER([action-if-true],[action-if-false])
dnl D*/
AC_DEFUN([PAC_PROG_FC_CRAY_POINTER],[
AC_CACHE_CHECK([whether Fortran 90 supports Cray-style pointer],
pac_cv_prog_fc_has_pointer,[
AC_LANG_PUSH([Fortran])
AC_LANG_CONFTEST([
    AC_LANG_PROGRAM([],[
        integer M
        pointer (MPTR,M)
        data MPTR/0/
    ])
])
saved_FCFLAGS="$FCFLAGS"
pac_cv_prog_fc_has_pointer=no
CRAYPTR_FCFLAGS=""
for ptrflag in '' '-fcray-pointer' ; do
    FCFLAGS="$saved_FCFLAGS $ptrflag"
    AC_COMPILE_IFELSE([],[
        pac_cv_prog_fc_has_pointer=yes
        CRAYPTR_FCFLAGS="$ptrflag"
        break
    ])
done
dnl Restore FCFLAGS first, since user may not want to modify FCFLAGS
FCFLAGS="$saved_FCFLAGS"
dnl remove conftest after ac_lang_conftest
rm -f conftest.$ac_ext
AC_LANG_POP([Fortran])
])
if test "$pac_cv_prog_fc_has_pointer" = "yes" ; then
    AC_MSG_CHECKING([for Fortran 90 compiler flag for Cray-style pointer])
    if test "X$CRAYPTR_FCFLAGS" != "X" ; then
        AC_MSG_RESULT([$CRAYPTR_FCFLAGS])
    else
        AC_MSG_RESULT([none])
    fi
    ifelse([$1],[],[:],[$1])
else
    ifelse([$2],[],[:],[$2])
fi
])