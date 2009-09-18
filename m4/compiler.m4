#######################################################################################
# Get libtool configuration variable
# Arguments:
#  libtool config tag (e.g CXX)
#  libtool variable name
#  variable in which to store result
#######################################################################################
AC_DEFUN([ITAPS_LIBTOOL_VAR], [
  echo "SED=$SED" > .tmp
  ./libtool --tag=$1 --config >>.tmp
  echo "echo \$$2" >> .tmp
  chmod +x .tmp
  $3=`./.tmp`
  rm .tmp
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
# Arguments:  three strings that msut be either "yes" or "no".
#             - test for C compiler
#             - test for C++ compiler
#             - test for Fortran compiler
#######################################################################################
AC_DEFUN([SNL_CHECK_COMPILERS], [

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
AC_ARG_WITH( [mpi], AC_HELP_STRING([[--with-mpi(=DIR)]], [Enable parallel support]),
             [WITH_MPI=$withval],[WITH_MPI=no] )
case "x$WITH_MPI" in
  xno)
    CC_LIST="cc gcc cl egcs pgcc"
    CXX_LIST="CC aCC cxx xlC_r xlC c++ g++ pgCC gpp cc++ cl FCC KCC RCC"
    FC_LIST="gfortran ifort pgf90"
    F77_LIST="gfortran ifort pgf77"
    ;;
  xyes)
    CC_LIST="mpicc mpcc"
    CXX_LIST="mpiCC mpCC mpicxx"
    FC_LIST="mpif90"
    F77_LIST="mpif77"
    ;;
  x*)
    if test "x" = "x$CC"; then
      for prog in mpicc mpcc; do
        if test -x ${WITH_MPI}/bin/$prog; then
          CC="${WITH_MPI}/bin/$prog"
          CC_LIST="$prog"
        fi
      done
    else
      CC_LIST="$CC"
    fi
    if test "x" = "x$CXX"; then
      for prog in mpicxx mpiCC mpCC mpicxx; do
        if test -x ${WITH_MPI}/bin/$prog; then
          CXX="${WITH_MPI}/bin/$prog"
          CXX_LIST="$prog"
        fi
      done
    else
      CXX_LIST="$CXX"
    fi
    if test "x" = "x$FC"; then
      for prog in mpif90; do
        if test -x ${WITH_MPI}/bin/$prog; then
          FC="${WITH_MPI}/bin/$prog"
          FC_LIST="$prog"
        fi
      done
    else
      FC_LIST="$FC"
    fi
    if test "x" = "x$F77";then
      for prog in mpif77; do
        if test -x ${WITH_MPI}/bin/$prog; then
          F77="${WITH_MPI}/bin/$prog"
          F77_LIST="$prog"
        fi
      done
    else
      F77_LIST="$F77"
    fi
    WITH_MPI=yes
    ;;
esac

if test "xno" != "x$CHECK_CC"; then
  AC_PROG_CC( [$CC_LIST] )
fi
AC_PROG_CPP
if test "xno" != "x$CHECK_CXX"; then
  AC_PROG_CXX( [$CXX_LIST] )
  AC_PROG_CXXCPP
fi
if test "xno" != "x$CHECK_FC"; then
  AC_PROG_FC( [$FC_LIST] )
  AC_PROG_F77( [$F77_LIST] )
fi

]) # SNL_CHECK_COMPILERS



#######################################################################################
# Implement checks for C and C++ compiler options
#
#  CFLAGS   - C compiler flags
#  CXXFLAGS - C++ compiler flags
#
#######################################################################################
AC_DEFUN([SNL_COMPILER_FLAGS], [


if test "xno" != "x$CHECK_CC"; then
  SNL_CC_FLAGS
fi
if test "xno" != "x$CHECK_CXX"; then
  SNL_CXX_FLAGS
fi


# Try to determine compiler-specific flags.  This must be done
# before setting up libtool so that it can override libtool settings.
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
                enable_cc_optimize=$enableval
		enable_fc_optimize=$enableval
		], 
               [enable_cxx_optimize=
                enable_cc_optimize=
		enable_fc_optimize=
		] )

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
  if test "x$enable_fc_optimize" = "x"; then
    if test "x$USER_FCFLAGS" = "x"; then
      enable_fc_optimize=yes
    fi
  fi
fi

# Choose compiler flags from CLI args
if test "xyes" = "x$enable_debug"; then
  CXXFLAGS="$CXXFLAGS -g"
  CFLAGS="$CFLAGS -g"
  FCFLAGS="$FCFLAGS -g"
fi
if test "xyes" = "x$enable_cxx_optimize"; then
  CXXFLAGS="$CXXFLAGS -O2 -DNDEBUG"
fi
if test "xyes" = "x$enable_cc_optimize"; then
  CFLAGS="$CFLAGS -O2 -DNDEBUG"
fi
if test "xyes" = "x$enable_fc_optimize"; then
  FCFLAGS="$FCFLAGS -O2 -DNDEBUG"
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

]) # SNL_COMPILER_FLAGS

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
AC_LANG_PUSH([C++])

# Detect compiler 
AC_MSG_CHECKING([for known c++ compilers])
# Autoconf does G++ for us
if test x$GXX = xyes; then
  cxx_compiler=GNU
  # Intel claims to be GCC, check for it here
  SNL_TRY_COMPILER_DEFINE([__INTEL_COMPILER],[cxx_compiler=Intel])
# Search for other compiler types
# For efficiency, limit checks to relevant OSs
else
  cxx_compiler=unknown
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
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    SNL_CXX_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  GNU:powerpc*)
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    SNL_CXX_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  GNU:i?86|GNU:x86_64)
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    SNL_CXX_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  GNU:mips*)
    SNL_CXX_32BIT="-mips32 -mabi=32"
    SNL_CXX_64BIT="-mips64 -mabi=64"
    SNL_CXX_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  GNU:*)
    SNL_CXX_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  Intel:*)
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    SNL_CXX_SPECIAL="$EXTRA_GNU_FLAGS -wd981 -wd383"
    ;;
  VisualAge:*)
    SNL_CXX_32BIT=-q32
    SNL_CXX_64BIT=-q64
    # Do V5.0 namemangling for compatibility with ACIS, and enable RTTI
    case "$target_vendor" in
      bgp)
        SNL_CXX_SPECIAL=""
        AR="ar"
        NM="nm -B"
        ;;
      *)
        SNL_CXX_SPECIAL="-qrtti=all -qnamemangling=v5"
        AR="ar"
        NM="nm -B -X 32_64"
        ;;
    esac
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
AC_LANG_PUSH([C])

# Detect compiler 

AC_MSG_CHECKING([for known C compilers])
# Autoconf does gcc for us
if test x$GCC = xyes; then
  cc_compiler=GNU
  # Intel claims to be GCC, check for it here
  SNL_TRY_COMPILER_DEFINE([__INTEL_COMPILER],[cc_compiler=Intel])
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
AC_LANG_POP([C])
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
    SNL_CC_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  GNU:powerpc*)
    SNL_CC_32BIT=-m32
    SNL_CC_64BIT=-m64
    SNL_CC_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  GNU:i?86|GNU:x86_64)
    SNL_CC_32BIT=-m32
    SNL_CC_64BIT=-m64
    SNL_CC_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  Intel:*)
    SNL_CC_32BIT=-m32
    SNL_CC_64BIT=-m64
    SNL_CC_SPECIAL="$EXTRA_GNU_FLAGS -wd981 -wd383"
    ;;
  GNU:mips*)
    SNL_CC_32BIT="-mips32 -mabi=32"
    SNL_CC_64BIT="-mips64 -mabi=64"
    SNL_CC_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  GNU:*)
    SNL_CC_SPECIAL="$EXTRA_GNU_FLAGS"
    ;;
  VisualAge:*)
    case "$target_vendor" in
      bgp)
        SNL_CC_32BIT=-q32
        SNL_CC_64BIT=-q64
        AR="ar"
        NM="nm -B"
        ;;
      *)
        SNL_CC_32BIT=-q32
        SNL_CC_64BIT=-q64
        AR="ar -X 32_64"
        NM="nm -B -X 32_64"
        ;;
    esac
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
