# Check for compiler-specific flags.
# Sets the following environmental variables:
#   SNL_CXX_SPECIAL : Any compiler-specific flags which must be specified
#   SNL_CXX_32BIT   : Flag to force compilation of 32-bit code
#   SNL_CXX_64BIT   : Flag to force compilation of 64-bit code
AC_DEFUN([SNL_CXX_FLAGS], [
AC_REQUIRE([AC_PROG_CXX])

# Detect compiler 
AC_MSG_CHECKING([for known c++ compilers])
# Autoconf does G++ for us
if test x$GXX = xyes; then
  cxx_compiler=GNU
  SNL_CXX_SPECIAL="-Wall -pipe"
# Detect IBM VisualAge compiler
elif $CXX 2>&1 | grep VisualAge > /dev/null; then
  cxx_compiler=VisualAge
# Sun Workshop/Forte
elif $CXX -V 2>&1 | grep Sun > /dev/null; then
  cxx_compiler=Sun
# Intel Compiler
elif $CXX -V 2>&1 | grep Intel > /dev/null; then
  cxx_compiler=Intel
# Compaq compiler
elif $CXX -V 2>&1 | grep Compaq > /dev/null; then
  cxx_compiler=Compaq
# IRIX
elif $CXX -version 2>&1 | grep MIPSpro > /dev/null; then
  cxx_compiler=MIPSpro
else
  # Try deciding based on compiler name
  case "$CXX" in
    xlC*)
      cxx_compiler=VisualAge
      ;;
    aCC)
      cxx_compiler=HP
      ;;
    icc)
      cxx_compiler=Intel
      ;;
    cl|cl.exe)
      cxx_compiler=VisualStudio
      ;;
    CC)
      case "$target_os" in
        solaris*)
          cxx_compiler=Sun
          ;;
        irix*)
          cxx_compiler=MIPSpro
          ;;
      esac
  esac
fi
if test x$cxx_compiler = x; then
  AC_MSG_RESULT([unknown])
  AC_MSG_WARN([Unrecognized c++ compiler: $CXX])
  cxx_compiler=unknown
else
  AC_MSG_RESULT([$cxx_compiler])
fi

# Now decide special compiler flags using compiler/cpu combination
AC_MSG_CHECKING([for known compiler/OS combinations])
case "$cxx_compiler:$host_cpu" in
  GNU:sparc*)
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    ;;
  GNU:powerpc*)
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    ;;
  GNU:i?86)
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    ;;
  GNU:mips*)
    SNL_CXX_32BIT="-mips32 -mabi=32"
    SNL_CXX_64BIT="-mips64 -mabi=64"
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
  Sun:sparc*)
    SNL_CXX_32BIT=-xarch=generic
    SNL_CXX_64BIT=-xarch=generic64
    ;;
  *)
    ;;
esac
AC_MSG_RESULT([$cxx_compiler:$host_cpu])
]) # end SNL_CXX_FLAGS

# Check for compiler-specific flags.
# Sets the following environmental variables:
#   SNL_CC_SPECIAL : Any compiler-specific flags which must be specified
#   SNL_CC_32BIT   : Flag to force compilation of 32-bit code
#   SNL_CC_64BIT   : Flag to force compilation of 64-bit code
AC_DEFUN([SNL_CC_FLAGS], [
AC_REQUIRE([AC_PROG_CC])

# Detect compiler 
AC_MSG_CHECKING([for known C compilers])
# Autoconf does gcc for us
if test x$GCC = xyes; then
  cc_compiler=GNU
  SNL_CC_SPECIAL="-Wall -pipe"
# Detect IBM VisualAge compiler
elif $CC 2>&1 | grep VisualAge > /dev/null; then
  cc_compiler=VisualAge
# Sun Workshop/Forte
elif $CC -V 2>&1 | grep Sun > /dev/null; then
  cc_compiler=Sun
# Intel Compiler
elif $CC -V 2>&1 | grep Intel > /dev/null; then
  cc_compiler=Intel
# Compaq compiler
elif $CC -V 2>&1 | grep Compaq > /dev/null; then
  cc_compiler=Compaq
# IRIX
elif $CC -version 2>&1 | grep MIPSpro > /dev/null; then
  cc_compiler=MIPSpro
else
  # Try deciding based on compiler name
  case "$CC" in
    xlC*)
      cc_compiler=VisualAge
      ;;
    aCC)
      cc_compiler=HP
      ;;
    icc)
      cc_compiler=Intel
      ;;
    cl|cl.exe)
      cc_compiler=VisualStudio
      ;;
    CC)
      case "$target_os" in
        solaris*)
          cc_compiler=Sun
          ;;
        irix*)
          cc_compiler=MIPSpro
          ;;
      esac
  esac
fi
if test x$cc_compiler = x; then
  AC_MSG_RESULT([unknown])
  AC_MSG_WARN([Unrecognized C compiler: $CC])
  cc_compiler=unknown
else
  AC_MSG_RESULT([$cc_compiler])
fi

# Now decide special compiler flags using compiler/cpu combination
AC_MSG_CHECKING([for known compiler/OS combinations])
case "$cc_compiler:$host_cpu" in
  GNU:sparc*)
    SNL_CC_32BIT=-m32
    SNL_CC_64BIT=-m64
    ;;
  GNU:powerpc*)
    SNL_CC_32BIT=-m32
    SNL_CC_64BIT=-m64
    ;;
  GNU:i?86)
    SNL_CC_32BIT=-m32
    SNL_CC_64BIT=-m64
    ;;
  GNU:mips*)
    SNL_CC_32BIT="-mips32 -mabi=32"
    SNL_CC_64BIT="-mips64 -mabi=64"
    ;;
  VisualAge:*)
    SNL_CC_32BIT=-q32
    SNL_CC_64BIT=-q64
    # Do V5.0 namemangling for compatibility with ACIS, and enable RTTI
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
  Sun:sparc*)
    SNL_CC_32BIT=-xarch=generic
    SNL_CC_64BIT=-xarch=generic64
    ;;
  *)
    ;;
esac
AC_MSG_RESULT([$cc_compiler:$host_cpu])
]) # end SNL_CC_FLAGS


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




# Check if template definitions (.cpp files) must be
# included (in the .hpp files).
# Sets TEMPLATE_DEFS_INCLUDED=1
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




