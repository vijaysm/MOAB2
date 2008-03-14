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
# Sets TEMPLATE_DEFS_INCLUDED=-DTEMPLATE_DEFS_INCLUDED
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
# Check if compiler supports template class specialization.
# Sets TEMPLATE_SPECIALIZATION=-DTEMPLATE_SPECIALIZATION
#######################################################################################
AC_DEFUN([SNL_TEMPLATE_SPECIALIZATION], [
AC_LANG_SAVE
AC_LANG_CPLUSPLUS

AC_MSG_CHECKING([if C++ compiler supports template specialization])
AC_TRY_COMPILE([
template <unsigned S> class MyTempl { public: char data[S]; };
template <> class MyTempl<0> { public: char value; };
],[
MyTempl<1> one;
MyTempl<0> zero;
one.data[0] = zero.value = '\0';
],
[TEMPLATE_SPECIALIZATION=-DTEMPLATE_SPECIALIZATION; AC_MSG_RESULT(yes)],
[TEMPLATE_SPECIALIZATION=; AC_MSG_RESULT(no)])

AC_LANG_RESTORE
]) # SNL_TEMPLATE_DEFS_INCLUDED
