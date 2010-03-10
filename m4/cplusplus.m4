#######################################################################################
# Check for existance unordered_map in either std:: or std::tr1:: namespace
# Executes first argument if found, second if not found.  If unordered_map
# is found, the variable 'result' will contain its namespace.
#######################################################################################
AC_DEFUN([MK_CHECK_UNORDERED_MAP],[

AC_CACHE_CHECK([for C++ unordered_map],
               [mk_cv_cxx_unordered_map],
               [AC_LANG_PUSH(C++)
                mk_cv_cxx_unordered_map=no
                AC_TRY_COMPILE([#include <unordered_map>],
                               [std::unordered_map<int,int> map],
                               [mk_cv_cxx_unordered_map="std"; incdir=])
                if test "xno" == "x$mk_cv_cxx_unordered_map"; then
                  AC_TRY_COMPILE([#include <tr1/unordered_map>],
                                 [std::tr1::unordered_map<int,int> map],
                                 [mk_cv_cxx_unordered_map="std::tr1"; incdir=tr1/])
                fi
                if test "xno" == "x$mk_cv_cxx_unordered_map"; then
                  AC_TRY_COMPILE([#include <boost/unordered_map>],
                                 [boost::unordered_map<int,int> map],
                                 [mk_cv_cxx_unordered_map="boost"; incdir=boost/])
                fi
                AC_LANG_POP(C++)])
if test "xno" != "x$mk_cv_cxx_unordered_map"; then
  result="$mk_cv_cxx_unordered_map"
  $1
else
  result=
  $2
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
AC_DEFUN([FATHOM_CANT_USE_STD], [

AC_CACHE_CHECK([for C++-standard library in std:: namespace],
               [snl_cv_cxx_stl_in_std],
               [AC_LANG_SAVE
                AC_LANG_CPLUSPLUS
                AC_TRY_COMPILE([#include <vector>
                                #include <string>
                                #include <map>
                                #include <algorithm>],
                               [std::vector<std::string> v;
                                std::map<std::string,std::string> m; ],
                               [snl_cv_cxx_stl_in_std=yes],
                               [snl_cv_cxx_stl_in_std=no])
                AC_LANG_RESTORE])
                
if test $snl_cv_cxx_stl_in_std = no; then
  CANT_USE_STD=-DCANT_USE_STD
else
  CANT_USE_STD=
fi

AC_CACHE_CHECK([for C++-standard I/O in std:: namespace],
               [snl_cv_cxx_io_in_std],
               [AC_LANG_SAVE
                AC_LANG_CPLUSPLUS
                AC_TRY_COMPILE([#include <iosfwd>
                                #include <iostream>
                                #include <ostream>
                                #include <sstream>],
                               [std::cout << std::endl;],
                               [snl_cv_cxx_io_in_std=yes],
                               [snl_cv_cxx_io_in_std=no])
                AC_LANG_RESTORE])
                
if test $snl_cv_cxx_io_in_std = no; then
  CANT_USE_STD_IO=-DCANT_USE_STD_IO
else
  CANT_USE_STD_IO=
fi

]) # FATHOM_CANT_USE_STD




#######################################################################################
# Check if template definitions (.cpp files) must be
# included (in the .hpp files).
# Sets TEMPLATE_DEFS_INCLUDED=-DTEMPLATE_DEFS_INCLUDED
#######################################################################################
AC_DEFUN([FATHOM_TEMPLATE_DEFS_INCLUDED], [

AC_CACHE_CHECK([if C++ template definitions should be included],
               [snl_cv_cxx_template_defs_included],
               [AC_LANG_SAVE
                AC_LANG_CPLUSPLUS

                src=conftest.cc
                templ=conftemp.cc
                exe=conftest

                echo "template <typename T> class MyTempl { public: T get() const; };" >$templ
                echo "template <typename T> T MyTempl<T>::get() const { return 0; }"  >>$templ
                echo "template <typename T> class MyTempl { public: T get() const; };" >$src
                echo "int main( ) { MyTempl<int> c; return c.get(); }"                >>$src
                if $CXX $CXXFLAGS $LDFLAGS $src $templ -o $exe >/dev/null 2>/dev/null; then
                  snl_cv_cxx_template_defs_included=no
                else
                  snl_cv_cxx_template_defs_included=yes
                fi

                rm -f $src $templ $exe
                AC_LANG_RESTORE])
                
if test $snl_cv_cxx_template_defs_included = yes; then
  TEMPLATE_DEFS_INCLUDED=-DTEMPLATE_DEFS_INCLUDED
else
  TEMPLATE_DEFS_INCLUDED=
fi

]) # FATHOM_TEMPLATE_DEFS_INCLUDED

#######################################################################################
# Check if compiler supports template class specialization.
# Sets TEMPLATE_SPECIALIZATION=-DTEMPLATE_SPECIALIZATION
#######################################################################################
AC_DEFUN([FATHOM_TEMPLATE_SPECIALIZATION], [
AC_CACHE_CHECK([if C++ compiler supports template class specialization],
               [snl_cv_template_specialization],
               [AC_LANG_SAVE
                AC_LANG_CPLUSPLUS
                AC_TRY_COMPILE([template <unsigned S> class MyTempl { public: char data[S]; };
                                template <> class MyTempl<0> { public: char value; };],
                               [MyTempl<1> one; MyTempl<0> zero; one.data[0] = zero.value = '\0';],
                               [snl_cv_template_specialization=yes],
                               [snl_cv_template_specialization=no])
                AC_LANG_RESTORE])

TEMPLATE_SPECIALIZATION=
if test $snl_cv_template_specialization = yes; then
  TEMPLATE_SPECIALIZATION=-DTEMPLATE_SPECIALIZATION
fi

]) # FATHOM_TEMPLATE_SPECIALIZATION

#######################################################################################
# Check if compiler supports template function specialization.
# Sets TEMPLATE_FUNC_SPECIALIZATION=-DTEMPLATE_FUNC_SPECIALIZATION
#######################################################################################
AC_DEFUN([FATHOM_TEMPLATE_FUNC_SPECIALIZATION], [
AC_CACHE_CHECK([if C++ compiler supports template function specialization],
               [snl_cv_template_func_specialization],
               [AC_LANG_SAVE
                AC_LANG_CPLUSPLUS
                AC_TRY_COMPILE([template <typename T> T templ_func( int i );
                                template <> int templ_func<int>( int j ) 
                                  { return j; }],
                               [return templ_func<int>( 0 );],
                               [snl_cv_template_func_specialization=yes],
                               [snl_cv_template_func_specialization=no])
                AC_LANG_RESTORE])

TEMPLATE_FUNC_SPECIALIZATION=
if test $snl_cv_template_func_specialization = yes; then
  TEMPLATE_FUNC_SPECIALIZATION=-DTEMPLATE_FUNC_SPECIALIZATION
fi

]) # FATHOM_TEMPLATE_FUNC_SPECIALIZATION
