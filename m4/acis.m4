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
