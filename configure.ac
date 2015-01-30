################################################################################
#                           Standard Stuff
################################################################################
AC_INIT([Lasso],[3.2.0])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([.])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE

AC_PROG_SED

AC_ARG_ENABLE([fortran],[AS_HELP_STRING([--disable-fortran],[No Fortran name mangling in ITAPS C headers])],
  [ENABLE_FORTRAN=$enableval],[ENABLE_FORTRAN=yes])

EXTRA_GNU_FLAGS='-Wall -pipe -pedantic -Wno-long-long'
EXTRA_INTEL_FLAGS='-Wall'
FATHOM_CHECK_COMPILERS([yes],[yes],[$ENABLE_FORTRAN])
AM_CONDITIONAL(PARALLEL,[test "x$WITH_MPI" == "xyes"])
LIBS="$LIBS -lm"

AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_DISABLE_SHARED
AC_PROG_FC
LT_INIT

FATHOM_CHECK_CXX_WORKS([], [AC_MSG_ERROR([Cannot build without C++ compiler])])

if test "xyes" = "x$ENABLE_FORTRAN" && test "x" != "x$FC"; then
  AC_FC_WRAPPERS
fi

FATHOM_COMPILER_FLAGS

################################################################################
#                              MPI OPTIONS
################################################################################

if test "x$WITH_MPI" == "xyes"; then
  DEFINES="$DEFINES -DUSE_MPI"
fi
AM_CONDITIONAL(USE_MPI, [test "xyes" == "x$WITH_MPI"])

################################################################################
#                           iGeom
################################################################################
AC_ARG_WITH(igeom, 
[AS_HELP_STRING([--with-igeom=DIR],[Specify directory containing iGeom libs/includes or dirs under which lib/ and include/ are located])],
[IGEOM_DIR=$withval
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-igeom=\"${withval}\""
 ], [IGEOM_DIR=])
case "x$IGEOM_DIR" in
  xno|x)
    IGEOM_MISSING=yes
    ;;
  xyes)
    AC_MSG_ERROR([--with-igeom requires an argument])
    ;;
  *)
    if ! test -d "$IGEOM_DIR"; then
      AC_MSG_WARN([$IGEOM_DIR : not a directory.])
    fi
    IGEOM_MISSING=no
    AC_CHECK_FILE([${IGEOM_DIR}/lib/iGeom-Defs.inc],
                  [IGEOM_CONFIG_OPTIONS="include ${IGEOM_DIR}/lib/iGeom-Defs.inc";
                   FATHOM_MAKE_INC_VAR([${IGEOM_DIR}/lib/iGeom-Defs.inc],[IGEOM_DEFINES], [IGEOM_DEFINES="$make_val"],
                      [IGEOM_DEFINES=])
                   DEFINES="$DEFINES -DIGEOM $IGEOM_DEFINES"],
                  [AC_CHECK_FILE([${IGEOM_DIR}/iGeom-Defs.inc],
                                 [IGEOM_CONFIG_OPTIONS="include ${IGEOM_DIR}/iGeom-Defs.inc";
                   DEFINES="$DEFINES -DIGEOM"],
                                 [AC_MSG_WARN([$IGEOM_DIR : not a configured iGeom]); 
                                  IGEOM_MISSING=yes])
                  ])
    ;;
esac
AC_SUBST(IGEOM_DIR)
AC_SUBST(IGEOM_CONFIG_OPTIONS)
AM_CONDITIONAL([ENABLE_IGEOM], [test x$IGEOM_MISSING = xno])

################################################################################
#                           FBiGeom
################################################################################
FBIGEOM_EXTRA_LIB= 
AC_ARG_WITH(fbigeom, 
[AS_HELP_STRING([--with-fbigeom=DIR],[Specify directory containing FBiGeom libs/includes or dirs under which lib/ and include/ are located])],
[FBIGEOM_DIR=$withval
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-fbigeom=\"${withval}\""
 ], [FBIGEOM_DIR=])
case "x$FBIGEOM_DIR" in
  xno|x)
    FBIGEOM_MISSING=yes
    ;;
  xyes)
    AC_MSG_ERROR([--with-fbigeom requires an argument])
    ;;
  *)
    if ! test -d "$FBIGEOM_DIR"; then
      AC_MSG_WARN([$FBIGEOM_DIR : not a directory.])
    fi
    FBIGEOM_MISSING=no
    AC_CHECK_FILE([${FBIGEOM_DIR}/lib/FBiGeom-Defs.inc],
                  [FBIGEOM_CONFIG_OPTIONS="include ${FBIGEOM_DIR}/lib/FBiGeom-Defs.inc";
                   DEFINES="$DEFINES -DFBIGEOM"; FBIGEOM_EXTRA_LIB="$FBIGEOM_EXTRA_LIB -lFBiGeomMOAB"],
                  [AC_CHECK_FILE([${FBIGEOM_DIR}/FBiGeom-Defs.inc],
                                 [FBIGEOM_CONFIG_OPTIONS="include ${FBIGEOM_DIR}/FBiGeom-Defs.inc";
                   DEFINES="$DEFINES -DFBIGEOM"],
                                 [AC_MSG_WARN([$FBIGEOM_DIR : not a configured iGeom]); 
                                  FBIGEOM_MISSING=yes])
                  ])
    ;;
esac
AC_SUBST(FBIGEOM_DIR)
AC_SUBST(FBIGEOM_CONFIG_OPTIONS)
AC_SUBST(FBIGEOM_EXTRA_LIB) 
AM_CONDITIONAL([ENABLE_FBIGEOM], [test x$FBIGEOM_MISSING = xno])

################################################################################
#                           iMesh
################################################################################
AC_ARG_WITH(imesh, 
[AS_HELP_STRING([--with-imesh=DIR],[Specify directory containing iMesh libs/includes or dirs under which lib/ and include/ are located])],
[IMESH_DIR=$withval
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-imesh=\"${withval}\""
 ], [IMESH_DIR=])
case "x$IMESH_DIR" in
  xno|x)
    IMESH_MISSING=yes
    ;;
  xyes)
    AC_MSG_ERROR([--with-imesh requires an argument])
    ;;
  *)
    if ! test -d "$IMESH_DIR"; then
      AC_MSG_WARN([$IMESH_DIR : not a directory.])
    fi
    IMESH_MISSING=no
    AC_CHECK_FILE([${IMESH_DIR}/lib/iMesh-Defs.inc],
                  [IMESH_CONFIG_OPTIONS="include ${IMESH_DIR}/lib/iMesh-Defs.inc";
                   DEFINES="$DEFINES -DIMESH"],
                  [AC_CHECK_FILE([${IMESH_DIR}/iMesh-Defs.inc],
                                 [IMESH_CONFIG_OPTIONS="include ${IMESH_DIR}/iMesh-Defs.inc";
                   DEFINES="$DEFINES -DIMESH"],
                                 [AC_MSG_WARN([$IMESH_DIR : not a configured iMesh]); 
                                  IMESH_MISSING=yes])
                  ])
    ;;
esac
AC_SUBST(IMESH_DIR)
AC_SUBST(IMESH_CONFIG_OPTIONS)
AM_CONDITIONAL([ENABLE_IMESH], [test x$IMESH_MISSING = xno])

################################################################################
#                           Output Files
################################################################################
AC_SUBST([DEFINES])
AC_SUBST([DISTCHECK_CONFIGURE_FLAGS])

AC_ARG_VAR([FC], [FORTRAN compiler command])
AC_CONFIG_HEADERS(config.h) # Don't want this, just working around autoheader issue
AC_CONFIG_HEADERS([itaps/iBase_FCDefs.h])
AC_CONFIG_FILES([Makefile 
                 test/Makefile
                 doc/lasso.dox
                 itaps/irel/iRel-Defs.inc
                 ])

# Generate iRel_FCDefs.h from iBase_FCDefs.h
AC_CONFIG_COMMANDS([itaps/irel/iRel_FCDefs.h],
  [sed -e "s/FC_FUNC/IREL_FC_FUNC/" itaps/iBase_FCDefs.h > itaps/irel/iRel_FCDefs.h])
# Remove iRel_protos.h if old format
AC_CONFIG_COMMANDS_POST([grep IREL_FC_FUNC itaps/irel/iRel_protos.h >/dev/null 2>&1 || rm -f itaps/irel/iRel_protos.h])
AC_OUTPUT
