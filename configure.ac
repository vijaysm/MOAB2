################################################################################
#                           Standard Stuff
################################################################################
AC_INIT(Lasso, 2.5)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([.])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE

AC_PROG_SED

AC_ARG_ENABLE([fortran],[AC_HELP_STRING([--disable-fortran],
  [No Fortran name mangling in ITAPS C headers])],
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
AC_PROG_LIBTOOL

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
[AC_HELP_STRING([--with-igeom=DIR], [Specify directory containing iGeom libs/includes or dirs under which lib/ and include/ are located])],
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
                   DEFINES="$DEFINES -DIGEOM"],
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

################################################################################
#                           iMesh
################################################################################
AC_ARG_WITH(imesh, 
[AC_HELP_STRING([--with-imesh=DIR], [Specify directory containing iMesh libs/includes or dirs under which lib/ and include/ are located])],
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

################################################################################
#                           Babel
################################################################################
AC_ARG_ENABLE( babel, AC_HELP_STRING([--enable-babel[=DIR]],[Enable Babel-based interface, with babel executable/libs/includes located under DIR]),
               [enable_babel=$enableval], [enable_babel=no] )  

USE_BABEL=no
if test "x$enable_babel" != "xno"; then
  if test "x$enable_babel" != "x"; then
     BABEL_DIR=/usr/local
  else
     BABEL_DIR=$enable_babel
  fi
  AC_CHECK_FILE([${BABEL_DIR}/bin/babel], [USE_BABEL=yes; 
         DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --enable-babel=\"${enableval}\""], [USE_BABEL=no])

  if test "x$IGEOM_MISSING" != "xno"; then
    AC_CHECK_FILE([${IGEOM_DIR}/lib/iGeom-SIDL-Defs.inc],
                [IGEOM_BABEL_CONFIG_OPTIONS="include ${IGEOM_DIR}/lib/iGeom-SIDL-Defs.inc"],
                [AC_CHECK_FILE([${IGEOM_DIR}/iGeom-SIDL-Defs.inc],
                               [IGEOM_BABEL_CONFIG_OPTIONS="include ${IGEOM_DIR}/iGeom-SIDL-Defs.inc"],
                               [AC_MSG_WARN([Can't enable Babel without a Babel-capable iGeom]); USE_BABEL=no])
                ])
  fi
  if test "x$IMESH_MISSING" != "xno"; then
    AC_CHECK_FILE([${IMESH_DIR}/lib/iMesh-SIDL-Defs.inc],
                [IMESH_BABEL_CONFIG_OPTIONS="include ${IMESH_DIR}/lib/iMesh-SIDL-Defs.inc"],
                [AC_CHECK_FILE([${IMESH_DIR}/iMesh-SIDL-Defs.inc],
                               [IMESH_BABEL_CONFIG_OPTIONS="include ${IMESH_DIR}/iMesh-SIDL-Defs.inc"],
                               [AC_MSG_WARN([Can't enable Babel without a Babel-capable iMesh]); USE_BABEL=no])
                ])
  fi
fi

AC_SUBST(BABEL_DIR)
AC_SUBST(USE_BABEL)
AC_SUBST(IGEOM_BABEL_CONFIG_OPTIONS)
AC_SUBST(IMESH_BABEL_CONFIG_OPTIONS)
AM_CONDITIONAL(USE_BABEL, [test "xno" != "x$USE_BABEL"])


################################################################################
#                           Output Files
################################################################################
AC_SUBST([INCLUDES])
AC_SUBST([DEFINES])
AC_SUBST([DISTCHECK_CONFIGURE_FLAGS])

AC_ARG_VAR([FC], [FORTRAN compiler command])
AC_CONFIG_HEADERS(config.h) # Don't want this, just working around autoheader issue
AC_CONFIG_HEADERS([iBase_FCDefs.h])
AC_CONFIG_FILES([Makefile 
                 iRel-Defs.inc
                 ])
#		 SIDL/iRel-SIDL-Defs.inc
#		 SIDL/Makefile

# Generate iRel_FCDefs.h from iBase_FCDefs.h
AC_CONFIG_COMMANDS([iRel_FCDefs.h],
  [sed -e "s/FC_FUNC/IREL_FC_FUNC/" iBase_FCDefs.h >iRel_FCDefs.h])
# Remove iRel_protos.h if old format
AC_CONFIG_COMMANDS_POST([grep IREL_FC_FUNC iRel_protos.h >/dev/null 2>&1 || rm -f iRel_protos.h])
AC_OUTPUT

if test "x$IGEOM_MISSING" = "xyes"; then
  AC_MSG_WARN([iRel build will fail!  iGeom not found!])
fi

if test "x$IMESH_MISSING" = "xyes"; then
  AC_MSG_WARN([iRel build will fail!  iMesh not found!])
fi

