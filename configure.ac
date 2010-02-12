################################################################################
#                           Standard Stuff
################################################################################
AC_INIT(MOAB, 3.99)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([.])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE

AC_CHECK_PROG( [SED], [sed], [sed], [true] )


################################################################################
#                           iMesh Part 1 of 3
################################################################################

CONFIGURE_FORTRAN=no
FORCE_ENABLE_IMESH=no
AC_ARG_ENABLE( [imesh],
[AC_HELP_STRING([--enable-imesh(=nofortran)],[Build iMesh interface.  If optional
'nofortran' argument is given, support for FORTRAN name mangling will be disabled.])
AC_HELP_STRING([--disable-imesh],[Don't build support for iMesh interface.])],
[FORCE_ENABLE_IMESH=yes
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --enable-imesh=\"${enableval}\""],
[enableval=yes] )

case $enableval in
  yes)
    ENABLE_imesh=yes
    CONFIGURE_FORTRAN=yes
    ;;
  nofortran|NoFortran|noFORTRAN|NoFORTRAN|Nofortran)
    ENABLE_imesh=yes
    ;;
  no)
    ENABLE_imesh=no
    ;;
  *)
    AC_MSG_ERROR([Invalid argument to --enable-imesh : $enableval])
    ;;
esac

################################################################################
#                           iGeom Part 1 of 3
################################################################################

CONFIGURE_FORTRAN=no
FORCE_ENABLE_IGEOM=no
AC_ARG_ENABLE( [igeom],
[AC_HELP_STRING([--enable-igeom(=nofortran)],[Build iGeom interface.  If optional
'nofortran' argument is given, support for FORTRAN name mangling will be disabled.])
AC_HELP_STRING([--disable-igeom],[Don't build support for iGeom interface (default).])],
[enableval="$enableval"],[enableval=no] )

case $enableval in
  yes)
    ENABLE_igeom=yes
    CONFIGURE_FORTRAN=yes
    ;;
  nofortran|NoFortran|noFORTRAN|NoFORTRAN|Nofortran)
    ENABLE_igeom=yes
    ;;
  no)
    ENABLE_igeom=no
    ;;
  *)
    AC_MSG_ERROR([Invalid argument to --enable-igeom : $enableval])
    ;;
esac

################################################################################
#                           Compilers
################################################################################

EXTRA_GNU_FLAGS='-Wall -pipe'
FATHOM_CHECK_COMPILERS([yes],[yes],[$CONFIGURE_FORTRAN])
AM_CONDITIONAL(PARALLEL,[test "x$WITH_MPI" != "xno"])
LIBS="-lm"

AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_DISABLE_SHARED
AC_PROG_LIBTOOL
FATHOM_COMPILER_FLAGS

AC_CHECK_PROG([ZCAT],[gunzip],[gunzip -c],[])
AC_CHECK_PROG([ZCAT],[gzip],[gzip -cd],[])
AC_CHECK_PROG([ZCAT],[zcat],[zcat],[])
if test "x" = "x$ZCAT"; then
  AC_MSG_WARN([[Cannot run tests: no way to uncompress input files.]])
fi
AC_SUBST([ZCAT])
AM_CONDITIONAL([HAVE_ZCAT],[test "x" = "x$ZCAT"])

FATHOM_TEMPLATE_SPECIALIZATION
FATHOM_TEMPLATE_FUNC_SPECIALIZATION
DEFINES="$DEFINES $TEMPLATE_SPECIALIZATION $TEMPLATE_FUNC_SPECIALIZATION"

if test "xno" != "x$CONFIGURE_FORTRAN"; then
  if test "x" != "x$FC"; then
    AC_FC_WRAPPERS
  fi
fi


################################################################################
#                              HANDLE SIZE
################################################################################
AC_CHECK_HEADER([inttypes.h], [AC_DEFINE(MOAB_HAVE_INTTYPES_H,  [], [MOAB qualified HAVE_INTTYPES_H]) ])
AC_CHECK_HEADER([stdint.h],   [AC_DEFINE(MOAB_HAVE_STDINT_H,    [], [MOAB qualified HAVE_STDINT_H])   ])
AC_CHECK_HEADER([sys/types.h],[AC_DEFINE(MOAB_HAVE_SYS_TYPES_H, [], [MOAB qualified HAVE_SYS_TYPES_H])])
AC_CHECK_HEADER([stddef.h],   [AC_DEFINE(MOAB_HAVE_STDDEF_H,    [], [MOAB qualified HAVE_STDDEF_H])   ])
AC_CHECK_HEADER([stdlib.h],   [AC_DEFINE(MOAB_HAVE_STDLIB_H,    [], [MOAB qualified HAVE_STDLIB_H])   ])
SIZEOF_MBENTITYHANDLE=0
AC_CHECK_SIZEOF(void *)
SIZEOF_VOID_P=$ac_cv_sizeof_void_p
AC_ARG_ENABLE([64bit-handles],
  [AC_HELP_STRING([--enable-64bit-handles],[Force MBEntityHandle to be 64 bits])],
  [AC_CHECK_TYPE([uint64_t],
    [AC_DEFINE(MOAB_FORCE_64_BIT_HANDLES,[],[Use int64_t for handles])
     SIZEOF_MBENTITYHANDLE=8
     DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --enable-64bit-handles"],
    [AC_MSG_ERROR([Platform does not have uint64_t. Cannot force 64-bit handles])])
])
AC_ARG_ENABLE([32bit-handles],
  [AC_HELP_STRING([--enable-32bit-handles],[Force MBEntityHandle to be 32 bits])],
  [AC_CHECK_TYPE([uint32_t],
    [AC_DEFINE(MOAB_FORCE_32_BIT_HANDLES,[],[Use int32_t for handles])
     SIZEOF_MBENTITYHANDLE=4
     DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --enable-32bit-handles"],
    [AC_MSG_ERROR([Platform does not have uint32_t. Cannot force 32-bit handles])])
])
if test "x$SIZEOF_MBENTITYHANDLE" = "x0"; then
  SIZEOF_MBENTITYHANDLE=$SIZEOF_VOID_P
  AC_CHECK_TYPE([size_t], 
   [AC_DEFINE(MOAB_HAVE_SIZE_T,[],[System provides size_t typedef])],
   [AC_CHECK_SIZEOF(unsigned long)
    if test "$SIZEOF_UNSIGNED_LONG" -ne "$SIZEOF_VOID_P"; then
      AC_MSG_ERROR([Cannot define handle with sizeof void*: size_t not found and long is of incorrect size.])
    fi]
  )
  AC_CHECK_TYPE([ptrdiff_t], 
   [AC_DEFINE(MOAB_HAVE_PTRDIFF_T,[],[System provides ptrdiff_t typedef])],
   [AC_CHECK_SIZEOF(long)
    AC_CHECK_SIZEOF(void *)
    if test "$SIZEOF_LONG" -ne "$SIZEOF_VOID_P"; then
      AC_MSG_ERROR([Cannot define handle with sizeof void*: ptrdiff_t not found and long is of incorrect size.])
    fi]
  )
fi


################################################################################
#                           iMesh Part 2 of 3
################################################################################

# The iMesh option must be checked before configuring the compilers
# so that we know if we need Fortran.  This stuff has to be done
# after configuring the compilers so that we know what size of 
# various types are.

AC_MSG_CHECKING([size of MBEntityHandle])
AC_MSG_RESULT([$SIZEOF_MBENTITYHANDLE])
IMESH_VALID_HANDLE_SIZE=yes
if test $SIZEOF_VOID_P -ne $SIZEOF_MBENTITYHANDLE; then
  IMESH_VALID_HANDLE_SIZE=no
  AC_MSG_WARN([Cannot build iMesh API because sizeof(MBEntityHandle) != sizeof(void*)])
  if test "x$FORCE_ENABLE_IMESH" = "xno"; then
    ENABLE_imesh=no
  fi
fi

AC_MSG_CHECKING([if iMesh support is to be built])
AC_MSG_RESULT([$ENABLE_imesh])
AM_CONDITIONAL([ENABLE_imesh],[test "xyes" = "x$ENABLE_imesh"])
if test "x$ENABLE_imesh" = "xyes"; then
  if test "x$IMESH_VALID_HANDLE_SIZE" != "xyes"; then
    AC_MSG_ERROR([Cannot build iMesh API with sizeof(MBEntityHandle) != sizeof(void*)])
  fi
fi

################################################################################
#                           iGeom Part 2 of 3
################################################################################

# The iGeom option must be checked before configuring the compilers
# so that we know if we need Fortran.  This stuff has to be done
# after configuring the compilers so that we know what size of 
# various types are.

AC_MSG_CHECKING([size of MBEntityHandle])
AC_MSG_RESULT([$SIZEOF_MBENTITYHANDLE])
IGEOM_VALID_HANDLE_SIZE=yes
if test $SIZEOF_VOID_P -ne $SIZEOF_MBENTITYHANDLE; then
  IMESH_VALID_HANDLE_SIZE=no
  AC_MSG_WARN([Cannot build iMesh API because sizeof(MBEntityHandle) != sizeof(void*)])
  if test "x$FORCE_ENABLE_IGEOM" = "xno"; then
    ENABLE_igeom=no
  fi
fi

AC_MSG_CHECKING([if iGeom support is to be built])
AC_MSG_RESULT([$ENABLE_igeom])
AM_CONDITIONAL([ENABLE_igeom],[test "xyes" = "x$ENABLE_igeom"])
if test "x$ENABLE_igeom" = "xyes"; then
  if test "x$IGEOM_VALID_HANDLE_SIZE" != "xyes"; then
    AC_MSG_ERROR([Cannot build iGeom API with sizeof(MBEntityHandle) != sizeof(void*)])
  fi
fi

################################################################################
#                              Version Stuff
################################################################################

AC_DEFINE(MB_VERSION,["AC_PACKAGE_VERSION"],[MOAB Version])
VERSION_MAJOR=`expr AC_PACKAGE_VERSION : '\([[0-9]]*\)'`
VERSION_MINOR=`expr AC_PACKAGE_VERSION : '[[0-9]]*\.\([[0-9]]*\)'`
VERSION_PATCH=`expr AC_PACKAGE_VERSION : '[[0-9]]*\.[[0-9]]*\.\(.*\)'`
test "x" != "x$VERSION_MAJOR" || AC_MSG_ERROR("Invalid version string: AC_PACKAGE_VERSION")
test "x" != "x$VERSION_MINOR" || AC_MSG_ERROR("Invalid version string: AC_PACKAGE_VERSION")
AC_DEFINE_UNQUOTED(MB_VERSION_MAJOR,$VERSION_MAJOR,[MOAB Major Version])
AC_DEFINE_UNQUOTED(MB_VERSION_MINOR,$VERSION_MINOR,[MOAB Minor Version])
if test "x" != "x$VERSION_PATCH"; then
  AC_DEFINE_UNQUOTED(MB_VERSION_PATCH,$VERSION_PATCH,[MOAB Patch Level])
  VERSION_STRING="AC_PACKAGE_STRING"
elif test $VERSION_MINOR -eq 99; then
  VERSION_STRING="AC_PACKAGE_STRING (Alpha)"
else
  VERSION_STRING="AC_PACKAGE_STRING (Beta)"
fi
AC_DEFINE_UNQUOTED(MB_VERSION_STRING,"${VERSION_STRING}",[MOAB Version String])


################################################################################
#                              Extract libtool config
################################################################################
ITAPS_LIBTOOL_VAR( [CXX], [compiler_lib_search_path], [MOAB_CXX_LINKFLAGS])
ITAPS_LIBTOOL_VAR( [CXX], [postdeps], [MOAB_CXX_LIBS])
AC_SUBST(MOAB_CXX_LINKFLAGS)
AC_SUBST(MOAB_CXX_LIBS)

################################################################################
#                             System Headers
################################################################################
AC_CHECK_HEADER([ieeefp.h],[DEFINES="$DEFINES -DHAVE_IEEEFP_H"] )
if test "xyes" = "x$enable_debug"; then
  AC_CHECK_HEADER( [valgrind/memcheck.h], [CPPFLAGS="$CPPFLAGS -DVALGRIND"] )
fi

################################################################################
#                             BOOST OPTIONS
################################################################################

dnl FATHOM_CHECK_BOOST( [boost/pool/object_pool.hpp] )
MK_CHECK_UNORDERED_MAP([CPPFLAGS="$CPPFLAGS -DUNORDERED_MAP_NS=$result -DHAVE_UNORDERED_MAP=${incdir}unordered_map"])

################################################################################
#                              MPI OPTIONS
################################################################################

FATHOM_CHECK_MPI
FATHOM_CONFIG_MPI_EXEC
test "xyes" != "x$WITH_MPI" || DEFINES="$DEFINES -DUSE_MPI"
AM_CONDITIONAL(USE_MPI, [test "xno" != "x$WITH_MPI"])
# Used to generate MBmpi.h
if test "xyes" = "x$MPI_CXX_HELP_NEEDED"; then
  AC_DEFINE( [MB_MPI_CXX_CONFLICT], [1], 
    [MPICH_IGNORE_CXX_SEEK is not sufficient to avoid conflicts] )
  AC_DEFINE_UNQUOTED([MB_SEEK_SET],[$SEEK_SET],["Value of C SEEK_SET"])
  AC_DEFINE_UNQUOTED([MB_SEEK_CUR],[$SEEK_CUR],["Value of C SEEK_CUR"])
  AC_DEFINE_UNQUOTED([MB_SEEK_END],[$SEEK_END],["Value of C SEEK_END"])
fi

################################################################################
#                              HDF5 OPTIONS
################################################################################

FATHOM_CHECK_HDF5
if test "xno" = "x$HAVE_HDF5"; then
  AC_MSG_WARN([Support for native HDF5 file format disabled])
else
  DEFINES="$DEFINES -DHDF5_FILE"
fi
AM_CONDITIONAL(HDF5_FILE, [test "xno" != "x$HAVE_HDF5"])
AM_CPPFLAGS="$HDF5_CPPFLAGS $AM_CPPFLAGS"
EXPORT_LDFLAGS="$EXPORT_LDFLAGS $HDF5_LDFLAGS"
AC_SUBST(HDF5_LIBS)

HAVE_HDF5_PARALLEL=no
if test "xno" != "x$HAVE_HDF5"; then
  if test "xno" != "x$WITH_MPI"; then
    old_LDFLAGS="$LDFLAGS"
    LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
    AC_CHECK_LIB( [hdf5], [H5Pset_fapl_mpio], [HAVE_HDF5_PARALLEL=yes],
      [AC_MSG_WARN("libhdf5 library does not include parallel support.  Parallel HDF5 I/O disabled")],
      [$HDF5_LIBS]
      )
    LDFLAGS="$old_LDFLAGS"
  fi
fi
AM_CONDITIONAL(PARALLEL_HDF5, [test "xno" != "x$HAVE_HDF5_PARALLEL"])
if test "xno" != "x$HAVE_HDF5_PARALLEL"; then
  DEFINES="$DEFINES -DHDF5_PARALLEL"
fi

################################################################################
#                             CCMIO OPTIONS
################################################################################

old_LDFLAGS="$LDFLAGS"
LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
FATHOM_CHECK_CCMIO
LDFLAGS="$old_LDFLAGS"
if test "xno" != "x$HAVE_CCMIO"; then
  DEFINES="$DEFINES -DCCMIO_FILE"
fi
AM_CONDITIONAL(CCMIO_FILE, [test "xno" != "x$HAVE_CCMIO"])
AM_CPPFLAGS="$CCMIO_CPPFLAGS $AM_CPPFLAGS"
EXPORT_LDFLAGS="$CCMIO_LDFLAGS $EXPORT_LDFLAGS"
AC_SUBST(CCMIO_LIBS)


################################################################################
#                             NetCDF OPTIONS
################################################################################

old_LDFLAGS="$LDFLAGS"
LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
FATHOM_CHECK_NETCDF
LDFLAGS="$old_LDFLAGS"
if test "xno" = "x$HAVE_NETCDF"; then
  AC_MSG_WARN([Support for ExodusII file format disabled])
else
  DEFINES="$DEFINES -DNETCDF_FILE"
fi
AM_CONDITIONAL(NETCDF_FILE, [test "xno" != "x$HAVE_NETCDF"])
AM_CPPFLAGS="$NETCDF_CPPFLAGS $AM_CPPFLAGS"
EXPORT_LDFLAGS="$NETCDF_LDFLAGS $EXPORT_LDFLAGS"
AC_SUBST(NETCDF_LIBS)


#################################################################################
#                             Documentation
#################################################################################
AC_ARG_ENABLE([doxygen],
[AC_HELP_STRING([[--enable-doxygen(=DIR)]],[Specify directory where Doxygen program is installed])
AC_HELP_STRING([--disable-doxygen],[Do not generate API documentation (default)])],
                        [ENABLE_DOXYGEN="$enableval"],[ENABLE_DOXYGEN=no] )
if test "x$ENABLE_DOXYGEN" = "xyes"; then
  AC_PATH_PROG( [DOXYGEN], [doxygen], [no] )
elif test "x$ENABLE_DOXYGEN" != "xno"; then
  AC_PATH_PROG( [DOXYGEN], [doxygen], [no], [$ENABLE_DOXYGEN] )
fi
if test "x$ENABLE_DOXYGEN" != "xno"; then
  if test "x$DOXYGEN" = "xno"; then
    AC_MSG_ERROR("Doxygen executable not found.")
  fi
fi
AC_SUBST([DOXYGEN])
AM_CONDITIONAL([ENABLE_DOXYGEN],[test "x$ENABLE_DOXYGEN" != "xno"])



###############################################################################
#                           Optional Tools
###############################################################################
  # Define a macro to avoid typing this for each individual tool
  # Usage: MB_OPTIONAL_TOOL( name, default )
  #  name      - name of option
  #  default   - yes/no, $ENABLE_TOOLS overrides this if set
  # Actions:
  #  sets ENABLE_${tool} to 'yes' or 'no'
  #  creates ENABLE_${tool} automake conditional
AC_DEFUN([MB_OPTIONAL_TOOL],[
  mb_default_$1=$2
  if test "x" != "x$ENABLE_TOOLS"; then
    mb_default_$1="$ENABLE_TOOLS"
  fi
  AC_ARG_ENABLE( [$1],
[AC_HELP_STRING([--enable-$1],[Build tool: $1])
AC_HELP_STRING([--disable-$1],[Don't build $1])],
                 [ENABLE_$1=$enableval
DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --enable-$1=\"${enableval}\""
],[ENABLE_$1=${mb_default_$1}] )
  AM_CONDITIONAL([ENABLE_$1],[test "x${ENABLE_$1}" != "xno"])
  AC_MSG_CHECKING([if $1 is to be built])
  AC_MSG_RESULT([${ENABLE_$1}])
])
  # The default for all tools
AC_ARG_ENABLE( [tools], [
AC_HELP_STRING([--enable-tools],[Build all tools by default])
AC_HELP_STRING([--disable-tools],[Disable all tools by default])],
                         [ENABLE_TOOLS=$enableval],[ENABLE_TOOLS=] )
  # Individual tools
MB_OPTIONAL_TOOL([mbconvert],    [yes])
MB_OPTIONAL_TOOL([hexmodops],    [yes])
MB_OPTIONAL_TOOL([mbchaco],      [no] )
MB_OPTIONAL_TOOL([mbperf],       [yes])
MB_OPTIONAL_TOOL([qvdual],       [no] )
MB_OPTIONAL_TOOL([mbsize],       [yes])
MB_OPTIONAL_TOOL([mbskin],       [yes])
MB_OPTIONAL_TOOL([mbtagprop],    [yes])
MB_OPTIONAL_TOOL([mcnpmit],      [no])
MB_OPTIONAL_TOOL([mbcoupler],    [$ENABLE_mcnpmit] )
MB_OPTIONAL_TOOL([spheredecomp], [yes])
MB_OPTIONAL_TOOL([mbsurfplot],   [yes])
MB_OPTIONAL_TOOL([mbzoltan],     [no] )
MB_OPTIONAL_TOOL([dagmc],        [yes])
MB_OPTIONAL_TOOL([gsets],        [yes])
MB_OPTIONAL_TOOL([mbdepth],      [yes])


################################################################################
#                           Zoltan
################################################################################
if test $ENABLE_mbzoltan = yes; then
  if test "x$WITH_MPI" = "xno"; then
     AC_MSG_ERROR([--enable-mbzoltan requires MPI, try specifying --with-mpi])
  fi
  

  ZOLTAN_LIB_FLAGS=
  ZOLTAN_INC_FLAGS=
  ZOLTAN_LIBS=
  ZOLTAN_ARCH=
  ZOLTAN_DIR=

  AC_ARG_WITH( [zoltan], 
               [AC_HELP_STRING([--with-zoltan=DIR],[Specify Zoltan build directory])],
               [if test "xyes" = "$withval"; then 
                  AC_MSG_WARN([Ingoring --with-zoltan flag w/out value])
                elif test -d $withval/$include; then
                  ZOLTAN_INC_FLAGS="$ZOLTAN_INC_FLAGS -I$withval/include"
                  ZOLTAN_DIR="$withval"
                else
                  AC_MSG_ERROR("$withval/include : not a directory")
                fi
                DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-zoltan=\"${withval}\""
               ])

  AC_ARG_WITH([zoltan-arch],
   [AC_HELP_STRING([--with-zoltan-arch=OS],[Specify zoltan architecture])],
   [if test "xyes" = "x$withval"; then
      AC_MSG_ERROR([--with-zoltan-arch flag requires a value])
    elif test "x" = "x$ZOLTAN_DIR"; then
      AC_MSG_WARN([Ingoring --with-zoltan-arch without --with-zoltan])
    elif test -d "$ZOLTAN_DIR/Obj_$withval"; then
      ZOLTAN_LIB_FLAGS="$ZOLTAN_LIB_FLAGS -L$ZOLTAN_DIR/Obj_$withval"
    else
      AC_MSG_ERROR("$ZOLTAN_DIR/Obj_$withval : not a directory")
    fi],
   [if test "x" != "x$ZOLTAN_DIR"; then
      AC_MSG_WARN([Expect --with-zoltan-arch with --with-zoltan])
      if test -d $ZOLTAN_DIR/lib; then
        AC_MSG_WARN("Guessing '$ZOLTAN_DIR/lib'")
        ZOLTAN_LIB_FLAGS="$ZOLTAN_LIB_FLAGS -L$ZOLTAN_DIR/lib"
      else
        # if only one Obj_$ARCH dir, then us it
        found=no
        for dir in $ZOLTAN_DIR/Obj_*; do
          if test -d $dir; then
            if test "xyes" = "x$found"; then
              AC_MSG_ERROR("Multiple Zoltan Obj_$ARCH dirs.  --with-zoltan-arch required.")
            fi
            ZOLTAN_LIB_FLAGS="$ZOLTAN_LIB_FLAGS -L$dir"
            found=yes
          fi
        done
        if test "xno" = "x$found"; then
          AC_MSG_ERROR('No Zoltan Obj_$ARCH dir found.  --with-zoltan-arch required.')
        fi
      fi
    fi
   ])

  AC_ARG_WITH([parmetis], 
              [AC_HELP_STRING([--with-parmetis=DIR], 
                [Directory containing Parmetis library for Zoltan])],
              [if test "xyes" = "$withval"; then 
                 AC_MSG_WARN([Ingoring --with-parmetis flag w/out value])
               elif test -d "$withval"; then
                 ZOLTAN_LIB_FLAGS="$ZOLTAN_LIB_FLAGS -L$withval"
               else
                 AC_MSG_ERROR("$withval: not a directory")
               fi]
               DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-parmetis=\"${withval}\""
             )

  old_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $ZOLTAN_INC_FLAGS -DMPICH_IGNORE_CXX_SEEK"
  AC_LANG_PUSH([C++])
  AC_CHECK_HEADER([zoltan_cpp.h],[],[AC_MSG_ERROR("mbzoltan requires zoltan")])
  AC_LANG_POP([C++])
  CPPFLAGS="$old_CPPFLAGS"
  
  old_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $ZOLTAN_LIB_FLAGS"
  AC_CHECK_LIB([zoltan],[Zoltan_LB_Partition],[ZOLTAN_LIBS='-lzoltan'],
    [unset ac_cv_lib_zoltan_Zoltan_LB_Partition
     AC_CHECK_LIB([zoltan],[Zoltan_LB_Partition],
                  [ZOLTAN_LIBS="-lzoltan -lparmetis -lmetis"],
                  [AC_MSG_ERROR("mbzoltan requires zoltan")],
                  [-lparmetis -lmetis])])
  LDFLAGS="$old_LDFLAGS"
  

  AC_SUBST(ZOLTAN_LIBS)
  AC_SUBST(ZOLTAN_LIB_FLAGS)
  AC_SUBST(ZOLTAN_INC_FLAGS)
fi

################################################################################
#                           Chaco
################################################################################
if test $ENABLE_mbchaco = yes; then
  AC_ARG_WITH(chaco, 
  [AC_HELP_STRING([--with-chaco=DIR], [Specify directory containing chaco source and library])],
  [CHACO_DIR=$withval
   DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-chaco=\"${withval}\""
   ], [CHACO_DIR=])
  case "x$CHACO_DIR" in
    xno|x)
      CHACO_MISSING=yes
      ;;
    xyes)
      AC_MSG_ERROR([--with-chaco requires and argument])
      ;;
  esac
  if ! test -d "$CHACO_DIR"; then
    AC_MSG_WARN([$CHACO_DIR : not a directory.])
  fi
  AC_CHECK_FILES([${CHACO_DIR}/code/main/defs.h ${CHACO_DIR}/code/main/params.h],
                [], [AC_MSG_WARN([$CHACO_DIR : not a valid chaco source]); CHACO_MISSING=yes])
  old_LIBS="$LIBS"
  LIBS="-L${CHACO_DIR}/code"
  AC_CHECK_LIB([chaco],[read_params],[],
               [AC_MSG_WARN([Could not find libchaco in ${CHACO_DIR}/code]); CHACO_MISSING=yes],
               [-lm])
  CHACO_INCLUDES="-I${CHACO_DIR}/code/main"
  LIBS="$old_LIBS"
  CHACO_LIBS="-L${CHACO_DIR}/code -lchaco -lm"
  AC_SUBST(CHACO_LIBS)
  AC_SUBST(CHACO_INCLUDES)
fi


################################################################################
#                               qvdual
################################################################################

GRAPHVIZ_MISSING=no
X_MISSING=no
VTK_MISSING=no
QT_MISSING=no

if test $ENABLE_qvdual = yes; then


################################################################################
#                               GraphViz
################################################################################
  AC_ARG_WITH([graphviz],
  [AC_HELP_STRING([--with-graphviz=DIR], [Specify directory containing GraphViz])],
  [GRAPH_DIR=$withval
   DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-graphviz=\"${withval}\""
  ],[GRAPH_DIR=])

  case "x$GRAPH_DIR" in
    xno)
      GRAPHVIZ_MISSING=yes
      ;;
    xyes|x)
      GRAPH_INC_DIR=
      GRAPH_DIR=
      for dir in /usr /usr/local; do
        if test "x" = "x$GRAPH_DIR"; then
          AC_CHECK_FILE([${dir}/include/graphviz/gvc.h],[GRAPH_DIR="$dir"])
        fi
      done
      AC_MSG_CHECKING([for GraphViz include dir])
      if test "x" != x"$GRAPH_DIR"; then
        GRAPH_INC_DIR="-I${GRAPH_DIR}/include/graphviz"
      fi
      AC_MSG_RESULT([$GRAPH_INC_DIR])
      AC_MSG_CHECKING([for GraphViz library dir])
      GRAPH_LIB_DIR=
      for dir in ${GRAPH_DIR}/lib64/graphviz ${GRAPH_DIR}/lib/graphviz /usr/lib64/graphviz /usr/lib/graphviz /usr/local/lib64/graphviz /usr/local/lib/graphviz /lib64/graphviz /lib/graphviz; do
        if test "x" = "x$GRAPH_LIB_DIR"; then
          if test -d $dir; then
            GRAPH_LIB_DIR=-L$dir
          fi
        fi
      done
      AC_MSG_RESULT([$GRAPH_LIB_DIR])
      ;;
    x?*)
      GRAPH_INC_DIR=
      for dir in ${GRAPH_DIR}/include/graphviz ${GRAPH_DIR}/include ${GRAPH_DIR}; do
        if test "x" = "x$GRAPH_INC_DIR"; then
          AC_CHECK_FILE([${dir}/gvc.h],[GRAPH_INC_DIR="-I$dir"])
        fi
      done
      AC_MSG_CHECKING([for GraphViz include dir])
      if test "x" = x"$GRAPH_INC_DIR"; then
        AC_MSG_RESULT([no])
        GRAPHVIZ_MISSING=yes
      else
        AC_MSG_RESULT([$GRAPH_INC_DIR])
      fi
      AC_MSG_CHECKING([for GraphViz library dir])
      GRAPH_LIB_DIR=
      for dir in ${GRAPH_DIR}/lib/graphviz ${GRAPH_DIR}/lib ${GRAPH_DIR}; do
        if test "x" = "x$GRAPH_LIB_DIR"; then
          if test -d $dir; then
            GRAPH_LIB_DIR=-L$dir
          fi
        fi
      done
      if test "x" = x"$GRAPH_LIB_DIR"; then
        AC_MSG_RESULT([no])
        GRAPHVIZ_MISSING=yes
      else
        AC_MSG_RESULT([$GRAPH_LIB_DIR])
      fi
      ;;
  esac

  GRAPH_LIBS="$GRAPH_LIB_DIR -lgvc -lgraph"
  old_CXXFLAGS="$CXXFLAGS"
  old_CPPFLAGS="$CPPFLAGS"
  CXXFLAGS="$GRAPH_INC_DIR"
  CPPFLAGS="$GRAPH_INC_DIR"
  old_LIBS="$LIBS"
  LIBS=$GRAPH_LIBS
  AC_CHECK_HEADERS([gvc.h],[],[AC_MSG_WARN([Missing GraphViz header]); GRAPHVIZ_MISSING=yes])
  AC_CHECK_LIB([gvc],[gvContext],[],[AC_MSG_WARN([Cannot find GraphViz library: -lgraph]); GRAPHVIZ_MISSING=yes])
  LIBS="$old_LIBS"
  CXXFLAGS="$old_CXXFLAGS"
  CPPFLAGS="$old_CPPFLAGS"


################################################################################
#                                   X
################################################################################
  AC_PATH_X
  AC_PATH_XTRA
  if test "xyes" = "x$no_x"; then
    X_MISSING=yes
  fi


################################################################################
#                                  VTK
################################################################################
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  
  AC_ARG_WITH([vtk],
  [AC_HELP_STRING([--with-vtk=DIR], [Specify directory containing VTK])],
  [VTK_DIR=$withval
   DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-vtk=\"${withval}\""
  ],[VTK_DIR=])

  VTK_INCLUDES=  
  case "x$VTK_DIR" in
    xno)
      VTK_MISSING=yes
      ;;
    xyes|x)
      VTK_LIB_DIR=
      VTK_INC_DIR=
      for dir in /usr/include/vtk /usr/local/include/vtk /include/vtk /usr/include/vtk-5.0; do
        if test "x" = "x$VTK_INC_DIR"; then
          AC_CHECK_FILE([${dir}/vtkRenderer.h],[VTK_INC_DIR="-I$dir"])
        fi
      done
      if test x = x"$VTK_INC_DIR"; then
        AC_CHECK_HEADER([vtkRenderer.h],[],
                        [VTK_MISSING=yes])
      else
        VTK_INCLUDES="$VTK_INC_DIR"
      fi
      old_LIBS="$LIBS"
      for dir in /usr/lib/vtk /usr/lib /usr/bin /usr; do
        if test "x" = "x$VTK_LIB_DIR"; then
          LIBS="$old_LIBS -L$dir"
          unset ac_cv_lib_vtkCommon
	  unset ac_cv_lib_vtkCommon_main
          AC_HAVE_LIBRARY([vtkCommon],[VTK_LIB_DIR=$dir])
        fi
      done
      if test "x" = "x$VTK_LIB_DIR"; then
        VTK_MISSING=yes
      fi
      LIBS="$old_LIBS"
      ;;
    x?*)
      VTK_LIB_DIR=
      VTK_INC_DIR=
      for dir in ${VTK_DIR}/include/vtk ${VTK_DIR}/include/vtk/Rendering ${VTK_DIR}/include ${VTK_DIR}; do
        if test "x" = "x$VTK_INC_DIR"; then
          AC_CHECK_FILE([${dir}/vtkRenderer.h],[VTK_INC_DIR="-I$dir"])
        fi
      done
      if test "x" = "x$VTK_INC_DIR"; then
        # Check for non-installed VTK source
        AC_CHECK_FILE([${VTK_DIR}/Rendering/vtkRenderer.h],[
          for subdir in . Common Filtering GenericFiltering Graphics Hybrid Imaging Parallel Patented Rendering Utilities Utilities/vtkexpat Utilities/expat Utilities/DICOMParser Utilties/vtkpng Utilities/vtkzlib Utilities/freetype/include/freetype Wrapping; do
  	  AC_MSG_CHECKING([for directory ${VTK_DIR}/$subdir])
  	  if test -d ${VTK_DIR}/$subdir; then
  	    VTK_INC_DIR="$VTK_INC_DIR -I${VTK_DIR}/$subdir"
  	    AC_MSG_RESULT([yes])
  	  else
  	    AC_MSG_RESULT([no])
  	  fi
  	done
  	],[VTK_MISSING=yes]
        )
      fi
      VTK_INCLUDES="$VTK_INC_DIR"
      old_LIBS="$LIBS"
      for dir in ${VTK_DIR}/lib/vtk ${VTK_DIR}/lib ${VTK_DIR}/bin ${VTK_DIR}; do
        if test "x" = "x$VTK_LIB_DIR"; then
          LIBS="$old_LIBS -L$dir"
          unset ac_cv_lib_vtkCommon
  	unset ac_cv_lib_vtkCommon_main
          AC_HAVE_LIBRARY([vtkCommon],[VTK_LIB_DIR=$dir])
        fi
      done
      if test "x" = "x$VTK_LIB_DIR"; then
        VTK_MISSING=yes
      fi
      LIBS="$old_LIBS"
      ;;
  esac
  
dnl Define a function that checks for the existence of a library
dnl and if it exists, prepends the corresponding link flag to a 
dnl sepecified variable
dnl Arguments:
dnl  - The library name
dnl  - The name of the variable to which to prepend the library link flag
dnl  - List of extra libraries required for link test
dnl  - Optional library path linker argument (e.g. -L/opt/foo)
  AC_DEFUN([FATHOM_CHECK_LIB_LIST],[
    old_LIBS="$LIBS"
    LIBS="$4 ${$2} $old_LIBS"
    AC_HAVE_LIBRARY([$1],[$2="-l$1 ${$2}"],[],[$3])
    LIBS="$old_LIBS"
  ])
  
   # These are ordered such that a given library depends only on those libraries
   # that preceed it in the list.
  FATHOM_CHECK_LIB_LIST([pthread],[vtk_sys_libs])
  FATHOM_CHECK_LIB_LIST([dl],     [vtk_sys_libs])
  FATHOM_CHECK_LIB_LIST([m],      [vtk_sys_libs])
  FATHOM_CHECK_LIB_LIST([gcc],    [vtk_sys_libs])
  x_libs_arg="$X_EXTRA_LIBS -lX11 $X_PRE_LIBS"
  FATHOM_CHECK_LIB_LIST([Xext],[vtk_x_libs],[$x_libs_arg],[-L$x_libraries]) 
  FATHOM_CHECK_LIB_LIST([ICE], [vtk_x_libs],[$x_libs_arg],[-L$x_libraries]) 
  FATHOM_CHECK_LIB_LIST([SM],  [vtk_x_libs],[$x_libs_arg],[-L$x_libraries]) 
  FATHOM_CHECK_LIB_LIST([Xt],  [vtk_x_libs],[$x_libs_arg],[-L$x_libraries]) 
  FATHOM_CHECK_LIB_LIST([GL],  [vtk_x_libs],[$x_libs_arg],[-L$x_libraries]) 
  FATHOM_CHECK_LIB_LIST([GLU], [vtk_x_libs],[$x_libs_arg],[-L$x_libraries]) 
  vtk_libs_arg="-L$x_libraries $vtk_x_libs $X_EXTRA_LIBS -lX11 $X_PRE_LIBS $vtk_sys_libs"
  FATHOM_CHECK_LIB_LIST([vtksys],        [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkCommon],     [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkDICOMParser],[vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkfreetype],   [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkexpat],      [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkzlib],       [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkpng],        [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkjpeg],       [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtktiff],       [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkftgl],       [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkFiltering],  [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkImaging],    [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([verdict],       [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkGraphics],   [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkNetCDF],     [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkexoIIc],     [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkmetaio],     [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkIO],         [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkRendering],  [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkParallel],   [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([vtkHybrid],     [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  FATHOM_CHECK_LIB_LIST([QVTK],          [vtk_libs],[$vtk_libs_arg],[-L$VTK_LIB_DIR])
  
   # Make sure we found at least vtkCommon and vtkRendering
  if ! echo "x $vtk_libs" | grep lvtkCommon >/dev/null; then
    AC_MSG_WARN("VTK library 'vtkCommon' not found.")
    VTK_MISSING=yes
  fi
  if ! echo "x $vtk_libs" | grep lvtkRendering >/dev/null; then
    AC_MSG_WARN("VTK library 'vtkRendering' not found.")
    VTK_MISSING=yes
  fi
  
  AC_LANG_RESTORE
  
  
################################################################################
#                                   QT
################################################################################
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_ARG_VAR([QTDIR],[QT installation directory])
  if test x"$QTDIR" = x; then
    search_path="$PATH"
    qlibdir=
    incdir=
    for dir in /usr/include/qt3 /usr/local/include/qt3 /include/qt3; do
      if test x = x${incdir}; then
        if test -d $dir; then
          incdir=-I$dir
        fi
      fi
    done
  else
    search_path="${QTDIR}/bin"
    qlibdir="-L${QTDIR}/lib"
    incdir="-I${QTDIR}/include"
  fi
  
  AC_PATH_PROG([MOC],[moc],[],[$search_path])
  if test x = x"$MOC";then
    AC_MSG_WARN([Cannot find moc (Qt meta object compiler)])
    QT_MISSING=yes
  fi
  AC_SUBST(MOC)
  AC_PATH_PROG([UIC],[uic],[],[$search_path])
  if test x = x"$UIC";then
    AC_MSG_WARN([Cannot find uic (Qt user interface compiler)])
    QT_MISSING=yes
  fi
  AC_SUBST(UIC)
  
  old_CXXFLAGS="$CXXFLAGS"
  old_CPPFLAGS="$CPPFLAGS"
  CXXFLAGS="$incdir"
  CPPFLAGS="$incdir"
  AC_CHECK_HEADERS([qlineedit.h qevent.h qapplication.h \
                    qwidget.h qtimer.h qpixmap.h qwidgetplugin.h \
                    qobject.h qmetaobject.h],[],
                   [QT_MISSING=yes])
  CXXFLAGS="$old_CXXFLAGS"
  CPPFLAGS="$old_CPPFLAGS"
  QT_INCLUDES="$incdir"
  
  old_LIBS="$LIBS"
  LIBS="$LIBS $qlibdir"
  # Cannot use AC_HAVE_LIBRARY for qt-mt because of
  # autoconf bug with libraries that have a '-' in their name.
  #AC_HAVE_LIBRARY([qt-mt],[QT_LIBS="$qlibdir -lqt-mt"],[
  AC_MSG_CHECKING([for main in -lqt-mt])
  QT_LIBS="$qlibdir -lqt-mt"
  AC_TRY_LINK_FUNC(main,[AC_MSG_RESULT([yes])],[
    AC_MSG_RESULT([no])
    LIBS="$old_LIBS $qlibdir"
    AC_HAVE_LIBRARY([qt],[QT_LIBS="$qlibdir -lqt"],[
      QT_MISSING=yes
    ])
  ])
  LIBS="$old_LIBS"
  
################################################################################
#                                 X/VTK/QT
################################################################################
  dnl Merge lists of libraries
  dnl Arguments
  dnl  - name of output variable
  dnl  - first list
  dnl  - second list
  AC_DEFUN([FATHOM_MERGE_LISTS],[
    $1="$2"
    # Reverse list two
    for item3 in $3; do 
      found=no
      for item2 in $2; do
        if test "x$item2" = "x$item3"; then found=yes; fi
      done
      if test "$found" = "no"; then $1="${$1} $item3"; fi
    done
    $1="${$1} $2"
  ])
  FATHOM_MERGE_LISTS([XLIBS_tmp],[$X_EXTRA_LIBS -lX11 $X_PRE_LIBS],[$vtk_x_libs])
  FATHOM_MERGE_LISTS([XLIBS],[$vtk_sys_libs],[$XLIBS_tmp])
  QVDUAL_LIBS="$LIBS $GRAPH_LIBS -L$VTK_LIB_DIR $vtk_libs $QT_LIBS $XLIBS"
  QVDUAL_INCLUDES="$GRAPH_INC_DIR $VTK_INCLUDES $QT_INCLUDES"
  AC_SUBST(VTK_LIB_DIR)
  AC_SUBST(QVDUAL_LIBS)
  AC_SUBST(QVDUAL_INCLUDES)
fi

################################################################################
#                           CGM
################################################################################
AC_ARG_WITH(cgm, 
[AC_HELP_STRING([--with-cgm=DIR], [Specify directory containing CGM])],
[CGM_DIR=$withval
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-cgm=\"${withval}\""
 ], [CGM_DIR=])

case "x$CGM_DIR" in
  xno|x)
    CGM_MISSING=yes
    ;;
  xyes)
    AC_MSG_ERROR([--with-cgm requires an argument])
    ;;
  *)
    if ! test -d "$CGM_DIR"; then
      AC_MSG_WARN([$CGM_DIR : not a directory.])
    fi
    CGM_MISSING=no
    CGM_MAKE=
    AC_CHECK_FILE([${CGM_DIR}/lib/cgm.make],
                  [CGM_MAKE="${CGM_DIR}/lib/cgm.make";
                   DEFINES="$DEFINES -DCGM"],
                  [AC_CHECK_FILE([${CGM_DIR}/cgm.make],
                                 [CGM_MAKE="${CGM_DIR}/cgm.make";
                   DEFINES="$DEFINES -DCGM"],
                                 [AC_MSG_WARN([$CGM_DIR : not a configured CGM]); 
                                  CGM_MISSING=yes])
                  ])
    ;;
esac

if test "x$CGM_MISSING" = "xno"; then
  CGM_CONFIG_OPTIONS="include $CGM_MAKE"

  AC_MSG_CHECKING([Checking for GeometryQueryTool::ray_fire(RefFace*,...)])
  FATHOM_MAKE_INC_VAR([$CGM_MAKE],[CGM_DEFINES], [CGM_CPPFLAGS="$make_val"])   
  FATHOM_MAKE_INC_VAR([$CGM_MAKE],[CGM_INCLUDES],[CGM_CPPFLAGS="$CGM_CPPFLAGS $make_val"])   
  FATHOM_MAKE_INC_VAR([$CGM_MAKE],[CGM_LTFLAGS], [CGM_LTFLAGS="$make_val"])   
  FATHOM_MAKE_INC_VAR([$CGM_MAKE],[CGM_LDFLAGS], [CGM_LDFLAGS="$make_val"])   
  FATHOM_MAKE_INC_VAR([$CGM_MAKE],[CGM_LIBS],    [CGM_LIBS="$make_val"])   
  AC_LANG_PUSH([C++])
  old_CPPFLAGS="$CPPFLAGS"
  old_LDFLAGS="$LDFLAGS"
  old_LIBS="$LIBS"
  CPPFLAGS="$CGM_CPPFLAGS $CPPFLAGS"
  LDFLAGS="$CGM_LTFLAGS $CGM_LDFLAGS $CPPFLAGS"
  LIBS="$CGM_LIBS $LIBS"
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include "GeometryQueryTool.hpp"
                      #include "CubitVector.hpp"],
                    [CubitVector a,b;
                     DLIList<double> c;
                     GeometryQueryTool::instance()->fire_ray((RefFace*)0, a, b, c );
                    ])],
    [MOAB_CGM_DEFINES="-DHAVE_CGM_FIRE_RAY"; AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)])
  CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$old_LDFLAGS"
  LIBS="$old_LIBS"
  AC_LANG_POP([C++])
  EXPORT_LDFLAGS="$EXPORT_LDFLAGS $CGM_LDFLAGS"
  EXPORT_LTFLAGS="$EXPORT_LTFLAGS $CGM_LTFLAGS"
fi

AC_SUBST(MOAB_CGM_DEFINES)
AC_SUBST(CGM_CPPFLAGS)
AC_SUBST(CGM_LDFLAGS)
AC_SUBST(CGM_LTFLAGS)
AC_SUBST(CGM_LIBS)

AM_CONDITIONAL( HAVE_CGM, [test "x$CGM_MISSING" = "xno"] )

################################################################################
#                           iMesh Part 3 of 3 (Babel)
################################################################################
if test "x$ENABLE_imesh" = "xyes"; then
  AC_ARG_WITH(babel, 
  [AC_HELP_STRING([--with-babel=DIR], [Specify directory containing babel])],
  [BABEL_DIR=$withval
   DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-babel=\"${withval}\""
   ], [BABEL_DIR=no])
  case "x$BABEL_DIR" in
    xno)
      BABEL_DIR=
      USE_BABEL=no
      ;;
    xyes|x)
      BABEL_DIR=/usr/local
      AC_CHECK_FILE([${BABEL_DIR}/bin/babel], [USE_BABEL=yes], [USE_BABEL=no])
      ;;
    *)
      AC_CHECK_FILE([${BABEL_DIR}/bin/babel], [USE_BABEL=yes], [USE_BABEL=no])
      ;;
  esac
  AC_SUBST(BABEL_DIR)
  AC_SUBST(USE_BABEL)
else
  USE_BABEL=no
fi
AM_CONDITIONAL(USE_BABEL, [test "xno" != "x$USE_BABEL"])

################################################################################
#                           Output Files
################################################################################
AC_SUBST([INCLUDES])
AC_SUBST([DEFINES])
AC_SUBST([AM_CPPFLAGS])
AM_LDFLAGS="$AM_LDFLAGS $EXPORT_LTFLAGS $EXPORT_LDFLAGS"
AC_SUBST([AM_LDFLAGS])
AC_SUBST([EXPORT_LTFLAGS])
AC_SUBST([EXPORT_LDFLAGS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([DISTCHECK_CONFIGURE_FLAGS])

AC_ARG_VAR([FC], [FORTRAN compiler command])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_HEADERS([MBVersion.h MBEntityHandle.h parallel/MBmpi_config.h])
AC_CONFIG_HEADERS([iBase_FCDefs.h])
AC_CONFIG_FILES([Makefile 
                 moab.make 
                 testdir.h
                 parallel/Makefile
                 mhdf/Makefile
                 test/Makefile
                 test/h5file/Makefile
                 test/dual/Makefile
                 test/obb/Makefile
                 test/perf/Makefile
                 tools/Makefile
                 tools/iMesh/Makefile
                 tools/iMesh/iMesh-Defs.inc
                 tools/iMesh/SIDL/iMesh-SIDL-Defs.inc
                 tools/iMesh/SIDL/Makefile
                 tools/iMesh/SIDL/common.make
		 tools/iGeom/Makefile
                 tools/iGeom/iGeom-Defs.inc
                 tools/mbperf/Makefile
                 tools/mbchaco/Makefile
                 tools/mbcoupler/Makefile
                 tools/mbzoltan/Config.moab
                 tools/mbzoltan/Makefile
                 tools/mcnpmit/Makefile
                 tools/qvdual/Makefile
                 tools/dagmc/Makefile
                 doc/Makefile
                 doxygen/Makefile
                 ])
AC_CONFIG_COMMANDS([MBCN_FCDefs.h],
  [sed -e "s/FC_FUNC/MBCN_FC_FUNC/" iBase_FCDefs.h >MBCN_FCDefs.h])

AC_OUTPUT

AC_MSG_RESULT([C:   $CC $CFLAGS $CPPFLAGS])
AC_MSG_RESULT([C++: $CXX $CXXFLAGS $CPPFLAGS])
if test "xno" != "x$CONFIGURE_FORTRAN"; then
  AC_MSG_RESULT([Fortran: $FC $FCFLAGS])
fi

if test "x$GRAPHVIZ_MISSING" = "xyes"; then
  AC_MSG_WARN([Could not find GraphViz library.  Build will FAIL!!! Try \"--with-graphviz=DIR\"])
fi
if test "x$X_MISSING" = "xyes"; then
  AC_MSG_WARN([Cound not find X client.  Build will FAIL!!!])
else
  if test "x$QT_MISSING" = "xyes"; then
    AC_MSG_WARN([Cound not find QT.  Build will FAIL!!! Try setting \"QTDIR\"])
  fi
  if test "x$VTK_MISSING" = "xyes"; then
    AC_MSG_WARN([Count not find VTK.  Build will FAIL!!! Try \"--with-VTK=DIR\"])
  fi
fi

if test "x$HAVE_HDF5" = "xno"; then
  AC_MSG_WARN([
*************************************************************************
*            MOAB has been configured w/out the HDF5 library.
*              Support for native file format is disabled!
*************************************************************************])
fi

if test "x$HAVE_NETCDF" = "xno"; then
  AC_MSG_WARN([
*************************************************************************
*           MOAB has been configured w/out the NetCDF library.
*     Support for ExodusII/Genesis and other file formats is disabled!
*************************************************************************])
fi

