AC_DEFUN([FATHOM_HDF5_LIBS_HELPER],[
if (test $HAVE_LIB_HDF5 == no); then
   unset "ac_cv_lib_${HDF5_LIBNAME}_H5Fopen"
   unset "ac_cv_lib_${HDF5_LIBNAME}___H5Fopen"
   AC_CHECK_LIB( [${HDF5_LIBNAME}], [H5Fopen], [HAVE_LIB_HDF5=yes; HDF5_LIBS="$HDF5_LIBS $1"], [], [$1] )
fi
])

#######################################################################################
# Helper function for FATHOM_CHECK_HDF5 and FATHOM_CHECK_NETCDF
# If HAVE_LIB_HDF5 == yes, then does nothing.
# Otherwise sets HAVE_LIB_HDF5 to yes or no depending on whether or not
# HDF5 libraries are detected and sets HDF5_LIBS
# Respects caller's LDFLAGS, CPPFLAGS, and LIBS.  Caller should set appropriately.
# Respections optional HDF5_LIBNAME uses to specify alternate library name.  If
# not specified, will be set to -lhdf5
#######################################################################################
AC_DEFUN([FATHOM_DETECT_HDF5_LIBS],[

 # if we've already done this check, then don't do it again
if test "xyes" != "x$HAVE_LIB_HDF5"; then
  test "x" != "x$HDF5_LIBNAME" || HDF5_LIBNAME=hdf5
  
  HAVE_LIB_HDF5=no
  FATHOM_HDF5_LIBS_HELPER([$LIBS])
fi
])

#######################################################################################
# Check for HDF5 library and related stuff
# Sets HAVE_HDF5 to 'yes' or 'no'
# If HAVE_HDF5 == yes, then sets:
#  HDF5_CPPFLAGS
#  HDF5_LDFLAGS
#  HDF5_LIBS
#######################################################################################
AC_DEFUN([FATHOM_CHECK_HDF5],[

  # Check whether user wants to autodownload HDF5
  # Call package Download/Configure/Installation procedures for MOAB, if requested by user
  PPREFIX=HDF5
  
  # VERSIONS: 1.8.4 (versions greater need CMake)
  m4_pushdef([HDF5_DOWNLOAD_VERSION],[1.8.12])dnl
 
  # Invoke the download-hdf5 command
  m4_case( HDF5_DOWNLOAD_VERSION, [1.8.15], [ CONFIGURE_EXTERNAL_PACKAGE([HDF5], [http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.15-patch1/src/hdf5-1.8.15-patch1.tar.gz], "hdf5-1.8.15.tar.gz") ],
                                  [1.8.14], [ CONFIGURE_EXTERNAL_PACKAGE([HDF5], [http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz], "hdf5-1.8.14.tar.gz") ],
                                  [1.8.12], [ CONFIGURE_EXTERNAL_PACKAGE([HDF5], [http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.12/src/hdf5-1.8.12.tar.gz], "hdf5-1.8.12.tar.gz") ],
                                  [1.8.10], [ CONFIGURE_EXTERNAL_PACKAGE([HDF5], [http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.10/src/hdf5-1.8.10.tar.gz], "hdf5-1.8.10.tar.gz") ],
                                  [1.8.4],  [ CONFIGURE_EXTERNAL_PACKAGE([HDF5], [http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.4/src/hdf5-1.8.4.tar.gz], "hdf5-1.8.4.tar.gz") ], 
                                  [ CONFIGURE_EXTERNAL_PACKAGE([HDF5], [http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.4/src/hdf5-1.8.14.tar.gz], "hdf5-1.8.14.tar.gz") ] )
  if (test "x$downloadhdf5" == "xyes") ; then
    # download the latest HDF5 sources, configure and install
    HDF5_SRCDIR="$hdf5_src_dir"
    AC_SUBST(HDF5_SRCDIR)
    # The default HDF5 installation is under libraries
    HDF5_DIR="$hdf5_install_dir"
    enablehdf5=yes
  fi  # if (test "$downloadhdf5" != no) ; then 

  # CLI option for linking zlib
AC_ARG_WITH(zlib,
  [AS_HELP_STRING([--with-zlib=DIR],[HDF5 requires zlib, and zlib can be found at...])],
  [WITH_ZLIB=$withval
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-zlib=\"${withval}\""
  ],[WITH_ZLIB=])
case "x$WITH_ZLIB" in
  xyes|xno|x)
    ;;
  *)
    if ! test -d  ${WITH_ZLIB}/lib; then
      AC_MSG_ERROR([Not a directory: ${WITH_ZLIB}/lib])
    fi
    HDF5_LDFLAGS="$HDF5_LDFLAGS -L${WITH_ZLIB}/lib"
    ;;
esac
HAVE_ZLIB=no
if test "x$WITH_ZLIB" != "xno"; then
  old_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
  AC_CHECK_LIB([z],[deflate],[HAVE_ZLIB=yes; HDF5_LIBS="$HDF5_LIBS -lz"],
    [if test "x$WITH_ZLIB" != "x"; then AC_MSG_ERROR([Could not find zlib]); fi])
  LDFLAGS="$old_LDFLAGS"
fi

  # CLI option for linking szip
AC_ARG_WITH(szip,
  [AS_HELP_STRING([--with-szip=DIR],[HDF5 requires szip, and szip an be found at...])],
  [WITH_SZIP=$withval
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-szip=\"${withval}\""
  ],[WITH_SZIP=])
case "x$WITH_SZIP" in
  xyes|xno|x)
    ;;
  *)
    if ! test -d  ${WITH_SZIP}/lib; then
      AC_MSG_ERROR([Not a directory: ${WITH_SZIP}/lib])
    fi
    HDF5_LDFLAGS="$HDF5_LDFLAGS -L${WITH_SZIP}/lib"
    ;;
esac
HAVE_SZIP=no
if test "x$WITH_SZIP" != "xno"; then
  old_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
  AC_CHECK_LIB([sz],[SZ_Decompress],[HAVE_SZIP=yes; HDF5_LIBS="$HDF5_LIBS -lsz"],
    [if test "x$WITH_SZIP" != "x"; then AC_MSG_ERROR([Could not find libsz]); fi])
  LDFLAGS="$old_LDFLAGS"
fi

  # CLI option for extra HDF5 link flags
AC_ARG_WITH([hdf5-ldflags],[AS_HELP_STRING([--with-hdf5-ldflags=...],
 [Extra LDFLAGS required for HDF5 library (e.g. parallel IO lib)])],
 [HDF5_LDFLAGS_WITHVAL="$withval"; 
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-hdf5-ldflags=\"${withval}\""
],[HDF5_LDFLAGS_WITHVAL=])
case "x$HDF5_LDFLAGS_WITHVAL" in
  xno)
    AC_MSG_ERROR("Invalid argument: --without-hdf5-ldflags")
    ;;
  xyes)
    AC_MSG_ERROR("Nonsensical argument:  --with-hdf5-ldflags without any specified flags")
    ;;
  *)
    HDF5_LDFLAGS="$HDF5_LDFLAGS $HDF5_LDFLAGS_WITHVAL"
    ;;
esac


  # CLI option for HDF5
AC_MSG_CHECKING([if HDF5 support is enabled])
AC_ARG_WITH(hdf5, 
[AS_HELP_STRING([--with-hdf5@<:@=DIR@:>@], [Specify HDF5 library to use for native file format])
AS_HELP_STRING([--without-hdf5], [Disable support for native HDF5 file format])],
[HDF5_ARG=$withval
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-hdf5=\"${withval}\""
], [HDF5_ARG=$HDF5_DIR])
if test "xno" = "x$HDF5_ARG"; then
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([yes])
fi

HAVE_HDF5=no
if test "xno" != "x$HDF5_ARG"; then
  HAVE_HDF5=yes

    # if a path is specified, update LIBS and INCLUDES accordingly
  if test "xyes" != "x$HDF5_ARG" && test "x" != "x$HDF5_ARG"; then
    if test -d "${HDF5_ARG}/dll"; then
      HDF5_LDFLAGS="$HDF5_LDFLAGS -L${HDF5_ARG}/dll"
      HDF5_LIBNAME=hdf5dll
    elif test -d "${HDF5_ARG}/lib"; then
      HDF5_LDFLAGS="$HDF5_LDFLAGS -L${HDF5_ARG}/lib"
    elif test -d "${HDF5_ARG}"; then
      HDF5_LDFLAGS="$HDF5_LDFLAGS -L${HDF5_ARG}"
    else
      AC_MSG_ERROR("$HDF5_ARG is not a directory.")
    fi
    if test -d "${HDF5_ARG}/include"; then
      HDF5_CPPFLAGS="$HDF5_CPPFLAGS -I${HDF5_ARG}/include"
      if test "x$GXX" = "xyes" && test "x$GCC" = "xyes"; then
        HDF5_CPPFLAGS="$HDF5_CPPFLAGS -isystem ${HDF5_ARG}/include"
      fi
    else
      HDF5_CPPFLAGS="$HDF5_CPPFLAGS -I${HDF5_ARG}"
    fi
  fi
 
  # Check for IBM parallel IO library
  if test "x$enablempi" != "xno"; then
    AC_CHECK_LIB([gpfs],[gpfs_stat],[LIBS="-lgpfs $LIBS"])
  fi

  # Add flag to defines
  old_CPPFLAGS="$CPPFLAGS"
  old_LDFLAGS="$LDFLAGS"
  old_LIBS="$LIBS"
  CPPFLAGS="$HDF5_CPPFLAGS $CPPFLAGS"
  LDFLAGS="$HDF5_LDFLAGS $LDFLAGS"
  LIBS="$HDF5_LIBS $LIBS"
  
    # check for libraries and headers
  AC_CHECK_HEADERS( [hdf5.h], [], [HAVE_HDF5=no] )

  HAVE_LIB_HDF5=no
  FATHOM_DETECT_HDF5_LIBS

  if test "x$HAVE_LIB_HDF5" = "xno"; then
    HAVE_HDF5=no
  fi
  
  if test "x$HAVE_HDF5" = "xno"; then
    if test "x" = "x$HDF5_ARG"; then
      AC_MSG_WARN([HDF5 library not found or not usable.])
    else
      AC_MSG_ERROR([HDF5 library not found or not usable.])
    fi
    HDF5_CPPFLAGS=
    HDF5_LDFLAGS=
    HDF5_LIBS=
  else
    HDF5_LIBS="-l$HDF5_LIBNAME $HDF5_LIBS"
  fi
  
  CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$old_LDFLAGS"
  LIBS="$old_LIBS"
fi

]) # FATHOM_CHECK_HDF5


dnl ---------------------------------------------------------------------------
dnl AUTOMATED SETUP PREPROCESS HDF5
dnl   Figure out what needs to be done to get a valid HDF5 installation.
dnl   Arguments: [PACKAGE, SRC_DIR, INSTALL_DIR, NEED_CONFIGURATION)
dnl ---------------------------------------------------------------------------
AC_DEFUN([AUTOMATED_SETUP_PREPROCESS_HDF5],[
  # uncompress and configure PACKAGE
  hdf5_src_dir="$2"
  hdf5_build_dir="$2/build"
  hdf5_install_dir="$3"
  hdf5_archive_name="$4"
  PPREFIX=HDF5

  if (test ! -d "$hdf5_src_dir" || test ! -f "$hdf5_src_dir/configure" ); then
    AC_MSG_ERROR([Invalid source configuration for HDF5. Source directory $hdf5_src_dir is invalid])
  fi

  # determine what steps we need to take to install hdf5
  hdf5_configured=false
  hdf5_made=false
  hdf5_installed=false
  if (test ! -d "$hdf5_build_dir" ); then
    AS_MKDIR_P( $hdf5_build_dir )
  else
    if (test -f "$hdf5_build_dir/src/H5config.h" ); then
      hdf5_configured=true
      if (test -f "$hdf5_build_dir/src/.libs/libhdf5.a" || test -f "$hdf5_build_dir/src/.libs/libhdf5.so" || test -f "$hdf5_build_dir/src/.libs/libhdf5.dylib" ); then
        hdf5_made=true
        if (test -f "$hdf5_install_dir/lib/libhdf5.settings"); then
          hdf5_installed=true
        fi
      fi
    fi
  fi
  # send the information back
  AS_IF([ ! $hdf5_configured || $need_configuration ], [need_configuration=true], [need_configuration=false])
  AS_IF([ ! $hdf5_made || $need_configuration ], [need_build=true], [need_build=false])
  AS_IF([ ! $hdf5_installed || $need_configuration ], [need_installation=true], [need_installation=false])
])


dnl ---------------------------------------------------------------------------
dnl AUTOMATED SETUP POSTPROCESS HDF5
dnl   Dummy macro to fit standard call pattern.  Tells MOAB we have HDF5.
dnl   Arguments: [PACKAGE, SRC_DIR, INSTALL_DIR, NEED_CONFIGURATION)
dnl ---------------------------------------------------------------------------
AC_DEFUN([AUTOMATED_SETUP_POSTPROCESS_HDF5],[
  # we have already checked configure/build/install logs for errors before getting here..
  enablehdf5=yes
])


dnl ---------------------------------------------------------------------------
dnl AUTOMATED CONFIGURE HDF5
dnl   Runs configure for HDF5 and looks for header files.
dnl   Arguments: [NEED_CONFIGURATION)
dnl ---------------------------------------------------------------------------
AC_DEFUN([AUTOMATED_CONFIGURE_HDF5],[
if [ $1 ]; then
  # configure HDF5
  if [ $need_configuration ]; then
    # configure PACKAGE with a minimal build: MPI
    export CC=$CC CXX=$CXX FC=$FC CFLAGS="$CFLAGS -fPIC -DPIC" CXXFLAGS="$CXXFLAGS -fPIC -DPIC" FCFLAGS="$FCFLAGS -fPIC" LDFLAGS=$LDFLAGS
    configure_command="$hdf5_src_dir/configure --prefix=$hdf5_install_dir --libdir=$hdf5_install_dir/lib --with-pic=1"
    # configure_command="$configure_command --enable-cxx --enable-unsupported"
    if (test "$enabledebug" != "no"); then
      configure_command="$configure_command --enable-debug=all"
    fi
    if (test "$enablefortran" != "no"); then
      configure_command="$configure_command --enable-fortran"
    fi
    if (test "$enableoptimize" != "no"); then
      configure_command="$configure_command --enable-production=yes"
    fi
    if (test "$enablempi" != "no"); then
      configure_command="$configure_command --enable-parallel"
    fi
    
    hdf5_configlog=`echo "Using configure command :==> cd $hdf5_build_dir ; $configure_command > $hdf5_src_dir/../config_hdf5.log ; cd \"\$OLDPWD\"" > "$hdf5_src_dir/../config_hdf5.log"`
    ##echo "Trying to use configure command:==> cd $hdf5_build_dir ; $configure_command > $hdf5_src_dir/../config_hdf5.log ; cd \"\$OLDPWD\""
    PREFIX_PRINT(Configuring with default options  {debug=$enabledebug production=$enableoptimize shared=$enable_shared parallel=$enablempi} )
    hdf5_configlog="`cd $hdf5_build_dir ; $configure_command >> $hdf5_src_dir/../config_hdf5.log 2>&1 ; cd \"\$OLDPWD\"`"
  fi

  # check if configuration - current or previous was successful
  if (test ! -f "$hdf5_build_dir/src/H5config.h" ); then
    AC_MSG_ERROR([HDF5 configuration was unsuccessful. Please refer to $hdf5_build_dir/config.log and $hdf5_src_dir/../config_hdf5.log for further details.])
  fi
  hdf5_configured=true
fi
])


dnl ---------------------------------------------------------------------------
dnl AUTOMATED BUILD HDF5
dnl   Calls make on HDF5 and looks for libraries.
dnl   Arguments: [NEED_BUILD)
dnl ---------------------------------------------------------------------------
AC_DEFUN([AUTOMATED_BUILD_HDF5],
[
  # if we need to build then call make all
  if [ $1 ]; then
    if [ $recompile_and_install || $need_build ]; then
      PREFIX_PRINT(Building the sources in parallel )
      hdf5_makelog="`cd $hdf5_build_dir ; make all -j4 > $hdf5_src_dir/../make_hdf5.log 2>&1 ; cd \"\$OLDPWD\"`"
    fi
  fi
  # check if it worked
  if (test -f "$hdf5_build_dir/src/.libs/libhdf5.a" || test -f "$hdf5_build_dir/src/.libs/libhdf5.so" || test -f "$hdf5_build_dir/src/.libs/libhdf5.dylib") ; then
    hdf5_made=true
  else
    AC_MSG_ERROR([HDF5 build was unsuccessful. Please refer to $hdf5_src_dir/../make_hdf5.log for further details.])
  fi
])


dnl ---------------------------------------------------------------------------
dnl AUTOMATED INSTALL HDF5
dnl   Calls make install on HDF5 and checks for libhdf5.settings
dnl   Arguments: [NEED_INSTALLATION)
dnl ---------------------------------------------------------------------------
AC_DEFUN([AUTOMATED_INSTALL_HDF5],
[
  # if we need to install then call make install
  if [ $1 ]; then
    if [ $recompile_and_install ]; then
      if [ $hdf5_installed ]; then
        hdf5_installlog="`cd $hdf5_build_dir ; make uninstall > $hdf5_src_dir/../uninstall_hdf5.log 2>&1 ; cd \"\$OLDPWD\"`"
      fi
      PREFIX_PRINT(Installing the headers and libraries in to directory ($hdf5_install_dir) )
      hdf5_installlog="`cd $hdf5_build_dir ; make install > $hdf5_src_dir/../install_hdf5.log 2>&1 ; cd \"\$OLDPWD\"`"
    fi
  fi
  # check if it worked
  if (test -f "$hdf5_install_dir/lib/libhdf5.settings"); then
    hdf5_installed=true
  else
    AC_MSG_ERROR([HDF5 installation was unsuccessful. Please refer to $hdf5_src_dir/../install_hdf5.log for further details.])
  fi
])
