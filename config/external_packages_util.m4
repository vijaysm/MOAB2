dnl -------------------------------------------------------------
dnl  external_packages_util.m4
dnl
dnl  This file contains general instructions for downloading,
dnl  unpacking, configuring, building, and installing 
dnl  dependent packages.  It also contains definitions of 
dnl  configure line download options.
dnl -------------------------------------------------------------

dnl Need the following
dnl  MOAB_ARCH

################
# DO NOT TOUCH #
################

m4_define([_m4_divert(HELP_BEGIN)],     100)
m4_define([_m4_divert(HELP_CANON)],     101)
m4_define([_m4_divert(HELP_ENABLE)],    102)
m4_define([_m4_divert(HELP_WITH)],      103)
m4_define([_m4_divert(HELP_DOWNLOAD)],  104)
m4_define([_m4_divert(HELP_VAR)],       105)
m4_define([_m4_divert(HELP_VAR_END)],   106)
m4_define([_m4_divert(HELP_END)],       107)

# --------------------------------------------------------------------
# AC_ARG_DOWNLOAD(PACKAGE, HELP-STRING, ACTION-IF-TRUE, [ACTION-IF-FALSE])
# --------------------------------------------------------------------
AC_DEFUN([AC_ARG_DOWNLOAD],
[AC_PROVIDE_IFELSE([AC_PRESERVE_HELP_ORDER],
[],
[m4_divert_once([HELP_DOWNLOAD], [[
Optional Downloads:
  --download-PACKAGE[=ARG]  download/configure PACKAGE with default options [ARG=yes,no,url]]])])dnl
m4_divert_once([HELP_DOWNLOAD], [$2])dnl
_AC_ENABLE_IF([download], [$1], [$3], [$4])dnl
])# AC_ARG_DOWNLOAD

AU_DEFUN([AC_DOWNLOAD],
[AC_ARG_DOWNLOAD([$1], [  --download-$1], [$2], [$3])])

# _AC_INIT_PARSE_ENABLE(OPTION-NAME)                                                                                                                                                                                             
# ----------------------------------
# A trivial front-end for _AC_INIT_PARSE_ENABLE2.
#
m4_define([_AC_INIT_PARSE_ENABLE_ORIG],
[m4_bmatch([$1], [^with],
     [_AC_INIT_PARSE_ENABLE2([$1], [with])],
     [m4_bmatch([$1], [^able],
        [_AC_INIT_PARSE_ENABLE2([$1], [enable])],
        [_AC_INIT_PARSE_ENABLE2([download], [download])])] )])

m4_define([_AC_INIT_PARSE_ENABLE],
[m4_bmatch([$1], [^with],
     [_AC_INIT_PARSE_ENABLE2([$1], [with])
      m4_ifdef([_AC_DEFINED_DOWNLOAD_MACRO], [], 
        [ _AC_INIT_PARSE_ENABLE2([download], [download]) 
          m4_define([_AC_DEFINED_DOWNLOAD_MACRO], [yes] )] )],
     [_AC_INIT_PARSE_ENABLE2([$1], [enable]) ]) ])


# _AC_INIT_PARSE_ENABLE2(OPTION-NAME, POSITIVE-NAME)
# --------------------------------------------------
# Handle an `--enable', a `--with' or a `--download` option.
#
# OPTION-NAME is `enable', `disable', `with', `without' or `download'.
# POSITIVE-NAME is the corresponding positive variant, i.e. `enable' or `with'.
#
# Positive variant of the option is recognized by the condition
# OPTION-NAME == POSITIVE-NAME .
#
m4_define([_AC_INIT_PARSE_ENABLE2],
[-$1-* | --$1-*)
    ac_useropt=`expr "x$ac_option" : 'x-*$1-\(m4_if([$1], [$2], [[[^=]]], [.])*\)'`
    # Reject names that are not valid shell variable names.
    expr "x$ac_useropt" : "[.*[^-+._$as_cr_alnum]]" >/dev/null &&
      AC_MSG_ERROR(
  [invalid ]m4_if([$2], [[with | download]], [package], [feature])[ name: $ac_useropt])                                                                                                                                                       
    ac_useropt_orig=$ac_useropt
    ac_useropt=`AS_ECHO(["$ac_useropt"]) | sed 's/[[-+.]]/_/g'`
    case $ac_user_opts in
      *"
"$2_$ac_useropt"
"*) ;;
      *) ac_unrecognized_opts="$ac_unrecognized_opts$ac_unrecognized_sep--$1-$ac_useropt_orig"
   ac_unrecognized_sep=', ';;
    esac
    eval $2_$ac_useropt=m4_if([$1], [$2], [\$ac_optarg], [no]) ;;dnl
])

####################
# RESTRICTION ENDS #
####################

AC_DEFUN([INITIALIZE_EXTERNAL_PACKAGES],
[
  # Check for command line utility/archive/download programs
  # that are essential for configure to work correctly.
  AC_CHECK_PROG(HAVE_READLINK, readlink, yes, no)
  AC_CHECK_PROG(HAVE_DIRNAME, dirname, yes, no)
  AC_CHECK_PROG(HAVE_BASENAME, basename, yes, no)
  AC_CHECK_PROG(HAVE_RSYNC, rsync, yes, no)
  # network download programs
  AC_CHECK_PROG(HAVE_WGET, wget, yes, no)
  AC_CHECK_PROG(HAVE_SCP, scp, yes, no)
  AC_CHECK_PROG(HAVE_CURL, curl, yes, no)
  # archive file inflation/deflation programs
  AC_CHECK_PROG(HAVE_TAR, tar, yes, no)
  AC_CHECK_PROG(HAVE_UNZIP, unzip, yes, no)
  AC_CHECK_PROG(HAVE_BZIP2, bzip2, yes, no)
  # file/directory hash computation programs
  AC_CHECK_PROG(HAVE_MD5SUM, md5sum, yes, no)
  if (test "$HAVE_MD5SUM" != "no"); then
    AC_PATH_PROG(MD5SUM_X, md5sum, "")
    HASHPRGM="$MD5SUM_X"
  else
    AC_CHECK_PROG(HAVE_SHASUM, shasum, yes, no)
    if (test "$HAVE_SHASUM" != "no"); then
      AC_PATH_PROG(SHASUM_X, shasum, "")
      HASHPRGM="$SHASUM_X -a1"
    else
      AC_ERROR([No file/directory hash computation program is available. Report to moab-dev@mcs.anl.gov])
    fi
  fi
  # other essential programs
  AC_PROG_LN_S
  AC_PATH_PROG(MKDIR_P, mkdir, "")
  if (test "x$MKDIR_P" != "x"); then
    MKDIR_P="$MKDIR_P -p"
  else
    AC_ERROR([Make directory command not found ? Seriously ? Report to moab-dev@mcs.anl.gov])
  fi
  AC_PROG_MKDIR_P
  AC_PROG_MAKE_SET

  # Some aliases for colors to pretty-print
  NORMAL=$(tput sgr0)
  GREEN=$(tput setaf 2)
  RED=$(tput setaf 1)

  MOAB_ARCH="$host"
  MOAB_ARCH_DIR="$PWD/sandbox/$MOAB_ARCH"
  MOAB_PACKAGES_DIR="$PWD/sandbox/archives"

  # The default PACKAGE installation is under libraries
  AS_MKDIR_P("$MOAB_ARCH_DIR")
  AS_MKDIR_P("$MOAB_PACKAGES_DIR")
])

# AC_PROG_MKDIR_P
# is a backport of autoconf-2.60's AC_PROG_MKDIR_P.
# Remove this macro when we can assume autoconf >= 2.60.
m4_ifdef([AC_PROG_MKDIR_P], [], [
  AC_DEFUN([AC_PROG_MKDIR_P],
    [AC_REQUIRE([AM_PROG_MKDIR_P])dnl defined by automake
     MKDIR_P='$(mkdir_p)'
     AC_SUBST([MKDIR_P])])])

dnl -------------------------------------------------------------
dnl Fetches an external package into MOAB_ARCH using:
dnl $1 = Package Name
dnl $2 = URL
dnl $3 = Storage location (Archive name)
dnl -------------------------------------------------------------
AC_DEFUN([DOWNLOAD_EXTERNAL_PACKAGE],
[
  PREFIX_PRINT(Downloading sources from URL: $2 )
  cdir=`dirname $3`
  if (test "x$cdir" != "x."); then
    op_dirname="$cdir"
  else
    op_dirname="$MOAB_PACKAGES_DIR"
  fi
  hashtarfile1="0"
  if (test -f "$3"); then
    hashtarfile1="`$HASHPRGM $3 | cut -d ' ' -f1`"
  fi
  filedownloaded=no
  remoteprotocol=yes
  
  # decipher protocol needed to download
  case $2 in
    @*) remoteprotocol=no ;;
    *)  remoteprotocol=yes ;;
  esac
  currdir="$PWD"
  if (test $remoteprotocol != no); then
    if (test "$HAVE_WGET" != "no" ); then
      PREFIX_PRINT([   WGET: $1 package downloading to $3 ])
      if (test -f "$3"); then
        # Check if the file requested exists in the remote directory -- inform user if there is a network error 
        op_checkifexists="`wget --spider -O/dev/null -q $2 && echo yes || echo no`"
        if (test "$op_checkifexists" != "yes"); then
          AC_ERROR([ --  Requested URL does not exist in remote host. Try again later. ($2)  -- ])
        fi
        #op_needdownload="`wget --spider -N -q $2 && echo yes || echo no; cd $currdir`"
        bnamerem="`basename $2`"
        bnameloc="`basename $3`"
        if (test "$bnamerem" != "$bnameloc"); then
          op_downloadlog$1="`wget -q -c -N --progress=bar $2 ; mv $bnamerem $bnameloc`"
        else
          op_downloadlog$1="`wget -q -c -N --progress=bar $2`"
        fi
      else
        # No need to check for time-stamping
        eval "wget -q -S --progress=bar $2 -O $3"
      fi
      filedownloaded=yes
    fi

    if (test "$HAVE_CURL" != "no" && test "$filedownloaded" != "yes"); then
      PREFIX_PRINT([   CURL: $1 package downloading to $3 ])
      op_downloadlog$1="`curl -R -s $2 -z $3 -o $3`"
      filedownloaded=yes
    fi
  else
    if (test "$HAVE_SCP" != "no" && test "$filedownloaded" != "yes"); then
      bnamerem="`echo $2 | cut -c 2-`"
      # op_downloadlog$1="`scp -q $bnamerem $3`"
      PREFIX_PRINT([   SCP: $1 package downloading to $3 ])
      op_downloadlog$1="`scp -q $2 $3`"
      filedownloaded=yes
    fi
  fi
  
  if (test "$filedownloaded" != "yes"); then
    AC_ERROR([ --  The archive URL ($2) specified cannot be handled by wget, curl or scp  -- ])
  fi

  hashtarfile2="`$HASHPRGM $3 | cut -d ' ' -f1`"
  if (test "$hashtarfile1" != "$hashtarfile2"); then
    new_download=true
  else
    new_download=false
  fi

])

dnl -------------------------------------------------------------
dnl Unpacks an external package using:
dnl $1 = Package Name
dnl $2 = Storage location (Archive name)
dnl $3 = Source location (to untar)
dnl -------------------------------------------------------------
AC_DEFUN([DEFLATE_EXTERNAL_PACKAGE],
[
  PREFIX_PRINT([Deflating archive ($2 => $3) ])

  currdir="$PWD"
  need_configuration=false
  op_pkg_subdir=""
  if (test "x`basename $2|grep -E '\.tar'`" != "x" && test "$HAVE_TAR" != "no" ); then
    ##op_pkg_subdir="`tar -tf $2 | head -n 1 | $SED -e 's/\/$//'`"
    op_pkg_subdir="`tar -tf $2 | $SED -e 's@/.*@@' | uniq`"
    PREFIX_PRINT([   Untar file: $2, and $1-SRCDIR=$pkg_srcdir ])
    if (test ! -d "$3/$op_pkg_subdir"); then
      op_unziplog$1="`cd $3 && tar -xf $2 && cd $currdir`"
      need_configuration=true
    fi
    if [ $new_download ]; then
      need_configuration=true
    fi
  elif (test "x`basename $2|grep -E '\.zip'`" != "x" && test "$HAVE_UNZIP" != "no" ); then
    PREFIX_PRINT([   Unzip file: $2, and $1-SRCDIR=$pkg_srcdir <<<])
    if ($new_download -eq true || test ! -d "$pkg_srcdir"); then
      op_unziplog$1="`cd $3 && unzip -q $2 -d $3 && cd $currdir`"
      need_configuration=true
    fi
  elif (test "x`basename $2|grep -E '\.bz'`" != "x" && test "$HAVE_BZIP2" != "no" ); then
    PREFIX_PRINT([   Bunzip file: $2, and $1-SRCDIR=$pkg_srcdir <<<])
    if ( $new_download -eq true || test ! -d "$pkg_srcdir"); then
      op_bziplog$1="`cd $3 && bzip2 -d -q $2 && cd $currdir`"
      need_configuration=true
    fi
  else
    AC_ERROR([ --  Unhandled file format for deflating package $1 -- Filename = $2 -- ])
  fi

  if (test "x$op_pkg_subdir" != "x"); then
    pkg_srcdir="$3/$op_pkg_subdir"
  else
    AC_ERROR([ --  Unhandled file format for getting the source tree name for package $1 -- Filename = $2 -- ])
  fi
])

dnl -------------------------------------------------------------
dnl $1 = Package Name
dnl $2 = Source tree location
dnl $3 = Tarball location
dnl -------------------------------------------------------------
AC_DEFUN([CHECK_SOURCE_RECOMPILATION_HASH],
[
  PREFIX_PRINTN([Checking whether $1 sources need compilation and installation... ])
  # Compute the hash of the source directory - Recompile only if sources have changed
  # ac_cv_sha_moabcpp="`find $moab_src_dir -name '*.cpp' \( -exec $HASHPRGM "$PWD"/{} \; -o -print \) | $HASHPRGM | cut -d ' ' -f1`"
  # defaultshasum="`find $2 -type f -regex '.*\(hpp\|cpp\|c\|h\|f\|f90\)$' \( -exec $HASHPRGM {} \; -o -print \) | $HASHPRGM | cut -d ' ' -f1`"
  # defaultshasum="`find $2/src $2/Source $2/SRC $2/include $2/inc $2/INC -type f -regex '.*\(hpp\|cpp\|c\|h\|f\|f90\)$' | xargs ls -al | $HASHPRGM | cut -d ' ' -f1`"
  defaultshasum="`cd $2/..; tar -tf $3 | xargs ls -l | $HASHPRGM | cut -d ' ' -f1`"
  AC_CACHE_VAL([ac_cv_sha_$1], [ac_cv_sha_$1="0"])
  if (test "$defaultshasum" != "$ac_cv_sha_$1" || test $need_configuration != false); then
    recompile_and_install=true
    ac_cv_sha_$1="$defaultshasum"
    AC_MSG_RESULT(yes)
  else
    recompile_and_install=false
    AC_MSG_RESULT(no)
  fi
])

dnl -------------------------------------------------------------
dnl Print out whether the configure, build, and install steps were succcessful
dnl -------------------------------------------------------------
AC_DEFUN([PRINT_AUTOMATION_STATUS],
[
  PREFIX_PRINT(Automation Status: )
  COLOR_PRINT([        xx  Configuration := ], m4_tolower($1)_configured)
  COLOR_PRINT([        xx  Build         := ], m4_tolower($1)_made)
  COLOR_PRINT([        xx  Installation  := ], m4_tolower($1)_installed)
])


dnl -------------------------------------------------------------
dnl AUSCM_CONFIGURE_EXTERNAL_PACKAGE(PACKAGE_NAME, DOWNLOAD_URL, PACKAGE_SRC_DIR, PACKAGE_INSTALL_DIR)
dnl Example: 
dnl AUSCM_CONFIGURE_EXTERNAL_PACKAGE(MOAB, "http://ftp.mcs.anl.gov/pub/fathom/moab-4.6-nightly.tar.gz", "moab-4.6-nightly.tar.gz")
dnl -------------------------------------------------------------
AC_DEFUN([AUSCM_CONFIGURE_EXTERNAL_PACKAGE],
[
  #m4_pushdef([pkg_short_name],[m4_tolower(m4_defn([CURRENT_PACKAGE]))])dnl
  pkg_short_name=m4_tolower($1)
  pkg_download_url="$2"
  current_build_dir=`pwd`
  pkg_archive_name="$3"
  pkg_srcdir="$current_build_dir/libraries/$pkg_short_name"
  download_ext_package=no

  # The default PACKAGE installation is under libraries
  pkg_install_dir="$current_build_dir/libraries/$MOAB_ARCH"
 
  AC_ARG_DOWNLOAD(m4_tolower($1),
    [AS_HELP_STRING([--download-m4_tolower($1)],[Download and configure $1 with default options (URL:$2)])],
    [case "x${downloadval}" in
		  xyes)  pkg_download_url="$2";   m4_tolower(enable$1)=yes;    download_ext_package=yes ;;
		  x)  pkg_download_url="$2";      m4_tolower(enable$1)=yes;    download_ext_package=yes ;;
		  *)  pkg_download_url="$downloadval"; m4_tolower(enable$1)=yes;   download_ext_package=yes ;;
		  xno)  pkg_download_url="none";  download_ext_package=no ;;
		esac],
    [pkg_download_url="$2"; download_ext_package=${m4_tolower(download$1)}])

  if (test "$download_ext_package" != "no") ; then
  
    # Check if the directory already exists
    # if found, we have already configured and linked the sources - do nothing
    # else, download, configure and make PACKAGE sources
    if (test ! -d "$pkg_install_dir" ); then
      AS_MKDIR_P($pkg_install_dir)
    fi
    
    if (test ! -d "$pkg_srcdir" ); then
      AS_MKDIR_P($pkg_srcdir)
    fi

    # Download the archive file containing the sources
    need_configuration=false
    need_build=false
    need_installation=false
    
    PPREFIX="$1"

	  MSG_ECHO_SEPARATOR

    pkg_archive_name="`basename $pkg_download_url`"

    # Check if we need to download an archive file
    DOWNLOAD_EXTERNAL_PACKAGE([$1], [$pkg_download_url], [$MOAB_PACKAGES_DIR/$pkg_archive_name])
    
    # Deflate the archive file containing the sources, if needed
    DEFLATE_EXTERNAL_PACKAGE([$1], [$MOAB_PACKAGES_DIR/$pkg_archive_name], [$pkg_srcdir])

    # Invoke the package specific configuration and build commands
    
    #m4_expand(m4_toupper([DEFAULT_CONFIGURE_MAKE_$1])([$1],"$pkg_srcdir","$pkg_install_dir", "$pkg_archive_name"))
    
    # Due to differences in autoconf we need to check if we should use m4_expand to call the package specific macros
    # Run the package preprocess and configure macros found in the package specific .m4 files
    m4_version_prereq(2.64, [ 
    	m4_expand(m4_toupper([AUSCM_AUTOMATED_SETUP_PREPROCESS_$1])([$1],"$pkg_srcdir","$pkg_install_dir", "$pkg_archive_name"))dnl
    	m4_expand(m4_toupper([AUSCM_AUTOMATED_CONFIGURE_$1])([$need_configuration]))dnl
    ],[
      	m4_toupper([AUSCM_AUTOMATED_SETUP_PREPROCESS_$1])([$1],"$pkg_srcdir","$pkg_install_dir", "$pkg_archive_name")dnl
    	  m4_toupper([AUSCM_AUTOMATED_CONFIGURE_$1])([$need_configuration])dnl
    ])

    CHECK_SOURCE_RECOMPILATION_HASH([$1],[$pkg_srcdir],[$MOAB_PACKAGES_DIR/$pkg_archive_name])
    
    # Run the build, install, and postprocess macros found in the package specific .m4 files.
    m4_version_prereq(2.64, [
    	m4_expand(m4_toupper([AUSCM_AUTOMATED_BUILD_$1])([$need_build]))dnl
    	m4_expand(m4_toupper([AUSCM_AUTOMATED_INSTALL_$1])([$need_installation]))dnl
    	m4_expand(m4_toupper([AUSCM_AUTOMATED_SETUP_POSTPROCESS_$1])([$1]))dnl
    ],[
	    m4_toupper([AUSCM_AUTOMATED_BUILD_$1])([$need_build])dnl
	    m4_toupper([AUSCM_AUTOMATED_INSTALL_$1])([$need_installation])dnl
	    m4_toupper([AUSCM_AUTOMATED_SETUP_POSTPROCESS_$1])([$1])dnl
    ])
    PRINT_AUTOMATION_STATUS([$1])

    # Determine if the installation process was successful
    if ( (test -f "$current_build_dir/libraries/$pkg_short_name/install_[]$pkg_short_name.log" ) &&
         ( "[$]m4_tolower($1)[]_configured" && "[$]m4_tolower($1)[]_made" && "[$]m4_tolower($1)[]_installed" ) ); then
      PREFIX_PRINT([Successful configure/build/install automated  (status=${GREEN}SUCCESS${NORMAL})])
      PREFIX_PRINT([Installed package under: $pkg_install_dir])
    else
      PREFIX_PRINT([Failed configure/build/install step  (status=${RED}FAILURE${NORMAL})])
    fi

	  MSG_ECHO_SEPARATOR

  fi  # if (test "$download_ext_package" != no) ; then 

  m4_tolower(download$1)="$download_ext_package"
  AC_SUBST(m4_tolower(download$1))

])

dnl ------------------------------------------------------------
dnl  Defines macros for printing colors, 
dnl  copying symlinks, and custom 
dnl  printing definitions.
dnl ------------------------------------------------------------
AC_DEFUN([COLOR_PRINT],
[
  if (test "x$PPREFIX" != "x"); then
    if [ ${$2} ]; then
      PREFIX_PRINT([$1 ${GREEN}SUCCESS${NORMAL}])
    else
      PREFIX_PRINT([$1 ${RED}FAIL${NORMAL}])
    fi
  else
    if [ ${$2} ]; then
      echo "$1 ${GREEN}SUCCESS${NORMAL}"
    else
      echo "$1 ${RED}FAIL${NORMAL}"
    fi
  fi
])


AC_DEFUN([RECURSIVE_COPY_DIR_SYMLINKS],
[
  if (test x$1 != x && test x$2 != x); then
    recs_srcdir="$1"
    recs_targetdir="$2"
    _AS_ECHO_LOG([--- Recursively copy directory symlinks ($PPREFIX) ---])
    if (test ! -d $recs_targetdir); then
      AS_MKDIR_P($recs_targetdir)
    fi
    _AS_ECHO_LOG([ Source directory: $recs_srcdir ])
    _AS_ECHO_LOG([ Target directory: $recs_targetdir ])
    recs_dirs="`rsync -ain $recs_srcdir/* $recs_targetdir --exclude '.svn' --exclude '.git' | grep 'cd+++' | cut -d ' ' -f2 `"
    if (test "x$recs_dirs" != "x"); then
      for recs_dname in $recs_dirs; do
        _AS_ECHO_LOG([Executing: $MKDIR_P $recs_targetdir/$recs_dname ])
        recs_mkdir_log="`mkdir -p $recs_targetdir/$recs_dname`"
      done
    fi
    recs_files="`rsync -ain $recs_srcdir/* $recs_targetdir --exclude '.svn' --exclude '.git' | grep 'f+++' | cut -d ' ' -f2`"
    for recs_fname in $recs_files; do
      _AS_ECHO_LOG([Executing: $LN_S -f $recs_srcdir/$recs_fname $recs_targetdir/$recs_fname ])
      recs_symlink_log="`$LN_S -f $recs_srcdir/$recs_fname $recs_targetdir/$recs_fname`"
    done
  fi
])


# Finds the parent path of a file or directory
# AC_FIND_ABSPATH(PATH TO A FILE OR DIR)
# ------------------------
m4_define([AC_FIND_ABSPATH], ["`perl -e 'use Cwd "abs_path";print abs_path(shift)' $1 | xargs dirname`"])


# Finds the parent path of a file or directory
# PREFIX_PRINT(PATH TO A FILE OR DIR)
# ------------------------
m4_define([PREFIX_PRINT], 
[_AS_ECHO_LOG([[[ $PPREFIX ]] --   $1 ]);
  AS_ECHO(["[[ $PPREFIX ]] --   $1 "])])


# Finds the parent path of a file or directory
# PREFIX_PRINTN(PATH TO A FILE OR DIR)
# ------------------------
m4_define([PREFIX_PRINTN], 
[_AS_ECHO_LOG([[[ $PPREFIX ]] --   $1 ]);
  AS_ECHO_N(["[[ $PPREFIX ]] --   $1 "])])


# MSG_ECHO_CUSTOM(FEATURE)
# ------------------------
m4_define([MSG_ECHO_CUSTOM],
[_AS_ECHO_LOG([$1]);
  AS_ECHO(["$1"])])

# MSG_ECHON_CUSTOM(FEATURE)
# ------------------------
m4_define([MSG_ECHON_CUSTOM],
[_AS_ECHO_LOG([$1]);
  AS_ECHO(["$1"])])

# MSG_ECHO_SEPARATOR
# ------------------------
m4_define([MSG_ECHO_SEPARATOR],
[ MSG_ECHO_CUSTOM([*** ================================================================================================================== ***]) ])

# Finds the parent path of a file or directory
# ECHO_EVAL(PATH TO A FILE, COMMAND)
# ------------------------
m4_define([ECHO_EVAL], 
[ echo "$2" >> $1;
  eval $2
])

