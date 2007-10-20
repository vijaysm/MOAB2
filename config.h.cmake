/* config.h.cmake.  Generated from thin air by dcthomp.  */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
#undef FC_DUMMY_MAIN

/* Define if F77 and FC dummy `main' functions are identical. */
#undef FC_DUMMY_MAIN_EQ_F77

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#undef FC_FUNC

/* As FC_FUNC, but for C identifiers containing underscores. */
#undef FC_FUNC_

/* Define to 1 if you have the <dlfcn.h> header file. */
#undef HAVE_DLFCN_H

/* Define to 1 if you have the <gvc.h> header file. */
#undef HAVE_GVC_H

/* Define to 1 if you have the <hdf5.h> header file. */
#cmakedefine HDF5_FOUND
#ifdef HDF5_FOUND
#  define HAVE_HDF5_H
#endif /* HDF5_FOUND */

/* Define to 1 if you have the <inttypes.h> header file. */
#cmakedefine HAVE_INTTYPES_H

/* Define to 1 if you have the `chaco' library (-lchaco). */
#cmakedefine HAVE_LIBCHACO

/* Define to 1 if you have the `gvc' library (-lgvc). */
#cmakedefine HAVE_LIBGVC

/* Define to 1 if you have the <memory.h> header file. */
#cmakedefine HAVE_MEMORY_H

/* Define to 1 if you have the <netcdf.h> header file. */
#cmakedefine HAVE_NETCDF_H

/* Define to 1 if you have the <qapplication.h> header file. */
#cmakedefine HAVE_QAPPLICATION_H

/* Define to 1 if you have the <qevent.h> header file. */
#cmakedefine HAVE_QEVENT_H

/* Define to 1 if you have the <qlineedit.h> header file. */
#cmakedefine HAVE_QLINEEDIT_H

/* Define to 1 if you have the <qmetaobject.h> header file. */
#cmakedefine HAVE_QMETAOBJECT_H

/* Define to 1 if you have the <qobject.h> header file. */
#cmakedefine HAVE_QOBJECT_H

/* Define to 1 if you have the <qpixmap.h> header file. */
#cmakedefine HAVE_QPIXMAP_H

/* Define to 1 if you have the <qtimer.h> header file. */
#cmakedefine HAVE_QTIMER_H

/* Define to 1 if you have the <qwidgetplugin.h> header file. */
#cmakedefine HAVE_QWIDGETPLUGIN_H

/* Define to 1 if you have the <qwidget.h> header file. */
#cmakedefine HAVE_QWIDGET_H

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
#cmakedefine HAVE_STDLIB_H

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine HAVE_STRINGS_H

/* Define to 1 if you have the <string.h> header file. */
#cmakedefine HAVE_STRING_H

/* Define to 1 if you have the <sys/stat.h> header file. */
#cmakedefine HAVE_SYS_STAT_H

/* Define to 1 if you have the <sys/types.h> header file. */
#cmakedefine HAVE_SYS_TYPES_H

/* Define to 1 if you have the <unistd.h> header file. */
#cmakedefine HAVE_UNISTD_H

/* MOAB Version */
#define MB_VERSION @MOAB_VERSION@

/* MOAB Major Version */
#define MB_VERSION_MAJOR @MOAB_VERSION_MAJOR@

/* MOAB Minor Version */
#define MB_VERSION_MINOR @MOAB_VERSION_MINOR@

/* MOAB Patch Level */
#define MB_VERSION_PATCH @MOAB_VERSION_PATCH@

/* MOAB Version String */
#define MB_VERSION_STRING "@MOAB_VERSION_STRING@"

/* Use int32_t for handles */
#cmakedefine MOAB_FORCE_32_BIT_HANDLES

/* Use int64_t for handles */
#cmakedefine MOAB_FORCE_64_BIT_HANDLES

/* MOAB qualified HAVE_INTTYPES_H */
#cmakedefine MOAB_HAVE_INTTYPES_H

/* System provides ptrdiff_t typedef */
#cmakedefine MOAB_HAVE_PTRDIFF_T

/* System provides size_t typedef */
#cmakedefine MOAB_HAVE_SIZE_T

/* MOAB qualified HAVE_STDDEF_H */
#cmakedefine MOAB_HAVE_STDDEF_H

/* MOAB qualified HAVE_STDINT_H */
#cmakedefine MOAB_HAVE_STDINT_H

/* MOAB qualified HAVE_STDLIB_H */
#cmakedefine MOAB_HAVE_STDLIB_H

/* MOAB qualified HAVE_SYS_TYPES_H */
#cmakedefine MOAB_HAVE_SYS_TYPES_H

/* Name of package */
#define PACKAGE "${PROJECT_NAME}"

/* Define to the address where bug reports for this package should be sent. */
#undef PACKAGE_BUGREPORT

/* Define to the full name of this package. */
#define PACKAGE_NAME "${PROJECT_NAME}"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "${PROJECT_NAME} ${${PROJECT_NAME}_VERSION_STRING}"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "${PROJECT_NAME}"

/* Define to the version of this package. */
#define PACKAGE_VERSION "${${PROJECT_NAME}_VERSION_STRING}"

/* Define to 1 if you have the ANSI C header files. */
#undef STDC_HEADERS

/* Version number of package */
#undef VERSION

/* Define to 1 if the X Window System is missing or not being used. */
#undef X_DISPLAY_MISSING
