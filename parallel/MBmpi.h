/* MBmpi.h.  Generated from MBmpi.h.in by configure.  */
#ifndef MB_MPI_H
#define MB_MPI_H
#include "MBmpi_config.h"

#ifndef __cplusplus
#  include <mpi.h>
#elif !defined(MB_MPI_CXX_CONFLICT)
#  ifndef MPICH_IGNORE_CXX_SEEK
#    define MPICH_IGNORE_CXX_SEEK
#  endif
#  include <mpi.h>
#else
#  include <stdio.h>
#  ifdef SEEK_SET
#    define MB_SEEK_SET SEEK_SET
#    define MB_SEEK_CUR SEEK_CUR
#    define MB_SEEK_END SEEK_END
#    undef SEEK_SET
#    undef SEEK_CUR
#    undef SEEK_END
#  endif
#  include <mpi.h>
#  ifdef MB_SEEK_SET
#    define SEEK_SET MB_SEEK_SET
#    define SEEK_CUR MB_SEEK_CUR
#    define SEEK_END MB_SEEK_END
#    undef MB_SEEK_SET
#    undef MB_SEEK_CUR
#    undef MB_SEEK_END
#  endif
#endif


#endif
