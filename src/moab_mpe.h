#ifndef MOAB_MPE_H
#define MOAB_MPE_H

#define MOAB_MPE_FIRST_EVENT 101

#ifdef __cplusplus
extern "C" {
#endif
  int MPE_Allocate_event(void);
#ifdef __cplusplus
 } /* extern "C" */
#endif


#ifdef USE_MPE
# include "moab_mpi.h"
# include "mpe.h"
#else

/* Define dummy logging functions */

/* mpe_misc.h */

#define MPE_Seq_begin( A, B )
#define MPE_Seq_end( A, B )

#define MPE_DelTag( A, B, C, D ) (MPI_SUCCESS)
#define MPE_GetTags( A, B, C, D ) (A = *C, *D = 0, MPI_SUCCESS)
#define MPE_ReturnTags( A, B, C ) (0)
#define MPE_TagsEnd() (MPI_SUCCESS)

#define MPE_IO_Stdout_to_file( A, B )

#define MPE_GetHostName( A, B )

#define MPI_Start_debugger()

/*
#if (defined(__STDC__) || defined(__cplusplus))
#define  MPE_Errors_to_dbx ( MPI_Comm *, int *, ... );
#else
void MPE_Errors_to_dbx ( MPI_Comm *, int *, char *, char *, int * );
#endif
void MPE_Errors_call_debugger ( char *, char *, char ** );
void MPE_Errors_call_xdbx     ( char *, char * );
void MPE_Errors_call_dbx_in_xterm ( char *, char * );
void MPE_Signals_call_debugger ( void );

int  MPE_Decomp1d ( int, int, int, int *, int * );

void MPE_Comm_global_rank ( MPI_Comm, int, int * );
*/

/* mpe_log.h */

#define MPE_LOG_OK                0
#define MPE_Log_OK                MPE_LOG_OK
  /* no problems */
#define MPE_LOG_LOCKED_OUT        1
#define MPE_Log_LOCKED_OUT        MPE_LOG_LOCKED_OUT
  /* logs are being worked on, cannot insert any new entries */
#define MPE_LOG_NO_MEMORY         2
#define MPE_Log_NO_MEMORY         MPE_LOG_NO_MEMORY
  /* could not allocate memory for logging data */
#define MPE_LOG_FILE_PROB         3
#define MPE_Log_FILE_PROB         MPE_LOG_FILE_PROB
  /* cound not open file for writing out the logged info */
#define MPE_LOG_NOT_INITIALIZED   4
#define MPE_Log_NOT_INITIALIZED   MPE_LOG_NOT_INITIALIZED
  /* logging not initialized */
#define MPE_LOG_PACK_FAIL         5
#define MPE_Log_PACK_FAIL         MPE_LOG_PACK_FAIL

#define MPE_Init_log() (MPI_SUCCESS)
#define MPE_Initialized_logging() 1

#define MPE_Describe_state( A, B, C, D )
#define MPE_Describe_event( A, B, C )
#define MPE_Log_get_event_number()
#define MPE_Log_send( A, B, C )
#define MPE_Log_receive( A, B, C )
#define MPE_Log_event( A, B, C )
#define MPE_Start_log()
#define MPE_Stop_log()
#define MPE_Finish_log( A )

#endif

#endif /* MOAB_MPE_H */
