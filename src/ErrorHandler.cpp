#include "moab/ErrorHandler.hpp"
#include "ErrorOutput.hpp"
#ifdef USE_MPI
#include "moab_mpi.h"
#endif

#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

namespace moab {

static ErrorOutput* errorOutput = NULL;
static bool hasError = false;

void MBErrorHandler_Init()
{
  if (NULL == errorOutput) {
    errorOutput = new (std::nothrow) ErrorOutput(stderr);
    assert(NULL != errorOutput);
    errorOutput->use_world_rank();
    hasError = false;
  }
}

void MBErrorHandler_Finalize()
{
  if (NULL != errorOutput) {
    delete errorOutput;
    errorOutput = NULL;
    hasError = false;
  }
}

bool MBErrorHandler_Initialized()
{
  return (NULL != errorOutput);
}

void MBTraceBackErrorHandler(int line, const char* func, const char* file, const char* dir, const char* err_msg, ErrorType err_type)
{
  if (NULL == errorOutput)
    return;

  // For a globally fatal error, get world rank of current processor, so that it is only printed from processor 0
  // For a per-processor relevant error, set rank of current processor to 0, so that it is always printed
  int rank = 0;
  if (MB_ERROR_TYPE_NEW_GLOBAL == err_type && errorOutput->have_rank())
    rank = errorOutput->get_rank();

  if (0 == rank) {
    // Print the error messages if it is a new error
    if (MB_ERROR_TYPE_EXISTING != err_type && NULL != err_msg) {
      errorOutput->print("--------------------- Error Message ------------------------------------\n");
      errorOutput->printf("%s!\n", err_msg);
      hasError = true;
    }

    // Print a line of stack trace for a new error or an existing one
    if (hasError)
      errorOutput->printf("%s() line %d in %s%s\n", func, line, dir, file);
  }
  else {
    // Do not print the error messages, since processor 0 will print them
    // Sleep 10 seconds before aborting so it will not accidently kill process 0
    sleep(10);
    abort();
  }
}

ErrorCode MBError(int line, const char* func, const char* file, const char* dir, ErrorCode err_code, const char* err_msg, ErrorType err_type)
{
  MBTraceBackErrorHandler(line, func, file, dir, err_msg, err_type);

#ifdef USE_MPI
  // If this is called from the main() routine we call MPI_Abort() to allow
  // the parallel program to be properly shutdown
  if (strncmp(func, "main", 4) == 0)
    MPI_Abort(MPI_COMM_WORLD, err_code);
#endif

  return err_code;
}

} // namespace moab
