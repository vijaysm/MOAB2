#ifndef MOAB_ERROR_HANDLER_HPP
#define MOAB_ERROR_HANDLER_HPP

#include "moab/Types.hpp"

#include <sstream>
#include <string.h>

namespace moab {

//! ErrorType - passed to the error handling routines indicating if this is a new error (globally fatal
//! or per-processor relevant) or an existing one
enum ErrorType {MB_ERROR_TYPE_NEW_GLOBAL = 0, MB_ERROR_TYPE_NEW_LOCAL = 1, MB_ERROR_TYPE_EXISTING = 2};

//! Initialize MOAB error handler (e.g. create a utility object for printing error output)
void MBErrorHandler_Init();

//! Finalize MOAB error handler (e.g. delete the utility object for printing error output)
void MBErrorHandler_Finalize();

//! Indicates whether MBErrorHandler_Init has been called
bool MBErrorHandler_Initialized();

//! Routine that is called when an error has been detected
ErrorCode MBError(int line, const char* func, const char* file, const char* dir,
                  ErrorCode err_code, const char* err_msg, ErrorType err_type);

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#define MBSTRINGIFY_(X) #X
#define MBSTRINGIFY(X) MBSTRINGIFY_(X)

#ifdef LOCDIR
#define __SDIR__ MBSTRINGIFY(LOCDIR)
#else
#define __SDIR__ ""
#endif

//! Set a new error with passed error code and passed error message string
#define SET_ERR(err_code, err_msg) \
  return MBError(__LINE__, __func__, __FILENAME__, __SDIR__, err_code, err_msg, MB_ERROR_TYPE_NEW_LOCAL)

//! Set a new error with passed error code and passed error message string stream
#define SET_ERR_STR(err_code, err_msg_str) \
  do { \
    std::ostringstream ostr; \
    ostr << err_msg_str; \
    return MBError(__LINE__, __func__, __FILENAME__, __SDIR__, err_code, ostr.str().c_str(), MB_ERROR_TYPE_NEW_LOCAL); \
  } while (false)

//! Set a new global error with passed error code and passed error message string
#define SET_GLB_ERR(err_code, err_msg) \
  return MBError(__LINE__, __func__, __FILENAME__, __SDIR__, err_code, err_msg, MB_ERROR_TYPE_NEW_GLOBAL)

//! Set a new global error with passed error code and passed error message string stream
#define SET_GLB_ERR_STR(err_code, err_msg_str) \
  do { \
    std::ostringstream ostr; \
    ostr << err_msg_str; \
    return MBError(__LINE__, __func__, __FILENAME__, __SDIR__, err_code, ostr.str().c_str(), MB_ERROR_TYPE_NEW_GLOBAL); \
  } while (false)

//! Check returned error code against MB_SUCCESS
#define CHK_ERR(err_code) \
  do { \
    if (MB_SUCCESS != err_code) \
      return MBError(__LINE__, __func__, __FILENAME__, __SDIR__, err_code, "", MB_ERROR_TYPE_EXISTING); \
  } while (false)

//! Check returned error code against MB_SUCCESS
//! Set a new error with the returned error code and passed error message string
#define CHK_ERR1(err_code, err_msg_to_set) \
  do { \
    if (MB_SUCCESS != err_code) \
      SET_ERR(err_code, err_msg_to_set); \
  } while (false)

//! Check returned error code against MB_SUCCESS
//! Set a new error with the returned error code and passed error message string stream
#define CHK_ERR1_STR(err_code, err_msg_str_to_set) \
  do { \
    if (MB_SUCCESS != err_code) \
      SET_ERR_STR(err_code, err_msg_str_to_set); \
  } while (false)

} // namespace moab

#endif
