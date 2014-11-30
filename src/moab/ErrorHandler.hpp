#ifndef MOAB_ERROR_HANDLER_HPP
#define MOAB_ERROR_HANDLER_HPP

#include "moab/Types.hpp"

#include <sstream>
#include <string.h>

namespace moab {

//! ErrorType - passed to the error handling routines indicating whether this is a new error (globally
//! fatal or per-processor relevant) to be created, or an existing one to be handled
enum ErrorType {MB_ERROR_TYPE_NEW_GLOBAL = 0, MB_ERROR_TYPE_NEW_LOCAL = 1, MB_ERROR_TYPE_EXISTING = 2};

//! Initialize MOAB error handler (e.g. create a utility object for printing error output)
void MBErrorHandler_Init();

//! Finalize MOAB error handler (e.g. delete the utility object for printing error output)
void MBErrorHandler_Finalize();

//! Indicates whether MBErrorHandler_Init has been called
bool MBErrorHandler_Initialized();

//! Get information about the last error
void MBErrorHandler_GetLastError(std::string& error);

//! Routine that is called to create a new error or handle an existing one
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

//! Set a new error with the given error message and return the given error code
//! Used in functions which return ErrorCode
#define SET_ERR(err_code, err_msg) \
  return MBError(__LINE__, __func__, __FILENAME__, __SDIR__, err_code, err_msg, MB_ERROR_TYPE_NEW_LOCAL)

//! Similar to SET_ERR except that the given error message is a stream instead of a string
#define SET_ERR_STR(err_code, err_msg_str) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg_str; \
    return MBError(__LINE__, __func__, __FILENAME__, __SDIR__, err_code, err_ostr.str().c_str(), MB_ERROR_TYPE_NEW_LOCAL); \
  } while (false)

//! Set a new error with the given error message and return (void)
//! Used in functions which return void
#define SET_ERR_RET_VOID(err_msg) \
  do { \
    MBError(__LINE__, __func__, __FILENAME__, __SDIR__, MB_FAILURE, err_msg, MB_ERROR_TYPE_NEW_LOCAL); \
    return; \
  } while (false)

//! Similar to SET_ERR_RET_VOID except that the given error message is a stream instead of a string
#define SET_ERR_STR_RET_VOID(err_msg_str) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg_str; \
    MBError(__LINE__, __func__, __FILENAME__, __SDIR__, MB_FAILURE, err_ostr.str().c_str(), MB_ERROR_TYPE_NEW_LOCAL); \
    return; \
  } while (false)

//! Set a new error with the given error message and return the given value
//! Used in functions which return any data type
#define SET_ERR_RET_VAL(err_msg, ret_val) \
  do { \
    MBError(__LINE__, __func__, __FILENAME__, __SDIR__, MB_FAILURE, err_msg, MB_ERROR_TYPE_NEW_LOCAL); \
    return ret_val; \
  } while (false)

//! Similar to SET_ERR_RET_VAL except that the given error message is a stream instead of a string
#define SET_ERR_STR_RET_VAL(err_msg_str, ret_val) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg_str; \
    MBError(__LINE__, __func__, __FILENAME__, __SDIR__, MB_FAILURE, err_ostr.str().c_str(), MB_ERROR_TYPE_NEW_LOCAL); \
    return ret_val; \
  } while (false)

//! Set a new error with the given error message and continue
//! Used in functions which return any data type
#define SET_ERR_CONT(err_msg) \
  MBError(__LINE__, __func__, __FILENAME__, __SDIR__, MB_FAILURE, err_msg, MB_ERROR_TYPE_NEW_LOCAL)

//! Similar to SET_ERR_CONT except that the given error message is a stream instead of a string
#define SET_ERR_STR_CONT(err_msg_str) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg_str; \
    MBError(__LINE__, __func__, __FILENAME__, __SDIR__, MB_FAILURE, err_ostr.str().c_str(), MB_ERROR_TYPE_NEW_LOCAL); \
  } while (false)

//! Similar to SET_ERR except that the error is considered globally fatal
#define SET_GLB_ERR(err_code, err_msg) \
  return MBError(__LINE__, __func__, __FILENAME__, __SDIR__, err_code, err_msg, MB_ERROR_TYPE_NEW_GLOBAL)

//! Similar to SET_ERR_STR except that the error is considered globally fatal
#define SET_GLB_ERR_STR(err_code, err_msg_str) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg_str; \
    return MBError(__LINE__, __func__, __FILENAME__, __SDIR__, err_code, err_ostr.str().c_str(), MB_ERROR_TYPE_NEW_GLOBAL); \
  } while (false)

//! Similar to SET_ERR_RET_VOID except that the error is considered globally fatal
#define SET_GLB_ERR_RET_VOID(err_msg) \
  do { \
    MBError(__LINE__, __func__, __FILENAME__, __SDIR__, MB_FAILURE, err_msg, MB_ERROR_TYPE_NEW_GLOBAL); \
    return; \
  } while (false)

//! Similar to SET_ERR_STR_RET_VOID except that the error is considered globally fatal
#define SET_GLB_ERR_STR_RET_VOID(err_msg_str) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg_str; \
    MBError(__LINE__, __func__, __FILENAME__, __SDIR__, MB_FAILURE, err_ostr.str().c_str(), MB_ERROR_TYPE_NEW_GLOBAL); \
    return; \
  } while (false)

//! Similar to SET_ERR_RET_VAL except that the error is considered globally fatal
#define SET_GLB_ERR_RET_VAL(err_msg, ret_val) \
  do { \
    MBError(__LINE__, __func__, __FILENAME__, __SDIR__, MB_FAILURE, err_msg, MB_ERROR_TYPE_NEW_GLOBAL); \
    return ret_val; \
  } while (false)

//! Similar to SET_ERR_STR_RET_VAL except that the error is considered globally fatal
#define SET_GLB_ERR_STR_RET_VAL(err_msg_str, ret_val) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg_str; \
    MBError(__LINE__, __func__, __FILENAME__, __SDIR__, MB_FAILURE, err_ostr.str().c_str(), MB_ERROR_TYPE_NEW_GLOBAL); \
    return ret_val; \
  } while (false)

//! Similar to SET_ERR_CONT except that the error is considered globally fatal
#define SET_GLB_ERR_CONT(err_msg) \
  MBError(__LINE__, __func__, __FILENAME__, __SDIR__, MB_FAILURE, err_msg, MB_ERROR_TYPE_NEW_GLOBAL)

//! Similar to SET_ERR_STR_CONT except that the error is considered globally fatal
#define SET_GLB_ERR_STR_CONT(err_msg_str) \
  do { \
    std::ostringstream err_ostr; \
    err_ostr << err_msg_str; \
    MBError(__LINE__, __func__, __FILENAME__, __SDIR__, MB_FAILURE, err_ostr.str().c_str(), MB_ERROR_TYPE_NEW_GLOBAL); \
  } while (false)

//! Check error code, if not MB_SUCCESS, call the error handler and return the given error code
//! Used in functions which return ErrorCode
#define CHK_ERR(err_code) \
  do { \
    if (MB_SUCCESS != err_code) \
      return MBError(__LINE__, __func__, __FILENAME__, __SDIR__, err_code, "", MB_ERROR_TYPE_EXISTING); \
  } while (false)

//! Check error code, if not MB_SUCCESS, call the error handler and return (void)
//! Used in functions which return void
#define CHK_ERR_RET_VOID(err_code) \
  do { \
    if (MB_SUCCESS != err_code) { \
      MBError(__LINE__, __func__, __FILENAME__, __SDIR__, err_code, "", MB_ERROR_TYPE_EXISTING); \
      return; \
    } \
  } while (false)

//! Check error code, if not MB_SUCCESS, call the error handler and return the given value
//! Used in functions which return any data type
#define CHK_ERR_RET_VAL(err_code, ret_val) \
  do { \
    if (MB_SUCCESS != err_code) { \
      MBError(__LINE__, __func__, __FILENAME__, __SDIR__, err_code, "", MB_ERROR_TYPE_EXISTING); \
      return ret_val; \
    } \
  } while (false)

//! Check error code, if not MB_SUCCESS, call the error handler and continue
//! Used in functions which return any data type
#define CHK_ERR_CONT(err_code) \
  do { \
    if (MB_SUCCESS != err_code) { \
      MBError(__LINE__, __func__, __FILENAME__, __SDIR__, err_code, "", MB_ERROR_TYPE_EXISTING); \
    } \
  } while (false)

//! Check error code, if not MB_SUCCESS, set a new error with the given error message and return the given error code
//! Used in functions which return ErrorCode
#define CHK_SET_ERR(err_code, err_msg) \
  do { \
    if (MB_SUCCESS != err_code) \
      SET_ERR(err_code, err_msg); \
  } while (false)

//! Similar to CHK_SET_ERR except that the given error message is a stream instead of a string
#define CHK_SET_ERR_STR(err_code, err_msg_str) \
  do { \
    if (MB_SUCCESS != err_code) \
      SET_ERR_STR(err_code, err_msg_str); \
  } while (false)

//! Check error code, if not MB_SUCCESS, set a new error with the given error message and return (void)
//! Used in functions which return void
#define CHK_SET_ERR_RET_VOID(err_code, err_msg) \
  do { \
    if (MB_SUCCESS != err_code) \
      SET_ERR_RET_VOID(err_msg); \
  } while (false)

//! Similar to CHK_SET_ERR_RET_VOID except that the given error message is a stream instead of a string
#define CHK_SET_ERR_STR_RET_VOID(err_code, err_msg_str) \
  do { \
    if (MB_SUCCESS != err_code) \
      SET_ERR_STR_RET_VOID(err_msg_str); \
  } while (false)

//! Check error code, if not MB_SUCCESS, set a new error with the given error message and return the given value
//! Used in functions which return any data type
#define CHK_SET_ERR_RET_VAL(err_code, err_msg, ret_val) \
  do { \
    if (MB_SUCCESS != err_code) \
      SET_ERR_RET_VAL(err_msg, ret_val); \
  } while (false)

//! Similar to CHK_SET_ERR_RET_VAL except that the given error message is a stream instead of a string
#define CHK_SET_ERR_STR_RET_VAL(err_code, err_msg_str, ret_val) \
  do { \
    if (MB_SUCCESS != err_code) \
      SET_ERR_STR_RET_VAL(err_msg_str, ret_val); \
  } while (false)

//! Check error code, if not MB_SUCCESS, set a new error with the given error message and continue
//! Used in functions which return any data type
#define CHK_SET_ERR_CONT(err_code, err_msg) \
  do { \
    if (MB_SUCCESS != err_code) \
      SET_ERR_CONT(err_msg); \
  } while (false)

//! Similar to CHK_SET_ERR_CONT except that the given error message is a stream instead of a string
#define CHK_SET_ERR_STR_CONT(err_code, err_msg_str) \
  do { \
    if (MB_SUCCESS != err_code) \
      SET_ERR_STR_CONT(err_msg_str); \
  } while (false)

} // namespace moab

#endif
