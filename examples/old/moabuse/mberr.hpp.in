/// \brief The MBERR macro is useful for checking the return value from MOAB
///        function calls
/// 
/// It takes an error message as a string and a moab return code. If
/// rval is not moab::MB_SUCCESS, the following will happen:
///
/// -# The current file and line number will be printed
/// -# The error code (rval) will be printed
/// -# The error message (msg) will be printed
/// -# rval will be thrown as an exception
///
/// If rval is moab::MB_SUCCESS then nothing will happen
///
/// \param msg [in] The error message
/// \param rval [in] The moab return value to check
#define MBERR(msg, rval) mberr(msg, rval, __LINE__, __FILE__)

void mberr(const std::string& msg, moab::ErrorCode rval, 
           unsigned line, const std::string &file)
{
  if(rval != moab::MB_SUCCESS) {
    std::cerr << file << ":" << line << ", MOAB Error: " << rval << "\n"
              << msg << std::endl;
    throw(rval);
  }
}
