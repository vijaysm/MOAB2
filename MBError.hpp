/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */


/*  
 *
 *  File:      MBError.hpp
 *
 *  Purpose:   To keep track of detail information about errors that occur
 *             in MB.
 *
 *  Date:      9-26-2002
 *
 *  Author:    Clinton Stimpson
 *
 */



#ifndef MB_ERROR_HPP
#define MB_ERROR_HPP

#ifndef IS_BUILDING_MB
#error "MBError.hpp isn't supposed to be included into an application"
#endif

#include <string>
#include <stdarg.h>
#include <stdio.h>

#ifdef WIN32
#define VSNPRINTF _vsnprintf
#else
#define VSNPRINTF vsnprintf
#endif

class MBError
{
  //! string to hold the last error that occurred in MB
  std::string mLastError;

public:

  MBError() {}
  ~MBError(){}

  MBErrorCode set_last_error(const std::string& error) 
  { 
    mLastError = error; 
    return MB_SUCCESS; 
  }

  MBErrorCode set_last_error(const char* fmt, ...)
  {
    if(!fmt)
      return MB_FAILURE;
    
    va_list ap;
    static char text[1024];

    va_start(ap, fmt);
    VSNPRINTF(text, 1024, fmt, ap);
    va_end(ap);

    mLastError = text;

    return MB_SUCCESS;
  }

  MBErrorCode get_last_error(std::string& error) const
  { 
    error = mLastError; 
    return MB_SUCCESS;
  }

};



#endif


