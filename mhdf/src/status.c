#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "status.h"


int mhdf_isError( mhdf_Status const* status )
{
  return !!status->message[0];
}

const char* mhdf_message( mhdf_Status const* status )
{
  return status->message;
}

void mhdf_setOkay( mhdf_Status* status )
{
  if (status) status->message[0] = '\0';
}

void mhdf_setFail( mhdf_Status* status, const char* fmt, ... )
{
  if (status)
  {
    va_list args;
    va_start( args, fmt );
    vsnprintf( status->message, MHDF_MESSAGE_BUFFER_LEN, fmt, args );
    va_end(args);
    if (!status->message[0])
      strncpy( status->message, "(Uknown error)", MHDF_MESSAGE_BUFFER_LEN );
  }
}

