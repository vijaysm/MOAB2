#ifndef MHDF_STATUS_INTERNAL_H
#define MHDF_STATUS_INTERNAL_H

#include "mhdf.h"

void mhdf_setOkay( mhdf_Status* );

void mhdf_setFail( mhdf_Status*, const char*, ... );

#endif
