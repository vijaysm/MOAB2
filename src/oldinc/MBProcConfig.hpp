#ifndef MBProcConfig_HEADER
#define MBProcConfig_HEADER

#include "MBTypes.h"
#include "MBRange.hpp"
#ifdef USE_MPI
#  include "MBmpi.h"
#endif

#include "moab/ProcConfig.hpp"
typedef moab::ProcConfig MBProcConfig;

#endif
