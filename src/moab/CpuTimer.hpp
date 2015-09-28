#ifndef CPUTIMER_HPP
#define CPUTIMER_HPP

#include "moab/MOABConfig.h"
#ifdef MOAB_HAVE_MPI
#  include "moab_mpi.h"
#endif

#include <time.h>

namespace moab 
{
    
class CpuTimer {
private:
  int mpi_initialized;
  double tAtBirth, tAtLast;
  double runtime();
public:
  CpuTimer() : mpi_initialized(0)
  {
#ifdef MOAB_HAVE_MPI
  	int flag=0;
  	if (MPI_SUCCESS==MPI_Initialized(&flag) && flag)
  	{
  		mpi_initialized=1;
  	}
#endif
  	tAtBirth = runtime();
  	tAtLast = tAtBirth;
  }
  double time_since_birth() { return (tAtLast = runtime()) - tAtBirth; };
  double time_elapsed() { double tmp = tAtLast; return (tAtLast = runtime()) - tmp; }
};

inline double CpuTimer::runtime()
{
#ifdef MOAB_HAVE_MPI
	if (mpi_initialized)
		return MPI_Wtime();
	else
		return (double)clock() / CLOCKS_PER_SEC;
#else
    return (double)clock() / CLOCKS_PER_SEC;
#endif
    }
}

#endif
