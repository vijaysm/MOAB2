#ifndef CPUTIMER_HPP
#define CPUTIMER_HPP

#include "moab/MOABConfig.h"
#ifdef MOAB_HAVE_MPI
#  include "moab_mpi.h"
/* it is needed if MPI is not initialized */
#  include <time.h>
#else
#  if defined(_MSC_VER)
#    include <time.h>
#  else
#    include <sys/resource.h>
#  endif
#endif

namespace moab 
{
    
class CpuTimer {
private:
  double tAtBirth, tAtLast;
  double mAtBirth, mAtLast;
  long rssAtBirth, rssAtLast;

  double runtime();
  long runmem();
  
public:
  CpuTimer() : tAtBirth(runtime()), tAtLast(tAtBirth) {}
  double time_since_birth() { return (tAtLast = runtime()) - tAtBirth; };
  double time_elapsed() { double tmp = tAtLast; return (tAtLast = runtime()) - tmp; }
  long mem_since_birth() {return (mAtLast=runmem()) - mAtBirth;}
  long mem_elapsed() {long tmp = mAtLast; return (mAtLast=runmem()) - tmp;}
};

    inline double CpuTimer::runtime() 
    {
#if defined(_MSC_VER) || defined(__MINGW32__)
      return (double)clock() / CLOCKS_PER_SEC;
#elif defined(MOAB_HAVE_MPI)
      int flag=0;
      if (MPI_SUCCESS==MPI_Initialized(&flag) && flag)
      {
        return MPI_Wtime();
      }
      else
      {
        return (double)clock() / CLOCKS_PER_SEC;
      }
#else      
      return (double)clock() / CLOCKS_PER_SEC;
#endif
    }

    inline long CpuTimer::runmem() 
    {
#if defined(_MSC_VER) || defined(__MINGW32__)
      return 0;
#elif defined(MOAB_HAVE_MPI)
      return 0;
#else
      struct rusage r_usage;
      getrusage(RUSAGE_SELF, &r_usage);
      mAtLast = r_usage.ru_maxrss;
      return mAtLast;
#endif
    }

}

#endif
