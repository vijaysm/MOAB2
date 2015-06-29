/*! \file MeasureTime.hpp
 * This class provides a timing routine that can be included with any code to measure performance of code snippets.
*/

#ifndef MEASURE_TIME_HPP
#define MEASURE_TIME_HPP

#include "moab/MOABConfig.h"
#ifdef MOAB_HAVE_MPI
#include <mpi.h>
#else
#include <sys/time.h>
#endif

namespace moab {

  class MeasureTime{
  public:
    MeasureTime() {};
    ~MeasureTime();

    double wtime();

  };

  double MeasureTime::wtime()
  {
    double y = -1;
#ifdef MOAB_HAVE_MPI
    y = MPI_Wtime();
#else
    struct timeval cur_time;
    gettimeofday(&cur_time, NULL);
    y = (double)(cur_time.tv_sec) + (double)(cur_time.tv_usec)*1.e-6;
#endif
    return (y);
  }
}

#endif
