#ifndef PARALLELMERGEMESH_HPP
#define PARALLELMERGEMESH_HPP

#include "moab/Types.hpp"
extern "C" 
{
#include "errmem.h"
#include "types.h"
#include "sort.h"
#include "tuple_list.h"
}

namespace moab {

class ParallelComm;

class ParallelMergeMesh {
public:
  ParallelMergeMesh(ParallelComm *pc, 
		    const double epsilon);

  ErrorCode merge();

private:
  ParallelComm *myPcomm;
  double myEps;

  static void temp_tuple_swap_real(tuple_list *tup, 
                                   unsigned long a, unsigned long b);

  static void temp_tuple_sort_real(tuple_list *tup,
				   double eps2=0);

  static void temp_perform_sort_real(tuple_list *tup, 
                                     unsigned long left, 
                                     unsigned long right,
                                     double eps2);

  static bool greaterThan(tuple_list *tup, unsigned long vrI, unsigned long vrJ, double eps2);
};

} // namespace moab

#endif
