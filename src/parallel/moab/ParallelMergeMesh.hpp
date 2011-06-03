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

/*
  Class to merge meshes in parallel
  Requires a ParallelComm and tolerance epsilon
  Currently uses a 1 dimensional partition of the global box
*/

namespace moab {

class ParallelComm;

class ParallelMergeMesh {
public:
  ParallelMergeMesh(ParallelComm *pc, 
		    const double epsilon);

  //Perform the merge
  ErrorCode merge();

private:
  ParallelComm *myPcomm;
  double myEps;

  //Swap 2 tuples
  static void tuple_swap_real(tuple_list *tup, 
			      unsigned long a, 
			      unsigned long b);
  
  //Sort a tuple list by its real values
  static void tuple_sort_real(tuple_list *tup,
			      double eps2=0);

  ////The recursive sorting function
  static void perform_sort_real(tuple_list *tup, 
				unsigned long left, 
				unsigned long right,
				double eps2);
  
  //Determines whether tuple i is greater than tuple j
  static bool greater_than(tuple_list *tup, unsigned long vrI, unsigned long vrJ, double eps2);
};

} // namespace moab

#endif
