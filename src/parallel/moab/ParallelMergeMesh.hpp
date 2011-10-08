#ifndef PARALLELMERGEMESH_HPP
#define PARALLELMERGEMESH_HPP

#include "moab/Types.hpp"
#include <vector>
#include "moab/Range.hpp"
#include "moab/ParallelComm.hpp"

extern "C" 
{
#include "errmem.h"
#include "types.h"
#include "sort.h"
#include "tuple_list.h"
#include "crystal.h"
#include "transfer.h"
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

  //Public Function to identified shared elements
  ErrorCode merge();

private:
  ParallelComm *myPcomm;
  Interface *myMB;
  std::vector<Range> mySkinEnts;
  double myEps;
  tuple_list myTup, myMatches;
  crystal_data myCD;
  bool cdAllocated;

  //Wrapper of merge() that performs the merge
  ErrorCode PerformMerge();
  //Determine the local skin entities (fills mySkinEnts)
  ErrorCode PopulateMySkinEnts(int dim);
  //Get the global bounding box
  ErrorCode GetGlobalBox(double *gbox);
  //Fill out the local myTup before the first gather-scatter
  ErrorCode PopulateMyTup(double * gbox);
  //Once myTup is filled and gather scattered, figure out the matches
  ErrorCode PopulateMyMatches();
  //Sort the matching tuples
  ErrorCode SortMyMatches();
  //Tag the shared elements once the myMatches has been filled
  ErrorCode TagSharedElements(int dim);
  //Cleanup any data allocated by class members
  void CleanUp();
  //Partition the global box by the number of procs
  //Lengths and parts needs to be of length 3
  ErrorCode PartitionGlobalBox(double *gbox, double *lengths, int *parts);
  //A function for determining how many parts a side should be split into
  static int PartitionSide(double sideLeng, double restLen, unsigned numProcs, bool altRatio);

  //Swap 2 tuples
  static void SwapTuples(tuple_list *tup, 
			 unsigned long a, 
			 unsigned long b);
  
  //Sort a tuple list by its real values
  static void SortTuplesByReal(tuple_list *tup,
			       double eps2=0);

  ////The recursive sorting function
  static void PerformRealSort(tuple_list *tup, 
			      unsigned long left, 
			      unsigned long right,
			      double eps2);
  
  //Determines whether tuple i is greater than tuple j
  static bool TupleGreaterThan(tuple_list *tup, unsigned long vrI, unsigned long vrJ, double eps2);
};

} // namespace moab

#endif
