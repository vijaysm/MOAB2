/**
 * ComputePartition: class to partition a mesh based on MOAB and Chaco
 *
 */

#ifndef COMPUTE_PARTITION_HPP
#define COMPUTE_PARTITION_HPP

#include <vector>

#include "MBInterface.hpp"

class ComputePartition 
{

public:
  ComputePartition(MBInterface *impl = NULL, const bool use_coords = false) 
      : mbImpl(impl), newMoab(false), useCoords(use_coords)
    {}

  ~ComputePartition() {}

    // compute a partition; NULL filename means MOAB already contains the mesh
  MBErrorCode compute_partition(const int nprocs, const char *filename = NULL, 
                                const bool write_file = false, const char *out_file = NULL);
  
private:

  MBInterface *mbImpl;

  bool newMoab;
  
  bool useCoords;

  MBErrorCode drive_chaco(const int nprocs,
                          std::vector<int> &adjacencies, 
                          std::vector<int> &start, 
                          short *&assignment);
  
    // given the dimension, assemble the graph and store in adjacencies and start
  MBErrorCode assemble_graph(const int dimension, 
                             std::vector<int> &adjacencies, 
                             std::vector<int> &start, 
                             MBRange &elems);
  
    // given a processor assignment returned from Chaco, write that as a processor
    // assignment to MOAB
  MBErrorCode write_partition(const int nprocs, MBRange &elems, 
                              const short *assignment);
};

#endif
