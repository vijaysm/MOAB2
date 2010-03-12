/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/**
 * ComputePartition: class to partition a mesh based on MOAB and Chaco
 *
 */

#ifndef COMPUTE_PARTITION_HPP
#define COMPUTE_PARTITION_HPP

#include <vector>

#include "moab/Interface.hpp"

using namespace moab;

class ComputePartition 
{

public:
  ComputePartition(Interface *impl = NULL, const bool use_coords = false) 
      : mbImpl(impl), newMoab(false), useCoords(use_coords)
    {}

  ~ComputePartition() {}

    // compute a partition; NULL filename means MOAB already contains the mesh
  ErrorCode compute_partition(const int nprocs, const char *filename = NULL, 
                                const bool write_file = false, const char *out_file = NULL);
  
private:

  Interface *mbImpl;

  bool newMoab;
  
  bool useCoords;

  ErrorCode drive_chaco(const int nprocs,
                          std::vector<int> &adjacencies, 
                          std::vector<int> &start, 
                          short *&assignment);
  
    // given the dimension, assemble the graph and store in adjacencies and start
  ErrorCode assemble_graph(const int dimension, 
                             std::vector<int> &adjacencies, 
                             std::vector<int> &start, 
                             Range &elems);
  
    // given a processor assignment returned from Chaco, write that as a processor
    // assignment to MOAB
  ErrorCode write_partition(const int nprocs, Range &elems, 
                              const short *assignment);
};

#endif
