/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 *
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */

#ifndef __partitioner_base_hpp__
#define __partitioner_base_hpp__

#include <stdlib.h>
#include "moab_mpi.h"
#include "moab/Range.hpp"
#include <vector>
#include "moab/Types.hpp"

#include "moab/ParallelComm.hpp"
namespace moab {

  class Interface;
  class Range;
}

using namespace moab;

  class PartitionerBase 
  {

  public:
    PartitionerBase( Interface *impl = NULL,
                      const bool use_coords = false);
    
    ~PartitionerBase();

    /*
    ErrorCode partition_mesh_geom(const int nparts,
                                  const char *method,
                                  const int part_dim = 3, 
                                  const bool write_as_sets = true,
                                  const bool write_as_tags = false,
                                  const bool partition_tagged_sets = false,
                                  const bool partition_tagged_ents = false,
                                  const char *aggregating_tag = NULL);
    
    int get_mesh(std::vector<double> &pts, std::vector<int> &ids,
                 std::vector<int> &adjs, std::vector<int> &length,
                 Range &elems);
    */

    virtual ErrorCode write_partition(const int nparts, Range &elems, 
                                const int *assignment,
                                const bool write_as_sets,
                                const bool write_as_tags) = 0;

    virtual ErrorCode write_file(const char *filename, const char *out_file) = 0;
 
    Range &part_sets() {return partSets;};
    
    const Range &part_sets() const {return partSets;};

  protected:

    Interface *mbImpl;

    ParallelComm *mbpc;
    
    bool write_output;
    bool useCoords;
    bool newComm;

    Range partSets;

  };

inline
PartitionerBase::PartitionerBase(Interface *impl,
                                  const bool use_coords)
    : mbImpl(impl), useCoords(use_coords), newComm(false)
{
  mbpc = ParallelComm::get_pcomm(mbImpl, 0);
  if (!mbpc) {
    mbpc = new ParallelComm(impl, MPI_COMM_WORLD, 0);
    newComm = true;
  }
}

inline
PartitionerBase::~PartitionerBase()
{
  if (newComm)
    delete mbpc;

  mbImpl = NULL;
}

#endif

