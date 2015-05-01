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

// Contributed by Lorenzo Alessio Botti (SpaFEDTe)
// This implementation is mostly borrowed from the mbzoltan MOAB partitioning tool

#ifndef __moabpartitioner_hpp__
#define __moabpartitioner_hpp__

#include <stdlib.h>
#include "moab_mpi.h"
#include "moab/Range.hpp"
#include "moab/Interface.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/Skinner.hpp"
#include "moab/WriteUtilIface.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/ParallelComm.hpp"
#include "MBTagConventions.hpp"
#include "moab/CN.hpp"
#include <vector>
#include "moab/Types.hpp"

#include "metis.h"

using namespace moab;

  class MOABPartitioner 
  {

  public:
    MOABPartitioner( Interface *impl = NULL,
                      const bool use_coords = false,
                      int argc = 0, 
                      char **argv = NULL);
    
    ~MOABPartitioner();

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

    ErrorCode write_partition(const int nparts, Range &elems, 
                                const int *assignment,
                                const bool write_as_sets,
                                const bool write_as_tags) = 0;

    ErrorCode write_file(const char *filename, const char *out_file) = 0;
 
    Range &part_sets() {return partSets;};
    
    const Range &part_sets() const {return partSets;};

  protected:

    bool write_output;
    bool useCoords;

    Range partSets;

  };

#endif

