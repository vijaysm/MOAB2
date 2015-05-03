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

#ifndef __metispartitioner_hpp__
#define __metispartitioner_hpp__

#include <stdlib.h>
#include "moab/PartitionerBase.hpp"
#include "metis.h"

namespace moab {

  class Interface;
  class Range;
}

using namespace moab;

  class MetisPartitioner : public PartitionerBase
  {

  public:
    MetisPartitioner( Interface *impl = NULL,
                          const bool use_coords = false,
                          int argc = 0, 
                          char **argv = NULL);
    
    virtual ~MetisPartitioner();

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

    virtual ErrorCode write_partition(const int nparts, Range &elems, 
                                const int *assignment,
                                const bool write_as_sets,
                                const bool write_as_tags);
    
    ErrorCode write_aggregationtag_partition(const int nparts, Range &elems, 
                                             const int *assignment,
                                             const bool write_as_sets,
                                             const bool write_as_tags);

    // virtual ErrorCode write_file(const char *filename, const char *out_file);
  
  private:

    int argcArg;
    
    char **argvArg;

    ErrorCode assemble_graph(const int dimension, 
                             std::vector<double> &coords,
                             std::vector<int> &moab_ids,
                             std::vector<int> &adjacencies, 
                             std::vector<int> &length,
                             Range &elems);
    
    ErrorCode assemble_taggedsets_graph(const int dimension, 
                                        std::vector<double> &coords,
                                        std::vector<int> &moab_ids,
                                        std::vector<int> &adjacencies, 
                                        std::vector<int> &length,
                                        Range &elems,
                                        const char *aggregating_tag);
    
    ErrorCode assemble_taggedents_graph(const int dimension, 
                                        std::vector<double> &coords,
                                        std::vector<int> &moab_ids,
                                        std::vector<int> &adjacencies, 
                                        std::vector<int> &length,
                                        Range &elems,
                                        const char *aggregating_tag);
  };

#endif

