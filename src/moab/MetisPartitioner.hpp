// SpaFEDTe, a Template based C++ library for creating 
// Discontinuous Finite Element Spaces,
// Copyright (C) 2012 Lorenzo Alessio Botti

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */ 
/* License as published by the Free Software Foundation either */ 
/* version 3.0 of the License, or (at your option) any later version. */

/* This software is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this software; if not, a copy of the full */
/* GNU Lesser General Public License can be found at */
/* http://www.gnu.org/licenses/ */

// This implementation is mostly borrowed from the mbzoltan MOAB partitioning tool

#ifndef __metis_moab_partitioner_hpp__
#define __metis_moab_partitioner_hpp__

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

  class MetisMOABPartitioner 
  {

  public:
    MetisMOABPartitioner( Interface *impl = NULL,
                          const bool use_coords = false,
                          int argc = 0, 
                          char **argv = NULL);
    
    ~MetisMOABPartitioner();

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
                                const bool write_as_tags);
    
    ErrorCode write_aggregationtag_partition(const int nparts, Range &elems, 
                                             const int *assignment,
                                             const bool write_as_sets,
                                             const bool write_as_tags);

    ErrorCode write_file(const char *filename, const char *out_file);
  
  private:

    Interface *mbImpl;

    ParallelComm *mbpc;

    Range partSets;
  
    bool useCoords;

    bool write_output;

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
