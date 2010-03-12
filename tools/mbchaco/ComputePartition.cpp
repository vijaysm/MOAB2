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

/** Compute a partition of a MOAB mesh using Chaco library
 * 
 * 
 */

#include <iostream>
#include <ostream>

#include "ComputePartition.hpp"
#define IS_BUILDING_MB
#include "moab/Core.hpp"
#undef IS_BUILDING_MB
#include "moab/Range.hpp"
#include "moab/WriteUtilIface.hpp"
#include "moab/MeshTopoUtil.hpp"

#include "defs.h"
#include "params.h"

using namespace moab;

extern char *PARAMS_FILENAME;	/* name of file with parameter updates */
extern double EIGEN_TOLERANCE;	/* tolerance for eigen calculations */
extern long RANDOM_SEED;	/* seed for random number generators */
extern int MATCH_TYPE;      /* matching routine to call */
extern int FREE_GRAPH;	/* free graph data structure after reformat? */
FILE     *params_file;	/* file with parameter value updates */
double   *goal;		/* desired set sizes */
float    *x, *y, *z;	/* coordinates for inertial method */
int       global_method;	/* global partitioning method */
int       local_method;	/* local partitioning method */
short    *assignment;	/* set number of each vtx (length nvtxs+1) */
double    eigtol;		/* tolerance in eigenvector calculation */
int       ndims;		/* dimension of recursive partitioning */
int       ndims_tot;	/* total number of cube dimensions to divide */
int       mesh_dims[3];	/* dimensions of mesh of processors */
long      seed;		/* for random graph mutations */
int       rqi_flag;		/* use RQI/Symmlq eigensolver? */
int       vmax;		/* if so, how many vertices to coarsen down to? */

extern "C" 
{
  int interface(int, int*, int*, int*,
                float*, float*, float*, float*,
                char*, char*, short*,
                int, int, int*, double*, int, int, int, int, int,
                double, long);
  void read_params(FILE*);
}

#define RR if (MB_SUCCESS != result) return result

ErrorCode ComputePartition::compute_partition(const int nprocs, 
                                                const char *filename, 
                                                const bool write_file,
                                                const char *out_file) 
{
    // check input
  if (NULL == filename && NULL == out_file) return MB_FAILURE;

  bool new_moab = false;
  
  if (NULL == mbImpl) {
      // instantiate MOAB & read the mesh
    mbImpl = new Core();
    new_moab = true;
  
    if (NULL == mbImpl) return MB_FAILURE;
  }
  
  ErrorCode result = mbImpl->load_mesh(filename); RR;
  
    // assemble the graph
  short *assignment;
  std::vector<int> adjacencies, start;
  Range elems;
  
  result = assemble_graph(3, adjacencies, start, elems); RR;

    // read input to Chaco and drive partition computation
  result = drive_chaco(nprocs, adjacencies, start, assignment); RR;

    // take results & write onto partition sets
  result = write_partition(nprocs, elems, assignment); RR;

  free((char *) assignment);

  if (write_file) {
    std::string outfile;
    if (out_file != NULL) {
      outfile = out_file;
        // check extension
    }
    else {
      outfile = filename;
      outfile.append(".h5m");
    }
    
    result = mbImpl->write_mesh(outfile.c_str()); RR;
  }
  
  if (new_moab) {
    delete mbImpl;
    mbImpl = NULL;
  }
  
  return MB_SUCCESS;
}

ErrorCode ComputePartition::drive_chaco(const int nprocs,
                                          std::vector<int> &adjacencies, 
                                          std::vector<int> &start, 
                                          short *&assignment) 
{
  params_file = fopen(PARAMS_FILENAME, "r");
  if (params_file == NULL) {
    printf("Parameter file `%s' not found; using default parameters.\n",
           PARAMS_FILENAME);
  }

  x = y = z = NULL;
  goal = NULL;
  assignment = NULL;
  global_method = 1;
  rqi_flag = 0;
  vmax = 75; // recommended to be between 50..100 in Chaco UG

    // always default to KL for local method
  local_method = 1;

  if (nprocs > 0) {
    ndims_tot = 1;
    mesh_dims[0] = nprocs;
  }
  else {
    std::cout << "Must input nprocs > 0 for now." << std::endl;
    return MB_FAILURE;
  }

    // choose partitioning dimension based on nprocs for now
  if (nprocs < 8) ndims = 1;
  else if (nprocs < 16) ndims = 2;
  else ndims = 3;
  
  read_params(params_file);

  if (7 == global_method ||
      3 == global_method) {
    std::cout << "Chaco global method #" << global_method 
              << "7 not supported by mbchaco at this time." << std::endl;
    return MB_NOT_IMPLEMENTED;
  }

  if (global_method == 2 && MATCH_TYPE == 5) {
    std::cout << "Chaco spectral method with geometric matching not supported at this time."
              << std::endl;
    return MB_NOT_IMPLEMENTED;
  }

  assignment = (short *) malloc((unsigned) start.size() * sizeof(short));

  if (global_method == 3 ||
      (MATCH_TYPE == 5 && (global_method == 1 || 
                           (global_method == 2 && rqi_flag)))) {
    useCoords = true;
  }

  eigtol = EIGEN_TOLERANCE;
  seed = RANDOM_SEED;

    // set free_graph to false
  FREE_GRAPH = 0;
  
  interface(start.size()-1, &start[0], &adjacencies[0], /*vwgts*/ NULL, 
            /*ewgts*/ NULL, /*x*/ NULL, /*y*/ NULL, /*z*/ NULL,
            /*outassignptr*/ NULL, /*outfileptr*/ NULL,
            assignment,
            /*ARCHITECTURE*/ 1, ndims_tot, mesh_dims, goal,
            global_method, local_method, rqi_flag, vmax, ndims,
            eigtol, seed);

  if (params_file != NULL)
    fclose(params_file);

  return MB_SUCCESS;
}

ErrorCode ComputePartition::assemble_graph(const int dimension, 
                                             std::vector<int> &adjacencies, 
                                             std::vector<int> &start, 
                                             Range &elems) 
{
    // assemble a graph with vertices equal to elements of specified dimension, edges
    // signified by list of other elements to which an element is connected

    // get the elements of that dimension
  ErrorCode result = mbImpl->get_entities_by_dimension(0, dimension, elems);
  if (MB_SUCCESS != result || elems.empty()) return result;
  
    // get a tag for graph vertex number
  Tag gvert_id;
  result = mbImpl->tag_create("__graph_vertex_id", 4, MB_TAG_DENSE, gvert_id, NULL);
  if (MB_SUCCESS != result) return result;
  
    // assign increasing graph vertex ids
  WriteUtilIface *iface = NULL;
  result = mbImpl->query_interface("WriteUtilIface", reinterpret_cast<void**>(&iface));
  if (MB_SUCCESS != result || NULL == iface) return (MB_SUCCESS != result ? result : MB_FAILURE);

#define START_ID 1
  
  result = iface->assign_ids(elems, gvert_id, START_ID); RR;
  
    // now assemble the graph, calling MeshTopoUtil to get bridge adjacencies through d-1 dimensional
    // neighbors
  MeshTopoUtil mtu(mbImpl);
  Range adjs;
    // can use a fixed-size array 'cuz the number of lower-dimensional neighbors is limited
    // by MBCN
  int neighbors[MB_MAX_SUB_ENTITIES];
  
  for (Range::iterator rit = elems.begin(); rit != elems.end(); rit++) {

      // get bridge adjacencies
    adjs.clear();
    result = mtu.get_bridge_adjacencies(*rit, dimension-1, dimension, adjs); RR;
    
    
      // get the graph vertex ids of those
    if (!adjs.empty()) {
      result = mbImpl->tag_get_data(gvert_id, adjs, neighbors); RR;
    }

      // copy those into adjacencies vector
    start.push_back((int)adjacencies.size());
    std::copy(neighbors, neighbors+adjs.size(), std::back_inserter(adjacencies));
  }
  
    // need to push an extra value onto start, due to indexing starting at 1
  start.push_back((int)adjacencies.size());

  std::cout << "Start vector: " << std::endl;
  std::copy(start.begin(), start.end(), std::ostream_iterator<int>(std::cout, ", "));
  std::cout << std::endl;
  std::cout << "Adjacencies vector: " << std::endl;
  std::copy(adjacencies.begin(), adjacencies.end(), std::ostream_iterator<int>(std::cout, ", "));
  std::cout << std::endl;

 
  return MB_SUCCESS;
}

ErrorCode ComputePartition::write_partition(const int nprocs, Range &elems, 
                                              const short *assignment) 
{
    // first, create partition sets and store in vector
  EntityHandle *part_sets = new EntityHandle[nprocs];
  if (NULL == part_sets) return MB_FAILURE;
  
  ErrorCode result;
  for (int i = 0; i < nprocs; i++) {
    result = mbImpl->create_meshset(MESHSET_SET, part_sets[i]); RR;
  }
  
    // write a tag to those sets denoting they're partition sets, with a value of the
    // proc number
  Tag part_set_tag;
  int dum_id = -1;
  result = mbImpl->tag_create("PARALLEL_PARTITION", 4, MB_TAG_SPARSE, part_set_tag, &dum_id); RR;
  
  int *dum_ids = new int[nprocs];
  for (int i = 0; i < nprocs; i++) dum_ids[i] = i;
  
  result = mbImpl->tag_set_data(part_set_tag, part_sets, nprocs, dum_ids); RR;
  

    // assign entities to the relevant sets
  Range::iterator rit;
  int i;
  int *tag_data = new int[elems.size()];
  for (rit = elems.begin(), i = 0; rit != elems.end(); rit++, i++) {
    int part = assignment[i];
    result = mbImpl->add_entities(part_sets[part], &(*rit), 1); RR;
    tag_data[i] = part;
  }

    // allocate integer-size partitions
  Tag elem_part_set_tag;
  dum_id = -1;
  result = mbImpl->tag_create("ELEM_PARALLEL_PARTITION", 4, MB_TAG_SPARSE, MB_TYPE_INTEGER,
                              elem_part_set_tag, &dum_id); RR;
  result = mbImpl->tag_set_data(part_set_tag, elems, (void*) tag_data);
  if (MB_SUCCESS != result) return result;

  return MB_SUCCESS;
}
