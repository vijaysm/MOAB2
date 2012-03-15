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
 * Zoltan: class to get a mesh from MOAB and write a Zoltan partition set for
 * that mesh back into MOAB and to a file
 *
 */

#ifndef MB_ZOLTAN_HPP
#define MB_ZOLTAN_HPP

#include <stdlib.h>
#include "moab_mpi.h"
#include "zoltan_cpp.h"
#include "moab/Range.hpp"

#ifdef CGM
#include <map>
#include "GeometryQueryTool.hpp"
#include "DLIList.hpp"
class RefEntity;
#endif

extern "C" 
{
  int mbGetNumberOfAssignedObjects(void *userDefinedData, int *err);
  
  void mbGetObjectList(void *userDefinedData, int numGlobalIds, int numLids,
                       ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int wgt_dim, float *obj_wgts,
                       int *err);
  
  int mbGetObjectSize(void *userDefinedData, int *err);
  
  void mbGetObject(void *userDefinedData, int numGlobalIds, int numLids, int numObjs,
                   ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int numDim, double *pts, int *err);
  
  void mbGetNumberOfEdges(void *userDefinedData, int numGlobalIds, int numLids,
                          int numObjs, 
                          ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,	int *numEdges,
                          int *err);
  
  void mbGetEdgeList(void *userDefinedData, int numGlobalIds, int numLids,
                     int numObjs,
                     ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int *numEdges,
                     ZOLTAN_ID_PTR nborGlobalIds, int *nborProcs, int wgt_dim,
                     float *edge_wgts, int *err);

  void mbGetPart(void *userDefinedData, int numGlobalIds, int numLids,
                 int numObjs, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                 int *part, int *err);
  
  void mbShowError(int val, const char *s, int me);
  
  void mbPrintGlobalResult(const char *s, 
                           int begin, int import, int exp, int change);
}

#include <vector>
#include "moab/Types.hpp"

namespace moab {

  class ParallelComm;
  class Interface;
  class Range;
}

using namespace moab;

  class MBZoltan 
  {

  public:
    MBZoltan( Interface *impl = NULL,
              const bool use_coords = false,
              int argc = 0, 
              char **argv = NULL
#ifdef CGM
              , GeometryQueryTool *gqt = NULL
#endif
              );

    
    ~MBZoltan();

    ErrorCode balance_mesh(const char *zmethod,
                           const char *other_method,
                           const bool write_as_sets = true,
                           const bool write_as_tags = false);

    ErrorCode partition_mesh_geom(const double part_geom_mesh_size,
                                  const int nparts,
                                  const char *zmethod,
                                  const char *other_method,
                                  double imbal_tol,
                                  const bool write_as_sets = true,
                                  const bool write_as_tags = false,
                                  const int part_dim = 3,
                                  const int obj_weight = 0,
                                  const int edge_weight = 0,
                                  const bool part_surf = false,
                                  const bool ghost = false);
    
    int get_mesh(std::vector<double> &pts, std::vector<int> &ids,
                 std::vector<int> &adjs, std::vector<int> &length,
                 Range &elems);

      // given a processor assignment returned from Zoltan, write that as a
      // processor assignment to MOAB
    ErrorCode write_partition(const int nparts, Range &elems, 
                              const int *assignment,
                              const bool write_as_sets,
                              const bool write_as_tags);

#ifdef CGM
    ErrorCode write_partition(const int nparts,
                              DLIList<RefEntity*> entities,
                              const int *assignment,
                              std::vector<double> &obj_weights,
                              const bool part_surf,
                              const bool ghost);

    ErrorCode partition_surface(const int nparts,
                                DLIList<RefEntity*> entities,
                                const int *assignment,
                                std::vector<double> &obj_weights);
#endif
    
      // put closure of entities in the part sets too
    ErrorCode include_closure();
    
    ErrorCode write_file(const char *filename, const char *out_file);
  
    void SetOCTPART_Parameters(const char *oct_method);
  
    void SetPARMETIS_Parameters(const char *parmetis_method);
  
    void SetHypergraph_Parameters(const char *phg_method);
  
    void SetHSFC_Parameters();
  
    void SetRIB_Parameters();
  
    void SetRCB_Parameters();

    Range &part_sets() {return partSets;};
    
    const Range &part_sets() const {return partSets;};
  
  private:

    Interface *mbImpl;

    ParallelComm *mbpc;

    Zoltan *myZZ;

    Range partSets;
  
    bool newMoab;
  
    bool useCoords;

    bool write_output;

    int myNumPts;

    int argcArg;
    
    char **argvArg;

    int mbGlobalSuccess(int rc);
  
    void mbPrintGlobalResult(const char *s,
                             int begin, int import, int exp, int change);
  
    void mbShowError(int val, const char *s);
  
      // given the dimension, assemble the vertices and store in coords and
      // moab_ids
    ErrorCode assemble_graph(const int dimension, 
                             std::vector<double> &coords,
                             std::vector<int> &moab_ids,
                             std::vector<int> &adjacencies, 
                             std::vector<int> &length,
                             Range &elems);
    
#ifdef CGM
    std::map<int, int> body_vertex_map, surf_vertex_map;

    std::vector<double> obj_weights;

    ErrorCode assemble_graph(const int dimension, 
                             std::vector<double> &coords,
                             std::vector<int> &moab_ids,
                             std::vector<int> &adjacencies, 
                             std::vector<int> &length,
                             std::vector<double> &obj_weights,
                             std::vector<double> &edge_weights,
                             std::vector<int> &parts,
                             DLIList<RefEntity*> &entities,
                             const double part_geom_mesh_size,
                             const int n_part);

    ErrorCode partition_round_robin(const int n_part);

    ErrorCode partition_child_entities(const int dim,
                                       const int n_part,
                                       const bool partition_surf,
                                       const bool ghost = false);

    double estimate_face_mesh_load(RefEntity* face, const double h);
    double estimate_face_comm_load(RefEntity* face, const double h);
#endif
    
    void mbFinalizePoints(int npts, int numExport,
                          ZOLTAN_ID_PTR exportLocalIDs, int *exportProcs,
                          int **assignment);
  
    int mbInitializePoints(int npts, double *pts, int *ids, 
                           int *adjs, int *length,
                           double *obj_weights = NULL,
                           double *edge_weights = NULL,
                           int *parts = NULL);

#ifdef CGM
    GeometryQueryTool *gti;
#endif
  };

#endif
