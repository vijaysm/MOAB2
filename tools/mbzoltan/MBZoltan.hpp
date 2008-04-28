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
 * MBZoltan: class to get a mesh from MOAB and write a Zoltan partition set for
 * that mesh back into MOAB and to a file
 *
 */

#ifndef MB_ZOLTAN_HPP
#define MB_ZOLTAN_HPP

#include "zoltan_cpp.h"

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
  
  void mbShowError(int val, const char *s, int me);
  
  void mbPrintGlobalResult(const char *s, 
                           int begin, int import, int exp, int change);
}

#include <vector>

#include "MBInterface.hpp"

  class MBZoltan 
  {

  public:
    MBZoltan(MBInterface *impl = NULL, const bool use_coords = false,
             int argc = 0, char **argv = NULL) 
        : mbImpl(impl), myZZ(NULL), newMoab(false), useCoords(use_coords),
          argcArg(argc), argvArg(argv)
      {}

    ~MBZoltan();

    MBErrorCode balance_mesh(const char *zmethod,
                             const char *other_method);
    
    MBErrorCode partition_mesh(const int nparts,
                               const char *zmethod,
                               const char *other_method);
    
    int get_mesh(std::vector<double> &pts, std::vector<int> &ids,
                 std::vector<int> &adjs, std::vector<int> &length,
                 MBRange &elems);

      // given a processor assignment returned from Zoltan, write that as a
      // processor assignment to MOAB
    MBErrorCode write_partition(const int nparts, MBRange &elems, const int *assignment);

    MBErrorCode write_file(const char *filename, const char *out_file);
  
    void SetOCTPART_Parameters(const char *oct_method);
  
    void SetPARMETIS_Parameters(const char *parmetis_method);
  
    void SetHypergraph_Parameters(const char *phg_method);
  
    void SetHSFC_Parameters();
  
    void SetRIB_Parameters();
  
    void SetRCB_Parameters();
  
  private:

    MBInterface *mbImpl;

    Zoltan *myZZ;
  
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
    MBErrorCode assemble_graph(const int dimension, 
                               std::vector<double> &coords,
                               std::vector<int> &moab_ids,
                               std::vector<int> &adjacencies, 
                               std::vector<int> &length,
                               MBRange &elems);

    void mbFinalizePoints(int npts, int numExport,
                          ZOLTAN_ID_PTR exportLocalIDs, int *exportProcs,
                          int **assignment);
  
    int mbInitializePoints(int npts, double *pts, int *ids, 
                           int *adjs, int *length);
  
  };

#endif
