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

#include "ReadVtk.hpp"
#include "MBReadUtilIface.hpp"
#include "MBRange.hpp"
#include "MBInternals.hpp"

#include <fstream>
#include <iostream>
#include <assert.h>

MBReaderIface* ReadVtk::factory( MBInterface* iface )
  { return new ReadVtk( iface ); }

ReadVtk::ReadVtk(MBInterface* impl)
    : mdbImpl(impl)
{
  assert(impl != NULL);
  std::string iface_name = "MBReadUtilIface";
  impl->query_interface(iface_name, reinterpret_cast<void**>(&readMeshIface));
}

MBErrorCode ReadVtk::load_file(const char *file_name,
                               const int*, const int) 
{
  MBErrorCode result;
  fileName = file_name;
  
  std::ifstream ifs(file_name);

  if (!ifs) {
    std::cerr << "file " << file_name << " not found\n";
    return MB_FILE_DOES_NOT_EXIST;
  }

  char line[256];
  int i;
  for (i=0; i<4; ++i)
    ifs.getline(line,256);

    // gets the number of vertices.
  char word[81];
  ifs >> word;
  if (strcmp(word, "POINTS") != 0) {
    std::cerr << "ReadVtk: expecting word POINTS, not " << word << " .\n";
    return MB_FAILURE;
  }
  
  int NbVertices;
  ifs >> NbVertices;
  ifs.getline(line,256);

  MBEntityHandle start_handle = 0;
  std::vector<double*> arrays;
  readMeshIface->get_node_arrays(3, NbVertices,
                                 MB_START_ID, start_handle, arrays);

  for(i=0;i<NbVertices;++i)
    ifs >> (arrays[0])[i] >> (arrays[1])[i] >> (arrays[2])[i];

    // make a meshset for this mesh
  result = mdbImpl->create_meshset(MESHSET_SET, mCurrentMeshHandle);
  if (MB_SUCCESS != result) return result;
  
    // add the vertices
  MBRange vert_range(start_handle, start_handle+NbVertices-1);
  result = mdbImpl->add_entities(mCurrentMeshHandle, vert_range);
  if (MB_SUCCESS != result) return result;

    // gets the number of regions
  ifs >> word;
  if (strcmp(word, "CELLS") != 0) {
    std::cerr << "ReadVtk: expecting word CELLS, not " << word << " .\n";
    return MB_FAILURE;
  }

  int NbRegions, tableSize;
  ifs >> NbRegions >> tableSize;

    // fills up the connectivity table
  std::vector<int> connectivity;
  for (i=0; i<tableSize; ++i) {
    int entry;
    ifs >> entry;
    connectivity.push_back(entry);
  }

    // gets the number of CELL_TYPES identifier
  ifs >> word;
  if (strcmp(word, "CELL_TYPES") != 0) {
    std::cerr << "ReadVtk: expecting word CELL_TYPES, not " << word << " .\n";
    return MB_FAILURE;
  }

  int NbCellTypes;
  ifs >> NbCellTypes;
  if (NbCellTypes != NbRegions) {
    std::cerr << "ReadVtk: nb of CELL_TYPES != nb of CELLS.\n";
    return MB_FAILURE;
  }

  int* cellType = new int[NbCellTypes];
  for (i=0; i < NbCellTypes; ++i)
    ifs >> cellType[i];

  MBEntityHandle vtx[8];
  std::vector<int>::const_iterator connectIter;
  connectIter = connectivity.begin();
  MBEntityHandle new_region;
  
  for (i=0; i<NbRegions; ++i) {

    switch (cellType[i]) {
      case 5:  // Triangle
        if (*connectIter++ != 3) {
          std::cerr << "ReadVtk: expecting 3 vtx for a triangle.\n";
          return MB_FAILURE;
        }
        vtx[0] = start_handle+*connectIter++;
        vtx[1] = start_handle+*connectIter++;
        vtx[2] = start_handle+*connectIter++;
        result = mdbImpl->create_element(MBTRI, vtx, 3, new_region);
        if (MB_SUCCESS != result) return result;
        break;
      case 9:  // Quad
        if (*connectIter++ != 4) {
          std::cerr << "ReadVtk: expecting 4 vtx for a Quad.\n";
          return MB_FAILURE;
        }
        vtx[0] = start_handle+*connectIter++;
        vtx[1] = start_handle+*connectIter++;
        vtx[2] = start_handle+*connectIter++;
        vtx[3] = start_handle+*connectIter++;
        result = mdbImpl->create_element(MBQUAD, vtx, 4, new_region);
        if (MB_SUCCESS != result) return result;
        break;
      case 10: // Tet
        if (*connectIter++ != 4) {
          std::cerr << "ReadVtk: expecting 4 vtx for a Tet.\n";
          return MB_FAILURE;
        }
        vtx[0] = start_handle+*connectIter++;
        vtx[1] = start_handle+*connectIter++;
        vtx[2] = start_handle+*connectIter++;
        vtx[3] = start_handle+*connectIter++;
        result = mdbImpl->create_element(MBTET, vtx, 4, new_region);
        if (MB_SUCCESS != result) return result;
        break;
      case 12: // Hex 
        if (*connectIter++ != 8) {
          std::cerr << "ReadVtk: expecting 8 vtx for an Hex.\n";
          return MB_FAILURE;
        }
        vtx[0] = start_handle+*connectIter++;
        vtx[1] = start_handle+*connectIter++;
        vtx[2] = start_handle+*connectIter++;
        vtx[3] = start_handle+*connectIter++;
        vtx[4] = start_handle+*connectIter++;
        vtx[5] = start_handle+*connectIter++;
        vtx[6] = start_handle+*connectIter++;
        vtx[7] = start_handle+*connectIter++;
        result = mdbImpl->create_element(MBHEX, vtx, 8, new_region);
        if (MB_SUCCESS != result) return result;
        break;
    }

      // add the new region to the set
    result = mdbImpl->add_entities(mCurrentMeshHandle, &new_region, 1);
    if (MB_SUCCESS != result) return result;
  }

    // initialize all vertices with non-boundary tags 
  MBTag bdy_tag;
  result = mdbImpl->tag_get_handle("boundary", bdy_tag);
  int zero = 0;
  if (MB_TAG_NOT_FOUND == result) {
    result = mdbImpl->tag_create("boundary", 4, MB_TAG_DENSE, bdy_tag, &zero);
  }
  if (MB_SUCCESS != result) return result;
  
    // gets the number of BOUNDARY_POINTS
  ifs >> word;
  if (strcmp(word, "POINT_DATA") != 0) {
    std::cerr << "ReadVtk: no boundary points specified" << word << " .\n";
  } else {
    int NbVtxWithValue;
    ifs >> NbVtxWithValue;
    for (int k=1; k<=5; k++) ifs >> word;

    MBEntityHandle dum_vert;
    for (int vtxNb=0; vtxNb<NbVtxWithValue; ++vtxNb)
    {
      int vtxValue;
      ifs >> vtxValue;
      if (vtxValue == 0 || vtxValue == 1) {
        dum_vert = start_handle + vtxNb;
        result = mdbImpl->tag_set_data(bdy_tag, &dum_vert, 1, &vtxValue);
        if (MB_SUCCESS != result) return result;
      }
      else
        std::cerr << "mVtk.cc ERROR: invalid boundary value for vertex.\n"; 
    }
  }
  ifs.close();

  delete[] cellType;

  return MB_SUCCESS;
}

