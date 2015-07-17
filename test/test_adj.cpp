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




#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif

#include <iostream>
#include "moab/Core.hpp"

#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif

#include "Internals.hpp"

using namespace moab;

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

int main()
{

  Interface* iface = new Core;

  const char* filename = 0;

#ifdef MESHDIR
#ifdef MOAB_HAVE_HDF5
  filename = STRINGIFY(MESHDIR) "/testquad-cyl.h5m";
#else
  filename = STRINGIFY(MESHDIR) "/testquad-cyl.g";
#endif
#else
#ifdef MOAB_HAVE_HDF5
  filename = "testquad-cyl.h5m";
#else
  filename = "testquad-cyl.g";
#endif
#endif

  ErrorCode err;
  err =  iface->load_mesh( filename);MB_CHK_ERR(err);

  Range quads, verts;
  err = iface->get_entities_by_dimension(0,0,verts);MB_CHK_ERR(err);
  err = iface->get_entities_by_dimension(0, 2, quads);MB_CHK_ERR(err);

  for (Range::iterator it = verts.begin(); it != verts.end(); it++)
    std::cout<<"verts["<<(*it-*verts.begin())<<"] = "<<*it<<std::endl;

  std::vector<EntityHandle> conn;
  for (Range::iterator it = quads.begin(); it != quads.end(); it++)
    {
      conn.clear();
      err = iface->get_connectivity(&*it, 1, conn);MB_CHK_ERR(err);
      std::cout<<"verts["<<(*it-*quads.begin())<<"] = "<<*it<<" :: conn = ["<<conn[0]<<", "<<conn[1]<<", "<<conn[2]<<", "<<conn[3]<<"]"<<std::endl;
    }

  std::vector<EntityHandle> nodes;
  int error;
  EntityHandle h = CREATE_HANDLE(MBQUAD, MB_START_ID, error);

  std::cout<<"h = "<<h<<std::endl;

  err = iface->get_adjacencies( &h, 1, 0, true, nodes);MB_CHK_ERR(err);

  for (int i=0; i<(int)nodes.size(); i++)
    std::cout<<"nodes["<<i<<"] = "<<nodes[i]<<" ";
  std::cout<<std::endl;

  std::vector<EntityHandle> edgs;
  err = iface->get_adjacencies( &h, 1, 1, true, edgs);MB_CHK_ERR(err);

  for (int i=0; i<(int)edgs.size(); i++)
    std::cout<<"edgs["<<i<<"] = "<<edgs[i]<<" ";
  std::cout<<std::endl;

  //faces to nodes

  for(unsigned int i=0; i<edgs.size(); i++)
  {
      nodes.clear();
      err = iface->get_adjacencies( &edgs[i], 1, 0, true, nodes );MB_CHK_ERR(err);
      std::cout << "edge " << ID_FROM_HANDLE(edgs[i]) << std::endl;
      std::cout << "nodes = ";
      for(unsigned int j=0; j<nodes.size(); j++)
        std::cout << nodes[j] << " ";
      std::cout << std::endl;
    }

  delete iface;

  return 0;
}



