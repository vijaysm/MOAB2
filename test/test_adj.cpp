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

int main()
{

  Interface* iface = new Core;

  iface->load_mesh( "../test/quad/quad1.gen" );

  std::vector<EntityHandle> nodes;
  int err;
  EntityHandle h = CREATE_HANDLE(MBQUAD, MB_START_ID, err);
  err = iface->get_adjacencies( &h, 1, 0, true, nodes);

  std::vector<EntityHandle> tris;
  err = iface->get_adjacencies( &h, 1, 1, true, tris);

  //faces to nodes

  for(unsigned int i=0; i<tris.size(); i++)
  {
    err = iface->get_adjacencies( &tris[i], 1, 0, true, nodes );
    std::cout << "edge " << ID_FROM_HANDLE(tris[i]) << std::endl;
    std::cout << "nodes = ";
    for(unsigned int j=0; j<nodes.size(); j++)
      std::cout << nodes[j] << " ";
    std::cout << std::endl;
  }

  delete iface;

  return 0;
}



