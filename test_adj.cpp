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
#include "MBCore.hpp"
#include "ReadWriteExoII.hpp"


int main()
{

  MBInterface* iface = new MBCore;

  ReadWriteExoII exo_reader( static_cast<MBCore*>(iface));
  exo_reader.load_file("../test/quad/quad1.gen", 0);

  std::vector<MBEntityHandle> nodes;
  int err;
  iface->get_adjacencies( CREATE_HANDLE(MBQUAD, 0, err), 0, nodes, &err);

  std::vector<MBEntityHandle> tris;
  iface->get_adjacencies( CREATE_HANDLE(MBQUAD, 0, err), 1, tris, &err);

  //faces to nodes

  for(int i=0; i<tris.size(); i++)
  {
    iface->get_adjacencies( tris[i], 0, nodes, &err );
    std::cout << "edge " << ID_FROM_HANDLE(tris[i]) << std::endl;
    std::cout << "nodes = ";
    for(int j=0; j<nodes.size(); j++)
      std::cout << nodes[j] << " ";
    std::cout << std::endl;
  }

  delete iface;

  return 0;
}



