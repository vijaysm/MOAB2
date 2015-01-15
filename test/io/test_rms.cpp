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

#include "TSTT_MB_QueryInterface.h"
#include "moab/_RMBSet.hpp"

using namespace moab;

int compare_coords(double *xval, double *yval, double *zval,
                   double *nodex, double *nodey, double *nodez, 
                   const int num_nodes);
int compare_connect(int *connect1, int *connect2, const int num_comps);

#include <iostream>

int main() 
{
    // very basic test of RMBSet
    // The test: instantiate a mesh with 2 tet elements sharing a face;
    // represent just the 5 nodes (1-5) and 2 elements (1234 & 3215)
    // 

    // nodal coordinates: shared tri in x-y plane with edges on x and y
    // axes and a node at the origin, plus nodes at z = +- 1
    //
    //  x4                
    //  . .              
    //  .  x1               z
    //   . .\               /\  /\y
    //   . .  \ E1            . |
    //    ..    \              .|
    //    2x------x3            .---->x
    //      .     .      
    //       .E2 .       
    //        . .        
    //         x5         
    //         . .               
    //         .  x6             
    //          . .\ 
    //          . .  \ E3
    //           ..    \         
    //           7x------x8      
    //             .     .       
    //              .   .  E4
    //               . .         
    //                x9
    //      
    //      
  double nodex[] = {0.0, 0.0, 1.0, 0.0, 0.0, 
                    0.0, 0.0, 1.0, 0.0};
  double nodey[] = {1.0, 0.0, 0.0, 0.0, 0.0, 
                    1.0, 0.0, 0.0, 0.0};
  double nodez[] = {0.0, 0.0, 0.0, 1.0, -1.0, 
                    -2.0, -2.0, -2.0, -3.0};
  int connect[] = {
    1, 2, 3, 4,
    3, 2, 1, 5,
    6, 7, 8, 5,
    8, 7, 6, 9
  };
  
  const int NUM_NODES = 9;

    // construct the RMS holding the nodes
  MB_RMBSet node_rms(NUM_NODES, nodex, nodey, nodez, 1, 0);
  
    // construct two RMS's holding each pair of elements
  MB_RMBSet elem_rms1(2, connect, 1, TSTT_REGION, TSTT_TETRAHEDRON);
  MB_RMBSet elem_rms2(2, &connect[8], 3, TSTT_REGION, TSTT_TETRAHEDRON);
  
  
    // now do some querying on this mesh

    // INFO FUNCTIONS
  int entity_type1 = node_rms.entity_type();
  int entity_type2 = elem_rms1.entity_type();
  if (entity_type1 != TSTT_VERTEX ||
      entity_type2 != TSTT_REGION)
    std::cout << "entity_type() function failed." << std::endl;
  
  int entity_topo = elem_rms1.entity_topology();
  if (entity_topo != TSTT_TETRAHEDRON)
    std::cout << "entity_topology() function failed." << std::endl;
  
  int num_ents1 = node_rms.num_entities();
  int num_ents2 = elem_rms1.num_entities();
  if (num_ents1 != NUM_NODES || num_ents2 != 2)
    std::cout << "num_entities() function failed for" << 
      (num_ents1 != NUM_NODES ? "(nodes)" : "") <<
      (num_ents2 != 2 ? "(elems)" : "") <<
      std::endl;
  
  int vpe = elem_rms1.vertices_per_element();
  if (vpe != 4)
    std::cout << "vertices_per_element() failed." << std::endl;
  
    // NODES
/*
    // get_coordinates
  double coords[15];
  coords[1] = node_rms.get_coordinates(1);
  coords[2] = node_rms.get_coordinates(2);
  coords[3] = node_rms.get_coordinates(3);
  coords[4] = node_rms.get_coordinates(4);
  coords[5] = node_rms.get_coordinates(5);
  int result = compare_coords(coords, nodex, nodey, nodez, 5);
  if (result == 0) std::cout << "get_coordinates works." << std::endl;
  else std::cout << "get_coordinates didn't work; result = " << result << "." << std::endl;
*/
    // node_x, node_y, node_z
  int num_nodes = NUM_NODES;
  double xval[NUM_NODES], yval[NUM_NODES], zval[NUM_NODES];
  double *xvalp = &xval[0], *yvalp = &yval[0], *zvalp = &zval[0];
  node_rms.node_x(1, NUM_NODES, &xvalp, &num_nodes);
  node_rms.node_y(1, NUM_NODES, &yvalp, &num_nodes);
  node_rms.node_z(1, NUM_NODES, &zvalp, &num_nodes);
  int result = compare_coords(xval, yval, zval, nodex, nodey, nodez, NUM_NODES);
  if (result != 0) 
    std::cout << "node_[xyz] didn't work; result = " << result << "." << std::endl;

    // set_node_x, set_node_y, set_node_z
  int i;
  for (i = 1; i <= NUM_NODES; i++) {
    nodex[i-1] = (double) i;
    nodey[i-1] = (double) i;
    nodez[i-1] = (double) i;
  }
  node_rms.set_node_x(1, NUM_NODES, nodex, NUM_NODES);
  node_rms.set_node_y(1, NUM_NODES, nodey, NUM_NODES);
  node_rms.set_node_z(1, NUM_NODES, nodez, NUM_NODES);
  node_rms.node_x(1, NUM_NODES, &xvalp, &num_nodes);
  node_rms.node_y(1, NUM_NODES, &yvalp, &num_nodes);
  node_rms.node_z(1, NUM_NODES, &zvalp, &num_nodes);
  result = compare_coords(xval, yval, zval, nodex, nodey, nodez, NUM_NODES);
  if (result != 0) 
    std::cout << "node_[xyz] didn't work; result = " << result << "." << std::endl;

    // ELEMENTS
    // elem_connectivity
  int *connect2 = NULL;
  int size_connect2 = 0;
  bool status = elem_rms1.elem_connectivity(1, 2, &connect2, &size_connect2);
  if (status != true)
    std::cout << "elem_connectivity() RETURN VALUE failed." << std::endl;
  if (8 != size_connect2)
    std::cout << "re-sizing of connect2 vector failed." << std::endl;
  
  result = compare_connect(connect2, connect, 8);
  if (result != 0)
    std::cout << "elem_connectivity() VALUES failed." << std::endl;

    // set_elem_connectivity
    // reverse the connectivity
  for (i = 1; i <= 8; i++)
    connect2[i-1] = connect[8-i];
  elem_rms1.set_elem_connectivity(1, 2, connect2, size_connect2);
  status = elem_rms1.elem_connectivity(1, 2, &connect2, &size_connect2);
  if (status != true)
    std::cout << "set_elem_connectivity() RETURN VALUE failed." << std::endl;
  
  result = compare_connect(connect2, connect, 8);
  if (result != 0)
    std::cout << "elem_connectivity() VALUES failed." << std::endl;

    // RMESHSET FIND FUNCTIONS
    // find the rmeshsets for a node, and an element in each set
  MB_RMBSet *new_rms1, *new_rms2, *new_rms3;
  new_rms1 = MB_RMBSet::find_rmeshset(TSTT_VERTEX, TSTT_LAST_TOPOLOGY, 
                                         reinterpret_cast<const void*>(2));
  new_rms2 = MB_RMBSet::find_rmeshset(TSTT_REGION, TSTT_TETRAHEDRON, 
                                         reinterpret_cast<const void*>(2));
  new_rms3 = MB_RMBSet::find_rmeshset(TSTT_REGION, TSTT_TETRAHEDRON, 
                                         reinterpret_cast<const void*>(4));
  if (new_rms1 != &node_rms || new_rms2 != &elem_rms1 || new_rms3 != &elem_rms2)
    std::cout << "find_rmeshset() function failed." << std::endl;

    // now test NULL returns
  new_rms1 = MB_RMBSet::find_rmeshset(TSTT_VERTEX, TSTT_LAST_TOPOLOGY, 
                                         reinterpret_cast<const void*>(10));
  new_rms2 = MB_RMBSet::find_rmeshset(TSTT_REGION, TSTT_TETRAHEDRON, 
                                         reinterpret_cast<const void*>(0));
  new_rms3 = MB_RMBSet::find_rmeshset(TSTT_REGION, TSTT_TETRAHEDRON, 
                                         reinterpret_cast<const void*>(5));
  if (NULL != new_rms1 || NULL != new_rms2 || NULL != new_rms3)
    std::cout << "find_rmeshset() for NULL RETURN failed." << std::endl;
  
    // test is_in_rmeshset
  bool result1, result2, result3;
  result1 = node_rms.is_in_rmeshset(reinterpret_cast<const void*>(6));
  result2 = elem_rms1.is_in_rmeshset(reinterpret_cast<const void*>(2));
  result3 = elem_rms2.is_in_rmeshset(reinterpret_cast<const void*>(4));
  if (false == result1 || false == result2 || false == result3)
    std::cout << "is_in_rmeshset() failed." << std::endl;
  
    // test is_in_rmeshset
  result1 = node_rms.is_in_rmeshset(reinterpret_cast<const void*>(10));
  result2 = elem_rms1.is_in_rmeshset(reinterpret_cast<const void*>(4));
  result3 = elem_rms2.is_in_rmeshset(reinterpret_cast<const void*>(2));
  if (true == result1 || true == result2 || true == result3)
    std::cout << "is_in_rmeshset() for NULL RETURN failed." << std::endl;

  return 1;
}


int compare_coords(double *xval, double *yval, double *zval, 
                   double *nodex, double *nodey, double *nodez, 
                   const int num_nodes) 
{
  int i, result = 0;
  for (i = 0; i < num_nodes; i++) {
    if (xval[i] != nodex[i] ||
        yval[i] != nodey[i] ||
        zval[i] != nodez[i]) result++;
    xval[i] = yval[i] = zval[i] = -2.0;
  }
  return result;
}

int compare_connect(int *connect1, int *connect2, const int num_comps) 
{
  int i, result = 0;
  for (i = 0; i < num_comps; i++) {
    if (connect1[i] != connect2[i]) result++;
    connect2[i] = -1;
  }

  return result;
}
