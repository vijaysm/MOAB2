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


  /*!
    MBTest.cpp
    Test Harness for MB mesh database system 
  */

#ifdef WIN32
#ifdef _DEBUG
  // turn off warnings that say they debugging identifier has been truncated
  // this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "MBInterface.hpp"
#include "MBTagConventions.hpp"
#include "MBRange.hpp"
#include "MBSkinner.hpp"
#include "MeshTopoUtil.hpp"
#include "MBCN.hpp"
#include "MBOrientedBox.hpp"
#include "MBCartVect.hpp"

#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "MBInternals.hpp"
#include "MBCore.hpp"
#include "SequenceManager.hpp"
#include "EntitySequence.hpp"
#include "MBRangeSeqIntersectIter.hpp"

using namespace std;

#include "testdir.h"
string TestDir( TEST_DIR );


  /*!
    @test 
    Vertex Coordinates
    @li Get coordinates of vertex 1 correctly
    @li Get coordinates of vertex 8 correctly
    @li Get coordinates of vertex 6 correctly
  */

MBErrorCode mb_vertex_coordinate_test(MBInterface *MB)
{
  double coords[3];
  MBEntityHandle handle;
  MBErrorCode error;
  int err;

    // coordinate 2 should be {1.5, -1.5, 3.5}

  handle = CREATE_HANDLE(MBVERTEX, 2, err);
  error = MB->get_coords(&handle, 1, coords );
  if (error != MB_SUCCESS)
    return error;

  if(coords[0] != 1.5 || coords[1] != -1.5 || coords[2] != 3.5)
    return MB_FAILURE;

  double xyz[3];
  error = MB->get_coords(&handle, 1, xyz); 
  if (error != MB_SUCCESS)
    return error;

  if(xyz[0] != 1.5 || xyz[1] != -1.5 || xyz[2] != 3.5)
    return MB_FAILURE;


    // coordinate 9 should be {1, -2, 3.5}
  handle = CREATE_HANDLE(MBVERTEX, 9, err);
  error = MB->get_coords(&handle, 1, coords );
  if (error != MB_SUCCESS)
    return error;

  if(coords[0] != 1.0 || coords[1] != -2.0 || coords[2] != 3.5)
    return MB_FAILURE;

    // coordinate 7 should be {0.5, -2, 3.5}
  handle = CREATE_HANDLE(MBVERTEX, 7, err);
  error = MB->get_coords(&handle, 1, coords);
  if (error != MB_SUCCESS)
    return error;

  if(coords[0] != 0.5 || coords[1] != -2.0 || coords[2] != 3.5)
    return MB_FAILURE;

    // Walk the entire list.  Stop when the error code indicates
    // failure (testing walking off the end of the list).  Count
    // the items and ensure a proper count.


    // tag names must be defined by convention

  MBRange vertices;
  error = MB->get_entities_by_type(0,  MBVERTEX, vertices);

  if (error != MB_SUCCESS)
    return error;

  int node_count = 0;
  double all_coords[3*47];
  double* coord_iter = all_coords;
  for ( MBRange::iterator iter = vertices.begin();
        iter != vertices.end(); ++iter)
  {
    error = MB->get_coords(&(*iter), 1, coord_iter );
    coord_iter += 3;
    if (error == MB_SUCCESS)
      node_count++;
  }
    // Number of vertices (node_count) should be 83 assuming no gaps in the handle space
  if ( node_count != 47 )
    return MB_FAILURE;
    
    
    // check blocked coordinates
  double x[48], y[48], z[48];
  error = MB->get_coords( vertices, x, y, z );
  int num_inequal = 0;
  for (int i = 0; i < 47; ++i) {
    if (x[i] != all_coords[3*i  ])
      ++num_inequal;
    if (y[i] != all_coords[3*i+1])
      ++num_inequal;
    if (z[i] != all_coords[3*i+2])
      ++num_inequal;
  }
  if (num_inequal)
    return MB_FAILURE;
    
    // add invalid handle to end of range and try query again
  vertices.insert( vertices.back() + 1 );
  error = MB->get_coords( vertices, x, y, z );
  if (MB_ENTITY_NOT_FOUND != error)
    return MB_FAILURE;  

    // Try getting coordinates for a hex (should fail)
  handle = CREATE_HANDLE(MBHEX, 0, err);
  error = MB->get_coords(&handle, 1, coords);
  if (error == MB_SUCCESS)
    return MB_FAILURE;

  return MB_SUCCESS;
}

  /*!
    @test
    MB Vertex Tag Test
    @li Add, Set and correctly get an int tag
    @li Add, Set and correctly get a boolean tag
    @li Add, Set and correctly get a double tag
    @li Add, Set and correctly get a struct tag
  */

MBErrorCode mb_vertex_tag_test(MBInterface *MB)
{
    // Add an int Vertex Tag to the database

  int tag_size = sizeof(int);
  MBTag tag_id;

    // Create a dense tag for all vertices
  MBErrorCode error = MB->tag_create("int_tag", tag_size, MB_TAG_SPARSE, tag_id, 0);
  if (error != MB_SUCCESS)
    return error;

    // put a value in vertex 1 and retrieve

  int err;
  MBEntityHandle handle = CREATE_HANDLE(MBVERTEX, 1, err);
  int input_value = 11;
  error = MB->tag_set_data(tag_id, &handle, 1, &input_value);
  if (MB_SUCCESS != error) return error;

  int output_value;
  error = MB->tag_get_data(tag_id, &handle, 1, &output_value);
  if(MB_SUCCESS != error) return error;
  else if (output_value != input_value)
    return MB_FAILURE;

    // put a value in vertex 5 and retrieve

  handle = CREATE_HANDLE(MBVERTEX, 5, err);
  input_value = 11;
  error = MB->tag_set_data(tag_id, &handle, 1, &input_value);
  if (MB_SUCCESS != error) return error;
  error = MB->tag_get_data(tag_id, &handle, 1, &output_value);
  if (MB_SUCCESS != error) return error;
  else if(output_value != input_value)
    return MB_FAILURE;

    // put a value in vertex 98088234 which doesn't exist

  handle = CREATE_HANDLE(MBVERTEX, 98088234, err);
  input_value = 11;

  error = MB->tag_set_data(tag_id, &handle, 1, &input_value);
  if (error == MB_SUCCESS)
    return error;

  error = MB->tag_get_data(tag_id, &handle, 1, &output_value);
  if (error == MB_SUCCESS)
    return error;

  if(output_value != input_value)
    return MB_FAILURE;

    // Add a bool Vertex Tag to the database

  tag_size = sizeof(bool);
  error = MB->tag_create("bool_tag", tag_size, MB_TAG_SPARSE, tag_id, 0);
  if (error != MB_SUCCESS)
    return error;

    // put a value in vertex 5 and retrieve

  handle = CREATE_HANDLE(MBVERTEX, 5, err);
  bool bool_input_value = true;
  bool bool_output_value = false;
  error = MB->tag_set_data(tag_id, &handle, 1, &bool_input_value);
  if (error != MB_SUCCESS) return error;
  error = MB->tag_get_data(tag_id, &handle, 1, &bool_output_value);
  if (error != MB_SUCCESS) return error;
  else if(bool_output_value != bool_input_value)
    return MB_FAILURE;

    // Add a double Vertex Tag to the database

  tag_size = sizeof(double);
  error = MB->tag_create("double_tag", tag_size, MB_TAG_SPARSE, tag_id, 0);
  if (error != MB_SUCCESS)
    return error;

    // put a value in vertex 8: and retrieve

  handle = CREATE_HANDLE(MBVERTEX, 8, err);
  double double_input_value = true;
  double double_output_value = false;
  error = MB->tag_set_data(tag_id, &handle, 1, &double_input_value);
  if (error != MB_SUCCESS) return error;

  error = MB->tag_get_data(tag_id, &handle, 1, &double_output_value);
  if (error != MB_SUCCESS) return error;
  else if(double_output_value != double_input_value)
    return MB_FAILURE;

    // Add a struct Vertex Tag to the database

  struct TagStruct {
    int test_int;
    double test_double;
  };
  tag_size = sizeof(TagStruct);
  error = MB->tag_create("struct_tag", tag_size, MB_TAG_SPARSE, tag_id, 0);
  if (error != MB_SUCCESS)
    return error;

    // put a value in vertex 7 and retrieve

  handle = CREATE_HANDLE(MBVERTEX, 7, err);
  TagStruct input_tag_struct;
  input_tag_struct.test_int = 55;
  input_tag_struct.test_double = -1.2345;
  TagStruct output_tag_struct;
  error = MB->tag_set_data(tag_id, &handle, 1, &input_tag_struct);
  if (error != MB_SUCCESS) return error;
  error = MB->tag_get_data(tag_id, &handle, 1, &output_tag_struct);
  if (error != MB_SUCCESS) return error;
  else if(output_tag_struct.test_int != input_tag_struct.test_int ||
          output_tag_struct.test_double != input_tag_struct.test_double)
    return MB_FAILURE;


    // Create sparse tags for 10 random entities including some outside the 
    // range of allowable entities.

  unsigned int node_ids[] = {1, 933, 5, 9327, 8, 13400, 11, 8344, 1, 0x00FFFFFF}; 

  error = MB->tag_create("sparse_int_tag", tag_size, MB_TAG_SPARSE, tag_id, 0);

  if (error != MB_SUCCESS )  
    return error;

    //print_yes = true;
  int i;
  for (i=0; i<10; i++)
  {
    int err=0;
    MBEntityHandle handle = CREATE_HANDLE(MBVERTEX, node_ids[i], err);

    if (err != 0)
      return MB_FAILURE;

    int input_value = 11;
    error = MB->tag_set_data(tag_id, &handle, 1, &input_value);

      // should fail on odd values of i
    if ( !(i % 2) )
    {
      if (error != MB_SUCCESS )  // even case and if failed
        return error;
    }
    else 
    {
      if ( error == MB_SUCCESS)  // odd case and it says it worked!
        return MB_FAILURE;
    }

    int output_value;
    error = MB->tag_get_data(tag_id, &handle, 1, &output_value);
    if ( (i % 2) && error != MB_FAILURE && error != MB_TAG_NOT_FOUND)
      return error;

    if( (i % 2) && output_value != input_value)
      return MB_FAILURE;
  }

    // get the tag_name of the last tag created above
  std::string int_tag_name;
  error = MB->tag_get_name (tag_id, int_tag_name);
  if (error != MB_SUCCESS)
    return error;

  if (int_tag_name != "sparse_int_tag")
    return MB_FAILURE;

    // get the tag handle of the last tag created above
  MBTag int_tag_handle;
  error = MB->tag_get_handle (int_tag_name.c_str(), int_tag_handle);
  if (MB_SUCCESS != error) return error;
    
  if (int_tag_handle != tag_id)
    return MB_FAILURE;

    // test tag_get_tags_on_entity and tag_delete_data
  std::vector<MBTag> all_tags;
  handle = CREATE_HANDLE(MBVERTEX, node_ids[0], err);
  error = MB->tag_get_tags_on_entity(handle, all_tags);
  if (MB_SUCCESS != error)
    return error;

  if (!all_tags.empty()) {
    error = MB->tag_delete_data(all_tags[0], &handle, 1);
    if (MB_SUCCESS != error)
      return error;

    error = MB->tag_delete(all_tags[0]);
    if (MB_SUCCESS != error)
      return error;
  }
    // delete tags test

    // delete 2 of the sparse tags that were created above.
  handle = CREATE_HANDLE(MBVERTEX, node_ids[2], err);
  error = MB->tag_delete_data(tag_id, &handle, 1);
  if (error != MB_SUCCESS )  
    return error;

  handle = CREATE_HANDLE(MBVERTEX, node_ids[6], err);
  error = MB->tag_delete_data(tag_id, &handle, 1);
  if (error != MB_SUCCESS )  
    return error;

    // delete all the rest of the sparse tags.

  error = MB->tag_delete(tag_id);
  if (error != MB_SUCCESS )  
    return error;

    // delete the dense tag named bool_tag 
  MBTag bool_tag_handle;
  error = MB->tag_get_handle ("bool_tag", bool_tag_handle);
  if (error != MB_SUCCESS) return error;

  error = MB->tag_delete(bool_tag_handle);
  if (error != MB_SUCCESS )  
    return error;

  return error;
}


  /*!
    @test
    MB Bar Element Connectivity Test
    @li Get coordinates for 2 node bar elements
  */

MBErrorCode mb_bar_connectivity_test(MBInterface *MB)
{

  std::vector<MBEntityHandle> conn;
  MBErrorCode error;

  MBRange bars;

  error = MB->get_entities_by_type(0, MBEDGE, bars);

  if (error != MB_SUCCESS)
    return error;

    // get the connectivity of the second bar
  MBEntityHandle handle = *(++bars.begin());

  error = MB->get_connectivity(&handle, 1, conn);
  if (error != MB_SUCCESS )  
    return error;

  if (conn.size() != 2)
    return MB_FAILURE;

    // from ncdump the connectivity of bar 2 (0 based) is
    //  14, 13 

  if ( conn[0] != 14)
    return MB_FAILURE;

  if ( conn[1] != 13)  
    return MB_FAILURE;

    // Now try getting the connectivity of one of the vertices for fun.
    // just return the vertex in the connectivity
  handle = conn[0];
  error = MB->get_connectivity(&handle, 1, conn);
  if (error != MB_SUCCESS && handle != conn[0] && conn.size() != 1)  
    return error;

  return MB_SUCCESS;
}

MBErrorCode mb_tri_connectivity_test(MBInterface *MB)
{

  std::vector<MBEntityHandle> conn; 
  MBErrorCode error;

  MBRange tris;
  error = MB->get_entities_by_type(0, MBTRI, tris);

  if (error != MB_SUCCESS)
    return error;

    // get the connectivity of the second tri
  MBEntityHandle handle = *(++tris.begin());

  error = MB->get_connectivity(&handle, 1, conn);
  if (error != MB_SUCCESS )  
    return error;

  if (conn.size() != 3)
    return MB_FAILURE;

    // from ncdump the connectivity of tri 2 (0 based) is
    //  45, 37, 38

  if (conn[0] != 45)
    return MB_FAILURE;

  if (conn[1] != 37)  
    return MB_FAILURE;

  if (conn[2] != 38) 
    return MB_FAILURE;

  return MB_SUCCESS;
}

MBErrorCode mb_quad_connectivity_test(MBInterface *MB)
{

  std::vector<MBEntityHandle> conn;

  MBRange quads;

  MBErrorCode error = MB->get_entities_by_type(0, MBQUAD, quads);

  if (error != MB_SUCCESS)
    return error;

    // get the connectivity of the second quad
  MBEntityHandle handle = *(++quads.begin());

  error = MB->get_connectivity(&handle, 1, conn);
  if (error != MB_SUCCESS )  
    return error;

  if (conn.size() != 4)
    return MB_FAILURE;

    // from ncdump the connectivity of quad 2 (0 based) is
    // 20, 11, 12, 26,

  if (conn[0] != 20)
    return MB_FAILURE;

  if (conn[1] != 11)  
    return MB_FAILURE;

  if (conn[2] != 12) 
    return MB_FAILURE;

  if (conn[3] != 26)
    return MB_FAILURE;

  return MB_SUCCESS;
}

MBErrorCode mb_hex_connectivity_test(MBInterface *MB)
{

  std::vector<MBEntityHandle> conn;

  MBRange hexes;

  MBErrorCode error = MB->get_entities_by_type(0,  MBHEX, hexes);

  if (error != MB_SUCCESS)
    return error;

    // get the connectivity of the second hex
  MBEntityHandle handle = *(++hexes.begin());

  error = MB->get_connectivity(&handle, 1, conn);
  if (error != MB_SUCCESS )  
    return error;

  if (conn.size() != 8)
    return MB_FAILURE;

    // from ncdump the connectivity of hex 1 (0 based) is
    //19, 13, 16, 23, 21, 14, 18, 27

  if (conn[0] != 19)
    return MB_FAILURE;

  if (conn[1] != 13)  
    return MB_FAILURE;

  if (conn[2] != 16) 
    return MB_FAILURE;

  if (conn[3] != 23)
    return MB_FAILURE;

  if (conn[4] != 21)
    return MB_FAILURE;

  if (conn[5] != 14)  
    return MB_FAILURE;

  if (conn[6] != 18) 
    return MB_FAILURE;

  if (conn[7] != 27)
    return MB_FAILURE;

  return MB_SUCCESS;
}

MBErrorCode mb_tet_connectivity_test(MBInterface *MB)
{
  std::vector<MBEntityHandle> conn; 

  MBRange tets;

  MBErrorCode error = MB->get_entities_by_type(0, MBTET, tets);

  if (error != MB_SUCCESS)
    return error;

    // get the connectivity of the second tet
  MBEntityHandle handle = *(++tets.begin());

  error = MB->get_connectivity(&handle, 1, conn);
  if (error != MB_SUCCESS )  
    return error;

  if (conn.size() != 4)
    return MB_FAILURE;

    // from ncdump the connectivity of tet 2 (0 based) is: 
    // 35, 34, 32, 43 

  if (conn[0] != 35)
    return MB_FAILURE;

  if (conn[1] != 34)  
    return MB_FAILURE;

  if (conn[2] != 32) 
    return MB_FAILURE;

  if (conn[3] != 43)
    return MB_FAILURE;

  return MB_SUCCESS;
}
MBErrorCode mb_temporary_test( MBInterface *gMB )
{

  double array[3] = {0.0, 0.0, 0.0};
  MBEntityHandle h_node1;
  MBErrorCode result = gMB->create_vertex(array, h_node1);
  if (MB_SUCCESS != result)
    return result;

  MBEntityHandle ordered_meshset1;
  result = gMB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, ordered_meshset1);
  if (MB_SUCCESS != result)
    return result;

  MBEntityHandle ordered_meshset2;
  result = gMB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, ordered_meshset2);
  if (MB_SUCCESS != result)
    return result;

  result = gMB->add_entities(ordered_meshset1, &h_node1, 1);
  if (MB_SUCCESS != result)
    return result;

  result = gMB->remove_entities(ordered_meshset1, &h_node1, 1);
  if (MB_SUCCESS != result)
    return result;

  result = gMB->add_entities(ordered_meshset2, &h_node1, 1);
  if (MB_SUCCESS != result)
    return result;

  bool create_if_missing = false;
  std::vector<MBEntityHandle> meshsets;
  result = gMB->get_adjacencies(&h_node1, 1, 4, create_if_missing, meshsets);
  if (MB_SUCCESS != result)
    return result;

  int num_adj = meshsets.size();
  assert(1 == num_adj);

  return MB_SUCCESS;
}

MBErrorCode mb_adjacent_vertex_test( MBInterface* mb )
{
  MBErrorCode rval;
  MBRange hexes, expt_vert, got_vert, some_hexes;
  MBRange::const_iterator i, j;
  int n;
  
    // get all hexes
  rval = mb->get_entities_by_type( 0, MBHEX, hexes );
  if (rval != MB_SUCCESS)
    return rval;
  if (hexes.empty())  // can't do test if no elements
    return MB_FAILURE;
  
    // get every third hex and its vertices
  n = 0;
  for (i = hexes.begin(); i != hexes.end(); ++i) {
    if (++n % 3)
      continue;
    some_hexes.insert( *i );
    const MBEntityHandle* conn;
    int len;
    rval = mb->get_connectivity( *i, conn, len );
    if (MB_SUCCESS != rval)
      return rval;
    for (int k = 0; k < len; ++k)
      expt_vert.insert( conn[k] );
  }
  
    // use get_adjacencies to get vertices
  rval = mb->get_adjacencies( some_hexes, 0, false, got_vert, MBInterface::UNION );
  if (MB_SUCCESS != rval) {
    std::cout << "get_adjacencies failed with error code " << rval << std::endl;
    return rval;
  }
  
  i = expt_vert.begin();
  j = got_vert.begin();
  while (i != expt_vert.end() && j != got_vert.end()) {
    if (*i < *j) {
      std::cout << "Result missing vertex: " << *i << std::endl;
      return MB_FAILURE;
    }
    else if (*j < *i) {
      std::cout << "Result contains extra vertex: " << *j << std::endl;
      return MB_FAILURE;
    }
    ++i;
    ++j;
  }
  
  if (i != expt_vert.end()) {
    std::cout << "Result missing vertex: " << *i << std::endl;
    return MB_FAILURE;
  }
  else if (j != got_vert.end()) {
    std::cout << "Result contains extra vertex: " << *j << std::endl;
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}
  
MBErrorCode mb_adjacencies_test(MBInterface *mb) 
{
    // this test does the following:
    // 1. For each element, creates vertex-element adjacencies (only for
    //    lowest-id vertex)
    // 2. Creates all lower-order ancillary entities
    // 3. Checks for proper number of same
    //
    // assume mesh has already been read

  MBEntityType seq_type;
  MBErrorCode result = MB_SUCCESS;
  MBRange handle_range;

    // Some code to make the test simpler locally (this will cause other tests to
    // fail though, so watch out!)
    //mb->delete_mesh();
    //result = mb->load_mesh("mbtest2.g", 0);
    //if (MB_SUCCESS != result)
  return result;

    // lets create a skin of the hexes
    // this may be far from being the most efficient, but
    // it certainly does exercise the adjacency routines

  MBRange::iterator iter;
  MBRange::reverse_iterator riter;

    // first get the hexes
  MBRange hexes;
  result = mb->get_entities_by_type(0, MBHEX, hexes);
  if (MB_SUCCESS != result)
    return result;


  unsigned int num_hexes = hexes.size();

    // make sure we got hexes
  for(riter = hexes.rbegin(); riter != hexes.rend(); ++riter)
  {
    if( TYPE_FROM_HANDLE(*riter) != MBHEX)
      return MB_FAILURE;
  }



    // get all the nodes that these hexes are connected to
  MBRange nodes;
  result = mb->get_adjacencies(hexes, 0, false, nodes, MBInterface::UNION);
  if (MB_SUCCESS != result)
    return result;

    // make sure we got nodes
  for(iter = nodes.begin(); iter != nodes.end(); ++iter)
  {
    if( TYPE_FROM_HANDLE(*iter) != MBVERTEX)
      return MB_FAILURE;
  }

    // find the interior nodes; assume a structured mesh
  MBRange interior_nodes;
  std::vector<MBEntityHandle> attached_hexes;
  for( iter = nodes.begin(); iter != nodes.end();)
  {
    attached_hexes.clear();
    result = mb->get_adjacencies(&(*iter), 1, 3, false, attached_hexes);
    if (MB_SUCCESS != result) 
      return result;
    attached_hexes.erase(std::remove_if(attached_hexes.begin(), 
                                        attached_hexes.end(), 
                                        type_not_equals(mb, MBHEX)), 
                         attached_hexes.end());



      // if this node has less than 8 hexes attached to it, it is not
      // an interior node
    if (attached_hexes.size() == 8)
    {
        // add to the interior nodes list and remove from the nodes list
      interior_nodes.insert(*iter);
      iter = nodes.erase(iter);
    }
    else
      iter++;

  }

    // get interior quads from interior nodes
  MBRange interior_quads;
  result = mb->get_adjacencies(interior_nodes, 2, true, interior_quads, MBInterface::UNION);
  if (MB_SUCCESS != result) 
    return result;

    // get a list of quads generated adjacent to the exterior nodes
  MBRange temp_quads, exterior_quads;
  result = mb->get_adjacencies(nodes, 2, true, temp_quads, MBInterface::UNION);
  if (MB_SUCCESS != result) 
    return result;

    // now remove any interior quads from the previous quads list
    // and we should be left with exterior quads
  std::set_difference(temp_quads.begin(), temp_quads.end(),
                      interior_quads.begin(), interior_quads.end(),
                      mb_range_inserter(exterior_quads));

    // check to make sure we have the right number of quads; for hexes, should be
    // .5(6*num_hexes - num_exterior_quads)
  unsigned int num_expected_int_quads = (6*num_hexes - exterior_quads.size())/2;
  if (num_expected_int_quads != interior_quads.size())
    return MB_FAILURE;

    // delete the interior quads
  result = mb->delete_entities(interior_quads);
  if (MB_SUCCESS != result) 
    return result;

  MBRange remaining_quads;
  result = mb->get_entities_by_type(0, MBQUAD, remaining_quads);
  if (MB_SUCCESS != result)
    return result;

  if(remaining_quads.size() != exterior_quads.size())
    return MB_FAILURE;


    // 3. Checks for proper number of same
    // re-retrieve and store the handle ranges for all element types

  int num_ents;
    
  for (seq_type = MBEDGE; seq_type != MBENTITYSET; seq_type++) 
  {
    handle_range.clear();

    result = mb->get_entities_by_type(0, seq_type, handle_range);
    if (MB_SUCCESS != result)
      return result;

    result = mb->get_number_entities_by_type(0, seq_type, num_ents);
    if (MB_SUCCESS != result)
      return result;

    if(handle_range.size() != (unsigned int) num_ents)
      return MB_FAILURE;
  }

  return result;

}

static MBErrorCode create_two_hex_full_mesh( MBInterface* mb,
                                             MBEntityHandle vertices[12],
                                             MBEntityHandle hexes[2],
                                             MBEntityHandle hex1_faces[6],
                                             MBEntityHandle hex2_faces[6],
                                             MBEntityHandle hex1_edges[12],
                                             MBEntityHandle hex2_edges[12] )
{
  MBErrorCode rval;
 // create a simple mesh containing 2 hexes
  const double coords[] = { 0, 0, 0, 
                            1, 0, 0,
                            2, 0, 0,
                            0, 1, 0,
                            1, 1, 0,
                            2, 1, 0,
                            0, 0, 1, 
                            1, 0, 1,
                            2, 0, 1,
                            0, 1, 1,
                            1, 1, 1,
                            2, 1, 1 };
  for (int i = 0; i < 12; ++i)
    if (MB_SUCCESS != mb->create_vertex( coords + 3*i, vertices[i] ))
      return MB_FAILURE;
  MBEntityHandle hex1_conn[] = { vertices[6], vertices[7], vertices[1], vertices[0],
                                 vertices[9], vertices[10],vertices[4], vertices[3] };
  MBEntityHandle hex2_conn[] = { vertices[7], vertices[8], vertices[2], vertices[1],
                                 vertices[10],vertices[11],vertices[5], vertices[4] };
  MBEntityHandle shared_quad_conn[] = { vertices[7], vertices[1], vertices[4], vertices[10] };
  MBEntityHandle hex1_face_conn[][4] = {
    { vertices[6], vertices[7], vertices[10],vertices[9] },
    { vertices[7], vertices[6], vertices[0], vertices[1] },
    { vertices[1], vertices[0], vertices[3], vertices[4] },
    { vertices[9], vertices[10],vertices[4], vertices[3] },
    { vertices[3], vertices[0], vertices[6], vertices[9] } };
  MBEntityHandle hex2_face_conn[][4] = {
    { vertices[7], vertices[8], vertices[11],vertices[10] },
    { vertices[8], vertices[7], vertices[1], vertices[2] },
    { vertices[2], vertices[1], vertices[4], vertices[5] },
    { vertices[10],vertices[11],vertices[5], vertices[4] },
    { vertices[5], vertices[2], vertices[8], vertices[11] } };
  MBEntityHandle shared_edge_conn[][2] = { { vertices[1], vertices[7] },
                                           { vertices[7], vertices[10]},
                                           { vertices[10],vertices[4] },
                                           { vertices[4], vertices[1] } };
  MBEntityHandle hex1_edge_conn[][2] = { { vertices[6], vertices[7] },
                                         { vertices[9], vertices[10] },
                                         { vertices[3], vertices[4] },
                                         { vertices[0], vertices[1] },
                                         { vertices[6], vertices[9] },
                                         { vertices[9], vertices[3] },
                                         { vertices[3], vertices[0] },
                                         { vertices[0], vertices[6] } };
  MBEntityHandle hex2_edge_conn[][2] = { { vertices[7], vertices[8] },
                                         { vertices[10], vertices[11] },
                                         { vertices[4], vertices[5] },
                                         { vertices[1], vertices[2] },
                                         { vertices[8], vertices[11] },
                                         { vertices[11], vertices[5] },
                                         { vertices[5], vertices[2] },
                                         { vertices[2], vertices[8] } };
  rval = mb->create_element( MBHEX, hex1_conn, 8, hexes[0] );
  if (MB_SUCCESS != rval)
    return rval;
  rval = mb->create_element( MBHEX, hex2_conn, 8, hexes[1] );
  if (MB_SUCCESS != rval)
    return rval;
  rval = mb->create_element( MBQUAD, shared_quad_conn, 4, hex1_faces[0] );
  if (MB_SUCCESS != rval)
    return rval;
  hex2_faces[0] = hex1_faces[0];
  for (int i = 0; i < 5; ++i) {
    rval = mb->create_element( MBQUAD, hex1_face_conn[i], 5, hex1_faces[i+1] );
    if (MB_SUCCESS != rval)
      return rval;
    rval = mb->create_element( MBQUAD, hex2_face_conn[i], 5, hex2_faces[i+1] );
    if (MB_SUCCESS != rval)
      return rval;
  }
  for (int i = 0; i < 4; ++i) {
    rval = mb->create_element( MBEDGE, shared_edge_conn[i], 2, hex1_edges[i] );
    if (MB_SUCCESS != rval)
      return rval;
    hex2_edges[i] = hex1_edges[i];
  }
  for (int i = 0; i < 8; ++i) {
    rval = mb->create_element( MBEDGE, hex1_edge_conn[i], 2, hex1_edges[i+4] );
    if (MB_SUCCESS != rval)
      return rval;
    rval = mb->create_element( MBEDGE, hex2_edge_conn[i], 2, hex2_edges[i+4] );
    if (MB_SUCCESS != rval)
      return rval;
  }
  return MB_SUCCESS;
}


MBErrorCode mb_upward_adjacencies_test(MBInterface *mb) 
{
  MBErrorCode rval;
  MBCore moab;
  mb = &moab;
  
  // create a simple mesh containing 2 hexes
  MBEntityHandle vertices[12], hexes[2], hex1_faces[6], hex2_faces[6], hex1_edges[12], hex2_edges[12];
  rval = create_two_hex_full_mesh( mb, vertices, hexes, hex1_faces, hex2_faces, hex1_edges, hex2_edges );
  if (MB_SUCCESS != rval)
    return rval;

    // test adjacences from dim to 3
  for (int dim = 0; dim < 3; ++dim) {
    std::vector<MBEntityHandle> hex1_ent, hex2_ent, shared;
    const MBEntityHandle *list1, *list2;
    int n;
    switch (dim) {
      case 0:
        rval = mb->get_connectivity( hexes[0], list1, n );
        if (MB_SUCCESS != rval)
          return rval;
        rval = mb->get_connectivity( hexes[1], list2, n );
        if (MB_SUCCESS != rval)
          return rval;
        break;
      case 1:
        list1 = hex1_edges;
        list2 = hex2_edges;
        n = 12;
        break;
      case 2:
        list1 = hex1_faces;
        list2 = hex2_faces;
        n = 6;
        break;
    }
    for (int i = 0; i < n; ++i) {
      if (std::find(list2, list2+n, list1[i]) - list2 == n)
        hex1_ent.push_back(list1[i]);
      else
        shared.push_back(list1[i]);
      if (std::find(list1, list1+n, list2[i]) - list1 == n)
        hex2_ent.push_back(list2[i]);
    }
    
    for (size_t j = 0; j < shared.size(); ++j) {
      std::vector<MBEntityHandle> adj;
      rval = mb->get_adjacencies( &shared[j], 1, 3, false, adj );
      if (MB_SUCCESS != rval)
        return rval;
      if (adj.size() != 2)
        return MB_FAILURE;
      if (!(adj[0] == hexes[0] && adj[1] == hexes[1]) &&
          !(adj[0] == hexes[1] && adj[1] == hexes[0]))
        return MB_FAILURE;
    }
    
    for (size_t j = 0; j < hex1_ent.size(); ++j) {
      std::vector<MBEntityHandle> adj;
      rval = mb->get_adjacencies( &hex1_ent[j], 1, 3, false, adj );
      if (MB_SUCCESS != rval)
        return rval;
      if (adj.size() != 1 || adj[0] != hexes[0])
        return MB_FAILURE;
    }
    
    for (size_t j = 0; j < hex2_ent.size(); ++j) {
      std::vector<MBEntityHandle> adj;
      rval = mb->get_adjacencies( &hex2_ent[j], 1, 3, false, adj );
      if (MB_SUCCESS != rval)
        return rval;
      if (adj.size() != 1 || adj[0] != hexes[1])
        return MB_FAILURE;
    }
  }
    
    // For each edge, get adjacent faces, and for each face
    // get adjacent hexes.  Result should be the same as
    // direct query from edges to hexes
  std::vector<MBEntityHandle> all_edges(24);
  std::copy( hex1_edges, hex1_edges+12, all_edges.begin() );
  std::copy( hex2_edges, hex2_edges+12, all_edges.begin()+12 );
  std::sort( all_edges.begin(), all_edges.end() );
  all_edges.erase( std::unique(all_edges.begin(), all_edges.end()), all_edges.end() );
  for (size_t j = 0; j < all_edges.size(); ++j) {
    std::vector<MBEntityHandle> edge_hexes, edge_faces, face_hexes;
    rval = mb->get_adjacencies( &all_edges[j], 1, 3, false, edge_hexes );
    if (MB_SUCCESS != rval)
      return rval;
    rval = mb->get_adjacencies( &all_edges[j], 1, 2, false, edge_faces );
    if (MB_SUCCESS != rval)
      return rval;
    rval = mb->get_adjacencies( &edge_faces[0], edge_faces.size(), 3,
                                false, face_hexes, MBInterface::UNION );
    if (MB_SUCCESS != rval)
      return rval;
    if (edge_hexes.size() != face_hexes.size())
      return MB_FAILURE;
    switch (edge_hexes.size()) {
      case 1:
        if (edge_hexes[0] != face_hexes[0])
          return MB_FAILURE;
        break;
      case 2:
        if (!(edge_hexes[0] == face_hexes[0] && edge_hexes[1] == face_hexes[1]) &&
            !(edge_hexes[0] == face_hexes[1] && edge_hexes[1] == face_hexes[0]))
        return MB_FAILURE;
        break;
      default:
        return MB_FAILURE;
    }
  }
  return MB_SUCCESS;
}

MBErrorCode nothing_but_type( MBRange& range, MBEntityType type )
{ 

    //make sure there's nothing but hexes in hex_ms
  MBRange::iterator iter, end_iter;
  iter = range.begin();
  end_iter = range.end();

  for(; iter != end_iter; iter++)
  {
    if( TYPE_FROM_HANDLE(*iter) != type )
    {
      return MB_FAILURE; 
    }
  }
  return MB_SUCCESS;
}

MBErrorCode check_esets(MBInterface * MB, const int num_sets) 
{
  int entity_sets_size;
  MBErrorCode result = MB->get_number_entities_by_type(0, MBENTITYSET, entity_sets_size);
  if (MB_SUCCESS != result ||
      entity_sets_size != num_sets) return MB_FAILURE;
  
  return MB_SUCCESS;
}

MBErrorCode mb_mesh_sets_test(MBInterface * MB, int flags)
{

  MBRange temp_range;
  std::vector<MBEntityHandle> temp_vector;
  MBEntityType ent_type;

  MBEntityHandle ms_array[MBENTITYSET] = {0};
  unsigned int number_array[MBENTITYSET] = {0};
  unsigned int num_dim_array[4] = { 0, 0, 0, 0 };
  int count, start_num_sets;

  MBErrorCode result = MB->get_number_entities_by_type(0, MBENTITYSET, start_num_sets);
  if (MB_SUCCESS != result) return result;

    //add entities to meshsets 
  for (ent_type = MBEDGE; ent_type != MBENTITYSET; ent_type++) 
  {
    result = MB->create_meshset( flags, ms_array[ent_type] );
    if( result != MB_SUCCESS )
      return result;

    temp_range.clear();
    result = MB->get_entities_by_type(0, ent_type, temp_range );
    if( result != MB_SUCCESS ) 
      return result;
    result = MB->get_number_entities_by_type(0, ent_type, count);
    if (result != MB_SUCCESS)
      return result;
    if ((unsigned)count != temp_range.size())
      return MB_FAILURE;
    result = MB->add_entities( ms_array[ent_type], temp_range);
    if( result != MB_SUCCESS )
      return result;

    number_array[ent_type] = temp_range.size(); //KGM
    num_dim_array[MBCN::Dimension(ent_type)] += count;

      //Check to make sure mesh set really has correct number of entities in it
    temp_range.clear();
    temp_vector.clear();
    result = MB->get_entities_by_handle(ms_array[ent_type], temp_range);
    if(result != MB_SUCCESS)
      return result;
    if(number_array[ent_type] != temp_range.size())
    {
      cout<<"Number of entities in meshset test is not correct"<<endl;
      return MB_FAILURE; 
    }
    if (!temp_range.all_of_type( ent_type ))
      return MB_FAILURE;

    result = MB->get_entities_by_handle(ms_array[ent_type], temp_vector);
    if(result != MB_SUCCESS)
      return result;
    if(number_array[ent_type] != temp_vector.size())
    {
      cout<<"Number of entities in meshset test is not correct"<<endl;
      return MB_FAILURE; 
    }
      
    temp_range.clear();
    result = MB->get_entities_by_type( ms_array[ent_type], ent_type, temp_range);
    if(result != MB_SUCCESS)
      return result;
    if(number_array[ent_type] != temp_range.size())
    {
      cout<<"Number of entities by type in meshset test is not correct"<<endl;
      return MB_FAILURE; 
    }
    if (!temp_range.all_of_type( ent_type ))
      return MB_FAILURE;
      
    temp_range.clear();
    result = MB->get_entities_by_type( ms_array[ent_type], MBVERTEX, temp_range);
    if(result != MB_SUCCESS)
      return result;
    if(0 != temp_range.size())
      return MB_FAILURE;
      
    temp_range.clear();
    result = MB->get_entities_by_dimension( ms_array[ent_type], MBCN::Dimension(ent_type), temp_range);
    if(result != MB_SUCCESS)
      return result;
    if(number_array[ent_type] != temp_range.size())
    {
      cout<<"Number of entities by dimension in meshset test is not correct"<<endl;
      return MB_FAILURE; 
    }
    if (!temp_range.all_of_type( ent_type ))
      return MB_FAILURE;
      
    temp_range.clear();
    result = MB->get_entities_by_dimension( ms_array[ent_type], 0, temp_range);
    if(result != MB_SUCCESS)
      return result;
    if(0 != temp_range.size())
    {
      cout<<"Number of entities by dimension in meshset test is not correct"<<endl;
      return MB_FAILURE; 
    }
      
    result = MB->get_number_entities_by_handle( ms_array[ent_type], count );
    if (result != MB_SUCCESS)
      return result;
    if ((unsigned)count != number_array[ent_type])
      return MB_FAILURE;
      
    result = MB->get_number_entities_by_type( ms_array[ent_type], ent_type, count );
    if (result != MB_SUCCESS)
      return result;
    if ((unsigned)count != number_array[ent_type])
      return MB_FAILURE;
      
    result = MB->get_number_entities_by_type( ms_array[ent_type], MBVERTEX, count );
    if (result != MB_SUCCESS)
      return result;
    if (count != 0)
      return MB_FAILURE;
        
    result = MB->get_number_entities_by_dimension( ms_array[ent_type], MBCN::Dimension(ent_type), count );
    if (result != MB_SUCCESS)
      return result;
    if ((unsigned)count != number_array[ent_type])
      return MB_FAILURE;
      
    result = MB->get_number_entities_by_dimension( ms_array[ent_type], 0, count );
    if (result != MB_SUCCESS)
      return result;
    if (count != 0)
      return MB_FAILURE;
  }

  result = check_esets(MB, start_num_sets + MBENTITYSET - MBEDGE);
  if (MB_SUCCESS != result) return result;
    
  for (int dim = 1; dim < 4; ++dim) {
    result = MB->get_number_entities_by_dimension( 0, dim, count );
    if (MB_SUCCESS != result)
      return MB_FAILURE;
    if ((unsigned)count != num_dim_array[dim])
      return MB_FAILURE;
  }

    //----------TEST RECURSIVE OPERATIONS----------------//
  MBEntityHandle recursive1, recursive2;
  result = MB->create_meshset( MESHSET_SET, recursive1 );
  if( result != MB_SUCCESS )
    return result;
  result = MB->create_meshset( 0, recursive2 );
  if( result != MB_SUCCESS )
    return result;
  unsigned num_sets = MBENTITYSET-MBEDGE;
  result = MB->add_entities( recursive2, ms_array+MBEDGE, num_sets );
  if (MB_SUCCESS != result)
    return result;
  result = MB->add_entities( recursive1, &recursive2, 1 );
  if (MB_SUCCESS != result)
    return result;

  temp_range.clear();
  result = MB->get_entities_by_type( recursive1, MBENTITYSET, temp_range );
  if (MB_SUCCESS != result)
    return result;
  if (temp_range.size() != 1 || *(temp_range.begin()) != recursive2)
    return MB_FAILURE;

  temp_range.clear();
  result = MB->get_entities_by_type( recursive2, MBENTITYSET, temp_range );
  if (MB_SUCCESS != result)
    return result;
  if (temp_range.size() != num_sets || !temp_range.all_of_type(MBENTITYSET))
    return MB_FAILURE;

  temp_range.clear();
  result = MB->get_entities_by_handle( recursive1, temp_range );
  if (MB_SUCCESS != result)
    return result;
  if (temp_range.size() != 1 || *(temp_range.begin()) != recursive2)
    return MB_FAILURE;

  temp_range.clear();
  result = MB->get_entities_by_handle( recursive2, temp_range );
  if (MB_SUCCESS != result)
    return result;
  if (temp_range.size() != num_sets || !temp_range.all_of_type(MBENTITYSET))
    return MB_FAILURE;
    
    
  unsigned total = 0;
  for (ent_type = MBEDGE; ent_type != MBENTITYSET; ent_type++) 
  {
    total += number_array[ent_type];
    temp_range.clear();
    result = MB->get_entities_by_type( recursive1, ent_type, temp_range, true );
    if(result != MB_SUCCESS)
      return result;
    if(number_array[ent_type] != temp_range.size())
    {
      cout<<"Recursive number of entities by type in meshset test is not correct"<<endl;
      return MB_FAILURE; 
    }
    if (!temp_range.all_of_type( ent_type ))
      return MB_FAILURE;
    result = MB->get_number_entities_by_type( recursive1, ent_type, count, true );
    if(result != MB_SUCCESS)
      return result;
    if(number_array[ent_type] != (unsigned)count)
    {
      cout<<"Recursive number of entities by type in meshset test is not correct"<<endl;
      return MB_FAILURE; 
    }
    if (!temp_range.all_of_type( ent_type ))
      return MB_FAILURE;
  }
  if (0 == total) {
    cout << "Invalid test input.  No entities!" << endl;
    return MB_FAILURE;
  }

  for (int dim = 1; dim < 4; ++dim) {
    temp_range.clear();
    result = MB->get_entities_by_dimension( recursive1, dim, temp_range, true );
    if (MB_SUCCESS != result)
      return MB_FAILURE;
    if (temp_range.size() != num_dim_array[dim])
      return MB_FAILURE;
    if (!temp_range.all_of_dimension(dim))
      return MB_FAILURE;
    result = MB->get_number_entities_by_dimension( recursive1, dim, count, true );
    if (MB_SUCCESS != result)
      return MB_FAILURE;
    if ((unsigned)count != num_dim_array[dim])
      return MB_FAILURE;
  }

  temp_range.clear();
  result = MB->get_entities_by_handle( recursive1, temp_range, true );
  if(result != MB_SUCCESS)
    return result;
  if(total != temp_range.size())
  {
    cout<<"Recursive number of entities in meshset test is not correct"<<endl;
    return MB_FAILURE; 
  }
    
    // try circular relation
  result = MB->add_entities( recursive2, &recursive1, 1 );
  if (MB_SUCCESS != result) {
    std::cout << "Failed to create circular set containment" << std::endl;
    return result;
  }
  temp_range.clear();
  result = MB->get_entities_by_handle( recursive1, temp_range, true );
  if(result != MB_SUCCESS)
    return result;
  if(total != temp_range.size())
  {
    cout<<"Recursive number of entities in meshset test is not correct"<<endl;
    return MB_FAILURE; 
  }
    
  result = check_esets(MB, start_num_sets + MBENTITYSET - MBEDGE + 2);
  if (MB_SUCCESS != result) return result;

    //----------TEST BOOLEAN OPERATIONS----------------//

  MBEntityHandle temp_ms1, temp_ms2; 
  result = MB->create_meshset(flags, temp_ms1);
  if(result  != MB_SUCCESS ) 
    return result;
  result = MB->create_meshset(flags, temp_ms2);
  if(result != MB_SUCCESS )
    return result;

    //Subtract
    //add all edges and hexes of ms_array[MBHEX] and ms_array[MBEDGE] to temp_ms1
    //get all Edge entities
  temp_range.clear();
  result = MB->get_entities_by_handle(ms_array[MBEDGE], temp_range );
  if(result != MB_SUCCESS )
    return result;

    //add Edges to ms1 
  result = MB->add_entities( temp_ms1, temp_range );
  if(result != MB_SUCCESS )
    return result;

  temp_range.clear();
  result = MB->get_entities_by_handle(ms_array[MBHEX], temp_range );
  if(result != MB_SUCCESS )
    return result;

    //add Hexes to ms1 
  result = MB->add_entities( temp_ms1, temp_range );
  if(result != MB_SUCCESS )
    return result;


    //subtract bars meshset out of hex meshset 
  result = MB->subtract_meshset( temp_ms1, ms_array[MBEDGE]);
  if(result != MB_SUCCESS )
    return result;

    //Perform the check
  temp_range.clear(); 
  result = MB->get_entities_by_handle(temp_ms1, temp_range );
  if(result != MB_SUCCESS )
    return result;

  if(number_array[MBHEX] != temp_range.size())
  {
    cout<<"MBset subtraction is bad"<<endl;
    return MB_FAILURE; 
  }
    //make sure there's nothing but hexes in hex_ms
  if( nothing_but_type( temp_range, MBHEX ) != MB_SUCCESS )
    return MB_FAILURE; 

  result = check_esets(MB, start_num_sets + MBENTITYSET - MBEDGE + 4);
  if (MB_SUCCESS != result) return result;

    //------------Intersect------------
    //
    //clean out the temp_ms1
  MB->clear_meshset(&temp_ms1, 1);

  temp_range.clear();
    //get all quad entities 
  temp_range.clear();
  result = MB->get_entities_by_handle(ms_array[MBQUAD], temp_range );
  if(result != MB_SUCCESS )
    return result;


    //add tets them to ms1 
  result = MB->add_entities( temp_ms1, temp_range ) ;
  if(result != MB_SUCCESS )
    return result;

    //get all tet entities 
  temp_range.clear();
  result = MB->get_entities_by_handle(ms_array[MBTET], temp_range );
  if(result != MB_SUCCESS )
    return result;


    //add tets them to ms1 
  result = MB->add_entities( temp_ms1, temp_range ) ;
  if(result != MB_SUCCESS )
    return result;

    //intersect temp_ms1 (which contains all quads and tets) with tet meshset 
    //temp_ms1 meshset is altered
  result = MB->intersect_meshset(temp_ms1, ms_array[MBTET]) ;
  if(result != MB_SUCCESS )
    return result;

    //Perform the check, only tets should be in temp_ms1 
  temp_range.clear(); 
  result = MB->get_entities_by_handle(temp_ms1, temp_range );
  if(result != MB_SUCCESS )
    return result;
  if(number_array[MBTET] != temp_range.size())
  {
    cout<<"MBset intersection is bad"<<endl;
    return MB_FAILURE; 
  }

    //make sure there's nothing but tet in range 
  if( nothing_but_type( temp_range, MBTET ) != MB_SUCCESS )
    return MB_FAILURE; 

  result = check_esets(MB, start_num_sets + MBENTITYSET - MBEDGE + 4);
  if (MB_SUCCESS != result) return result;

    //-------------Unite--------------
    //fill temp_ms1 with tets and tris
    //get all tris 
  temp_range.clear();
  result = MB->get_entities_by_handle(ms_array[MBTRI], temp_range );
  if(result != MB_SUCCESS )
    return result;

    //add tets them to ms1 
  result = MB->add_entities( temp_ms1, temp_range ) ;
  if(result != MB_SUCCESS )
    return result;

  temp_range.clear();
  result = MB->get_entities_by_handle(ms_array[MBTET], temp_range );
  if(result != MB_SUCCESS )
    return result;

    //add tets them to ms1 
  result = MB->add_entities( temp_ms1, temp_range ) ;
  if(result != MB_SUCCESS )
    return result;


    //fill temp_ms2 with tris and hexes
  temp_range.clear();
  result = MB->get_entities_by_handle(ms_array[MBTRI], temp_range );
  if(result != MB_SUCCESS )
    return result;

    //add tets them to ms1 
  result = MB->add_entities( temp_ms2, temp_range ) ;
  if(result != MB_SUCCESS )
    return result;

  temp_range.clear();
  result = MB->get_entities_by_handle(ms_array[MBQUAD], temp_range );
  if(result != MB_SUCCESS )
    return result;

    //add tets them to ms1 
  result = MB->add_entities( temp_ms2, temp_range ) ;
  if(result != MB_SUCCESS )
    return result;


    //temp_ms1 is now filled with tets and tris 
    //temp_ms2 is now filled with quads and tris 
  int size1 = 0, size2 = 0;
  result = MB->get_number_entities_by_handle(ms_array[MBTRI], size1 ) ;
  if(result != MB_SUCCESS )
    return result;

  result = MB->intersect_meshset(temp_ms1, temp_ms2 ) ;
  if(result != MB_SUCCESS )
    return result;
  result = MB->get_number_entities_by_handle(temp_ms1, size2 ) ;
  if(result != MB_SUCCESS )
    return result;

  if(size1 != size2)
  {
    return MB_FAILURE;
  }

  temp_range.clear(); 
  result = MB->get_entities_by_handle(temp_ms1, temp_range );
  if(result != MB_SUCCESS )
    return result;

    //make sure there's nothing but tris in temp_range 
  if( nothing_but_type( temp_range, MBTRI ) != MB_SUCCESS) 
    return MB_FAILURE; 


  result = check_esets(MB, start_num_sets + MBENTITYSET - MBEDGE + 4);
  if (MB_SUCCESS != result) return result;

    //-------------Misc tests--------------
  MBEntityHandle temp_ms3;
  result = MB->create_meshset(flags, temp_ms3);
  if(result  != MB_SUCCESS ) 
    return result;

  MBEntityHandle handle_array[] = {1, 2, 3, 4, 5, 7, 8, 9, 10};
  const int num_handle = sizeof(handle_array)/sizeof(handle_array[0]);
    //add ents to set
  result = MB->add_entities( temp_ms3, handle_array, num_handle );
  if(result  != MB_SUCCESS ) 
    return result;

    // try adding again
  result = MB->add_entities( temp_ms3, handle_array, num_handle );
  if(result  != MB_SUCCESS ) 
    return result;

  int num_ents;
  result = MB->get_number_entities_by_handle(temp_ms3, num_ents);
  if(result  != MB_SUCCESS ) 
    return result;
    
  int num_expected = (flags & MESHSET_SET) ? num_handle : 2*num_handle;
  if (num_ents != num_expected)
    return MB_FAILURE;

  return MB_SUCCESS;
}

static bool compare_lists( std::vector<MBEntityHandle> vect,
                           const MBEntityHandle* array, 
                           int count,
                           bool ordered = true )
{
  if (vect.size() != (size_t)count)
    return false;
  for (int i = 0; i < count; ++i) {
    if (ordered) { 
      if (vect[i] != array[i])
        return false;
    }
    else if (std::find(vect.begin(),vect.end(),array[i]) == vect.end())
      return false;
  }
  return true;
}

    //Test parent/child stuff in meshsets
MBErrorCode mb_mesh_set_parent_child_test(MBInterface *MB)
{
  MBErrorCode rval;
  std::vector<MBEntityHandle> list;
  MBRange range;
  MBRange::iterator iter;
  int count;

    // create a few mesh sets
  const int num_sets = 10;
  MBEntityHandle sets[num_sets];
  for (int i = 0; i < num_sets; ++i) {
    rval = MB->create_meshset( i % 2 ? MESHSET_SET : 0, sets[i] );
    if (MB_SUCCESS != rval) 
      return rval;
  }
  
    // test adding child meshsets
    
    // add first child
  rval = MB->add_child_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 1 ))
    return MB_FAILURE;
    // try to add child again
  rval = MB->add_child_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 1 ))
    return MB_FAILURE;
    
    // add second child
  rval = MB->add_child_meshset( sets[0], sets[2] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 2 ))
    return MB_FAILURE;
    // try adding child again
  rval = MB->add_child_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 2 ))
    return MB_FAILURE;
  
    // add third child
  rval = MB->add_child_meshset( sets[0], sets[3] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 3 ))
    return MB_FAILURE;
    // try adding child again
  rval = MB->add_child_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 3 ))
    return MB_FAILURE;
  
    // add fourth child
  rval = MB->add_child_meshset( sets[0], sets[4] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 4 ))
    return MB_FAILURE;
  
    // make sure range query returns same result
  std::sort( list.begin(), list.end() );
  rval = MB->get_child_meshsets( sets[0], range );
  iter = range.begin();
  for (unsigned i = 0; i < 4; ++i, ++iter)
    if (*iter != list[i])
      return MB_FAILURE;
  
    // remove first child
  rval = MB->remove_child_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+2, 3 ))
    return MB_FAILURE;
    // try removing child again
  rval = MB->remove_child_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  
    // remove second child
  rval = MB->remove_child_meshset( sets[0], sets[2] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list );
  if (!compare_lists( list, sets+3, 2 ))
    return MB_FAILURE;
    // try removing child again
  rval = MB->remove_child_meshset( sets[0], sets[2] );
  if (MB_SUCCESS != rval) return rval;
  
    // remove third child
  rval = MB->remove_child_meshset( sets[0], sets[3] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list );
  if (list.size() != 1 || list[0] != sets[4])
    return MB_FAILURE;
  
    // remove fourth child
  rval = MB->remove_child_meshset( sets[0], sets[4] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list );
  if (!list.empty())
    return MB_FAILURE;

  
    // test adding parent meshsets
    
    // add first parent
  rval = MB->add_parent_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 1 ))
    return MB_FAILURE;
    // try to add parent again
  rval = MB->add_parent_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 1 ))
    return MB_FAILURE;
    
    // add second parent
  rval = MB->add_parent_meshset( sets[0], sets[2] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 2 ))
    return MB_FAILURE;
    // try adding parent again
  rval = MB->add_parent_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 2 ))
    return MB_FAILURE;
  
    // add third parent
  rval = MB->add_parent_meshset( sets[0], sets[3] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 3 ))
    return MB_FAILURE;
    // try adding parent again
  rval = MB->add_parent_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 3 ))
    return MB_FAILURE;
  
    // add fourth parent
  rval = MB->add_parent_meshset( sets[0], sets[4] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 4 ))
    return MB_FAILURE;
  
    // make sure range query returns same result
  std::sort( list.begin(), list.end() );
  rval = MB->get_parent_meshsets( sets[0], range );
  iter = range.begin();
  for (unsigned i = 0; i < 4; ++i, ++iter)
    if (*iter != list[i])
      return MB_FAILURE;
  
    // remove first parent
  rval = MB->remove_parent_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+2, 3 ))
    return MB_FAILURE;
    // try removing parent again
  rval = MB->remove_parent_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  
    // remove second parent
  rval = MB->remove_parent_meshset( sets[0], sets[2] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list );
  if (!compare_lists( list, sets+3, 2 ))
    return MB_FAILURE;
    // try removing parent again
  rval = MB->remove_parent_meshset( sets[0], sets[2] );
  if (MB_SUCCESS != rval) return rval;
  
    // remove third parent
  rval = MB->remove_parent_meshset( sets[0], sets[3] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list );
  if (list.size() != 1 || list[0] != sets[4])
    return MB_FAILURE;
  
    // remove fourth parent
  rval = MB->remove_parent_meshset( sets[0], sets[4] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list );
  if (!list.empty())
    return MB_FAILURE;
    
    
    // setup tests of recursive child query
    //              0
    //            /   \         .
    //          1       2
    //        /   \   /   \       .
    //      3       4       5
    //        \   /   \   /
    //          6       7
  rval = MB->add_child_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_child_meshset( sets[0], sets[2] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_child_meshset( sets[1], sets[3] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_child_meshset( sets[1], sets[4] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_child_meshset( sets[2], sets[4] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_child_meshset( sets[2], sets[5] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_child_meshset( sets[3], sets[6] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_child_meshset( sets[4], sets[6] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_child_meshset( sets[4], sets[7] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_child_meshset( sets[5], sets[7] );
  if (MB_SUCCESS != rval) return rval;
  
    // test query at depth of 1
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list, 1 );
  if (MB_SUCCESS != rval) return rval;
  if (list.size() != 2 || list[0] != sets[1] || list[1] != sets[2])
    return MB_FAILURE;
  rval = MB->num_child_meshsets( sets[0], &count, 1 );
  if (MB_SUCCESS != rval) return rval;
  if (count != 2)
    return MB_FAILURE;
  
    // test query at depth of 2
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list, 2 );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 5, false ))
    return MB_FAILURE;  
  rval = MB->num_child_meshsets( sets[0], &count, 2 );
  if (MB_SUCCESS != rval) return rval;
  if (count != 5)
    return MB_FAILURE;
  
    // test query at depth of 3
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list, 3 );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 7, false ))
    return MB_FAILURE;  
  rval = MB->num_child_meshsets( sets[0], &count, 3 );
  if (MB_SUCCESS != rval) return rval;
  if (count != 7)
    return MB_FAILURE;
  
    // test query at depth of 4
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list, 4 );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 7, false ))
    return MB_FAILURE;  
  rval = MB->num_child_meshsets( sets[0], &count, 4 );
  if (MB_SUCCESS != rval) return rval;
  if (count != 7)
    return MB_FAILURE;
  
    // test query at all
  list.clear();
  rval = MB->get_child_meshsets( sets[0], list, 0 );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 7, false ))
    return MB_FAILURE;  
  rval = MB->num_child_meshsets( sets[0], &count, 0 );
  if (MB_SUCCESS != rval) return rval;
  if (count != 7)
    return MB_FAILURE;
    
    // clean up child links
  rval = MB->remove_child_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_child_meshset( sets[0], sets[2] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_child_meshset( sets[1], sets[3] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_child_meshset( sets[1], sets[4] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_child_meshset( sets[2], sets[4] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_child_meshset( sets[2], sets[5] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_child_meshset( sets[3], sets[6] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_child_meshset( sets[4], sets[6] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_child_meshset( sets[4], sets[7] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_child_meshset( sets[5], sets[7] );
  if (MB_SUCCESS != rval) return rval;
  for (int i = 0; i < 5; ++i)
    if (MB_SUCCESS != MB->num_child_meshsets(sets[i], &count) || count)
      return MB_FAILURE;
     
    // setup tests of recursive parent query
    //          6       7
    //        /   \   /   \       .
    //      3       4       5
    //        \   /   \   /
    //          1       2
    //            \   / 
    //              0
  rval = MB->add_parent_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_parent_meshset( sets[0], sets[2] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_parent_meshset( sets[1], sets[3] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_parent_meshset( sets[1], sets[4] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_parent_meshset( sets[2], sets[4] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_parent_meshset( sets[2], sets[5] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_parent_meshset( sets[3], sets[6] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_parent_meshset( sets[4], sets[6] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_parent_meshset( sets[4], sets[7] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->add_parent_meshset( sets[5], sets[7] );
  if (MB_SUCCESS != rval) return rval;
  
    // test query at depth of 1
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list, 1 );
  if (MB_SUCCESS != rval) return rval;
  if (list.size() != 2 || list[0] != sets[1] || list[1] != sets[2])
    return MB_FAILURE;
  rval = MB->num_parent_meshsets( sets[0], &count, 1 );
  if (MB_SUCCESS != rval) return rval;
  if (count != 2)
    return MB_FAILURE;
  
    // test query at depth of 2
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list, 2 );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 5, false ))
    return MB_FAILURE;  
  rval = MB->num_parent_meshsets( sets[0], &count, 2 );
  if (MB_SUCCESS != rval) return rval;
  if (count != 5)
    return MB_FAILURE;
  
    // test query at depth of 3
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list, 3 );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 7, false ))
    return MB_FAILURE;  
  rval = MB->num_parent_meshsets( sets[0], &count, 3 );
  if (MB_SUCCESS != rval) return rval;
  if (count != 7)
    return MB_FAILURE;
  
    // test query at depth of 4
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list, 4 );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 7, false ))
    return MB_FAILURE;  
  rval = MB->num_parent_meshsets( sets[0], &count, 4 );
  if (MB_SUCCESS != rval) return rval;
  if (count != 7)
    return MB_FAILURE;
  
    // test query at all
  list.clear();
  rval = MB->get_parent_meshsets( sets[0], list, 0 );
  if (MB_SUCCESS != rval) return rval;
  if (!compare_lists( list, sets+1, 7, false ))
    return MB_FAILURE;  
  rval = MB->num_parent_meshsets( sets[0], &count, 0 );
  if (MB_SUCCESS != rval) return rval;
  if (count != 7)
    return MB_FAILURE;
    
    // clean up parent links
  rval = MB->remove_parent_meshset( sets[0], sets[1] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_parent_meshset( sets[0], sets[2] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_parent_meshset( sets[1], sets[3] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_parent_meshset( sets[1], sets[4] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_parent_meshset( sets[2], sets[4] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_parent_meshset( sets[2], sets[5] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_parent_meshset( sets[3], sets[6] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_parent_meshset( sets[4], sets[6] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_parent_meshset( sets[4], sets[7] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->remove_parent_meshset( sets[5], sets[7] );
  if (MB_SUCCESS != rval) return rval;
  for (int i = 0; i < 5; ++i)
    if (MB_SUCCESS != MB->num_parent_meshsets(sets[i], &count) || count)
      return MB_FAILURE;
   
    
    // test combined parent/child links

    // test creation
  rval = MB->add_parent_child( sets[9], sets[8] );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_child_meshsets( sets[9], list );
  if (MB_SUCCESS != rval) return rval;
  if (list.size() != 1 || list[0] != sets[8])
    return MB_FAILURE;  
  list.clear();
  rval = MB->get_parent_meshsets( sets[8], list );
  if (MB_SUCCESS != rval) return rval;
  if (list.size() != 1 || list[0] != sets[9])
    return MB_FAILURE;  

    // test deletion of parent/child
  rval = MB->add_parent_child( sets[7], sets[9] );
  if (MB_SUCCESS != rval) return rval;
  rval = MB->delete_entities( &sets[9], 1 );
  if (MB_SUCCESS != rval) return rval;
  list.clear();
  rval = MB->get_parent_meshsets( sets[8], list );
  if (!list.empty())
    return MB_FAILURE;
  list.clear();
  rval = MB->get_child_meshsets( sets[7], list );
  if (!list.empty())
    return MB_FAILURE;
  
    // clean up remaining sets
  return MB->delete_entities( sets, 9 );
}
  
MBErrorCode mb_mesh_sets_set_test( MBInterface* mb )
{
  return mb_mesh_sets_test( mb, MESHSET_SET );
}
  
MBErrorCode mb_mesh_sets_list_test( MBInterface* mb )
{
  return mb_mesh_sets_test( mb, MESHSET_ORDERED );
} 

// Verify that all query functions *append* to output MBRange
MBErrorCode mb_mesh_set_appends( MBInterface* mb, int flags )
{
  MBErrorCode rval;
  
    // get all handles and subdivide into vertex and non-vertex ents
  MBRange all_ents, verts, elems, results;
  rval = mb->get_entities_by_handle( 0, all_ents );
  if (MB_SUCCESS != rval)
    return rval;
  MBRange::iterator ve = all_ents.upper_bound( MBVERTEX );
  verts.merge( all_ents.begin(), ve );
  elems.merge( ve, all_ents.end() );

    // If we're not testing queries from the root set, 
    // create a set containing all the vertices
  MBEntityHandle set = 0;
  if (flags != -1) {
    rval = mb->create_meshset( flags, set );
    if (MB_SUCCESS != rval)
      return rval;
    rval = mb->add_entities( set, verts );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
    // Test get_entities_by_handle.  This one doesn't make
    // much sense if we're testing the root set, but it doesn't
    // hurt to test it anyway.
  results = elems;
  rval = mb->get_entities_by_handle( set, results );
  if (MB_SUCCESS != rval)
    return rval;
  if (results != all_ents)
    return MB_FAILURE;
  
    // Test get_entities_by_dimension
  results = elems;
  rval = mb->get_entities_by_dimension( set, 0, results );
  if (MB_SUCCESS != rval)
    return rval;
  if (results != all_ents)
    return MB_FAILURE;
  
    // Test get_entities_by_type
  results = elems;
  rval = mb->get_entities_by_type( set, MBVERTEX, results );
  if (MB_SUCCESS != rval)
    return rval;
  if (results != all_ents)
    return MB_FAILURE;
  
    // choose a single entity for testing tag queries
  MBEntityHandle entity = verts.front();
  MBRange expected( elems );
  expected.insert( entity );
  
  MBTag sparse, dense;
  const int zero = 0, one = 1;
  const void* vals[] = {&one};

    // Test get_entities_by_type_and_tag w/ sparse tag and no value
  rval = mb->tag_create( "mb_mesh_set_appends_sparse", 
                         sizeof(int), 
                         MB_TAG_SPARSE, 
                         MB_TYPE_INTEGER,
                         sparse, 0 );
  if (MB_SUCCESS != rval)
    return rval;
  rval = mb->tag_set_data( sparse, &entity, 1, &one );
  if (MB_SUCCESS != rval)
    return rval;
  results = elems;
  rval = mb->get_entities_by_type_and_tag( set, 
                                           TYPE_FROM_HANDLE(entity),
                                           &sparse, 0, 1, 
                                           results, MBInterface::UNION );
  if (MB_SUCCESS != rval)
    return rval;
  if (results != expected)
    return MB_FAILURE;
    // Test again, but specify tag value
  results = elems;
  rval = mb->get_entities_by_type_and_tag( set, 
                                           TYPE_FROM_HANDLE(entity),
                                           &sparse, vals, 1, 
                                           results, MBInterface::UNION );
  if (MB_SUCCESS != rval)
    return rval;
  if (results != expected)
    return MB_FAILURE;
   
    // Test get_entities_by_type_and_tag w/ dense tag
  rval = mb->tag_create( "mb_mesh_set_appends_dense", 
                         sizeof(int), 
                         MB_TAG_DENSE, 
                         MB_TYPE_INTEGER,
                         dense, &zero );
  if (MB_SUCCESS != rval)
    return rval;
  rval = mb->tag_set_data( dense, &entity, 1, &one );
  if (MB_SUCCESS != rval)
    return rval;
  results = elems;
  rval = mb->get_entities_by_type_and_tag( set, 
                                           TYPE_FROM_HANDLE(entity),
                                           &dense, vals, 1, 
                                           results, MBInterface::UNION );
  if (MB_SUCCESS != rval)
    return rval;
  if (results != expected)
    return MB_FAILURE;
 
  return MB_SUCCESS;
}

MBErrorCode mb_mesh_set_set_appends( MBInterface* mb )
{
  return mb_mesh_set_appends( mb, MESHSET_SET );
}

MBErrorCode mb_mesh_set_list_appends( MBInterface* mb )
{
  return mb_mesh_set_appends( mb, MESHSET_ORDERED );
}

MBErrorCode mb_mesh_set_root_appends( MBInterface* mb )
{
  return mb_mesh_set_appends( mb, -1 );
}

  // number of entities of type MBVERTEX, MBEDGE, MBDTri, MBQUAD, MBTET, and MBHEX
  // in mbtest1.g  (all other values are 0.
static const unsigned int num_entities[MBMAXTYPE] = {47,12,18,8,22,8};

MBErrorCode mb_delete_mesh_test(MBInterface *gMB)
{
  MBErrorCode error = MB_SUCCESS;

    // Lets also test the global MB pointer (gMB) here.
  error = gMB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;

    // load the mesh again 
  std::string file_name = TestDir + "/mbtest1.g";
  error = gMB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;


  MBRange entities;
  error = gMB->get_entities_by_type(0,  MBVERTEX, entities);
  if (error != MB_SUCCESS)
    return error;

    // As before the number of vertices had better be 83
  if ( entities.size() != num_entities[MBVERTEX] )
    return MB_FAILURE;

  MBTag tag_handle = 0;
  MBEntityType type;
    // step through each element type
  for (type = MBEDGE; type != MBENTITYSET; type++)
  {
      // There should be entities
    error = gMB->tag_get_handle("connectivity", tag_handle);
    if (error == MB_SUCCESS)
    {
      entities.clear();
      error = gMB->get_entities_by_type_and_tag(0,  type, &tag_handle, NULL, 1, entities);
      if (error != MB_SUCCESS)
        return error;

      if (!entities.empty())
      {
        if ( entities.size() != num_entities[type] )
          return MB_FAILURE;
      }
    }
  }

  return MB_SUCCESS;
}


MBErrorCode mb_meshset_tracking_test( MBInterface *MB )
{

    //read in a file so you have some data in the database
  std::string file_name = TestDir + "/mbtest1.g";
  MBErrorCode error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;

  MBEntityHandle ms1, ms2, ms3;

    //create meshsets 
  MBErrorCode result = MB->create_meshset( MESHSET_TRACK_OWNER | MESHSET_ORDERED, ms1 ) ;
  if(result != MB_SUCCESS )
    return result;
  result = MB->create_meshset( MESHSET_SET | MESHSET_TRACK_OWNER, ms2 ) ;
  if(result != MB_SUCCESS )
    return result;
  result = MB->create_meshset( MESHSET_SET | MESHSET_TRACK_OWNER, ms3 ) ;
  if(result != MB_SUCCESS )
    return result;

    // get all hexes 
  MBRange hexes;
  result = MB->get_entities_by_type(0, MBHEX, hexes);
  if(result != MB_SUCCESS )
    return result;

    // get all tris 
  MBRange tris;
  result = MB->get_entities_by_type(0, MBTRI, tris );
  if(result != MB_SUCCESS )
    return result;

    // get all tets 
  MBRange temp_range;
  result = MB->get_entities_by_type(0, MBTET, temp_range);
  if(result != MB_SUCCESS )
    return result;

    //copy all the tets from the range into a vector 'tets'
  std::vector<MBEntityHandle> tets( temp_range.size() );
  std::copy(temp_range.begin(), temp_range.end(), tets.begin() );

    //Quick check on 'get_entities_by_dimension()'
  MBRange dim_3_range;
  result = MB->get_entities_by_dimension(0, 3, dim_3_range) ;
  if(result != MB_SUCCESS )
    return result;

    //hexes and tets should be only dimension 3 entities
  if( hexes.size() + tets.size() != dim_3_range.size() )
    return MB_FAILURE;

    //put all hexes in ms1, ms2, ms3
  result = MB->add_entities(ms1, hexes); //add ents in a range 
  if(result != MB_SUCCESS )
    return result;
    // to ordered meshset 

  result = MB->add_entities(ms2, hexes); //add ents in a range 
  if(result != MB_SUCCESS )
    return result;
    //to unordered meshset

  result = MB->add_entities(ms3, hexes);
  if(result != MB_SUCCESS )
    return result;

    //put all tets in ms1, ms2
  if(MB->add_entities(ms1, &tets[0], tets.size()) != MB_SUCCESS )  //add ents in a vector
    return MB_FAILURE;                             //to ordered meshset 

  if(MB->add_entities(ms2, &tets[0], tets.size()) != MB_SUCCESS )  //add ents in a vector
    return MB_FAILURE;                             //to unordered meshset

    //put all tris in ms1
  result = MB->add_entities(ms1, tris) ;
  if(result != MB_SUCCESS )
    return result;

  MBRange::iterator iter;
  iter = tris.begin();

  std::vector< MBEntityHandle > temp_vec;

    //ask a tri which meshsets it is in
  result = MB->get_adjacencies( &(*iter), 1, 4, false, temp_vec ) ;
  if(result != MB_SUCCESS )
    return result;

    //cout<<"tris temp_vec.size() = "<<temp_vec.size()<<endl; 
  if( temp_vec.size() != 2 )
    return MB_FAILURE;

    //ask a tet which meshsets it is in
  temp_vec.clear();
  std::vector<MBEntityHandle>::iterator vec_iter = tets.begin();
  result = MB->get_adjacencies( &(*vec_iter), 1, 4, false, temp_vec ) ;
  if(result != MB_SUCCESS )
    return result;

    //cout<<"tet temp_vec.size() = "<<temp_vec.size()<<endl; 
  if( temp_vec.size() != 3 )
    return MB_FAILURE;

    //ask a hex which meshsets it is in  
  temp_vec.clear();
  iter = hexes.begin();
  result = MB->get_adjacencies( &(*iter), 1, 4, false, temp_vec ) ;
  if(result != MB_SUCCESS )
    return result;

    //should be in 4 temp_vec
  if( temp_vec.size() != 4 )
    return MB_FAILURE;

    //take this hex out of the ms1, ms2, ms3
  if(MB->remove_entities(ms1, &(*iter), 1) != MB_SUCCESS ) //remove ents in a vector  
    return MB_FAILURE;                                   //from ordered meshset

  if(MB->remove_entities(ms2, &(*iter), 1) != MB_SUCCESS ) //remove ents in a vector
    return MB_FAILURE;                                   //from unordered meshset

  temp_range.clear();
  temp_range.insert(*iter);
  if(MB->remove_entities(ms3, temp_range) != MB_SUCCESS ) //remove ents in a range
    return MB_FAILURE;                                     //from unordered meshset


    //ask the hex how many meshsets it is in
  temp_vec.clear();
  result = MB->get_adjacencies( &(*iter), 1, 4, false, temp_vec ) ;
  if(result != MB_SUCCESS )
    return result;
  if( temp_vec.size() != 1 )
    return MB_FAILURE;

    //add the hex back into ms1
  result = MB->add_entities(ms1, temp_range) ;
  if(result != MB_SUCCESS )
    return result;

    //ask the hex how many meshsets it is in
  temp_vec.clear();
  result = MB->get_adjacencies( &(*iter), 1, 4, false, temp_vec ) ;
  if(result != MB_SUCCESS )
    return result;
  if( temp_vec.size() != 2 )
    return MB_FAILURE;

  temp_range.clear();
  temp_range.insert(*iter);
  if(MB->remove_entities(ms1, temp_range) != MB_SUCCESS ) //remove ents in a range 
    return MB_FAILURE;                                     //from an ordered meshset 

    //ask the hex how many meshsets it is in
  temp_vec.clear();
  result = MB->get_adjacencies( &(*iter), 1, 4, false, temp_vec ) ;
  if(result != MB_SUCCESS )
    return result;
  if( temp_vec.size() != 1 )
    return MB_FAILURE;


    //Deleting a meshset

  iter = tris.begin();
  temp_vec.clear();
    //ask a tri which meshsets it is in
  result = MB->get_adjacencies( &(*iter), 1, 4, false, temp_vec ) ;
  if(result != MB_SUCCESS )
    return result;

  if( temp_vec.size() != 2 )
    return MB_FAILURE;

    //Try deleting a meshset
  result = MB->delete_entities(&ms1, 1);
  if(result != MB_SUCCESS )
    return result;

  temp_vec.clear();
    //Query tri again for meshsets it's in
  result = MB->get_adjacencies( &(*iter), 1, 4, false, temp_vec ) ;
  if(result != MB_SUCCESS )
    return result;

  if( temp_vec.size() != 1 )
    return MB_FAILURE;

    //Delete an entitiy from ms2....make sure it's removed out of ms2 
  int num_before = 0;
  MB->get_number_entities_by_handle(ms2, num_before);
  vec_iter = tets.begin();
  result = MB->delete_entities( &(*vec_iter), 1);
  if(result != MB_SUCCESS )
    return result;

  int num_after = 0;
  MB->get_number_entities_by_handle(ms2, num_after);
  if( num_before != num_after + 1)
    return MB_FAILURE;


  return MB_SUCCESS;

}

MBErrorCode mb_write_mesh_test(MBInterface *MB)
{
  std::string file_name = "test/mb_write.g";

    // no need to get lists, write out the whole mesh
  MBErrorCode result = MB->write_mesh(file_name.c_str());
  if(result != MB_SUCCESS )
    return result;

    //---------The following tests outputting meshsets that are in meshsets of blocks ---/

    //lets create a block meshset and put some entities and meshsets into it
  MBEntityHandle block_ms;
  result = MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, block_ms );
  if(result != MB_SUCCESS )
    return result;

    //make another meshset to put quads in, so SHELLs can be written out
  MBEntityHandle block_of_shells;
  result = MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, block_of_shells); 
  if(result != MB_SUCCESS )
    return result;

    //tag the meshset so it's a block, with id 100
  int id = 100;
  MBTag tag_handle;
  result = MB->tag_get_handle( MATERIAL_SET_TAG_NAME, tag_handle ) ;
  if(result != MB_SUCCESS )
    return result;
  result = MB->tag_set_data( tag_handle, &block_ms, 1, &id ) ;
  if(result != MB_SUCCESS )
    return result;
  id = 101;
  result = MB->tag_set_data( tag_handle, &block_of_shells, 1, &id ) ;
  if(result != MB_SUCCESS )
    return result;

    // set dimension tag on this to ensure shells get output; reuse id variable
  result = MB->tag_get_handle( GEOM_DIMENSION_TAG_NAME, tag_handle) ;
  if(result != MB_SUCCESS )
    return result;
  id = 3;
  result = MB->tag_set_data( tag_handle, &block_of_shells, 1, &id ) ;
  if(result != MB_SUCCESS )
    return result;

    //get some entities (tets) 
  MBRange temp_range;
  result = MB->get_entities_by_type(0,  MBHEX, temp_range ) ;
  if(result != MB_SUCCESS )
    return result;

  MBRange::iterator iter, end_iter;
  iter = temp_range.begin();
  end_iter = temp_range.end();

    //add evens to 'block_ms'
  std::vector<MBEntityHandle> temp_vec; 
  for(; iter != end_iter; iter++)
  {
    if( ID_FROM_HANDLE( *iter ) % 2 == 0 ) 
      temp_vec.push_back( *iter );
  }
  result = MB->add_entities( block_ms, &temp_vec[0], temp_vec.size()); 
  if(result != MB_SUCCESS )
    return result;


    //make another meshset
  MBEntityHandle ms_of_block_ms;
  result = MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, ms_of_block_ms);
  if(result != MB_SUCCESS )
    return result;

    //add some entities to it
  temp_vec.clear();
  iter = temp_range.begin();
  for(; iter != end_iter; iter++)
  {
    if( ID_FROM_HANDLE( *iter ) % 2 )  //add all odds
      temp_vec.push_back( *iter );
  }
  result = MB->add_entities( ms_of_block_ms, &temp_vec[0], temp_vec.size() ); 
  if(result != MB_SUCCESS )
    return result;

    //add the other meshset to the block's meshset
  result = MB->add_entities( block_ms, &ms_of_block_ms, 1);
  if(result != MB_SUCCESS )
    return result;


    //---------------testing sidesets----------------/

    //lets create a sideset meshset and put some entities and meshsets into it
  MBEntityHandle sideset_ms;
  result = MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, sideset_ms );
  if(result != MB_SUCCESS )
    return result;

    //tag the meshset so it's a sideset, with id 104
  id = 104;
  result = MB->tag_get_handle( NEUMANN_SET_TAG_NAME, tag_handle ) ;
  if(result != MB_SUCCESS )
    return result;

  result = MB->tag_set_data( tag_handle, &sideset_ms, 1, &id ) ;
  if(result != MB_SUCCESS )
    return result;

    //get some entities (tris) 
  temp_range.clear();
  result = MB->get_entities_by_type(0,  MBQUAD, temp_range ) ;
  if(result != MB_SUCCESS )
    return result;

  iter = temp_range.begin();
  end_iter = temp_range.end();

    //add evens to 'sideset_ms'
  temp_vec.clear(); 
  for(; iter != end_iter; iter++)
  {
    if( ID_FROM_HANDLE( *iter ) % 2 == 0 ) 
      temp_vec.push_back( *iter );
  }
  result = MB->add_entities( sideset_ms, &temp_vec[0], temp_vec.size() ); 
  if(result != MB_SUCCESS )
    return result;

    //make another meshset
  MBEntityHandle ms_of_sideset_ms;
  result = MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, ms_of_sideset_ms);
  if(result != MB_SUCCESS )
    return result;

    //add some entities to it
  temp_vec.clear();
  iter = temp_range.begin();
  for(; iter != end_iter; iter++)
  {
    if( ID_FROM_HANDLE( *iter ) % 2 )  //add all odds
      temp_vec.push_back( *iter );
  }
  result = MB->add_entities( ms_of_sideset_ms, &temp_vec[0], temp_vec.size() ); 
  if(result != MB_SUCCESS )
    return result;

    //add the other meshset to the sideset's meshset
  result = MB->add_entities( sideset_ms, &ms_of_sideset_ms, 1);
  if(result != MB_SUCCESS )
    return result;

    //---------test sense on meshsets (reverse/foward)-------//

    //get all quads whose x-coord = 2.5 and put them into a meshset_a 
  MBEntityHandle meshset_a;
  result = MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, meshset_a );
  if(result != MB_SUCCESS )
    return result;

  temp_range.clear();
  result = MB->get_entities_by_type(0,  MBQUAD, temp_range ) ;
  if(result != MB_SUCCESS )
    return result;

  std::vector<MBEntityHandle> nodes, entity_vec;
  std::copy(temp_range.begin(), temp_range.end(), std::back_inserter(entity_vec));
  result = MB->get_connectivity(&entity_vec[0], entity_vec.size(), nodes);
  if(result != MB_SUCCESS ) 
    return result;
  assert( nodes.size() == 4 * temp_range.size() );
  temp_vec.clear(); 
  std::vector<double> coords(3*nodes.size());
  result = MB->get_coords(&nodes[0], nodes.size(), &coords[0]);
  if(result != MB_SUCCESS ) 
    return result;
  
  unsigned int k = 0;
  for(MBRange::iterator it = temp_range.begin(); it != temp_range.end(); it++) {
    if( coords[12*k] == 2.5 && coords[12*k+3] == 2.5 &&
        coords[12*k+6] == 2.5 && coords[12*k+9] == 2.5 )
      temp_vec.push_back(*it);
    k++;
  }
  result = MB->add_entities( meshset_a, &temp_vec[0], temp_vec.size() );
  if(result != MB_SUCCESS ) 
    return result;
  result = MB->add_entities( block_of_shells, &temp_vec[0], temp_vec.size());
  if(result != MB_SUCCESS ) 
    return result;

    //put these quads into a different meshset_b and tag them with a reverse sense tag
  MBEntityHandle meshset_b;
  result = MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, meshset_b );
  if(result != MB_SUCCESS ) 
    return result;

  result = MB->add_entities( meshset_b, &meshset_a, 1);
  if(result != MB_SUCCESS ) 
    return result;


  result = MB->tag_get_handle( "SENSE", tag_handle );

  if(result != MB_SUCCESS ) 
  {
      //create the tag
    int default_value = 0;
    result = MB->tag_create( "SENSE", sizeof(int), MB_TAG_SPARSE, tag_handle, 
                             &default_value );
    if(result != MB_SUCCESS)
      return result;
  }

  int reverse_value = -1;
  result = MB->tag_set_data( tag_handle, &meshset_b, 1, &reverse_value ) ; 
  if(result != MB_SUCCESS)
    return result;


    //get some random quad, whose x-coord != 2.5, and put it into a different meshset_c
    //and tag it with a reverse sense tag

  iter = temp_range.begin();
  end_iter = temp_range.end();

  temp_vec.clear();
  for(; iter != end_iter; iter++ )
  {
    std::vector<MBEntityHandle> nodes;
    result = MB->get_connectivity( &(*iter), 1, nodes );
    if(result != MB_SUCCESS)
      return result;

    bool not_equal_2_5 = true; 
    for(unsigned int k=0; k<nodes.size(); k++ )
    {
      double coords[3] = {0};

      result = MB->get_coords( &(nodes[k]), 1, coords );
      if(result != MB_SUCCESS)
        return result;

      if( coords[0] == 2.5 )
      {
        not_equal_2_5 = false;
        break;
      }
    }

    if( not_equal_2_5 && nodes.size()> 0)
    {
      temp_vec.push_back( *iter );
      break;
    }
  }

  MBEntityHandle meshset_c;
  MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, meshset_c );
    
  
  result = MB->tag_get_handle( "SENSE", tag_handle ); 
  if(result != MB_SUCCESS)
    return result;

  reverse_value = -1;
  result = MB->tag_set_data( tag_handle, &meshset_c, 1, &reverse_value ) ; 
  if(result != MB_SUCCESS)
    return result;

  MB->add_entities( meshset_c, &temp_vec[0], temp_vec.size() );
  MB->add_entities( block_of_shells, &temp_vec[0], temp_vec.size());


    //create another meshset_abc, adding meshset_a, meshset_b, meshset_c 
  MBEntityHandle meshset_abc;
  MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, meshset_abc );

  temp_vec.clear();
  temp_vec.push_back( meshset_a );
  temp_vec.push_back( meshset_b );
  temp_vec.push_back( meshset_c );

  MB->add_entities( meshset_abc, &temp_vec[0], temp_vec.size());


    //tag it so it's a sideset
  id = 444;
  result = MB->tag_get_handle( "NEUMANN_SET", tag_handle ) ;
  if(result != MB_SUCCESS)
    return result;

  result = MB->tag_set_data( tag_handle, &meshset_abc, 1, &id ) ;
  if(result != MB_SUCCESS)
    return result;



    //---------------do nodesets now -----------------//


    //lets create a nodeset meshset and put some entities and meshsets into it
  MBEntityHandle nodeset_ms;
  MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, nodeset_ms );

    //tag the meshset so it's a nodeset, with id 119
  id = 119;
  result = MB->tag_get_handle( DIRICHLET_SET_TAG_NAME, tag_handle ) ;
  if(result != MB_SUCCESS)
    return result;

  result = MB->tag_set_data( tag_handle, &nodeset_ms, 1, &id ) ;
  if(result != MB_SUCCESS)
    return result;

    //get all Quads 
  temp_range.clear();
  result = MB->get_entities_by_type(0,  MBQUAD, temp_range ) ;
  if(result != MB_SUCCESS)
    return result;


    //get all the nodes of the tris
  MBRange nodes_of_quads;
  iter = temp_range.begin();
  end_iter = temp_range.end();


  for(; iter != end_iter; iter++ )
  {
    std::vector<MBEntityHandle> nodes;
    result = MB->get_connectivity( &(*iter), 1, nodes);
    if(result != MB_SUCCESS)
      return result;

    for(unsigned int k=0; k<nodes.size(); k++ )
      nodes_of_quads.insert( nodes[k] ); 

  }

  iter = nodes_of_quads.begin();
  end_iter = nodes_of_quads.end();

    //add evens to 'nodeset_ms'
  temp_vec.clear(); 
  for(; iter != end_iter; iter++)
  {
    if( ID_FROM_HANDLE( *iter ) % 2 == 0 ) 
      temp_vec.push_back( *iter );
  }
  MB->add_entities( nodeset_ms, &temp_vec[0], temp_vec.size() ); 


    //make another meshset
  MBEntityHandle ms_of_nodeset_ms;
  MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, ms_of_nodeset_ms);

    //add some entities to it
  temp_vec.clear();
  iter = nodes_of_quads.begin();
  end_iter = nodes_of_quads.end();
  for(; iter != end_iter; iter++)
  {
    if( ID_FROM_HANDLE( *iter ) % 2 )  //add all odds
      temp_vec.push_back( *iter );
  }
  MB->add_entities( ms_of_nodeset_ms, &temp_vec[0], temp_vec.size() ); 

    //add the other meshset to the nodeset's meshset
  MB->add_entities( nodeset_ms, &ms_of_nodeset_ms, 1);


    // no need to get lists, write out the whole mesh
  file_name = "test/mb_write2.g";
  std::vector<MBEntityHandle> output_list;
  output_list.push_back( block_ms );
  output_list.push_back( sideset_ms );
  output_list.push_back( meshset_abc );
  output_list.push_back( nodeset_ms );
  output_list.push_back( block_of_shells );
  MBErrorCode error = MB->write_mesh(file_name.c_str(), &output_list[0], output_list.size());

  return error;
}


MBErrorCode mb_higher_order_test(MBInterface *MB)
{
  double nodes_array [7][3];

  nodes_array[0][0] = 0;
  nodes_array[0][1] = 0;
  nodes_array[0][2] = 0;
  nodes_array[1][0] = 2;
  nodes_array[1][1] = 0; 
  nodes_array[1][2] = 0; 
  nodes_array[2][0] = 1; 
  nodes_array[2][1] = 2; 
  nodes_array[2][2] = 1; 
  nodes_array[3][0] = 1; 
  nodes_array[3][1] = 0; 
  nodes_array[3][2] = 0; 
  nodes_array[4][0] = 1.5; 
  nodes_array[4][1] = 0.5; 
  nodes_array[4][2] = 0.5; 
  nodes_array[5][0] = 0.5; 
  nodes_array[5][1] = 0.5; 
  nodes_array[5][2] = 0.5; 
  nodes_array[6][0] = 1; 
  nodes_array[6][1] = 1; 
  nodes_array[6][2] = 0.5; 


    //create the nodes
  std::vector<MBEntityHandle> connectivity(8);
  MBEntityHandle node_handle;
  int i;
  for( i=0; i<7; i++)
  {
    if(MB->create_vertex( nodes_array[i], node_handle ) != MB_SUCCESS )
      return MB_FAILURE;
    connectivity[i] = node_handle;
  }

    //create the higher order tri
  MBEntityHandle tri_handle;
  MBErrorCode result = MB->create_element(MBTRI, &connectivity[0], 6, tri_handle);
  if(result != MB_SUCCESS)
    return result;

    //create the higher order tri
  std::vector<MBEntityHandle> other_conn(3);

  double other_nodes[3][3];
  other_nodes[0][0] = 1.999;
  other_nodes[0][1] = 1.999;
  other_nodes[0][2] = 1.999;
  other_nodes[1][0] = 2.999; 
  other_nodes[1][1] = 2.999;
  other_nodes[1][2] = 2.999;
  other_nodes[2][0] = 3.999;
  other_nodes[2][1] = 3.999;
  other_nodes[2][2] = 3.999;

  for( i=0; i<3; i++)
  {
    if(MB->create_vertex( other_nodes[i], node_handle ) != MB_SUCCESS )
      return MB_FAILURE;
    other_conn[i] = node_handle;
  }

  MBEntityHandle other_tri_handle;
  result = MB->create_element(MBTRI, &other_conn[0], 3, other_tri_handle);
  if(result != MB_SUCCESS)
    return result;

    //get the connectivity now
  std::vector<MBEntityHandle> retrieved_conn; 

  result = MB->get_connectivity(&tri_handle, 1, retrieved_conn) ;
  if(result != MB_SUCCESS)
    return result;

  unsigned int k;
  for( k=0; k< retrieved_conn.size(); k++)
    if( connectivity[k] != retrieved_conn[k] )
      return MB_FAILURE;

  result = MB->get_connectivity(&other_tri_handle, 1, retrieved_conn) ;
  if(result != MB_SUCCESS)
    return result;

  for( k=0; k< other_conn.size(); k++)
    if( other_conn[k] != retrieved_conn[k] )
      return MB_FAILURE;

    // now let's just try getting the topological connectivity (the 3 corner vertices)
  std::vector<MBEntityHandle> topo_conn; 
  result = MB->get_connectivity(&other_tri_handle, 1, topo_conn, true) ;
  if(result != MB_SUCCESS)
    return result;

  if (topo_conn.size() != 3)
    return MB_FAILURE;

  for ( k=0; k<3; k++)
    if (topo_conn[k] != retrieved_conn[k] )
      return MB_FAILURE;

    // short check to make sure that MBCore::handle_from_id() works
  unsigned long handle_id = MB->id_from_handle( node_handle); 

  MBEntityHandle test_handle; 
  result = MB->handle_from_id( MBVERTEX, handle_id, test_handle ); 
  if(result != MB_SUCCESS)
    return result;

  if( test_handle != node_handle )
    return MB_FAILURE; 


  handle_id = MB->id_from_handle( tri_handle); 

  result = MB->handle_from_id( MBTRI, handle_id, test_handle ); 
  if(result != MB_SUCCESS)
    return result;

  if( test_handle != tri_handle )
    return MB_FAILURE; 


    //make up some bogus id 
  handle_id = 2140824; 

  result = MB->handle_from_id( MBTRI, handle_id, test_handle );
  if (result != MB_ENTITY_NOT_FOUND ) 
    return MB_FAILURE; 

  result = MB->handle_from_id( MBVERTEX, handle_id, test_handle );
  if (result != MB_ENTITY_NOT_FOUND ) 
    return MB_FAILURE; 


  return MB_SUCCESS;

}

MBErrorCode mb_bit_tags_test(MBInterface* MB)
{

  MBTag bit_tag;
  MBRange entities;
  MB->get_entities_by_type(0, MBVERTEX, entities);
  MBErrorCode success = MB_SUCCESS;

  if(MB->tag_create("bit on vertex", 3, MB_TAG_BIT, bit_tag, NULL) != MB_SUCCESS)
  {
    cout << "couldn't create bit tag" << endl;
    return MB_FAILURE;
  }

  MBRange::iterator iter;
  unsigned char bits;
  for(iter = entities.begin();
      iter != entities.end(); ++iter)
  {
      // tag each vertex with the low 3 bits of the entity handle
    bits = ((*iter) & 0x7);
    success = MB->tag_set_data(bit_tag, &(*iter), 1, &bits);
    if(success != MB_SUCCESS)
      return MB_FAILURE;
  }

  bits = 0;
  for(iter = entities.begin();
      iter != entities.end(); ++iter)
  {
      // tag each vertex with the low 3 bits of the entity handle
    success = MB->tag_get_data(bit_tag, &(*iter), 1, &bits);
    if(success != MB_SUCCESS)
      return MB_FAILURE;

    if(bits != ((*iter) & 0x7))
      return MB_FAILURE;
  }

  return MB_SUCCESS;
}

MBErrorCode mb_tags_test(MBInterface *MB)
{

  MBTag stale_bits, stale_dense, stale_sparse;
  MBErrorCode result = MB->tag_create("stale data", 5, MB_TAG_BIT, stale_bits, NULL);
  if (MB_SUCCESS != result)
    return result;
     
  int def_data = 9;
  result = MB->tag_create("dense stale_data", sizeof(int), MB_TAG_DENSE, stale_dense, &def_data);
  if (MB_SUCCESS != result)
    return result;
  result = MB->tag_create("sparse stale data", sizeof(int), MB_TAG_SPARSE, stale_sparse, NULL);
  if (MB_SUCCESS != result)
    return result;

  double coords[3] = { 0,0,0 };
  MBEntityHandle stale_handle1, stale_handle2;
  result = MB->create_vertex( coords, stale_handle1 );
  if (MB_SUCCESS != result)
    return result;

  unsigned char bits = 0x5;
  result = MB->tag_set_data(stale_bits, &stale_handle1, 1, &bits);
  if (MB_SUCCESS != result)
    return result;
  bits = 0;
  result = MB->tag_get_data(stale_bits, &stale_handle1, 1, &bits);
  if (MB_SUCCESS != result)
    return result;
  if (bits != 0x5)
    return MB_FAILURE;
    
  def_data = 1;
  result = MB->tag_set_data(stale_dense, &stale_handle1, 1, &def_data);
  if (MB_SUCCESS != result)
    return result;
  def_data = 0;
  result = MB->tag_get_data(stale_dense, &stale_handle1, 1, &def_data);
  if (MB_SUCCESS != result)
    return result;
  if (def_data != 1)
    return MB_FAILURE;
    
  def_data = 100;
  result = MB->tag_set_data(stale_sparse, &stale_handle1, 1, &def_data);
  if (MB_SUCCESS != result)
    return result;
  def_data = 0;
  result = MB->tag_get_data(stale_sparse, &stale_handle1, 1, &def_data);
  if (MB_SUCCESS != result)
    return result;
  if (def_data != 100)
    return MB_FAILURE;

  result = MB->delete_entities(&stale_handle1, 1);
  if (MB_SUCCESS != result)
    return result;
  result = MB->create_vertex(coords, stale_handle2);
  if (MB_SUCCESS != result)
    return result;

  if(stale_handle1 != stale_handle2)
    cout<< "Tag test could test stale data" << endl;
  else
  {
    bits=0;
    result = MB->tag_get_data(stale_bits, &stale_handle2, 1, &bits);
    if (MB_SUCCESS != result)
      return result;
    if(bits != 0)
      return MB_FAILURE;

    def_data = 3;
    result = MB->tag_get_data(stale_dense, &stale_handle2, 1, &def_data);
    if (MB_SUCCESS != result)
      return result;
    if(def_data != 9)
      return MB_FAILURE;

    def_data = 3;
    MBErrorCode stale_result = MB->tag_get_data(stale_sparse, &stale_handle2, 1, &def_data);
      // we are supposed to fail here
    if(stale_result != MB_TAG_NOT_FOUND)
      return MB_FAILURE;
  }

  result = MB->tag_delete(stale_dense);
  if (MB_SUCCESS != result)
    return result;
  
  result = MB->delete_entities(&stale_handle2, 1);
  if (MB_SUCCESS != result)
    return result;


    //get all blocks with material tag and with tag_value of 1 (should only be 1)
  MBRange entities;
  int value = 1;
  const void *dum_ptr = &value;
  MBTag material_tag;
  result = MB->tag_get_handle( MATERIAL_SET_TAG_NAME, material_tag);
  if (MB_SUCCESS != result)
    return result;
  if(MB->get_entities_by_type_and_tag( 0, MBENTITYSET, &material_tag, 
                                       &dum_ptr, 
                                       1, entities) != MB_SUCCESS)
    return MB_FAILURE;

  if( entities.size() != 1)
    return MB_FAILURE;
 
    //add a dense tag to hexes
  MBTag junk_tag;
  if(MB->tag_create( "junk_tag", sizeof(int), MB_TAG_DENSE, junk_tag, 0) 
     != MB_SUCCESS)
    return MB_FAILURE;    

    //Set the dense tag on 5 hexes to 3489 
  MBRange test_range;
  result = MB->get_entities_by_type(0,  MBHEX, test_range ) ;
  if(result != MB_SUCCESS)
    return result;

  MBRange::iterator iter, end_iter;
  iter = test_range.begin();
  end_iter = test_range.end();

  int data = 3489;
  const void *ptr_data = &data;

    //mark approxiamtely the first 20% of the hex entities; also put a bit tag on them
  unsigned int times = test_range.size()/5; 
  bits = 0x5;
  const void *ptr_bits = &bits;
  
  for(unsigned int i=0; i<times; i++) 
  {
    if(MB->tag_set_data( junk_tag, &(*iter), 1, &data ) != MB_SUCCESS ) 
      return MB_FAILURE;
    if(MB->tag_set_data( stale_bits, &(*iter), 1, &bits ) != MB_SUCCESS ) 
      return MB_FAILURE;
    iter++;
  }

  entities.clear();
    //fetch the hex entities of type--MBHEX, tag-"junk_tag", and tag value -- 3489
  if(MB->get_entities_by_type_and_tag(0, MBHEX, &junk_tag, 
                                      &ptr_data, 
                                      1, entities ) != MB_SUCCESS)
    return MB_FAILURE;
 
  if( entities.size() != times)  //should get as many hexes as you perviously marked
    return MB_FAILURE;

    //fetch the hex entities of type--MBHEX, tag-"junk_tag", and tag value -- 3489
  entities.clear();
  if(MB->get_entities_by_type_and_tag(0, MBHEX, &stale_bits, 
                                      &ptr_bits, 1, entities ) != MB_SUCCESS)
    return MB_FAILURE;
 
  if( entities.size() != times)  //should get as many hexes as you perviously marked
    return MB_FAILURE;

    // test fetch by tag value again, this time limiting the results
    // to the contents of an entity set
  MBEntityHandle meshset;
  result = MB->create_meshset( MESHSET_SET, meshset );
  if (MB_SUCCESS != result)
    return result;
  result = MB->add_entities( meshset, test_range );
  if (MB_SUCCESS != result)
    return result;
  
  
    //fetch the hex entities of type--MBHEX, tag-"junk_tag", and tag value -- 3489
  entities.clear();
  result = MB->get_entities_by_type_and_tag(meshset, MBHEX, &junk_tag, 
                                      &ptr_data, 1, entities );
  if (MB_SUCCESS != result)
    return result;
 
  if( entities.size() != times)  //should get as many hexes as you perviously marked
    return MB_FAILURE;

    //fetch the hex entities of type--MBHEX, tag-"stale_bits", and tag value -- 0x5
  entities.clear();
  result = MB->get_entities_by_type_and_tag(meshset, MBHEX, &stale_bits, 
                                      &ptr_bits, 1, entities );
  if (MB_SUCCESS != result)
    return result;
 
  if( entities.size() != times)  //should get as many hexes as you perviously marked
    return MB_FAILURE;


    // now try the query with an empty meshset, expecting to get back
    // an empty MBRange

  result = MB->create_meshset( MESHSET_SET, meshset );
  if (MB_SUCCESS != result)
    return result;
    
  entities.clear();
  result = MB->get_entities_by_type_and_tag(meshset, MBHEX, &junk_tag, 
                                      &ptr_data, 1, entities );
  if (MB_SUCCESS != result)
    return result;
    
  if(!entities.empty())
    return MB_FAILURE;
  


  result = MB->tag_delete(stale_bits);
  if (MB_SUCCESS != result)
    return result;

  return MB_SUCCESS;
}

MBErrorCode mb_common_tag_test( MBTagType storage, MBInterface* mb )
{
  char tagname[64];
  sprintf( tagname, "t%d", rand() );
  
  MBTag tag;
  const MBEntityHandle def_val = ~(MBEntityHandle)0;
  MBErrorCode rval = mb->tag_create( tagname, 
                                     sizeof(MBEntityHandle), 
                                     storage, 
                                     MB_TYPE_HANDLE, 
                                     tag, 
                                     &def_val );
  if (MB_SUCCESS != rval)
    return rval;
  
  
  MBRange entities;
  mb->get_entities_by_handle( 0, entities );
  if (entities.empty())
    return MB_FAILURE;
  
    // set tag on every other entity to be the entities handle
  MBRange::const_iterator i;
  bool odd = true;
  for (i = entities.begin(); i != entities.end(); ++i, odd = !odd) {
    if (odd) {
      const MBEntityHandle h = *i;
      rval = mb->tag_set_data( tag, &h, 1, &h );
      if (MB_SUCCESS != rval)
        return rval;
    }
  }
  
    // check values on every entity -- expect default for every other entity
  odd = true;
  for (i = entities.begin(); i != entities.end(); ++i, odd = !odd) {
    MBEntityHandle val = 0;
    rval = mb->tag_get_data( tag, &*i, 1, &val );
    if (MB_SUCCESS != rval)
      return rval;
    
    if (odd) {
      if (val != *i)
        return MB_FAILURE;
    }
    else {
      if (val != def_val)
        return MB_FAILURE;
    }
  }
  
    // set tag values on all entities
  std::vector<MBEntityHandle> values( entities.size() );
  std::copy( entities.begin(), entities.end(), values.begin() );
  rval = mb->tag_set_data( tag, entities, &values[0] );
  if (MB_SUCCESS != rval)
    return rval;
  
    // check values on every entity -- expect default for every other entity
  for (i = entities.begin(); i != entities.end(); ++i) {
    MBEntityHandle val = 0;
    rval = mb->tag_get_data( tag, &*i, 1, &val );
    if (MB_SUCCESS != rval)
      return rval;
    if (val != *i)
      return MB_FAILURE;
  }
  
    // find each entity by tag value
  for (i = entities.begin(); i != entities.end(); ++i) {
    const MBEntityHandle h = *i;
    const MBEntityType type = mb->type_from_handle( h );
    const void* const tag_vals[] = { &h };
    MBRange result;
    rval = mb->get_entities_by_type_and_tag( 0, type, &tag, tag_vals, 1, result );
    if (MB_SUCCESS != rval)
      return rval;
    if (result.size() != 1)
      return MB_FAILURE;
    if (result.front() != h)
      return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}


MBErrorCode mb_dense_tag_test( MBInterface* mb )
{
  return mb_common_tag_test( MB_TAG_DENSE, mb );
}

MBErrorCode mb_sparse_tag_test( MBInterface* mb )
{
  return mb_common_tag_test( MB_TAG_SPARSE, mb );
}
  
// class to offset hex center nodes
class OffsetHexCenterNodes : public MBInterface::HONodeAddedRemoved
{
public:
  OffsetHexCenterNodes(MBInterface* mb, double x, double y, double z)
      : gMB(mb)
    { 
      mOffset[0] = x; mOffset[1] = y; mOffset[2] = z; 
    }
    
  ~OffsetHexCenterNodes(){}

  void node_added(MBEntityHandle node, MBEntityHandle)
    {
      gMB->get_coords(&node, 1, mCoords);
      mCoords[0] += mOffset[0];
      mCoords[1] += mOffset[1];
      mCoords[2] += mOffset[2];
      gMB->set_coords(&node, 1, mCoords);
    }

    //do nothing
  void node_removed( MBEntityHandle /*node*/) {}

private:
  MBInterface* gMB;
  double mCoords[3];
  double mOffset[3];
};

MBErrorCode mb_entity_conversion_test(MBInterface *MB)
{
  MBErrorCode error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;

    //read in a file so you have some data in the database
  std::string file_name = TestDir + "/mbtest3.g";
  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;

  MBRange entities;
  MBEntityHandle meshset;
  MB->create_meshset(MESHSET_SET, meshset);
  
  MB->get_entities_by_type(0, MBHEX, entities);
  MB->add_entities(meshset, entities);


  OffsetHexCenterNodes function_object(MB,0.07, 0.15, 0);

  MB->convert_entities(meshset, false, false, true, &function_object);

  file_name = "test/hex_mid_volume_nodes.g";
  error = MB->write_mesh(file_name.c_str());
  if (error != MB_SUCCESS)
    return error;

  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;

  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;
  
  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;





  file_name = TestDir + "/mbtest3.g";
  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;

  entities.clear();
  MB->get_entities_by_type(0, MBHEX, entities);

  MB->create_meshset(MESHSET_SET, meshset);
  MB->add_entities(meshset, entities);
  MB->convert_entities(meshset, true, true, true);

  file_name = "test/hex_mid_edge_face_vol_nodes.g";
  error = MB->write_mesh(file_name.c_str());
  if (error != MB_SUCCESS)
    return error;

  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;

  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;
  
  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;





  file_name = TestDir + "/mbtest3.g";
  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;

  entities.clear();
  MB->get_entities_by_type(0, MBVERTEX, entities);
  unsigned int original_num_nodes = entities.size();
  entities.clear();
  MB->get_entities_by_type(0, MBHEX, entities);

  MB->create_meshset(MESHSET_SET, meshset);
  MB->add_entities(meshset, entities);
  MB->convert_entities(meshset, true, false, false);

  file_name = "test/hex_mid_edge_nodes.g";
  error = MB->write_mesh(file_name.c_str());
  if (error != MB_SUCCESS)
    return error;
  
    // convert them back to hex8's
  MB->convert_entities(meshset, false, false, false);

  entities.clear();
  MB->get_entities_by_type(0, MBVERTEX, entities);
    // make sure the higher order nodes really were deleted
  if(entities.size() != original_num_nodes)
    return MB_FAILURE;
  
  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;
  
  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;
  
  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;




  file_name = TestDir + "/mbtest1.g";
  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;

  entities.clear();
  MB->get_entities_by_type(0, MBTET, entities);
  
  MB->create_meshset(MESHSET_SET, meshset);
  MB->add_entities(meshset, entities);
  MB->convert_entities(meshset, true, false, false);

  file_name = "test/tet_mid_edge_nodes.g";
  error = MB->write_mesh(file_name.c_str());
  if (error != MB_SUCCESS)
    return error;
  
  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;
  
  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;
  
  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;

  
  
  file_name = TestDir + "/mbtest1.g";
  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;

  entities.clear();
  MB->get_entities_by_type(0, MBTET, entities);

  MB->create_meshset(MESHSET_SET, meshset);
  MB->add_entities(meshset, entities);
  MB->convert_entities(meshset, false, true, false);

  file_name = "test/tet_mid_face_nodes.g";
  error = MB->write_mesh(file_name.c_str());
  if (error != MB_SUCCESS)
    return error;
  
  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;
  
  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;
  
  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;







  file_name = TestDir + "/mbtest1.g";
  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;

  entities.clear();
  MB->get_entities_by_type(0, MBTET, entities);

  MB->create_meshset(MESHSET_SET, meshset);
  MB->add_entities(meshset, entities);
  MB->convert_entities(meshset, true, true, false);

  file_name = "test/tet_mid_edge_face_nodes.g";
  error = MB->write_mesh(file_name.c_str());
  if (error != MB_SUCCESS)
    return error;
  
  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;
  
  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;
  
  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;

  
  
  
  file_name = TestDir + "/mbtest1.g";
  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;

    // delete all MBTRI's
  entities.clear();
  error = MB->get_entities_by_type(0, MBTRI, entities);
  if (MB_SUCCESS != error)
    return error;
  error = MB->delete_entities(entities);
  if (MB_SUCCESS != error)
    return error;

  entities.clear();
  error = MB->get_entities_by_type(0, MBTET, entities);
  if (MB_SUCCESS != error)
    return error;
  
    // skin the model
  for(MBRange::iterator tet_iter = entities.begin(); tet_iter != entities.end(); ++tet_iter)
  {
    std::vector<MBEntityHandle> adj(32);
    error = MB->get_adjacencies(&(*tet_iter), 1, 2, true, adj);
    if (MB_SUCCESS != error)
      return error;
    for(std::vector<MBEntityHandle>::iterator tri_iter = adj.begin();
        tri_iter != adj.end(); ++tri_iter)
    {
      std::vector<MBEntityHandle> up_adj(2);
      MB->get_adjacencies(&(*tri_iter), 1, 3, false, up_adj);
      if(up_adj.size() > 1) {
        error = MB->delete_entities(&(*tri_iter), 1);
        if (MB_SUCCESS != error)
          return error;
      }
    }
  }

    // create a meshset of the skin
  MBEntityHandle export_meshset;
  MB->create_meshset( MESHSET_SET, export_meshset);
  MBTag material_tag;
  MB->tag_get_handle(MATERIAL_SET_TAG_NAME, material_tag);
  int block_id = 100;
  MB->tag_set_data(material_tag, &export_meshset, 1, &block_id);
  entities.clear();
  MB->get_entities_by_type(0, MBTRI, entities);
    // remove the first few tri's for fun
  MBRange tmp_ents;
  tmp_ents.insert(*entities.begin());
  entities.erase(entities.begin());
  tmp_ents.insert(*entities.begin());
  entities.erase(entities.begin());
  tmp_ents.insert(*entities.begin());
  entities.erase(entities.begin());
  tmp_ents.insert(*entities.begin());
  entities.erase(entities.begin());

  MB->add_entities(export_meshset, entities);

    // convert the skin
  MB->convert_entities(export_meshset, true, true, false);

    // make sure our first few tri's were untouched
  std::vector<MBEntityHandle> conn(3);
  for(MBRange::iterator kter=tmp_ents.begin(); kter != tmp_ents.end(); ++kter)
  {
    MB->get_connectivity(&(*kter), 1, conn);
    if(conn.size() != 3)
      return MB_FAILURE;
  }

    // output the skin
  file_name = "test/tri_mid_edge_face_nodes.g";
  error = MB->write_mesh(file_name.c_str(), &export_meshset, 1);
  if (error != MB_SUCCESS)
    return error;

  MB->delete_entities(&export_meshset, 1);
  
  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;
  
    //read the skin back in
  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;

  entities.clear();
  MB->get_entities_by_type(0, MBVERTEX, entities);
    // must have 101 nodes
  if(entities.size() != 101)
    return MB_FAILURE;

  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;


  return MB_SUCCESS;
}

//! Build two Quads with two edges shared between them.  The
//! edges share the same nodes.  We should be able to get
//! adjacencies on the edges and get one (correct) quad.  We
//! should be able to get the edge adjacencies of the quads
//! and only get 4 (not 5) edges.
//!

MBErrorCode mb_forced_adjacencies_test(MBInterface *MB)
{
    //! first clean up any existing mesh.
  MBErrorCode error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;


    //! create 6 nodes, 2 quads and 8 edges.  Edge 4 is adjacent
    //! to quad 1 and edge 5 is adjacent to quad 2.
    //!
    //!  4       5       6
    //!  o--e7---o--e8---o
    //!  |       |       |
    //!  e3 q1 e4 e5 q2  e6
    //!  |       |       |
    //!  o--e1---o--e2---o
    //!  1       2       3
    //!
  double node_coord1[3] = {0., 0., 0.};
  double node_coord2[3] = {1., 0., 0.};
  double node_coord3[3] = {2., 0., 0.};
  double node_coord4[3] = {0., 1., 0.};
  double node_coord5[3] = {1., 1., 0.};
  double node_coord6[3] = {2., 1., 0.};

  MBEntityHandle node1, node2, node3, node4, node5, node6;
  error = MB->create_vertex(node_coord1, node1);
  if (error != MB_SUCCESS)
    return error;

  error = MB->create_vertex(node_coord2, node2);
  if (error != MB_SUCCESS)
    return error;

  error = MB->create_vertex(node_coord3, node3);
  if (error != MB_SUCCESS)
    return error;

  error = MB->create_vertex(node_coord4, node4);
  if (error != MB_SUCCESS)
    return error;

  error = MB->create_vertex(node_coord5, node5);
  if (error != MB_SUCCESS)
    return error;

  error = MB->create_vertex(node_coord6, node6);
  if (error != MB_SUCCESS)
    return error;

  std::vector<MBEntityHandle> conn(4);
    //! create the first quad
  MBEntityHandle              quad1;
  conn[0] = node1;
  conn[1] = node2;
  conn[2] = node5;
  conn[3] = node4;
  error = MB->create_element(MBQUAD, &conn[0], 4, quad1);
  if (error != MB_SUCCESS)
    return error;

    //! create the second quad
  MBEntityHandle              quad2;
  conn[0] = node2;
  conn[1] = node3;
  conn[2] = node6;
  conn[3] = node5;
  error = MB->create_element(MBQUAD, &conn[0], 4, quad2);
  if (error != MB_SUCCESS)
    return error;

    //! create the edges
  MBEntityHandle edge1;
  conn.resize(2);
  conn[0] = node1;
  conn[1] = node2;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge1);
  if (error != MB_SUCCESS)
    return error;

  MBEntityHandle edge2;
  conn[0] = node2;
  conn[1] = node3;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge2);
  if (error != MB_SUCCESS)
    return error;

  MBEntityHandle edge3;
  conn[0] = node1;
  conn[1] = node4;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge3);
  if (error != MB_SUCCESS)
    return error;

  MBEntityHandle edge4;
  conn[0] = node2;
  conn[1] = node5;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge4);
  if (error != MB_SUCCESS)
    return error;

  MBEntityHandle edge5;
  conn[0] = node2;
  conn[1] = node5;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge5);
  if (error != MB_SUCCESS)
    return error;

  MBEntityHandle edge6;
  conn[0] = node3;
  conn[1] = node6;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge6);
  if (error != MB_SUCCESS)
    return error;

  MBEntityHandle edge7;
  conn[0] = node4;
  conn[1] = node5;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge7);
  if (error != MB_SUCCESS)
    return error;

  MBEntityHandle edge8;
  conn[0] = node5;
  conn[1] = node6;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge8);
  if (error != MB_SUCCESS)
    return error;


    //! Edge 4 and 5 share the same nodes, but should be different entities
  if (edge4 == edge5)
    return MB_FAILURE;

    //! Now that the geometry is created start adding the adjacency information
  std::vector<MBEntityHandle> edge_adjacencies1(4);
  edge_adjacencies1[0] = edge1;
  edge_adjacencies1[1] = edge4;
  edge_adjacencies1[2] = edge7;
  edge_adjacencies1[3] = edge3;

    //! does this (should this) say anything about order of the edges around the
    //! quad?  Also, does this also build the edge->quad adjacency, or is that
    //! does with a separate call?
  error = MB->add_adjacencies(quad1, &edge_adjacencies1[0], edge_adjacencies1.size(), true);
  if (error != MB_SUCCESS)
    return error;

  std::vector<MBEntityHandle> edge_adjacencies2(4);
  edge_adjacencies2[0] = edge2;
  edge_adjacencies2[1] = edge6;
  edge_adjacencies2[2] = edge8;
  edge_adjacencies2[3] = edge5;
  error = MB->add_adjacencies(quad2, &edge_adjacencies2[0], edge_adjacencies2.size(), true);
  if (error != MB_SUCCESS)
    return error;

    //! now get the adjacencies of each quad.
  std::vector<MBEntityHandle> quad1_adjacencies;
  error = MB->get_adjacencies(&(quad1), 1, 1, false, quad1_adjacencies);
  if (error != MB_SUCCESS)
    return error;


  std::sort(quad1_adjacencies.begin(), quad1_adjacencies.end());
  std::sort(edge_adjacencies1.begin(), edge_adjacencies1.end());

  if (quad1_adjacencies != edge_adjacencies1)
    return MB_FAILURE;
  
  std::vector<MBEntityHandle> quad2_adjacencies;
  error = MB->get_adjacencies(&(quad2), 1, 1, false, quad2_adjacencies);
  if (error != MB_SUCCESS)
    return error;

  std::sort(quad2_adjacencies.begin(), quad2_adjacencies.end());
  std::sort(edge_adjacencies2.begin(), edge_adjacencies2.end());

  if (quad2_adjacencies != edge_adjacencies2)
    return MB_FAILURE;

    //! try getting the adjacency of edge1 (should be quad1)
  std::vector<MBEntityHandle> edge1_adjacencies;
  error = MB->get_adjacencies(&(edge1), 1, 2, false, edge1_adjacencies);
  if (error != MB_SUCCESS)
    return error;

    //! there should be only 1 entity adjacent to edge1
  if (edge1_adjacencies.size() != 1)
    return error;

    //! and that entity should be quad1
  if (edge1_adjacencies[0] != quad1)
    return error;

    //! try getting the adjacency of edge6 (should be none)
  std::vector<MBEntityHandle> edge6_adjacencies;
  error = MB->get_adjacencies(&(edge6), 1, 2, false, edge6_adjacencies);
  if (error != MB_SUCCESS)
    return error;

    //! there should be only 1 entity adjacent to edge1
  if (edge1_adjacencies.size() != 0)
    return error;

    //! Now seal up the "gap" caused by edges 4 and 5.  Remove edge5
    //! from the adjacencies of quad2 and add edge 4 to quad2.

  std::vector<MBEntityHandle> edge5_adjacencies(1, edge5);
  error = MB->remove_adjacencies(quad2, &edge5_adjacencies[0], edge5_adjacencies.size());
  if (error != MB_SUCCESS)
    return error;


  std::vector<MBEntityHandle> edge4_adjacencies(1, edge4);
  error = MB->add_adjacencies(quad2, &edge4_adjacencies[0], edge4_adjacencies.size(), true);
  
    //! get the adjacencies of edge4 and it should return both quads.
  std::vector<MBEntityHandle> quad_adjacencies;
  error = MB->get_adjacencies(&(edge4), 1, 2, false, quad_adjacencies);
  if (error != MB_SUCCESS)
    return error;

    //! there should be 2 entities adjacent to edge4
  if (quad_adjacencies.size() != 2)
    return error;

    //! and they should be quad1 and quad2.  Note that we are not saying anything
    //! about order in the array.
  if ( (quad_adjacencies[0] != quad1 || quad_adjacencies[1] != quad1) &&
       (quad_adjacencies[0] != quad2 || quad_adjacencies[1] != quad2) ) 
    return error;
  
    //! clean up on exit
  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;


  return MB_SUCCESS;
}

/*bool lessnodesZ(const MBEntityHandle entity_handle1, const MBEntityHandle entity_handle2) 
  {
  double coords1[3], coords2[3];
  gMB->get_coords(entity_handle1, coords1);
  gMB->get_coords(entity_handle2, coords2);

  return coords2[2] < coords1[2];
  }*/

/*void sort_verts(MBRange vertices) 
  {
  std::vector<MBEntityHandle> vert_vec(vertices.size());
  MBRange::const_iterator iter;
  for (iter = vertices.begin(); iter != vertices.end(); iter++) 
  vert_vec.push_back(*iter);
  vert_vec.sort(lessnodesZ);
  }*/
bool points_are_coincident(const double *first, const double *second)
{
  double diff[3];
  diff[2] = first[2] - second[2];
    //  if (diff[2] > 0.001) return false;

  diff[0] = first[0] - second[0];
  diff[1] = first[1] - second[1];

  double length = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
  if(fabs(length) < .001)
    return true;

  return false;
}
MBErrorCode find_coincident_nodes(MBInterface* gMB, MBRange vertices,
                                  std::vector< std::pair<MBEntityHandle,MBEntityHandle> > &coin_nodes)
{
  double first_coords[3], second_coords[3];
  MBRange::const_iterator iter, jter;
  std::pair<MBEntityHandle, MBEntityHandle> coincident_pair;
  MBErrorCode result;

  for (iter = vertices.begin(); iter != vertices.end(); iter++)
  {
    result = gMB->get_coords(&(*iter),1, first_coords);
    if (result != MB_SUCCESS)
      return result;

    for (jter = iter; jter != vertices.end(); jter++)
    {
      if (*iter != *jter)
      {
        result = gMB->get_coords(&(*jter), 1, second_coords);
        if (result != MB_SUCCESS)
          return result;

        if(points_are_coincident(first_coords, second_coords))
        {
          coincident_pair.first  = *iter;
          coincident_pair.second = *jter;
          coin_nodes.push_back(coincident_pair);
        }
      }
    }
  }
  return MB_SUCCESS;
}

MBErrorCode find_coincident_elements(MBInterface* gMB, MBRange entities, int num_nodes,
                                     std::vector< std::pair<MBEntityHandle,MBEntityHandle> > &coin)
{
  double coords1[8][3], coords2[8][3];
  MBRange::iterator iter, jter;
  std::vector<MBEntityHandle> conn(8);
  std::pair<MBEntityHandle, MBEntityHandle> coincident_pair;
  int i = 0,/* j = 0,*/ ii = 0;

  for(iter = entities.begin(); iter != entities.end(); iter++)
  {
      // Get the coordinates for the element corners.
    if(gMB->get_connectivity(&(*iter), 1, conn) != MB_SUCCESS)
      return MB_FAILURE;
    for(ii=0; ii<num_nodes; ii++)
    {
      if(gMB->get_coords(&(conn[ii]), 1, coords1[ii]) != MB_SUCCESS)
        return MB_FAILURE;
    }

    for(jter = iter; jter != entities.end(); jter++)
    {
      if(*iter != *jter)
      {
          // Elements should be the same sense to merge.
        if(gMB->get_connectivity(&(*jter), 1, conn) != MB_SUCCESS)
          return MB_FAILURE;

        for (int tq = 0; tq < num_nodes; tq++)
          if(gMB->get_coords(&(conn[tq]), 1,  coords2[tq]) != MB_SUCCESS)
            return MB_FAILURE;
          //        if(gMB->get_coords(&(conn[0]), 1,  coords2[0]) != MB_SUCCESS)
          //          return MB_FAILURE;

          // Find if first node is coincident before testing the rest.
        bool first = false;
        for(i=0; i<num_nodes; i++)
        {
          if(points_are_coincident(coords1[i], coords2[0])) {
	    /*	    cout <<"first("<<i<<",0) - ";
                        cout <<" coords1["<<i<<"] = ("
                        <<coords1[i][0]<<","
                        <<coords1[i][1]<<","
                        <<coords1[i][2]<<") ";
                        cout <<" coords2["<<0<<"] = ("
                        <<coords2[0][0]<<","
                        <<coords2[0][1]<<","
                        <<coords2[0][2]<<")\n";*/
            first = true;
            break;
          }
        }
          // TEST -- Find if second node is coincident before testing the rest.
        bool second = false;
        for(int t2=0; t2<num_nodes; t2++)
        {
          if(points_are_coincident(coords1[t2], coords2[1])) {
	    /*	    cout <<"second("<<t2<<",1) - ";
                        cout <<" coords1["<<t2<<"] = ("
                        <<coords1[t2][0]<<","
                        <<coords1[t2][1]<<","
                        <<coords1[t2][2]<<") ";
                        cout <<" coords2["<<1<<"] = ("
                        <<coords2[1][0]<<","
                        <<coords2[1][1]<<","
                        <<coords2[1][2]<<")\n";*/
            second = true;
            break;
          }
        }
          // TEST -- Find if second node is coincident before testing the rest.
        bool third = false;
        for(int ti=0; ti<num_nodes; ti++)
        {
          if(points_are_coincident(coords1[ti], coords2[2])) {
	    /*	    cout <<"third("<<ti<<",2) - ";
                        cout <<" coords1["<<ti<<"] = ("
                        <<coords1[ti][0]<<","
                        <<coords1[ti][1]<<","
                        <<coords1[ti][2]<<") ";
                        cout <<" coords2["<<1<<"] = ("
                        <<coords2[2][0]<<","
                        <<coords2[2][1]<<","
                        <<coords2[2][2]<<")\n";*/
            third = true;
            break;
          }
        }
        if ((first)&&(second)&&(third)) {
          cout <<"i = "<<i<<"\n";
          for (int tii = 0; tii < num_nodes; tii++) {
            cout <<" coords1["<<tii<<"] = ("
                 <<coords1[tii][0]<<","
                 <<coords1[tii][1]<<","
                 <<coords1[tii][2]<<") ";
            cout <<" coords2["<<tii<<"] = ("
                 <<coords2[tii][0]<<","
                 <<coords2[tii][1]<<","
                 <<coords2[tii][2]<<")\n";
          }
        }

        if(i < num_nodes)
        {
          for(ii=1; ii<num_nodes; ii++)
          {
            if(gMB->get_coords(&(conn[ii]), 1,  coords2[ii]) != MB_SUCCESS)
              return MB_FAILURE;
          }
	  /*
              for(j=1; j<num_nodes; j++)
              {
	    
              if(!points_are_coincident(coords1[j], coords2[(j+i)%num_nodes]))
              break;
              }
              if(j == num_nodes)*/
          if ((first)&&(second)&&(third))
          {
            coincident_pair.first  = *iter;
            coincident_pair.second = *jter;
            coin.push_back(coincident_pair);
          }
        }
      }
    }
  }

  return MB_SUCCESS;
}


MBErrorCode mb_merge_test(MBInterface *MB)
{ 
  MB->delete_mesh();
  time_t begin_time = clock();
  unsigned int i;
  MBErrorCode result;
  MBSkinner MBSkinner_Obj(MB);

  std::string test_files[] = {std::string("test/cell1.gen"),
                              std::string("test/cell2.gen")};
    /*  			      std::string("cell3.gen"),
    			      std::string("cell4.gen"),
			      std::string("cell5.gen"),
			      std::string("cell6.gen"),
			      std::string("cell7.gen"),
			      std::string("cell8.gen"),
			      std::string("cell9.gen"),
			      std::string("cell10.gen"),
			      std::string("cell11.gen"),
			      std::string("cell12.gen"),
			      std::string("cell13.gen"),
			      std::string("cell14.gen"),
			      std::string("cell15.gen"),
			      std::string("cell16.gen"),
			      std::string("cell17.gen"),
			      std::string("cell18.gen"),
			      std::string("cell19.gen"),
			      std::string("cell20.gen"),
			      std::string("cell21.gen"),
			      std::string("cell22.gen"),
			      std::string("cell23.gen"),
			      std::string("cell24.gen")};*/

    /*std::vector<MBRange> entities(sizeof(test_files));
      std::vector<MBRange> forward_lower(sizeof(test_files));
      std::vector<MBRange> reverse_lower(sizeof(test_files));
      std::vector<MBRange> nodes(sizeof(test_files));*/
  MBRange entities;
  MBRange forward_lower;
  MBRange reverse_lower;
  MBRange faces;
  MBRange nodes;

  cout << "---Starting Merge Tests---" << endl << endl;
  for(i=0; i<(sizeof(test_files)/sizeof(std::string)); i++)
  {

    cout << "---Testing:\"" << test_files[i] << "\"---" << endl;
    result = MB->load_mesh(test_files[i].c_str(), NULL, 0);
    if (result == MB_SUCCESS) 
      cout <<"Loaded "<<test_files[i]<<"\n";
      //get Hexes from model
  }
  result = MB->get_entities_by_type(0, MBHEX, entities);
  MBSkinner_Obj.find_skin(entities,forward_lower,reverse_lower);
  cout <<"num hexes = "<<entities.size()<<"\n";
  cout <<"fl = "<<forward_lower.size()<<" rl = "<<reverse_lower.size()<<"\n";
  
    //  MBRange::const_iterator iter;
  int dim = 0;
    //  int num_ents = 1;
  result = MB->get_adjacencies(forward_lower, dim, true, nodes, MBInterface::UNION);
  cout <<"nodes.size() = "<<nodes.size() <<"\n";
    
  if(result == MB_SUCCESS)
    cout << "---Success---";
  else
    cout << "---Failure---";
  cout << endl << endl;

    //  result = MB->get_entities_by_type(0, MBQUAD, faces);
    //  cout <<"num faces = "<<faces.size() <<"\n";

  std::vector<std::pair<MBEntityHandle, MBEntityHandle> > coin_nodes;
    //  cout <<"Begining sort...\n";
    //  std::sort(nodes.begin(),nodes.end(),lessnodesZ);
    //  cout <<"Ending sort...\n";
  result = find_coincident_nodes(MB,nodes, coin_nodes);
  cout <<"coin_nodes.size() = "<<coin_nodes.size() <<"\n";
  std::vector< std::pair<MBEntityHandle, MBEntityHandle> >::iterator n_iter;
  for (n_iter=coin_nodes.begin(); n_iter != coin_nodes.end(); n_iter++) {
    result = MB->merge_entities((*n_iter).first, (*n_iter).second, false, true);
    if (MB_SUCCESS != result)
      return result;
  }
    /*  std::vector<std::pair<MBEntityHandle, MBEntityHandle> > coin_faces;
        int nodes_per_elt = 4;
        result = find_coincident_elements(forward_lower, nodes_per_elt, coin_faces);
        if (result != MB_SUCCESS) cout <<"find_coincident_elements fail!\n";
        cout <<"coin_faces.size() = "<<coin_faces.size() <<"\n";
        std::vector< std::pair<MBEntityHandle, MBEntityHandle> >::iterator f_iter;
        for (f_iter=coin_faces.begin(); f_iter != coin_faces.end(); f_iter++)
        MB->merge_entities((*f_iter).first, (*f_iter).second, true, true);*/
    /*
      std::vector<std::pair<MBEntityHandle, MBEntityHandle> > coin_fl;
      nodes_per_elt = 4;
      result = find_coincident_elements(entities, nodes_per_elt, coin_fl);
      cout <<"coin_fl.size() = "<<coin_fl.size() <<"\n";
    */
  int num_ents;
  if ((MB_SUCCESS == MB->get_number_entities_by_dimension(0, 3, num_ents) &&
       0 != num_ents) ||
      (MB_SUCCESS == MB->get_number_entities_by_dimension(0, 2, num_ents) &&
       0 != num_ents))
    result = MB->write_mesh("test/merge_test.ncdf");
  ;
  

  double clocks_per_sec = (double) CLOCKS_PER_SEC;
  double real_time =  difftime(time(NULL), begin_time);
  cout <<"TIME: "<<(real_time/clocks_per_sec)<<" seconds.\n";
  return result;
}

MBErrorCode mb_merge_update_test(MBInterface*)
{
  MBCore moab;
  MBInterface* mb = &moab;
  MBErrorCode rval;
  
    // create two quads with a coincident edge pair
  double coords[] = { 0, 0, 0,
                      1, 0, 0,
                      1, 1, 0,
                      0, 1, 0,
                      1, 1, 0,
                      1, 0, 0,
                      2, 0, 0,
                      2, 1, 0 };
  MBEntityHandle verts[8];
  for (int i = 0; i < 8; ++i) 
    mb->create_vertex( coords + 3*i, verts[i] );
  MBEntityHandle quad1, quad2, edge1, edge2;
  mb->create_element( MBQUAD, verts, 4, quad1 );
  mb->create_element( MBQUAD, verts+4, 4, quad2 );
  mb->create_element( MBEDGE, verts+1, 2, edge1 );
  mb->create_element( MBEDGE, verts+4, 2, edge2 );
  
    // create two tracking sets containing the vertices
    // and edge of each quad
  MBEntityHandle set1, set2;
  mb->create_meshset( MESHSET_TRACK_OWNER, set1 );
  mb->create_meshset( MESHSET_TRACK_OWNER, set2 );
  mb->add_entities( set1, verts, 4 );
  mb->add_entities( set2, verts+4, 4 );
  mb->add_entities( set1, &edge1, 1 );
  mb->add_entities( set2, &edge2, 1 );
  
    // now merge the coincident edges
  rval = mb->merge_entities( verts[1], verts[5], false, true );
  if (MB_SUCCESS != rval) {
    std::cerr << "Merge failed at " << __FILE__ << ":" << __LINE__ << std::endl;
    return rval;
  }
  rval = mb->merge_entities( verts[2], verts[4], false, true );
  if (MB_SUCCESS != rval) {
    std::cerr << "Merge failed at " << __FILE__ << ":" << __LINE__ << std::endl;
    return rval;
  }
  rval = mb->merge_entities( edge1, edge2, false, true );
  if (MB_SUCCESS != rval) {
    std::cerr << "Merge failed at " << __FILE__ << ":" << __LINE__ << std::endl;
    return rval;
  }
  
    // check that there is only one edge and that it has the correct connectivity
  MBRange r;
  mb->get_entities_by_type( 0, MBEDGE, r );
  if (r.size() != 1 || r.front() != edge1) {
    std::cerr << "Edge merge failed at " << __FILE__ << ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  std::vector<MBEntityHandle> exp(verts+1, verts+3), act;
  mb->get_connectivity( &edge1, 1, act );
  if (exp != act) {
    std::cerr << "Incorrect conn for edge at " << __FILE__ << ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  
    // check that quad connectivity is as expected
  exp = std::vector<MBEntityHandle>(verts, verts+4);
  act.clear();
  mb->get_connectivity( &quad1, 1, act );
  if (exp != act) {
    std::cerr << "Incorrect conn for quad at " << __FILE__ << ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  exp.resize(4);
  exp[0] = verts[2];
  exp[1] = verts[1];
  exp[2] = verts[6];
  exp[3] = verts[7];
  act.clear();
  mb->get_connectivity( &quad2, 1, act );
  if (exp != act) {
    std::cerr << "Incorrect conn for quad at " << __FILE__ << ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  
    // check that set contents are correctly updated
  exp = std::vector<MBEntityHandle>(verts, verts+4);
  exp.push_back( edge1 );
  act.clear();
  mb->get_entities_by_handle( set1, act );
  std::sort( exp.begin(), exp.end() );
  std::sort( act.begin(), act.end() );
  if (exp != act) {
    std::cerr << "Incorrect set contents at " << __FILE__ << ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  
  exp.resize(5);
  exp[0] = verts[2];
  exp[1] = verts[1];
  exp[2] = verts[6];
  exp[3] = verts[7];
  exp[4] = edge1;
  act.clear();
  mb->get_entities_by_handle( set2, act );
  std::sort( exp.begin(), exp.end() );
  std::sort( act.begin(), act.end() );
  if (exp != act) {
    std::cerr << "Incorrect set contents at " << __FILE__ << ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}
  
  

MBErrorCode mb_stress_test(MBInterface *MB)
{
    // delete the existing mesh
  MBErrorCode error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;

  cout << "    Beginning Stress Test . . ." << endl;
  cout << "\n        Reading elements" << endl;
  clock_t start = clock();
  clock_t total_start = clock();

    //read in a file so you have some data in the database
  std::string file_name = "test/mb_big_test.g";
  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;

  clock_t stop = clock();

  int num_entities;
  error = MB->get_number_entities_by_type(0, MBHEX, num_entities);
  if (error != MB_SUCCESS)
    return error;

  if (num_entities != 256000)
    return error;

  float time = static_cast<float>(stop - start)/CLOCKS_PER_SEC;
  float speed = num_entities/time;
  cout << "        Read " << num_entities << " entities"                
       << " in "   << time << " seconds" << endl; 
  cout << "        at " << speed << " elements per second." << endl;


  cout << "\n        Transforming and copying elements" << endl;
  start = clock();

  MBRange hexes;
  error = MB->get_entities_by_type(0, MBHEX, hexes);
  if (error != MB_SUCCESS)
    return error;


  std::vector<MBEntityHandle> conn;
  MBRange::iterator iter;
  for (iter = hexes.begin(); iter != hexes.end(); iter++)
  {
    error = MB->get_connectivity(&(*iter), 1, conn);
    if (error != MB_SUCCESS)
      return error;

    double coords[3];
    int i = 0;
    std::vector<MBEntityHandle> vertex_handle(8);
    MBEntityHandle element_handle;
    std::vector<MBEntityHandle>::iterator jter;

    for (jter = conn.begin(); jter != conn.end(); jter++)
    {
      error = MB->get_coords(&(*jter), 1, coords);
      if (error != MB_SUCCESS)
        return error;
      coords[2] += 20.0;
      error = MB->create_vertex(coords, vertex_handle[i++]);
      if (error != MB_SUCCESS)
        return error;
    }
    error = MB->create_element(MBHEX, &vertex_handle[0], 8, element_handle);
    if (error != MB_SUCCESS)
      return error;
  }


  stop = clock();
  time = static_cast<float>(stop - start)/CLOCKS_PER_SEC;

  cout << "        Transformed and created " << num_entities << " entities" 
       << " in "   << time << " seconds" << endl; 

    // Create mesh set
  cout << "\n        Creating meshset" << endl;
  start = clock();
  error = MB->get_entities_by_type(0, MBHEX, hexes);
  if (error != MB_SUCCESS)
    return error;

  if (hexes.size() != 512000)
    return MB_FAILURE;

  MBEntityHandle mesh_set;
  error = MB->create_meshset( MESHSET_SET, mesh_set );
  if (error != MB_SUCCESS)
    return error;

  error = MB->add_entities(mesh_set, hexes);
  if (error != MB_SUCCESS)
    return error;

  stop = clock();
  time = static_cast<float>(stop - start)/CLOCKS_PER_SEC;

  cout << "        Created meshset with " << hexes.size() << " entities" 
       << " in "   << time << " seconds" << endl; 

  cout << "\n        Writing 512K element file . . ." << endl;
  start = clock();

    // set the block tag
  MBTag tag_handle;
  MBErrorCode result = MB->tag_get_handle( MATERIAL_SET_TAG_NAME, tag_handle ) ;
  if(result != MB_SUCCESS)
    return result;

  int id = 1;
  result = MB->tag_set_data( tag_handle, &mesh_set, 1, &id ) ;
  if(result != MB_SUCCESS)
    return result;

  std::vector<MBEntityHandle> output_list;
  output_list.push_back(mesh_set);

  file_name = "test/mb_stress_out.g";
  error = MB->write_mesh(file_name.c_str(), &output_list[0], output_list.size());
  if (error != MB_SUCCESS)
    return error;

  stop = clock();
  time = static_cast<float>(stop - start)/CLOCKS_PER_SEC;

  cout << "        Wrote file with " << hexes.size() << " entities" 
       << " in "   << time << " seconds" << endl; 

  clock_t total_stop = clock();
  time = static_cast<float>(total_stop - total_start)/CLOCKS_PER_SEC;

  cout << "        Total time: " << time << " seconds." << endl;

  MB->delete_mesh();
  
  return MB_SUCCESS;
}

MBErrorCode mb_canon_number_test(MBInterface *MB) 
{
    // various tests for canonical ordering

    // MBCN::AdjacentSubEntities
  std::vector<int> vec1, vec2;
  MBErrorCode result;

  MBEntityType this_type;

  for (this_type = MBEDGE; this_type != MBKNIFE; this_type++) {
    
    for (int i = 0; i < MBCN::VerticesPerEntity(this_type); i++) {
        // test for edges and faces
      for (int dim = 1; dim <= MBCN::Dimension(this_type); dim++) {
          // get the sides adjacent to this vertex
        vec1.clear();
        int temp_result = MBCN::AdjacentSubEntities(this_type, &i, 1, 0, dim, vec1);
                                     
        if (0 != temp_result ||
            vec1.size() > (unsigned int) MBCN::NumSubEntities(this_type, dim)) {
          cout << "failed getting sides for type " << MBCN::EntityTypeName(this_type)
               << " dimension" << dim << endl;
          return MB_FAILURE;
        }
        

          // now get the vertices shared by these sides
        vec2.clear();
        temp_result = 
          MBCN::AdjacentSubEntities(this_type, &vec1[0], vec1.size(), dim, 0,
                                    vec2);

          // vertex side recovered should be i
        if (0 != temp_result || 
              // if dimension is same as DIMENSION(this_type), will get back all the
              // vertices in the entity
            (dim == MBCN::Dimension(this_type) && 
             vec2.size() != (unsigned int) MBCN::VerticesPerEntity(this_type)) ||
              // otherwise, we should get back only one vertex, and it should be the one
              // we started with
            (dim != MBCN::Dimension(this_type) &&
             (vec2.size() != 1 || vec2[0] != i))) {
          cout << "traversal from verts to sides to verts failed for " << endl
               << "vertex " << i << " type " << MBCN::EntityTypeName(this_type) 
               << " dimension " << dim << endl;
          return MB_FAILURE;
        }
      }
    }
  }
  
    // MBCN::side_number

    // create vertices to use later
  double xyz[3] = {0.0, 0.0, 0.0};
  MBEntityHandle vertex_handles[8];
  for (int i = 0; i < 8; i++) {
    result = MB->create_vertex(xyz, vertex_handles[i]);
    assert(result == MB_SUCCESS);
  }
  int side, sense, offset;
  
  MBEntityHandle this_entity;
  
  for (this_type = MBEDGE; this_type != MBKNIFE; this_type++) {
    
      // skip remainder of the loop for POLYGONS and POLYHEDRA, which don't follow
      // the standard canonical ordering 
    if (this_type == MBPOLYGON || this_type == MBPOLYHEDRON)
      continue;
    
      // make an entity of this type
    result = MB->create_element(this_type, vertex_handles, 
                                MBCN::VerticesPerEntity(this_type),
                                this_entity);
    if (MB_SUCCESS != result || 0 == this_entity) {
      cout << "failed to create entity of type " 
           << MBCN::EntityTypeName(this_type) << endl;
      return MB_FAILURE;
    }

      // now get the connectivity vector *
    const MBEntityHandle *entity_vertices;
    int num_verts;
    result = MB->get_connectivity(this_entity, entity_vertices, num_verts);
    if (MB_SUCCESS != result || 
        num_verts != MBCN::VerticesPerEntity(this_type)) {
      cout << "failed to get connectivity for entity type " 
           << MBCN::EntityTypeName(this_type) << endl;
      return MB_FAILURE;
    }
    
      // for each dimension
    for (int dim = 1; dim <= MBCN::Dimension(this_type); dim++) {
        // for each side of this dimension
      const MBCN::ConnMap &cm = MBCN::mConnectivityMap[this_type][dim-1];
      int tmp_conn[MB_MAX_SUB_ENTITY_VERTICES];

      for (int side_no = 0; side_no < MBCN::NumSubEntities(this_type, dim); side_no++) {

        for (int j = 0; j < MB_MAX_SUB_ENTITY_VERTICES; j++) tmp_conn[j] = cm.conn[side_no][j];
        int temp_result = 
          MBCN::SideNumber(this_type,
                           tmp_conn, 
                           MBCN::VerticesPerEntity(MBCN::SubEntityType(this_type, dim, side_no)),
                           dim, side, sense, offset);
        if (0 != temp_result) {
          cout << "call to MBCN::side_number failed with non-success result"
               << " for type " 
               << MBCN::EntityTypeName(this_type) << " dimension " << dim 
               << " side no " << side_no << endl;
          return MB_FAILURE;
        }

          // side number should be the same as side_no, sense should be forward, offset should
          // be zero
        if (side != side_no || sense != 1 || offset != 0) {
          cout << "call to MBCN::side_number failed for type " 
               << MBCN::EntityTypeName(this_type) << " dimension " << dim << " side no " 
               << side_no << endl
               << "side, sense, offset = " << side << " " << sense << " " << offset << endl;
          return MB_FAILURE;
        }
      }
    }

      // destroy the entity of this_type
    result = MB->delete_entities(&this_entity, 1);
    if (MB_SUCCESS != result)
      return result;
  }

  return MB_SUCCESS;
}

MBErrorCode mb_poly_test(MBInterface *mb) 
{
    // test polygon and polyhedron representation
    // create a couple of polygons; vertices first
  const double vert_pos[48] = {
    -1, 0, 0,
    1, 0, 0, 
    2, 0, 0,
    2, 1, 0,
    1, 1, 0,
    0, 2, 0,
    -1, 1, 0,
    -2, 1, 0,
    -2, 0, 0,
    -2, -1, 0,
    -1, -1, 0,
    -1, -2, 0,
    1, -2, 0,
    1, -1, 0,
    2, -1, 0,
    1.5, .5, 1};

  MBEntityHandle verts[16];
  MBErrorCode result;
  int i;
  for (i = 0; i < 16; i++) {
    result = mb->create_vertex(&vert_pos[3*i], verts[i]);
    if (MB_SUCCESS != result) {
      std::cout << "Failed to create vertex " << i << " in mb_poly_test." << std::endl;
      return result;
    }
  }

    // then polygons
  const int connect_idx[] = 
    {
      1, 2, 5, 6, 7,
      2, 3, 4, 5,
      5, 4, 6,
      7, 6, 8,
      0, 1, 7, 8,
      0, 1, 2, 3, 14, 13, 12, 11, 10, 9,
      0, 9, 10, 11, 12, 13, 14, 3, 4, 6, 8,
      2, 3, 15, 5,
      3, 4, 5, 15
    };
  
  int num_connect_idx[] = {5, 4, 3, 3, 4, 10, 11, 4, 4};
  
  MBEntityHandle polygons[9], temp_connect[12];
  int idx = 0, nump = 0;
  MBErrorCode tmp_result;
  while (nump < 9) {
    for (i = 0; i < num_connect_idx[nump]; i++)
      temp_connect[i] = verts[connect_idx[idx+i]];
      
    tmp_result = mb->create_element(MBPOLYGON, temp_connect, num_connect_idx[nump], 
                                    polygons[nump]);
    if (MB_SUCCESS != tmp_result) {
      std::cout << "mb_poly_test: create_element failed for polygon " << i << "." << std::endl;
      result = tmp_result;
      nump++;
      continue;
    }
    
    idx += num_connect_idx[nump];
    nump++;
  }

  if (MB_SUCCESS != result) return result;
  
    // ok, made 'em; now get all the vertices and make sure they're the same
  const MBEntityHandle *connect;
  int num_connect;
  idx = 0;
  int j;
  for (i = 0; i < 9; i++) {
    tmp_result = mb->get_connectivity(polygons[i], connect, num_connect);
    if (MB_SUCCESS != tmp_result || num_connect != num_connect_idx[i]) {
      std::cout << "mb_poly_test: get_connectivity test failed for polygon " << i << "." << std::endl;
      result = (tmp_result != MB_SUCCESS ? tmp_result : MB_FAILURE);
      continue;
    }

    for (j = 0; j < num_connect; j++) {
      if (connect[j] != verts[connect_idx[idx+j]]) {
        std::cout << "mb_poly_test: get_connectivity test returned wrong vertices for polygon " 
                  << i << "." << std::endl;
        result = MB_FAILURE;
        continue;
      }
    }
    
    idx += num_connect;
  }
  
  if (MB_SUCCESS != result) return result;

    // check a different way, with ranges
  MBRange vert_range, poly_range;
  for (i = 0; i < 9; i++) poly_range.insert(polygons[i]);
  result = mb->get_adjacencies(poly_range, 0, false, vert_range, 
                               MBInterface::UNION);
  if (MB_SUCCESS != result) {
    std::cout << "mb_poly_test: get_adjacencies failed for polygon " 
              << i << "." << std::endl;
    return result;
  }
  
  else if (vert_range.size() != 16) {
    std::cout << "mb_poly_test: get_adjacencies returned wrong # of vertices for polygon " 
              << i << "." << std::endl;
    return MB_FAILURE;
  }
  
    // make a couple polyhedra
  MBEntityHandle polyhedra[2];
  result = mb->create_element(MBPOLYHEDRON, polygons, 7, polyhedra[0]);
  if (MB_SUCCESS != result) {
    std::cout << "mb_poly_test: create_element failed for polyhedron 1." << std::endl;
    return result;
  }

  temp_connect[0] = polygons[1];
  temp_connect[1] = polygons[7];
  temp_connect[2] = polygons[8];
  result = mb->create_element(MBPOLYHEDRON, temp_connect, 3, polyhedra[1]);
  if (MB_SUCCESS != result) {
    std::cout << "mb_poly_test: create_element failed for polyhedron 2." << std::endl;
    return result;
  }
  
    // now look for vertices common to both
  std::vector<MBEntityHandle> temp_verts;
  result = mb->get_adjacencies(polyhedra, 2, 0, false, temp_verts);
  if (MB_SUCCESS != result) {
    std::cout << "mb_poly_test: get_adjacencies failed for polyhedra." << std::endl;
    return result;
  }

  if (4 != temp_verts.size()) {
    std::cout << "mb_poly_test: get_adjacencies for polyhedra returned " << temp_verts.size()
              << " vertices, should be 4." << std::endl;
    return MB_FAILURE;
  }
  
    // ok, we're probably fine
  return MB_SUCCESS;
}

MBErrorCode mb_range_test(MBInterface *) 
{

  MBRange r1, r2, rhs;
  MBErrorCode result = MB_SUCCESS;
  int err;

  MBEntityHandle h1 = CREATE_HANDLE(MBVERTEX, 1, err);
  MBEntityHandle h4 = CREATE_HANDLE(MBVERTEX, 4, err);
  MBEntityHandle h5 = CREATE_HANDLE(MBVERTEX, 5, err);
  MBEntityHandle h9 = CREATE_HANDLE(MBVERTEX, 9, err);
  MBEntityHandle h10 = CREATE_HANDLE(MBVERTEX, 10, err);
  MBEntityHandle h15 = CREATE_HANDLE(MBVERTEX, 15, err);
  MBEntityHandle h16 = CREATE_HANDLE(MBVERTEX, 16, err);
  MBEntityHandle h20 = CREATE_HANDLE(MBVERTEX, 20, err);

    // equal start/end
  r1.insert(h1, h5);
  r1.insert(h10, h15);
  r2.insert(h5, h10);
  r2.insert(h15, h20);
  
  rhs = intersect( r1, r2);
  if (rhs.size() != 3) {
    result = MB_FAILURE;
    std::cout << "Range test equal start/end failed." << std::endl;
  }
  
    // identical ranges test
  r1.clear(); r2.clear(); rhs.clear();
  r1.insert(h1, h5); r1.insert(h10, h20); 
  r2.insert(h1, h5); r2.insert(h10, h20); 
  rhs = intersect( r1, r2);
  if (rhs.size() != r1.size() || rhs.size() != r2.size()) {
    result = MB_FAILURE;
    std::cout << "Range test identical ranges failed." << std::endl;
  }
  
    // off by one test
  r1.clear(); r2.clear(); rhs.clear();
  r1.insert(h1, h4); r1.insert(h10, h15); 
  r2.insert(h5, h9); r2.insert(h16, h20);
  rhs = intersect( r1, r2);
  if (rhs.size() != 0) {
    result = MB_FAILURE;
    std::cout << "Range test off by one failed." << std::endl;
  }
  
    // interior test
  r1.clear(); r2.clear(); rhs.clear();
  r1.insert(h1, h20); 
  r2.insert(h5, h10); 
  rhs = intersect( r1, r2);
  if (rhs.size() != r2.size()) {
    result = MB_FAILURE;
    std::cout << "Range test interior failed." << std::endl;
  }
  
    // half-above test
  r1.clear(); r2.clear(); rhs.clear();
  r1.insert(h1, h10);
  r2.insert(h5, h20);
  rhs = intersect( r1, r2);
  if (rhs.size() != h10-h5+1) {
    result = MB_FAILURE;
    std::cout << "Range test half-above failed." << std::endl;
  }
  
    // half-below test
  r1.clear(); r2.clear(); rhs.clear();
  r1.insert(h5, h20);
  r2.insert(h1, h10);
  rhs = intersect( r1, r2);
  if (rhs.size() != h10-h5+1) {
    result = MB_FAILURE;
    std::cout << "Range test half-below failed." << std::endl;
  }
  
    // merge all test
  r1.clear();
  r2.clear();
  r1.insert(h1, h5);
  r1.insert(h10, h15);
  r2 = r1;
  r1.merge( r2 );
  if (r1.size() != r2.size()) {
    result = MB_FAILURE;
    std::cout << "Range merge test failed." << std::endl;
  }
  
    // merge subset test
  r1.clear();
  r2.clear();
  r1.insert(h1, h5);
  r1.insert(h10, h15);
  MBRange::const_iterator i1 = r1.begin();
  MBRange::const_iterator i2 = r1.end();
  i1 += 2;
  i2 -= 2;
  r2.merge( i1, i2 );
  if (r1.size() - 4 != r2.size()) {
    result = MB_FAILURE;
    std::cout << "Range merge subset test failed." << std::endl;
  }
  
    // const_pair_iterator test
  r1.clear();
  r1.insert( h1 );
  r1.insert( h4 );
  r1.insert( h5 );
  MBRange::const_pair_iterator pair_iter = r1.const_pair_begin();
  MBEntityHandle cpi_h1 = (*pair_iter).first;
  MBEntityHandle cpi_h2 = (*pair_iter).second;
  ++pair_iter;
  MBEntityHandle cpi_h3 = (*pair_iter).first;
  MBEntityHandle cpi_h4 = (*pair_iter).second;
  ++pair_iter;
  if (cpi_h1 != h1 || cpi_h2 != h1 || 
      cpi_h3 != h4 || cpi_h4 != h5 ||
      pair_iter != r1.const_pair_end()) {
    result = MB_FAILURE;
    std::cout << "Range::const_pair_iterator test failed.\n" << std::endl;
  }

    // swap tests
    // both non-empty
  r1.clear();
  r2.clear();
  r1.insert(h1);
  r1.insert(h4);
  r1.insert(h5);
  r2.insert(h9);
  r1.swap(r2);
  if (r1.size() != 1 || r2.size() != 3 ||
      r1.find(h9) == r1.end() ||
      r2.find(h1) == r2.end() ||
      r2.find(h4) == r2.end() ||
      r2.find(h5) == r2.end()) 
  {
    std::cout << "Range test: swap both non-empty failed." << std::endl;
    result = MB_FAILURE;
  }
  
    // destination empty
  r1.clear();
  r2.clear();
  r2.insert(h1);
  r2.insert(h4);
  r2.insert(h5);
  r1.swap(r2);
  if (r1.size() != 3 || r2.size() != 0 ||
      r1.find(h1) == r1.end() ||
      r1.find(h4) == r1.end() ||
      r1.find(h5) == r1.end()) 
  {
    std::cout << "Range test: swap destination empty failed." << std::endl;
    result = MB_FAILURE;
  }
  
    // source empty
  r1.clear();
  r2.clear();
  r1.insert(h1);
  r1.insert(h4);
  r1.insert(h5);
  r1.swap(r2);
  if (r1.size() != 0 || r2.size() != 3 ||
      r2.find(h1) == r2.end() ||
      r2.find(h4) == r2.end() ||
      r2.find(h5) == r2.end()) 
  {
    std::cout << "Range test: swap source empty failed." << std::endl;
    result = MB_FAILURE;
  }
  
    // both empty
  r1.clear();
  r2.clear();
  if (r1.size() != 0 || r2.size() != 0) 
  {
    std::cout << "Range test: swap both empty failed." << std::endl;
    result = MB_FAILURE;
  }

    // subtract tests
    // both non-empty
  r1.clear();
  r2.clear();
  r1.insert(h1);
  r1.insert(h4);
  r1.insert(h5);
  r2.insert(h4);
  MBRange r3 = subtract( r1, r2);
  if (r3.size() != 2 ||
      r3.find(h1) == r3.end() ||
      r3.find(h5) == r3.end() ||
      r3.find(h4) != r3.end()) 
  {
    std::cout << "Range test: subtract both non-empty failed." << std::endl;
    result = MB_FAILURE;
  }
  
    // destination empty
  r1.clear();
  r2.clear();
  r1.insert(h1);
  r1.insert(h4);
  r1.insert(h5);
  r3 = subtract( r1, r2);
  if (r3.size() != 3 ||
      r3.find(h1) == r3.end() ||
      r3.find(h4) == r3.end() ||
      r3.find(h5) == r3.end()) 
  {
    std::cout << "Range test: subtract destination empty failed." << std::endl;
    result = MB_FAILURE;
  }

  return result;
}

MBErrorCode mb_topo_util_test(MBInterface *gMB) 
{
  MeshTopoUtil mtu(gMB);

    // construct a four-hex mesh for testing purposes
  double grid_vert_pos[] = 
    {
      0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0,
      0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 2.0, 1.0, 0.0,
      0.0, 2.0, 0.0, 1.0, 2.0, 0.0, 2.0, 2.0, 0.0,
//
      0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 2.0, 0.0, 1.0,
      0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0,
      0.0, 2.0, 1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 1.0
    };
      
  MBEntityHandle grid_verts[18], grid_elems[4];
  MBErrorCode result;
#define RR if (result != MB_SUCCESS) return result
  int init_edges, init_faces;
  result = gMB->get_number_entities_by_dimension(0, 1, init_edges); RR;
  result = gMB->get_number_entities_by_dimension(0, 2, init_faces); RR;

    // make vertices
  for (int i = 0; i < 18; i++) {
    result = gMB->create_vertex(&grid_vert_pos[3*i], grid_verts[i]); RR;
  }

    // make hexes
  int numv = 3, numv_sq = 9;
#define VINDEX(i,j,k) (i + (j*numv) + (k*numv_sq))
  MBEntityHandle connect[8];
  for (int j = 0; j < 2; j++) {
    for (int i = 0; i < 2; i++) {
      int vijk = VINDEX(i,j,0);
      connect[0] = grid_verts[vijk];
      connect[1] = grid_verts[vijk+1];
      connect[2] = grid_verts[vijk+1+numv];
      connect[3] = grid_verts[vijk+numv];
      connect[4] = grid_verts[vijk+numv*numv];
      connect[5] = grid_verts[vijk+1+numv*numv];
      connect[6] = grid_verts[vijk+1+numv+numv*numv];
      connect[7] = grid_verts[vijk+numv+numv*numv];
      result = gMB->create_element(MBHEX, connect, 8, grid_elems[2*j+i]); RR;
    }
  }

  MBRange vert_range;
  std::copy(grid_verts, grid_verts+18, mb_range_inserter(vert_range));
  
    // generate aentities
  result = mtu.construct_aentities(vert_range); RR;
  
  int this_edges, this_faces;
  result = gMB->get_number_entities_by_dimension(0, 1, this_edges); RR;
  result = gMB->get_number_entities_by_dimension(0, 2, this_faces); RR;

  if (this_edges != init_edges+33 || this_faces != init_faces+20) {
    std::cout << "Wrong number of edges or faces in mb_topo_util test." << std::endl;
//    return MB_FAILURE;
  }
  
    // get average position
  double pos[3];
  for (int j = 0; j < 2; j++) {
    for (int i = 0; i < 2; i++) {
      result = mtu.get_average_position(grid_elems[2*j+i], pos);
      RR;
      if (pos[0] != .5+i || pos[1] != .5+j || pos[2] != .5) {
        std::cout << "Wrong position at i = " << i << ", j = " << j << std::endl;
        result = MB_FAILURE;
      }
    }
  }
  RR;

    // get star faces
  MBRange all_hexes, middle_edge;
  std::copy(grid_elems, grid_elems+4, mb_range_inserter(all_hexes));
    // get the shared edge
  result = gMB->get_adjacencies(all_hexes, 1, false, middle_edge);
  if (MB_SUCCESS != result || 1 != middle_edge.size()) {
    std::cout << "Bad result getting single shared edge." << std::endl;
    return MB_FAILURE;
  }
  
  std::vector<MBEntityHandle> star_faces, star_hexes;
  bool bdy_edge;
  result = mtu.star_entities(*middle_edge.begin(), star_faces, bdy_edge, 0, &star_hexes);
  if (MB_SUCCESS != result || bdy_edge || star_faces.size() != 4 || star_hexes.size() != 4) {
    std::cout << "Bad result from star_faces for non-bdy edge." << std::endl;
    return MB_FAILURE;
  }
  
    // now try for a different edge, which has to be on the bdy
  MBRange other_edges;
  all_hexes.clear(); all_hexes.insert(grid_elems[0]);
  result = gMB->get_adjacencies(all_hexes, 1, false, other_edges); RR;
  other_edges.erase(*middle_edge.begin());
  if (11 != other_edges.size()) {
    std::cout << "Wrong number of edges in hex." << std::endl;
    return MB_FAILURE;
  }
  star_faces.clear();
  star_hexes.clear();
  result = mtu.star_entities(*other_edges.begin(), star_faces, bdy_edge, 0, &star_hexes);
  if (MB_SUCCESS != result || !bdy_edge || 
      (star_faces.size() != 2 && star_faces.size() != 3) ||
      (star_hexes.size() != 1 && star_hexes.size() != 2)) {
    std::cout << "Bad result from star_faces for bdy edge." << std::endl;
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

MBErrorCode mb_split_test(MBInterface *gMB) 
{
  MeshTopoUtil mtu(gMB);

    // construct a four-hex mesh for testing purposes
  double grid_vert_pos[] = 
    {
      0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0,
      0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 2.0, 1.0, 0.0,
      0.0, 2.0, 0.0, 1.0, 2.0, 0.0, 2.0, 2.0, 0.0,
//
      0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 2.0, 0.0, 1.0,
      0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0,
      0.0, 2.0, 1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 1.0,
//
      0.0, 0.0, 2.0, 1.0, 0.0, 2.0, 2.0, 0.0, 2.0,
      0.0, 1.0, 2.0, 1.0, 1.0, 2.0, 2.0, 1.0, 2.0,
      0.0, 2.0, 2.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0
    };
      
  MBEntityHandle grid_verts[27], grid_elems[8];
  MBErrorCode result;
#define RR if (result != MB_SUCCESS) return result
  int init_edges, init_faces, init_regions;
  result = gMB->get_number_entities_by_dimension(0, 1, init_edges); RR;
  result = gMB->get_number_entities_by_dimension(0, 2, init_faces); RR;
  result = gMB->get_number_entities_by_dimension(0, 3, init_regions); RR;

    // make vertices
  for (int i = 0; i < 27; i++) {
    result = gMB->create_vertex(&grid_vert_pos[3*i], grid_verts[i]); RR;
  }

    // make hexes
  int numv = 3, numv_sq = 9;
#define VINDEX(i,j,k) (i + (j*numv) + (k*numv_sq))
  MBEntityHandle connect[8];
  for (int k = 0; k < 2; k++) {
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 2; i++) {
        int vijk = VINDEX(i,j,k);
        connect[0] = grid_verts[vijk];
        connect[1] = grid_verts[vijk+1];
        connect[2] = grid_verts[vijk+1+numv];
        connect[3] = grid_verts[vijk+numv];
        connect[4] = grid_verts[vijk+numv*numv];
        connect[5] = grid_verts[vijk+1+numv*numv];
        connect[6] = grid_verts[vijk+1+numv+numv*numv];
        connect[7] = grid_verts[vijk+numv+numv*numv];
        result = gMB->create_element(MBHEX, connect, 8, grid_elems[4*k+2*j+i]); RR;
      }
    }
  }
  
  MBRange vert_range;
  std::copy(grid_verts, grid_verts+27, mb_range_inserter(vert_range));
  
    // generate aentities
  result = mtu.construct_aentities(vert_range); RR;
  
  int this_edges, this_faces;
  result = gMB->get_number_entities_by_dimension(0, 1, this_edges); RR;
  result = gMB->get_number_entities_by_dimension(0, 2, this_faces); RR;

  if (this_edges != init_edges+54 || this_faces != init_faces+36) {
    std::cout << "Wrong number of edges or faces in mb_topo_util test." << std::endl;
    return MB_FAILURE;
  }

    // split the faces between the 2 layers
    // first get the faces
  MBRange split_faces, tmp_ents, tmp_faces;
  for (int i = 0; i < 4; i++) {
    tmp_ents.clear();
    tmp_ents.insert(grid_elems[i]);
    tmp_ents.insert(grid_elems[i+4]);
    tmp_faces.clear();
    result = gMB->get_adjacencies(tmp_ents, 2, false, tmp_faces);
    if (MB_SUCCESS != result || tmp_faces.size() != 1) {
      std::cout << "mb_split_test failed to get shared quad." << std::endl;
      return MB_FAILURE;
    }
    split_faces.insert(*tmp_faces.begin());
  }

  MBRange new_faces, new_regions;

    // NOTE: passing in non-NULL pointer for new_regions requests that region between the
    // split entities be filled with an element; in this case, since we're splitting faces,
    // the new entities are polyhedra
  result = mtu.split_entities_manifold(split_faces, new_faces, &new_regions);
  if (MB_SUCCESS != result || new_faces.size() != 4 ||
      new_regions.size() != 4) {
    std::cout << "mb_split_test failed to split quads." << std::endl;
    return MB_FAILURE;
  }

  int this_regions;
  result = gMB->get_number_entities_by_dimension(0, 1, this_edges); RR;
  result = gMB->get_number_entities_by_dimension(0, 2, this_faces); RR;
  result = gMB->get_number_entities_by_dimension(0, 3, this_regions); RR;

  if (this_edges != init_edges+54 || this_faces != init_faces+40 ||
      this_regions != init_regions+12) {
    std::cout << "Wrong number of edges or faces or regions after splitting in mb_topo_util test." 
              << std::endl;
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

MBErrorCode mb_range_seq_intersect_test( MBInterface*) 
{
  MBErrorCode rval;
  SequenceManager sequences;
  MBRangeSeqIntersectIter iter( &sequences );
  MBRange range;

    // create some entity sequences
  EntitySequence *ts1, *ts2, *ts3, *qs1;
  MBEntityHandle th1, th2, th3, qh1;
  const int nt1 = 100, nt2 = 10, nt3 = 1, nq1 = 20;
  rval = sequences.create_entity_sequence( MBTRI, nt1, 3, 5, th1, ts1 );
  if (MB_SUCCESS != rval) return rval;
  rval = sequences.create_entity_sequence( MBTRI, nt2, 6, 0, th2, ts2 );
  if (MB_SUCCESS != rval) return rval;
  rval = sequences.create_entity_sequence( MBTRI, nt3, 3, MB_END_ID, th3, ts3 );
  if (MB_SUCCESS != rval) return rval;
  rval = sequences.create_entity_sequence( MBQUAD, nq1, 4, 0, qh1, qs1 );
  if (MB_SUCCESS != rval) return rval;
  
    // we're going to assume this below, so verify it now
  if (th1 > th2 || th2 > th3 || th3 > qh1)
    return MB_FAILURE;
  
    // construct an MBRange containing all valid handles;
  range.clear();
  range.insert( ts1->start_handle(), ts1->end_handle() );
  range.insert( ts2->start_handle(), ts2->end_handle() );
  range.insert( ts3->start_handle(), ts3->end_handle() );
  range.insert( qs1->start_handle(), qs1->end_handle() );
  
    // iterate over all and check results

  rval = iter.init( range.begin(), range.end() );
  if (MB_SUCCESS != rval) 
    return rval;
  if (ts1 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != ts1->start_handle())
    return MB_FAILURE;
  if (iter.get_end_handle() != ts1->end_handle())
    return MB_FAILURE;

  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (ts2 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != ts2->start_handle())
    return MB_FAILURE;
  if (iter.get_end_handle() != ts2->end_handle())
    return MB_FAILURE;

  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (ts3 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != ts3->start_handle())
    return MB_FAILURE;
  if (iter.get_end_handle() != ts3->end_handle())
    return MB_FAILURE;

  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (qs1 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != qs1->start_handle())
    return MB_FAILURE;
  if (iter.get_end_handle() != qs1->end_handle())
    return MB_FAILURE;
  
  if (!iter.is_at_end())
    return MB_FAILURE;
  rval = iter.step();
  if (MB_FAILURE != rval)
    return MB_FAILURE;
  
    // iterate over just the quads

  rval = iter.init( range.lower_bound(MBQUAD), range.end() );
  if (MB_SUCCESS != rval) 
    return rval;
  if (qs1 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != qs1->start_handle())
    return MB_FAILURE;
  if (iter.get_end_handle() != qs1->end_handle())
    return MB_FAILURE;
  
  if (!iter.is_at_end())
    return MB_FAILURE;
  rval = iter.step();
  if (MB_FAILURE != rval)
    return MB_FAILURE;

    // iterate starting one past the beginning of the
    // triangles and stopping one before the end.  The last
    // sequence contains only one tri, so should stop and end
    // of second-to-last sequence

  rval = iter.init( ++(range.begin()), --(range.lower_bound(MBQUAD)) );
  if (MB_SUCCESS != rval) 
    return rval;
  if (ts1 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != ts1->start_handle() + 1)
    return MB_FAILURE;
  if (iter.get_end_handle() != ts1->end_handle())
    return MB_FAILURE;

  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (ts2 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != ts2->start_handle())
    return MB_FAILURE;
  if (iter.get_end_handle() != ts2->end_handle())
    return MB_FAILURE;
  
  if (!iter.is_at_end())
    return MB_FAILURE;
  rval = iter.step();
  if (MB_FAILURE != rval)
    return MB_FAILURE;
   
    // Iterate over the quad sequence, starting two past the
    // beginning and stopping one before the end

  rval = iter.init( ++(range.lower_bound(MBQUAD)), --(range.end()));
  if (MB_SUCCESS != rval) 
    return rval;
  if (qs1 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != qs1->start_handle() + 1)
    return MB_FAILURE;
  if (iter.get_end_handle() != qs1->end_handle() - 1)
    return MB_FAILURE;
  
  if (!iter.is_at_end())
    return MB_FAILURE;
  rval = iter.step();
  if (MB_FAILURE != rval)
    return MB_FAILURE;
   
    // Iterate over two subsets of the quad sequence
    
  MBRange quads = range.subset_by_type( MBQUAD );
  MBEntityHandle removed = qs1->start_handle() + nq1/2;
  if (quads.erase( removed ) == quads.end())
    return MB_FAILURE;

  rval = iter.init( quads.begin(), quads.end());
  if (MB_SUCCESS != rval) 
    return rval;
  if (qs1 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != qs1->start_handle())
    return MB_FAILURE;
  if (iter.get_end_handle() != removed - 1)
    return MB_FAILURE;

  rval = iter.step();
  if (MB_SUCCESS != rval) 
    return rval;
  if (qs1 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != removed + 1)
    return MB_FAILURE;
  if (iter.get_end_handle() != qs1->end_handle())
    return MB_FAILURE;
  
  if (!iter.is_at_end())
    return MB_FAILURE;
  rval = iter.step();
  if (MB_FAILURE != rval)
    return MB_FAILURE;
    
    // Iterate over everything, including a bunch of
    // invalid handles

  MBRange big;
  int junk;
  MBEntityHandle last = CREATE_HANDLE(MBQUAD+1, 0, junk);
  big.insert( ts1->start_handle() - 1, last );

    // first some invalid handles in the beginning of the range
  rval = iter.init( big.begin(), big.end() );
  if (MB_ENTITY_NOT_FOUND != rval) 
    return MB_FAILURE;
  if (NULL != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != *big.begin())
    return MB_FAILURE;
  if (iter.get_end_handle() != ts1->start_handle() - 1)
    return MB_FAILURE;

    // next the first triangle sequence
  rval = iter.step();
  if (MB_SUCCESS != rval) 
    return rval;
  if (ts1 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != ts1->start_handle())
    return MB_FAILURE;
  if (iter.get_end_handle() != ts1->end_handle())
    return MB_FAILURE;

    // next the the invalid handles between the first two tri sequences
  if (ts1->end_handle() + 1 != ts2->start_handle()) {
    rval = iter.step();
    if (MB_ENTITY_NOT_FOUND != rval) 
      return MB_FAILURE;
    if (NULL != iter.get_sequence())
      return MB_FAILURE;
    if (iter.get_start_handle() != ts1->end_handle()+1)
      return MB_FAILURE;
    if (iter.get_end_handle() != ts2->start_handle()-1)
      return MB_FAILURE;
  }

    // next the second triangle sequence
  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (ts2 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != ts2->start_handle())
    return MB_FAILURE;
  if (iter.get_end_handle() != ts2->end_handle())
    return MB_FAILURE;

    // next the the invalid handles between the 2nd and 3rd tri sequences
  if (ts2->end_handle() + 1 != ts3->start_handle()) {
    rval = iter.step();
    if (MB_ENTITY_NOT_FOUND != rval) 
      return MB_FAILURE;
    if (NULL != iter.get_sequence())
      return MB_FAILURE;
    if (iter.get_start_handle() != ts2->end_handle()+1)
      return MB_FAILURE;
    if (iter.get_end_handle() != ts3->start_handle()-1)
      return MB_FAILURE;
  }

    // next the third triangle sequence
  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (ts3 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != ts3->start_handle())
    return MB_FAILURE;
  if (iter.get_end_handle() != ts3->end_handle())
    return MB_FAILURE;

    // third tri sequence contains the MAX tri handle, so no more
    // invalid triangles.
    // next 1 invalid quad at the before MB_START_ID
  if (ts3->end_handle() + 1 != qs1->start_handle()) {
    rval = iter.step();
    if (MB_ENTITY_NOT_FOUND != rval) 
      return MB_FAILURE;
    if (NULL != iter.get_sequence())
      return MB_FAILURE;
    if (iter.get_start_handle() != ts3->end_handle()+1)
      return MB_FAILURE;
    if (iter.get_end_handle() != qs1->start_handle()-1)
      return MB_FAILURE;
  }

    // next the quad sequence
  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (qs1 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != qs1->start_handle())
    return MB_FAILURE;
  if (iter.get_end_handle() != qs1->end_handle())
    return MB_FAILURE;

    // next remaining invalid quad handles in the range
  rval = iter.step();
  if (MB_ENTITY_NOT_FOUND != rval)
    return MB_FAILURE;
  if (0 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != qs1->end_handle() + 1)
    return MB_FAILURE;
  if (iter.get_end_handle() != last - 1)
    return MB_FAILURE;

    // next invalid entity after the last quad in the range
  rval = iter.step();
  if (MB_ENTITY_NOT_FOUND != rval)
    return MB_FAILURE;
  if (0 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != last)
    return MB_FAILURE;
  if (iter.get_end_handle() != last)
    return MB_FAILURE;
  
    // now at the end
  if (!iter.is_at_end())
    return MB_FAILURE;
  rval = iter.step();
  if (MB_FAILURE != rval)
    return MB_FAILURE;
  
  
  
    // Create some holes
  MBEntityHandle ts1s  = ts1->start_handle();
  MBEntityHandle dead1 = ts1->start_handle() + 1;
  MBEntityHandle dead2 = ts1->start_handle() + 2;
  MBEntityHandle ts1e  = ts1->end_handle();
  MBEntityHandle dead3 = ts2->start_handle();
  MBEntityHandle dead4 = ts2->end_handle();
  MBEntityHandle qs1s  = qs1->start_handle();
  MBEntityHandle qs1e  = qs1->end_handle();
  MBEntityHandle dead5 = qs1->start_handle() + nq1/2;
  MBEntityHandle dead6 = dead5+1;
  rval = sequences.delete_entity( dead1 );
  if (MB_SUCCESS != rval) return rval;
  rval = sequences.delete_entity( dead2 );
  if (MB_SUCCESS != rval) return rval;
  rval = sequences.delete_entity( dead3 );
  if (MB_SUCCESS != rval) return rval;
  rval = sequences.delete_entity( dead4 );
  if (MB_SUCCESS != rval) return rval;
  rval = sequences.delete_entity( dead5 );
  if (MB_SUCCESS != rval) return rval;
  rval = sequences.delete_entity( dead6 );
  if (MB_SUCCESS != rval) return rval;
  
    // Iterate over sequences w/out removing deleted entities
    // from range.
  
    // first sequence should have one valid handle at beginning
  rval = iter.init( range.begin(), range.end() );
  if (MB_SUCCESS != rval) 
    return rval;
  if (0 == iter.get_sequence() ||
      ts1s != iter.get_sequence()->start_handle() ||
      dead1-1 != iter.get_sequence()->end_handle())
    return MB_FAILURE;
  if (iter.get_start_handle() != ts1->start_handle())
    return MB_FAILURE;
  if (iter.get_end_handle() != dead1 - 1)
    return MB_FAILURE;
  
    // next two invalid handles in sequence in first sequence
  rval = iter.step( );
  if (MB_ENTITY_NOT_FOUND != rval) 
    return MB_FAILURE;
  if (0 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != dead1)
    return MB_FAILURE;
  if (iter.get_end_handle() != dead2)
    return MB_FAILURE;
  
    // next the remainder of the fist sequence
  rval = iter.step();
  if (MB_SUCCESS != rval) 
    return rval;
  if (0 == iter.get_sequence() ||
      dead2+1 != iter.get_sequence()->start_handle() ||
      ts1e != iter.get_sequence()->end_handle())
    return MB_FAILURE;
  if (iter.get_start_handle() != dead2+1)
    return MB_FAILURE;
  if (iter.get_end_handle() != ts1e)
    return MB_FAILURE;

    // next an invalid handle at the start of the second sequence
  rval = iter.step();
  if (MB_ENTITY_NOT_FOUND != rval)
    return MB_FAILURE;
  if (0 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != dead3)
    return MB_FAILURE;
  if (iter.get_end_handle() != dead3)
    return MB_FAILURE;

    // next the second sequence up to the invalid handle at the end
  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (0 == iter.get_sequence() ||
      dead3+1 != iter.get_sequence()->start_handle() ||
      dead4-1 != iter.get_sequence()->end_handle())
    return MB_FAILURE;
  if (ts2 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != dead3+1)
    return MB_FAILURE;
  if (iter.get_end_handle() != dead4-1)
    return MB_FAILURE;

    // next invaild handle at the end of the second sequence
  rval = iter.step();
  if (MB_ENTITY_NOT_FOUND != rval)
    return MB_FAILURE;
  if (0 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != dead4)
    return MB_FAILURE;
  if (iter.get_end_handle() != dead4)
    return MB_FAILURE;

    // next the third sequence
  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (ts3 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != ts3->start_handle())
    return MB_FAILURE;
  if (iter.get_end_handle() != ts3->end_handle())
    return MB_FAILURE;

    // next the quad sequence up to the invalid handle in the middle
  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (0 == iter.get_sequence() ||
      qs1s != iter.get_sequence()->start_handle() ||
      dead5-1 != iter.get_sequence()->end_handle())
    return MB_FAILURE;
  if (qs1 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != qs1s)
    return MB_FAILURE;
  if (iter.get_end_handle() != dead5-1)
    return MB_FAILURE;

    // next the two invalid handles in the middle
  rval = iter.step();
  if (MB_ENTITY_NOT_FOUND != rval)
    return MB_FAILURE;
  if (0 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != dead5)
    return MB_FAILURE;
  if (iter.get_end_handle() != dead6)
    return MB_FAILURE;

    // next the remainder of the quad sequence
  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (0 == iter.get_sequence() ||
      dead6+1 != iter.get_sequence()->start_handle() ||
      qs1e != iter.get_sequence()->end_handle())
    return MB_FAILURE;
  if (iter.get_start_handle() != dead6+1)
    return MB_FAILURE;
  if (iter.get_end_handle() != qs1e)
    return MB_FAILURE;
  
    // now at the end
  if (!iter.is_at_end())
    return MB_FAILURE;
  rval = iter.step();
  if (MB_FAILURE != rval)
    return MB_FAILURE;

  
    // now remove the dead entities from the range and iterate again
    
  if (range.erase( dead1 ) == quads.end())
    return MB_FAILURE;
  if (range.erase( dead2 ) == quads.end())
    return MB_FAILURE;
  if (range.erase( dead3 ) == quads.end())
    return MB_FAILURE;
  if (range.erase( dead4 ) == quads.end())
    return MB_FAILURE;
  if (range.erase( dead5 ) == quads.end())
    return MB_FAILURE;
  if (range.erase( dead6 ) == quads.end())
    return MB_FAILURE;
  
  
    // first sequence should have one valid handle at beginning
  rval = iter.init( range.begin(), range.end() );
  if (MB_SUCCESS != rval) 
    return rval;
  if (0 == iter.get_sequence() ||
      ts1s != iter.get_sequence()->start_handle() ||
      dead1-1 != iter.get_sequence()->end_handle())
    return MB_FAILURE;
  if (iter.get_start_handle() != ts1s)
    return MB_FAILURE;
  if (iter.get_end_handle() != dead1 - 1)
    return MB_FAILURE;
  
    // next the remainder of the fist sequence after the hole
  rval = iter.step();
  if (MB_SUCCESS != rval) 
    return rval;
  if (0 == iter.get_sequence() ||
      dead2+1 != iter.get_sequence()->start_handle() ||
      ts1e != iter.get_sequence()->end_handle())
    return MB_FAILURE;
  if (iter.get_start_handle() != dead2+1)
    return MB_FAILURE;
  if (iter.get_end_handle() != ts1e)
    return MB_FAILURE;

    // next the second sequence between deleted start and end handles
  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (0 == iter.get_sequence() ||
      dead3+1 != iter.get_sequence()->start_handle() ||
      dead4-1 != iter.get_sequence()->end_handle())
    return MB_FAILURE;
  if (iter.get_start_handle() != dead3+1)
    return MB_FAILURE;
  if (iter.get_end_handle() != dead4-1)
    return MB_FAILURE;

    // next the third sequence
  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (ts3 != iter.get_sequence())
    return MB_FAILURE;
  if (iter.get_start_handle() != ts3->start_handle())
    return MB_FAILURE;
  if (iter.get_end_handle() != ts3->end_handle())
    return MB_FAILURE;

    // next the quad sequence up to the hole in the middle
  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (0 == iter.get_sequence() ||
      qs1s != iter.get_sequence()->start_handle() ||
      dead5-1 != iter.get_sequence()->end_handle())
    return MB_FAILURE;
  if (iter.get_start_handle() != qs1s)
    return MB_FAILURE;
  if (iter.get_end_handle() != dead5-1)
    return MB_FAILURE;

    // next the remainder of the quad sequence after the hole
  rval = iter.step();
  if (MB_SUCCESS != rval)
    return rval;
  if (0 == iter.get_sequence() ||
      dead6+1 != iter.get_sequence()->start_handle() ||
      qs1e != iter.get_sequence()->end_handle())
    return MB_FAILURE;
  if (iter.get_start_handle() != dead6+1)
    return MB_FAILURE;
  if (iter.get_end_handle() != qs1e)
    return MB_FAILURE;
  
    // now at the end
  if (!iter.is_at_end())
    return MB_FAILURE;
  rval = iter.step();
  if (MB_FAILURE != rval)
    return MB_FAILURE;
  
  return MB_SUCCESS;
}

#define ASSERT_EQUAL( A, B ) \
  do { if (!_assert_equal( (A), (B), #A, #B, __LINE__)) \
    return MB_FAILURE; } while(false)

#define ASSERT_NOT_EQUAL( A, B ) \
  do { if (!_assert_not_equal( (A), (B), #A, #B, __LINE__)) \
    return MB_FAILURE; } while(false)

template <typename T1, typename T2>
bool _assert_equal( T1 a, T2 b, const char* as, const char* bs, int line )
{
  if (a == b)
    return true;
  
  std::cout << "Assertion failed at line " << line << std::endl
            << "\t" << as << " == " << bs << std::endl
            << "\t" << as << " = " << a << std::endl
            << "\t" << bs << " = " << b << std::endl;
  return false;
}

template <typename T1, typename T2>
bool _assert_not_equal( T1 a, T2 b, const char* as, const char* bs, int line )
{
  if (a != b)
    return true;
  
  std::cout << "Assertion failed at line " << line << std::endl
            << "\t" << as << " != " << bs << std::endl
            << "\t" << as << " = " << a << std::endl
            << "\t" << bs << " = " << b << std::endl;
  return false;
}
    
std::ostream& operator<<( std::ostream& s, MBRange::const_iterator i ) {
  return s << *i;
}

MBErrorCode mb_poly_adjacency_test( MBInterface* )
{
  MBErrorCode rval;
  MBCore moab;
  MBInterface *mbImpl = &moab;
  
    // make a couple polygons and a polyhedron
  double coords[3] = {0,1,2};
  MBEntityHandle verts[10], polygons[2], polyhedron;
  
  for (int i = 0; i < 10; i++) {
    rval = mbImpl->create_vertex(coords, verts[i]);
    if (MB_SUCCESS != rval)
      return rval;
  }

  for (int i = 0; i < 2; i++) {
    rval = mbImpl->create_element(MBPOLYGON, verts, 5, polygons[i]);
    if (MB_SUCCESS != rval)
      return rval;
  }
  rval = mbImpl->create_element(MBPOLYHEDRON, polygons, 2, polyhedron);
  if (MB_SUCCESS != rval)
    return rval;
  
    // create the aentities
  MBRange dum_range;
  for (int dim = 0; dim < 3; dim++) {
    dum_range.clear();
    rval = mbImpl->get_adjacencies(&polyhedron, 1, dim, true, dum_range);
    if (MB_SUCCESS != rval)
      return rval;
  }
  
    // delete the polyhedron
  rval = mbImpl->delete_entities(&polyhedron, 1);
  if (MB_SUCCESS != rval)
    return rval;
  
    // check adjacencies
  return moab.check_adjacencies();
}

MBErrorCode mb_memory_use_test( MBInterface* ) 
{
  MBCore mb;
  unsigned long init_total, total_with_elem, total_with_tag, total_with_tag_data;
  mb.estimated_memory_use(0,0,0,&init_total);
  
  double coords[12] = { 1, 2, 0, 3, 4, 0, 5, 6, 0, 7, 8, 0 };
  MBEntityHandle verts[4];
  for (int i = 0; i < 4; ++i)
    if (MB_SUCCESS != mb.create_vertex( coords + 3*i, verts[i] ))
      return MB_FAILURE;
  
  MBEntityHandle elem;
  mb.create_element( MBQUAD, verts, 4, elem );
  
  mb.estimated_memory_use(0,0,0,&total_with_elem);
  if (total_with_elem <= init_total)
    return MB_FAILURE;
  
  unsigned long min, am;
  MBRange r;
  r.insert( elem );
  mb.estimated_memory_use( r, &min, &am );
  if (min != 4*sizeof(MBEntityHandle))
    return MB_FAILURE;
  
  r.clear();
  r.insert( verts[0] );
  r.insert( verts[1] );
  mb.estimated_memory_use( r, &min, &am );
  if (min != 6*sizeof(double))
    return MB_FAILURE;
  
  MBTag tag;
  if (MB_SUCCESS != mb.tag_create( "TMP_TAG", sizeof(int), MB_TAG_SPARSE, tag, 0 ))
    return MB_FAILURE;
  mb.estimated_memory_use( r, &min, &am );
  if (min != 6*sizeof(double))
    return MB_FAILURE;
    
  mb.estimated_memory_use(0,0,0,&total_with_tag);
  if (total_with_tag <= total_with_elem)
    return MB_FAILURE;

  int tag_data[] = { 0xA, 0xB };
  if (MB_SUCCESS != mb.tag_set_data( tag, r, &tag_data ))
    return MB_FAILURE;
  mb.estimated_memory_use( r, &min, &am );
  if (min <= 6*sizeof(double))
    return MB_FAILURE;
    
  mb.estimated_memory_use(0,0,0,&total_with_tag_data);
  if (total_with_tag_data <= total_with_tag)
    return MB_FAILURE;
  
  return MB_SUCCESS;
}
  


int number_tests = 0;
int number_tests_failed = 0;
#define RUN_TEST( A ) _run_test( (A), #A )

typedef MBErrorCode (*TestFunc)(MBInterface*);
static void _run_test( TestFunc func, const char* func_str ) 
{
  MBErrorCode error;
  MBCore moab;
  MBInterface* iface = &moab;
  
  std::string file_name = TestDir + "/mbtest1.g";
  error = iface->load_mesh( file_name.c_str() );
  if (MB_SUCCESS != error) {
    std::cout << "Failed to load input file: " << file_name << std::endl;
    std::string error_reason;
    iface->get_last_error(error_reason);
    cout << error_reason << std::endl;
    abort(); // going to try this for every test, so if it doesn't work, just abort
  }
    
  ++number_tests;
  cout << "   " << func_str << ": ";
  cout.flush();
  error = func( iface );
  
  if (MB_SUCCESS == error)
    std::cout << "Success" << std::endl;
  else if (MB_FAILURE == error)
    std::cout << "Failure" << std::endl;
  else {
    std::cout << "Failed: " << moab.get_error_string( error ) << std::endl;
  }
  
  if (MB_SUCCESS != error) {
    ++number_tests_failed;
    
    std::string error_reason;
    iface->get_last_error(error_reason);
    cout << error_reason << std::endl;
  }
}


static void usage(const char* exe) {
  cerr << "Usage: " << exe << " [-nostress] [-d input_file_dir]\n";
  exit (1);
}


/*!
  main routine for test harness 
*/


int main(int argc, char* argv[])
{

    // Check command line arg to see if we should avoid doing the stress test
  bool stress_test = true;

  std::cout << "Size of mConnMap = " << sizeof(MBCN::mConnectivityMap)
            << std::endl;
  std::cout << "Size of mUpConnMap = " << sizeof(MBCN::mUpConnMap)
            << std::endl;
  std::cout << "Size of MBCN = " << sizeof(MBCN)
            << std::endl;
    
  for (int i = 1; i < argc; ++i) {
    if (string(argv[i]) == "-nostress")
      stress_test = false;
    else if (string(argv[i]) == "-d" && (i+1) < argc)
      TestDir = argv[++i];
    else if (string(argv[i]) == "-h" || string(argv[i]) == "--help")
      usage( argv[0] );
    else {
      cerr << "Invalid argument: " << argv[i] << endl;
      usage( argv[0] );
    }
  }

    // Print out Header information

  cout << "\n\nMB TEST PROGRAM:\n\n";

  RUN_TEST( mb_adjacent_vertex_test );
  RUN_TEST( mb_adjacencies_test );
  RUN_TEST( mb_upward_adjacencies_test );
  RUN_TEST( mb_vertex_coordinate_test );
  RUN_TEST( mb_vertex_tag_test );
  RUN_TEST( mb_bar_connectivity_test );
  RUN_TEST( mb_tri_connectivity_test );
  RUN_TEST( mb_quad_connectivity_test );
  RUN_TEST( mb_tet_connectivity_test );
  RUN_TEST( mb_temporary_test );
  RUN_TEST( mb_hex_connectivity_test );
  RUN_TEST( mb_mesh_sets_set_test );
  RUN_TEST( mb_mesh_sets_list_test );
  RUN_TEST( mb_mesh_set_parent_child_test );
  RUN_TEST( mb_mesh_set_set_appends );
  RUN_TEST( mb_mesh_set_list_appends );
  RUN_TEST( mb_mesh_set_root_appends );
  RUN_TEST( mb_tags_test );
  RUN_TEST( mb_dense_tag_test );
  RUN_TEST( mb_sparse_tag_test );
  RUN_TEST( mb_write_mesh_test );
  RUN_TEST( mb_delete_mesh_test );
  RUN_TEST( mb_meshset_tracking_test );
  RUN_TEST( mb_higher_order_test );
  RUN_TEST( mb_bit_tags_test );
  RUN_TEST( mb_entity_conversion_test );
  RUN_TEST( mb_forced_adjacencies_test );
  RUN_TEST( mb_canon_number_test );
  RUN_TEST( mb_poly_test );
  RUN_TEST( mb_range_test );
  RUN_TEST( mb_topo_util_test );
  RUN_TEST( mb_split_test );
  RUN_TEST( mb_range_seq_intersect_test );
  RUN_TEST( mb_poly_adjacency_test );
  RUN_TEST( mb_memory_use_test );
  RUN_TEST( mb_merge_test );
  RUN_TEST( mb_merge_update_test );
  if (stress_test) RUN_TEST( mb_stress_test );

    // summary

  cout << "\nMB TEST SUMMARY: \n"
       << "   Number Tests:           " << number_tests << "\n"
       << "   Number Successful:      " << number_tests - number_tests_failed << "\n"
       << "   Number Failed:          " << number_tests_failed 
       << "\n\n";

  return number_tests_failed;
}
