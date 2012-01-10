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
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "moab/Interface.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/CN.hpp"
#include "OrientedBox.hpp"
#include "moab/CartVect.hpp"
#include "moab/WriteUtilIface.hpp"

#ifdef USE_MPI
#include "moab_mpi.h"
#endif

#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "Internals.hpp"
#include "moab/Core.hpp"
#include "SequenceManager.hpp"
#include "EntitySequence.hpp"
#include "RangeSeqIntersectIter.hpp"
#include "moab/Error.hpp"
#include "moab/ScdInterface.hpp"

using namespace std;
using namespace moab;

#define STRINGIFY_(A) #A
#define STRINGIFY(A) STRINGIFY_(A)
#ifdef MESHDIR
string TestDir( STRINGIFY(MESHDIR) );
#else
string TestDir( "." );
#endif

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)


#define CHECK_EQUAL( A, B ) do { if ((A) != (B)) { \
            std::cerr << "Equality Test failed at " __FILE__ ":" << __LINE__ << std::endl; \
            std::cerr << "  Expected: " << (A) << std::endl; \
            std::cerr << "  Actual:   " << (B) << std::endl; \
            return MB_FAILURE; } } while(false)
#define CHECK( A ) do { if (!(A)) { \
            std::cerr << "Test failed at " __FILE__ ":" << __LINE__ << std::endl; \
            return MB_FAILURE; } } while(false)


ErrorCode load_file_one( Interface* iface )
{
  std::string file_name = TestDir + "/mbtest1.g";
  ErrorCode error = iface->load_mesh( file_name.c_str() );
  if (MB_SUCCESS != error) {
    std::cout << "Failed to load input file: " << file_name << std::endl;
    std::string error_reason;
    iface->get_last_error(error_reason);
    cout << error_reason << std::endl;
  }
  return error;
}

/* Create a regular 2x2x2 hex mesh */
ErrorCode create_some_mesh( Interface* iface );

ErrorCode check_valid_connectivity( Interface* iface );


  /*!
    @test 
    Vertex Coordinates
    @li Get coordinates of vertex 1 correctly
    @li Get coordinates of vertex 8 correctly
    @li Get coordinates of vertex 6 correctly
  */
ErrorCode mb_vertex_coordinate_test()
{
  Core moab;
  Interface* MB = &moab;
  ErrorCode error = create_some_mesh( MB );
  CHKERR(error);

  Range vertices;
  error = MB->get_entities_by_type(0,  MBVERTEX, vertices);
  CHKERR(error);

  std::vector<double> all_coords(3*vertices.size());
  double* coord_iter = &all_coords[0];
  for ( Range::iterator iter = vertices.begin();
        iter != vertices.end(); ++iter)
  {
    error = MB->get_coords(&(*iter), 1, coord_iter );
    CHKERR(error);
    coord_iter += 3;
  }
    
    // check blocked coordinates
  const size_t N = vertices.size()+1;
  std::vector<double> x(N), y(N), z(N);
    // set last value so we can check later that nothing wrote past the
    // intended end of an array
  x[vertices.size()] = y[vertices.size()] = z[vertices.size()] = -3.14159;
  error = MB->get_coords( vertices, &x[0], &y[0], &z[0] );
  for (size_t i = 0; i < vertices.size(); ++i) {
    CHECK_EQUAL( all_coords[3*i  ], x[i] );
    CHECK_EQUAL( all_coords[3*i+1], y[i] );
    CHECK_EQUAL( all_coords[3*i+2], z[i] );
  }
    // checkthat get_coords did not write past intended end of arrays
  CHECK_EQUAL( -3.14159, x[vertices.size()] );
  CHECK_EQUAL( -3.14159, y[vertices.size()] );
  CHECK_EQUAL( -3.14159, z[vertices.size()] );
    
    // add invalid handle to end of range and try query again
  vertices.insert( vertices.back() + 1 );
  error = MB->get_coords( vertices, &x[0], &y[0], &z[0] );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, error );

    // Try getting coordinates for a hex (should fail)
  Range hexes;
  error = MB->get_entities_by_type( 0, MBHEX, hexes );
  EntityHandle handle = hexes.front();
  error = MB->get_coords(&handle, 1, &x[0]);
  CHECK_EQUAL( MB_TYPE_OUT_OF_RANGE, error );
  
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
ErrorCode mb_vertex_tag_test()
{
  Core moab;
  Interface* MB = &moab;
  ErrorCode error = create_some_mesh( MB );
  if (MB_SUCCESS != error)
    return error;

    // Add an int Vertex Tag to the database

  Tag tag_id;

    // Create a dense tag for all vertices
  error = MB->tag_get_handle("int_tag", 1, MB_TYPE_INTEGER, tag_id, MB_TAG_SPARSE|MB_TAG_EXCL);
  if (error != MB_SUCCESS)
    return error;

    // put a value in vertex 1 and retrieve
  std::vector<EntityHandle> verts;
  error = MB->get_entities_by_type( 0, MBVERTEX, verts );
  EntityHandle handle = verts[0];
  int input_value = 11;
  error = MB->tag_set_data(tag_id, &handle, 1, &input_value);
  if (MB_SUCCESS != error) return error;

  int output_value;
  error = MB->tag_get_data(tag_id, &handle, 1, &output_value);
  if(MB_SUCCESS != error) return error;
  else if (output_value != input_value)
    return MB_FAILURE;

    // put a value in vertex 5 and retrieve
  
  handle = verts[5];
  input_value = 11;
  error = MB->tag_set_data(tag_id, &handle, 1, &input_value);
  if (MB_SUCCESS != error) return error;
  error = MB->tag_get_data(tag_id, &handle, 1, &output_value);
  if (MB_SUCCESS != error) return error;
  else if(output_value != input_value)
    return MB_FAILURE;

    // put a value in vertex 98088234 which doesn't exist
  handle = *std::max_element( verts.begin(), verts.end() ) + 1;
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

  error = MB->tag_get_handle("bool_tag", sizeof(bool), MB_TYPE_OPAQUE, tag_id, MB_TAG_SPARSE|MB_TAG_EXCL);
  if (error != MB_SUCCESS)
    return error;

    // put a value in vertex 5 and retrieve

  handle = verts[5];
  bool bool_input_value = true;
  bool bool_output_value = false;
  error = MB->tag_set_data(tag_id, &handle, 1, &bool_input_value);
  if (error != MB_SUCCESS) return error;
  error = MB->tag_get_data(tag_id, &handle, 1, &bool_output_value);
  if (error != MB_SUCCESS) return error;
  else if(bool_output_value != bool_input_value)
    return MB_FAILURE;

    // Add a double Vertex Tag to the database

  error = MB->tag_get_handle("double_tag", 1, MB_TYPE_DOUBLE, tag_id, MB_TAG_SPARSE|MB_TAG_EXCL);
  if (error != MB_SUCCESS)
    return error;

    // put a value in vertex 8: and retrieve

  handle = verts[8];
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
  error = MB->tag_get_handle("struct_tag", sizeof(TagStruct), MB_TYPE_OPAQUE, tag_id, MB_TAG_SPARSE|MB_TAG_EXCL);
  if (error != MB_SUCCESS)
    return error;

    // put a value in vertex 7 and retrieve

  handle = verts[7];
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

  error = MB->tag_get_handle("sparse_int_tag", 1, MB_TYPE_INTEGER, tag_id, MB_TAG_SPARSE|MB_TAG_EXCL);

  if (error != MB_SUCCESS )  
    return error;

    //print_yes = true;
  int i;
  for (i=0; i<10; i++)
  {
      // use invalid handles for odd values
    EntityHandle handle;
    if (i % 2) 
      handle = verts[i] + *std::max_element( verts.begin(), verts.end() );
    else 
      handle = verts[i];
      
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
  Tag int_tag_handle;
  error = MB->tag_get_handle (int_tag_name.c_str(), 1, MB_TYPE_INTEGER, int_tag_handle);
  if (MB_SUCCESS != error) return error;
    
  if (int_tag_handle != tag_id)
    return MB_FAILURE;

    // test tag_get_tags_on_entity and tag_delete_data
  std::vector<Tag> all_tags;
  handle = verts[0];
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
  handle = verts[2];
  error = MB->tag_delete_data(tag_id, &handle, 1);
  if (error != MB_SUCCESS )  
    return error;

  handle = verts[6];
  error = MB->tag_delete_data(tag_id, &handle, 1);
  if (error != MB_SUCCESS )  
    return error;

    // delete all the rest of the sparse tags.

  error = MB->tag_delete(tag_id);
  if (error != MB_SUCCESS )  
    return error;

    // delete the dense tag named bool_tag 
  Tag bool_tag_handle;
  error = MB->tag_get_handle ("bool_tag", sizeof(bool), MB_TYPE_OPAQUE, bool_tag_handle);
  if (error != MB_SUCCESS) return error;

  error = MB->tag_delete(bool_tag_handle);
  if (error != MB_SUCCESS )  
    return error;

  return error;
}


ErrorCode mb_temporary_test()
{
  Core moab;
  Interface* gMB = &moab;

  double array[3] = {0.0, 0.0, 0.0};
  EntityHandle h_node1;
  ErrorCode result = gMB->create_vertex(array, h_node1);
  if (MB_SUCCESS != result)
    return result;

  EntityHandle ordered_meshset1;
  result = gMB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, ordered_meshset1);
  if (MB_SUCCESS != result)
    return result;

  EntityHandle ordered_meshset2;
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
  std::vector<EntityHandle> meshsets;
  result = gMB->get_adjacencies(&h_node1, 1, 4, create_if_missing, meshsets);
  if (MB_SUCCESS != result)
    return result;

  if (1u != meshsets.size())
    return MB_FAILURE;

  return MB_SUCCESS;
}

ErrorCode mb_adjacent_vertex_test()
{
  Core moab;
  Interface* mb = &moab;
  ErrorCode rval = create_some_mesh( mb );
  if (MB_SUCCESS != rval)
    return rval;

  Range hexes, expt_vert, got_vert, some_hexes;
  Range::const_iterator i, j;
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
    const EntityHandle* conn;
    int len;
    rval = mb->get_connectivity( *i, conn, len );
    if (MB_SUCCESS != rval)
      return rval;
    for (int k = 0; k < len; ++k)
      expt_vert.insert( conn[k] );
  }
  
    // use get_adjacencies to get vertices
  rval = mb->get_adjacencies( some_hexes, 0, false, got_vert, Interface::UNION );
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
  
ErrorCode mb_adjacencies_test() 
{
  Core moab;
  Interface* mb = &moab;
  ErrorCode result = load_file_one( mb );
  if (MB_SUCCESS != result)
    return result;

    // this test does the following:
    // 1. For each element, creates vertex-element adjacencies (only for
    //    lowest-id vertex)
    // 2. Creates all lower-order ancillary entities
    // 3. Checks for proper number of same
    //
    // assume mesh has already been read

  EntityType seq_type;
  Range handle_range;

    // lets create a skin of the hexes
    // this may be far from being the most efficient, but
    // it certainly does exercise the adjacency routines

  Range::iterator iter;
  Range::reverse_iterator riter;

    // first get the hexes
  Range hexes;
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
  Range nodes;
  result = mb->get_adjacencies(hexes, 0, false, nodes, Interface::UNION);
  if (MB_SUCCESS != result)
    return result;

    // make sure we got nodes
  for(iter = nodes.begin(); iter != nodes.end(); ++iter)
  {
    if( TYPE_FROM_HANDLE(*iter) != MBVERTEX)
      return MB_FAILURE;
  }

    // find the interior nodes; assume a structured mesh
  Range interior_nodes;
  std::vector<EntityHandle> attached_hexes;
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
  Range interior_quads;
  result = mb->get_adjacencies(interior_nodes, 2, true, interior_quads, Interface::UNION);
  if (MB_SUCCESS != result) 
    return result;

    // get a list of quads generated adjacent to the exterior nodes
  Range temp_quads, exterior_quads;
  result = mb->get_adjacencies(nodes, 2, true, temp_quads, Interface::UNION);
  if (MB_SUCCESS != result) 
    return result;

    // now remove any interior quads from the previous quads list
    // and we should be left with exterior quads
  std::set_difference(temp_quads.begin(), temp_quads.end(),
                      interior_quads.begin(), interior_quads.end(),
                      range_inserter(exterior_quads));

    // check to make sure we have the right number of quads; for hexes, should be
    // .5(6*num_hexes - num_exterior_quads)
  unsigned int num_expected_int_quads = (6*num_hexes - exterior_quads.size())/2;
  if (num_expected_int_quads != interior_quads.size())
    return MB_FAILURE;

    // delete the interior quads
  result = mb->delete_entities(interior_quads);
  if (MB_SUCCESS != result) 
    return result;

  Range remaining_quads;
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
  
ErrorCode mb_adjacencies_create_delete_test() 
{
  Core moab;
  Interface* mb = &moab;
  ErrorCode rval;

  double coords[] = { 0, 0, 0, 2, 0, 0, 1, 2, 0 };
  Range verts;
  rval = mb->create_vertices( coords, 3, verts );
  if (MB_SUCCESS != rval)
    return rval;
  if (verts.size() != 3)
    return MB_FAILURE;
  EntityHandle vert_arr[3];

  EntityHandle tri;
  std::copy( verts.begin(), verts.end(), vert_arr );
  rval = mb->create_element( MBTRI, vert_arr, 3, tri );
  if (MB_SUCCESS != rval)
    return rval;
  
  vert_arr[2] = vert_arr[0];
  EntityHandle forward_edge, reverse_edge;
  rval = mb->create_element( MBEDGE, vert_arr, 2, forward_edge );
  if (MB_SUCCESS != rval)
    return rval;
    
  std::vector<EntityHandle> results;
  rval = mb->get_adjacencies( &forward_edge, 1, 2, false, results );
  if (results.size() != 1 || results.front() != tri) {
    std::cerr << "Adjacency query from forward edge to tri failed at "
              << __FILE__ <<  ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  results.clear();
  rval = mb->get_adjacencies( &tri, 1, 1, false, results );
  if (results.size() != 1 || results.front() != forward_edge) {
    std::cerr << "Adjacency query from tri to forward edge failed at "
              << __FILE__ <<  ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  
  rval = mb->delete_entities( &forward_edge, 1 );
  if (MB_SUCCESS != rval)
    return rval;
  
  results.clear();
  rval = mb->get_adjacencies( &tri, 1, 1, false, results );
  if (!results.empty()) {
    std::cerr << "Adjacency query from tri returned non-existent edge at "
              << __FILE__ <<  ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  
  rval = mb->create_element( MBEDGE, vert_arr+1, 2, reverse_edge );
  if (MB_SUCCESS != rval)
    return rval;
    
  results.clear();
  rval = mb->get_adjacencies( &reverse_edge, 1, 2, false, results );
  if (results.size() != 1 || results.front() != tri) {
    std::cerr << "Adjacency query from reverse edge to tri failed at "
              << __FILE__ <<  ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  results.clear();
  rval = mb->get_adjacencies( &tri, 1, 1, false, results );
  if (results.size() != 1 || results.front() != reverse_edge) {
    std::cerr << "Adjacency query from tri to reverse edge failed at "
              << __FILE__ <<  ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

static ErrorCode create_two_hex_full_mesh( Interface* mb,
                                             EntityHandle vertices[12],
                                             EntityHandle hexes[2],
                                             EntityHandle hex1_faces[6],
                                             EntityHandle hex2_faces[6],
                                             EntityHandle hex1_edges[12],
                                             EntityHandle hex2_edges[12] )
{
  ErrorCode rval;
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
  EntityHandle hex1_conn[] = { vertices[6], vertices[7], vertices[1], vertices[0],
                                 vertices[9], vertices[10],vertices[4], vertices[3] };
  EntityHandle hex2_conn[] = { vertices[7], vertices[8], vertices[2], vertices[1],
                                 vertices[10],vertices[11],vertices[5], vertices[4] };
  EntityHandle shared_quad_conn[] = { vertices[7], vertices[1], vertices[4], vertices[10] };
  EntityHandle hex1_face_conn[][4] = {
    { vertices[6], vertices[7], vertices[10],vertices[9] },
    { vertices[7], vertices[6], vertices[0], vertices[1] },
    { vertices[1], vertices[0], vertices[3], vertices[4] },
    { vertices[9], vertices[10],vertices[4], vertices[3] },
    { vertices[3], vertices[0], vertices[6], vertices[9] } };
  EntityHandle hex2_face_conn[][4] = {
    { vertices[7], vertices[8], vertices[11],vertices[10] },
    { vertices[8], vertices[7], vertices[1], vertices[2] },
    { vertices[2], vertices[1], vertices[4], vertices[5] },
    { vertices[10],vertices[11],vertices[5], vertices[4] },
    { vertices[5], vertices[2], vertices[8], vertices[11] } };
  EntityHandle shared_edge_conn[][2] = { { vertices[1], vertices[7] },
                                           { vertices[7], vertices[10]},
                                           { vertices[10],vertices[4] },
                                           { vertices[4], vertices[1] } };
  EntityHandle hex1_edge_conn[][2] = { { vertices[6], vertices[7] },
                                         { vertices[9], vertices[10] },
                                         { vertices[3], vertices[4] },
                                         { vertices[0], vertices[1] },
                                         { vertices[6], vertices[9] },
                                         { vertices[9], vertices[3] },
                                         { vertices[3], vertices[0] },
                                         { vertices[0], vertices[6] } };
  EntityHandle hex2_edge_conn[][2] = { { vertices[7], vertices[8] },
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
    rval = mb->create_element( MBQUAD, hex1_face_conn[i], 4, hex1_faces[i+1] );
    if (MB_SUCCESS != rval)
      return rval;
    rval = mb->create_element( MBQUAD, hex2_face_conn[i], 4, hex2_faces[i+1] );
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

ErrorCode mb_upward_adjacencies_test() 
{
  ErrorCode rval;
  Core moab;
  Interface* mb = &moab;
  
  // create a simple mesh containing 2 hexes
  EntityHandle vertices[12], hexes[2], hex1_faces[6], hex2_faces[6], hex1_edges[12], hex2_edges[12];
  rval = create_two_hex_full_mesh( mb, vertices, hexes, hex1_faces, hex2_faces, hex1_edges, hex2_edges );
  CHKERR(rval);

    // test adjacences from dim to 3
  for (int dim = 0; dim < 3; ++dim) {
    std::vector<EntityHandle> hex1_ent, hex2_ent, shared;
    const EntityHandle *list1, *list2;
    int n;
    switch (dim) {
      case 0:
        rval = mb->get_connectivity( hexes[0], list1, n );
        CHKERR(rval);
        rval = mb->get_connectivity( hexes[1], list2, n );
        CHKERR(rval);
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
      // split entities into those uniquely in hex1, those uniquely in hex2
      // and those shared between the two
    for (int i = 0; i < n; ++i) {
      if (std::find(list2, list2+n, list1[i]) - list2 == n)
        hex1_ent.push_back(list1[i]);
      else
        shared.push_back(list1[i]);
      if (std::find(list1, list1+n, list2[i]) - list1 == n)
        hex2_ent.push_back(list2[i]);
    }
      // for each shared entity check that get_adjacencies returns both hexes
    for (size_t j = 0; j < shared.size(); ++j) {
      std::vector<EntityHandle> adj;
      rval = mb->get_adjacencies( &shared[j], 1, 3, false, adj );
      CHKERR(rval);
      if (adj.size() != 2) {
        std::cout << "Expected 2 hexes adjacent to " << dim << "D entity " << j
                  << ". Got " << adj.size() << " hexes." << std::endl;
        return MB_FAILURE;
      }
      if (!(adj[0] == hexes[0] && adj[1] == hexes[1]) &&
          !(adj[0] == hexes[1] && adj[1] == hexes[0])) {
        std::cout << "Got incorrect hexes adjacent to " << dim << "D entity " << j << std::endl;
        return MB_FAILURE;
      }
    }
    
    for (size_t j = 0; j < hex1_ent.size(); ++j) {
      std::vector<EntityHandle> adj;
      rval = mb->get_adjacencies( &hex1_ent[j], 1, 3, false, adj );
      CHKERR(rval);
      CHECK(adj.size() == 1 && adj[0] == hexes[0]);
    }
    
    for (size_t j = 0; j < hex2_ent.size(); ++j) {
      std::vector<EntityHandle> adj;
      rval = mb->get_adjacencies( &hex2_ent[j], 1, 3, false, adj );
      CHKERR(rval);
      CHECK(adj.size() == 1 && adj[0] == hexes[1]);
    }
  }
    
    // For each edge, get adjacent faces, and for each face
    // get adjacent hexes.  Result should be the same as
    // direct query from edges to hexes
  std::vector<EntityHandle> all_edges(24);
  std::copy( hex1_edges, hex1_edges+12, all_edges.begin() );
  std::copy( hex2_edges, hex2_edges+12, all_edges.begin()+12 );
  std::sort( all_edges.begin(), all_edges.end() );
  all_edges.erase( std::unique(all_edges.begin(), all_edges.end()), all_edges.end() );
  for (size_t j = 0; j < all_edges.size(); ++j) {
    std::vector<EntityHandle> edge_hexes, edge_faces, face_hexes;
    rval = mb->get_adjacencies( &all_edges[j], 1, 3, false, edge_hexes );
    CHKERR(rval);
    rval = mb->get_adjacencies( &all_edges[j], 1, 2, false, edge_faces );
    CHKERR(rval);
    rval = mb->get_adjacencies( &edge_faces[0], edge_faces.size(), 3,
                                false, face_hexes, Interface::UNION );
    CHKERR(rval);
    if (edge_hexes.size() != face_hexes.size()) {
      std::cout << "Inconsistent adjacency data for edge " << j
                << ". edge->face->hex resulted in " << face_hexes.size() 
                << "hexes while edge->hex resulted in " << edge_hexes.size()
                << std::endl;
      return MB_FAILURE;
    }
    switch (edge_hexes.size()) {
      case 1:
        CHECK(edge_hexes[0] == face_hexes[0]);
        break;
      case 2:
        CHECK((edge_hexes[0] == face_hexes[0] && edge_hexes[1] == face_hexes[1]) ||
              (edge_hexes[0] == face_hexes[1] && edge_hexes[1] == face_hexes[0]));
        break;
      default:
        std::cout << "Got " << edge_hexes.size() << " hexes adjacent to edge " << j << std::endl;
        return MB_FAILURE;
    }
  }
  return MB_SUCCESS;
}

ErrorCode mb_adjacent_create_test() 
{
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  
    // create vertices
  const double coords[][3] = 
    { {-0.5, -0.5,  0.5 },
      {-0.5, -0.5, -0.5 },
      {-0.5,  0.5, -0.5 },
      {-0.5,  0.5,  0.5 },
      { 0.5, -0.5,  0.5 },
      { 0.5, -0.5, -0.5 },
      { 0.5,  0.5, -0.5 },
      { 0.5,  0.5,  0.5 } };
  EntityHandle verts[8] = {0};
  for (int i = 0; i < 8; ++i) {
    rval = mb.create_vertex( coords[i], verts[i] );
    CHKERR(rval);
  }
    // create a single hex
  const EntityHandle hconn[8] = { verts[0], verts[1], verts[2], verts[3], 
                                  verts[4], verts[5], verts[6], verts[7] };
  EntityHandle hex;
  rval = mb.create_element( MBHEX, hconn, 8, hex );
  CHKERR(rval);
    // create hex faces
  std::vector<EntityHandle> quads;
  rval = mb.get_adjacencies( &hex, 1, 2, true, quads, Interface::UNION );
  CHKERR(rval);
  CHECK_EQUAL( (size_t)6, quads.size() );
    // check that we got each of the 6 expected faces, with outwards
    // normals assuming CCW order and correct connectivity
  const EntityHandle faces[6][4] = {
    { verts[0], verts[1], verts[5], verts[4] },
    { verts[1], verts[2], verts[6], verts[5] },
    { verts[2], verts[3], verts[7], verts[6] },
    { verts[3], verts[0], verts[4], verts[7] },
    { verts[3], verts[2], verts[1], verts[0] },
    { verts[4], verts[5], verts[6], verts[7] } };
  for (int i = 0; i < 6; ++i) { // for each expected face
    // get sorted list of verts first for easy comparison
    std::vector<EntityHandle> exp_sorted(4);
    std::copy( faces[i], faces[i]+4, exp_sorted.begin() );
    std::sort( exp_sorted.begin(), exp_sorted.end() );
    // search for matching face in output
    int j = 0;
    std::vector<EntityHandle> conn;
    for (; j < 6; ++j) {
      conn.clear();
      rval = mb.get_connectivity( &quads[j], 1, conn );
      CHKERR(rval);
      CHECK_EQUAL( (size_t)4, conn.size() );
      std::vector<EntityHandle> sorted(conn);
      std::sort( sorted.begin(), sorted.end() );
      if (sorted == exp_sorted)
        break;
    }
    if (j == 6) {
      std::cerr << "No adjacent face matching hex face " << i << std::endl;
      CHECK(j<6);
    }
    // check order
    int k = std::find( conn.begin(), conn.end(), faces[i][0] ) - conn.begin();
    for (j = 1; j < 4; ++j)
      if (faces[i][j] != conn[(j+k)%4])
        break;
    if (j != 4) {
      std::cerr << "Incorrect vertex order for hex face " << i << std::endl;
      std::cerr << "Expected: " << faces[i][0] << ", " << faces[i][1] << ", "
                                << faces[i][2] << ", " << faces[i][3] << std::endl;
      std::cerr << "Actual:   " << conn[0] << ", " << conn[1] << ", "
                                << conn[2] << ", " << conn[3] << std::endl;
      CHECK(false);
    }
  }
  return MB_SUCCESS;   
}

ErrorCode nothing_but_type( Range& range, EntityType type )
{ 

    //make sure there's nothing but hexes in hex_ms
  Range::iterator iter, end_iter;
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

ErrorCode check_esets(Interface * MB, const int num_sets) 
{
  int entity_sets_size;
  ErrorCode result = MB->get_number_entities_by_type(0, MBENTITYSET, entity_sets_size);
  if (MB_SUCCESS != result ||
      entity_sets_size != num_sets) return MB_FAILURE;
  
  return MB_SUCCESS;
}

ErrorCode mb_mesh_sets_test(int flags)
{
  Core moab;
  Interface* MB = &moab;
  ErrorCode result = create_some_mesh( MB );
  if (MB_SUCCESS != result)
    return result;

  Range temp_range;
  std::vector<EntityHandle> temp_vector;
  EntityType ent_type;

  EntityHandle ms_array[MBENTITYSET] = {0};
  unsigned int number_array[MBENTITYSET] = {0};
  unsigned int num_dim_array[4] = { 0, 0, 0, 0 };
  int count, start_num_sets;

  result = MB->get_number_entities_by_type(0, MBENTITYSET, start_num_sets);
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
    num_dim_array[CN::Dimension(ent_type)] += count;

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
    result = MB->get_entities_by_dimension( ms_array[ent_type], CN::Dimension(ent_type), temp_range);
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
        
    result = MB->get_number_entities_by_dimension( ms_array[ent_type], CN::Dimension(ent_type), count );
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
  EntityHandle recursive1, recursive2;
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

  EntityHandle temp_ms1, temp_ms2; 
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
  EntityHandle temp_ms3;
  result = MB->create_meshset(flags, temp_ms3);
  if(result  != MB_SUCCESS ) 
    return result;

  EntityHandle handle_array[] = {1, 2, 3, 4, 5, 7, 8, 9, 10};
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

static bool compare_lists( std::vector<EntityHandle> vect,
                           const EntityHandle* array, 
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
ErrorCode mb_mesh_set_parent_child_test()
{
  Core moab;
  Interface *MB = &moab;

  ErrorCode rval;
  std::vector<EntityHandle> list;
  Range range;
  Range::iterator iter;
  int count;

    // create a few mesh sets
  const int num_sets = 10;
  EntityHandle sets[num_sets];
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
  
ErrorCode mb_mesh_sets_set_test()
{
  return mb_mesh_sets_test( MESHSET_SET );
}
  
ErrorCode mb_mesh_sets_list_test()
{
  return mb_mesh_sets_test( MESHSET_ORDERED );
} 

// Verify that all query functions *append* to output Range
ErrorCode mb_mesh_set_appends( int flags )
{
  Core moab;
  Interface* mb = &moab;
  ErrorCode rval = create_some_mesh( mb );
  if (MB_SUCCESS != rval)
    return rval;
  
    // get all handles and subdivide into vertex and non-vertex ents
  Range all_ents, verts, elems, results;
  rval = mb->get_entities_by_handle( 0, all_ents );
  if (MB_SUCCESS != rval)
    return rval;
  Range::iterator ve = all_ents.upper_bound( MBVERTEX );
  verts.merge( all_ents.begin(), ve );
  elems.merge( ve, all_ents.end() );

    // If we're not testing queries from the root set, 
    // create a set containing all the vertices
  EntityHandle set = 0;
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
  EntityHandle entity = verts.front();
  Range expected( elems );
  expected.insert( entity );
  
  Tag sparse, dense;
  const int zero = 0, one = 1;
  const void* vals[] = {&one};

    // Test get_entities_by_type_and_tag w/ sparse tag and no value
  rval = mb->tag_get_handle( "mb_mesh_set_appends_sparse", 
                         1, MB_TYPE_INTEGER, sparse, MB_TAG_SPARSE|MB_TAG_EXCL );
  if (MB_SUCCESS != rval)
    return rval;
  rval = mb->tag_set_data( sparse, &entity, 1, &one );
  if (MB_SUCCESS != rval)
    return rval;
  results = elems;
  rval = mb->get_entities_by_type_and_tag( set, 
                                           TYPE_FROM_HANDLE(entity),
                                           &sparse, 0, 1, 
                                           results, Interface::UNION );
  if (MB_SUCCESS != rval)
    return rval;
  if (results != expected)
    return MB_FAILURE;
    // Test again, but specify tag value
  results = elems;
  rval = mb->get_entities_by_type_and_tag( set, 
                                           TYPE_FROM_HANDLE(entity),
                                           &sparse, vals, 1, 
                                           results, Interface::UNION );
  if (MB_SUCCESS != rval)
    return rval;
  if (results != expected)
    return MB_FAILURE;
   
    // Test get_entities_by_type_and_tag w/ dense tag
  rval = mb->tag_get_handle( "mb_mesh_set_appends_dense", 
                         1, MB_TYPE_INTEGER, dense,
                         MB_TAG_DENSE|MB_TAG_EXCL, &zero );
  if (MB_SUCCESS != rval)
    return rval;
  rval = mb->tag_set_data( dense, &entity, 1, &one );
  if (MB_SUCCESS != rval)
    return rval;
  results = elems;
  rval = mb->get_entities_by_type_and_tag( set, 
                                           TYPE_FROM_HANDLE(entity),
                                           &dense, vals, 1, 
                                           results, Interface::UNION );
  if (MB_SUCCESS != rval)
    return rval;
  if (results != expected)
    return MB_FAILURE;
 
  return MB_SUCCESS;
}

ErrorCode mb_mesh_set_set_appends()
{
  return mb_mesh_set_appends( MESHSET_SET );
}

ErrorCode mb_mesh_set_list_appends()
{
  return mb_mesh_set_appends( MESHSET_ORDERED );
}

ErrorCode mb_mesh_set_root_appends()
{
  return mb_mesh_set_appends( -1 );
}

ErrorCode mb_mesh_set_set_replace_test()
{
  Core moab;
  Interface* mb = &moab;
  ErrorCode rval;
  Range r;
    // create 10 vertices to put in set
  std::vector<double> coords(30);
  rval = mb->create_vertices( &coords[0], 10, r );
  CHKERR(rval);
  std::vector<EntityHandle> verts(r.size());
  std::copy( r.begin(), r.end(), verts.begin() );
  r.clear();
    // create a set
  EntityHandle set;
  rval = mb->create_meshset( MESHSET_SET, set );
  CHKERR(rval);
    // put every other vertex in set
  for (size_t i = 0; i < 10; i += 2)
    r.insert( verts[i] );
  rval = mb->add_entities( set, r );
  CHKERR(rval);
  r.clear();
    // swap 3 of the vertices
  EntityHandle old_ents[3] = { verts[2], verts[4], verts[6] };
  EntityHandle new_ents[3] = { verts[1], verts[9], verts[5] };
  rval = mb->replace_entities( set, old_ents, new_ents, 3 );
  CHKERR(rval);
    // check new set contents
  rval = mb->get_entities_by_handle( set, r );
  CHKERR(rval);
  Range r2;
  r2.insert( verts[0] );
  r2.insert( verts[1] );
  r2.insert( verts[9] );
  r2.insert( verts[5] );
  r2.insert( verts[8] );
  if (r != r2) {
    std::cerr << "Range does not contain expected values." << std::endl
              << "  Expected: " << r2 << std::endl
              << "  Actual  : " << r  << std::endl;
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

ErrorCode mb_mesh_set_list_replace_test()
{
  Core moab;
  Interface* mb = &moab;
  ErrorCode rval;
    // create 10 vertices to put in set
  Range r;
  std::vector<double> coords(30);
  rval = mb->create_vertices( &coords[0], 10, r );
  CHKERR(rval);
  std::vector<EntityHandle> verts(r.size());
  std::copy( r.begin(), r.end(), verts.begin() );
  r.clear();
    // create a set
  EntityHandle set;
  rval = mb->create_meshset( MESHSET_ORDERED, set );
  CHKERR(rval);
    // put all vertices in set, but add the first one a second time
  std::vector<EntityHandle> list( verts );
  list.push_back( verts.front() );
  rval = mb->add_entities( set, &list[0], list.size() );
    // swap 3 of the vertices
  EntityHandle old_ents[3] = { verts[2], verts[4], verts[6] };
  EntityHandle new_ents[3] = { verts[1], verts[9], verts[5] };
  rval = mb->replace_entities( set, old_ents, new_ents, 3 );
  CHKERR(rval);
    // check new set contents
  std::vector<EntityHandle> list2;
  rval = mb->get_entities_by_handle( set, list2 );
  CHKERR(rval);
  list[2] = verts[1];
  list[4] = verts[9];
  list[6] = verts[5];
  if (list != list2) {
    std::cerr << "Range does not contain expected values." << std::endl;
    std::cerr << "  Expected: ";
    std::copy( list.begin(), list.end(), std::ostream_iterator<EntityHandle>(std::cerr, " ") );
    std::cerr << std::endl << "  Actual  : ";
    std::copy( list2.begin(), list2.end(), std::ostream_iterator<EntityHandle>(std::cerr, " ") );
    std::cerr << std::endl;
    return MB_FAILURE;
  }
    // now try replacing a repeated value
  rval = mb->replace_entities( set, &verts[0], &verts[3], 1 );
  CHKERR(rval);
  list[0] = list[10] = verts[3];
  list2.clear();
  rval = mb->get_entities_by_handle( set, list2 );
  CHKERR(rval);
  if (list != list2) {
    std::cerr << "Range does not contain expected values." << std::endl;
    std::cerr << "  Expected: ";
    std::copy( list.begin(), list.end(), std::ostream_iterator<EntityHandle>(std::cerr, " ") );
    std::cerr << std::endl << "  Actual  : ";
    std::copy( list2.begin(), list2.end(), std::ostream_iterator<EntityHandle>(std::cerr, " ") );
    std::cerr << std::endl;
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

/* Test the following changes to a meshset:
  set       MB-> tracking
  tracking  MB-> set
  unordered MB-> ordered
  ordered   MB-> unordered
*/
ErrorCode mb_mesh_set_flag_test() 
{
  Core moab;
  Interface* mb = &moab;

  ErrorCode rval;
    // create 10 vertices to put in set
  Range verts;
  std::vector<double> coords(30);
  rval = mb->create_vertices( &coords[0], 10, verts );
  CHKERR(rval);
  
  // CHECK SET->TRACKING
  // create a set and add the verts
  EntityHandle set;
  rval = mb->create_meshset( MESHSET_SET, set );
  CHKERR(rval);
  rval = mb->add_entities( set, verts);
  CHKERR(rval);
  // the verts should not be tracking adjacencies
  Range adj_sets;
  rval = mb->get_adjacencies( verts, 4, false, adj_sets);
  CHKERR(rval);
  if(!adj_sets.empty()) {
    std::cerr << "range should be empty but contains:" << std::endl;
    rval = mb->list_entities( adj_sets );
    return MB_FAILURE;
  }
  // check to make sure the flags on MESHSET_SET
  unsigned int flags;
  rval = mb->get_meshset_options( set, flags );
  CHKERR(rval);
  if(!MESHSET_SET&flags || MESHSET_TRACK_OWNER&flags || MESHSET_ORDERED&flags){
    std::cerr << "set should be MESHSET_SET only, flags=" << flags << std::endl;
    return MB_FAILURE;
  }
  // change to a tracking set and check flags
  rval = mb->set_meshset_options( set, MESHSET_TRACK_OWNER);
  CHKERR(rval);
  rval = mb->get_meshset_options( set, flags );
  CHKERR(rval);
  if(MESHSET_SET&flags || !MESHSET_TRACK_OWNER&flags || MESHSET_ORDERED&flags){
    std::cerr << "set should be MESHSET_TRACK_OWNER only, flags=" << flags 
              << std::endl;
    return MB_FAILURE;
  }
  // check adjacencies
  rval = mb->get_adjacencies( verts, 4, false, adj_sets);
  CHKERR(rval);
  if(1 != adj_sets.size()) {
    std::cerr << "range should contain a set, adj_sets.size()=" 
              << adj_sets.size() << std::endl;
    rval = mb->list_entities( adj_sets );
    return MB_FAILURE;
  }

  // CHECK TRACKING->SET
  // change to a standard set and check flags
  rval = mb->set_meshset_options( set, MESHSET_SET);
  CHKERR(rval);
  rval = mb->get_meshset_options( set, flags );
  CHKERR(rval);
  if(!MESHSET_SET&flags || MESHSET_TRACK_OWNER&flags || MESHSET_ORDERED&flags){
    std::cerr << "set should be MESHSET_SET only, flags=" << flags 
              << std::endl;
    return MB_FAILURE;
  }
  // the set should no longer be adjacent to the vertices
  adj_sets.clear();
  rval = mb->get_adjacencies( verts, 4, false, adj_sets);
  CHKERR(rval);
  if(!adj_sets.empty()) {
    std::cerr << "range should be empty but contains:" << std::endl;
    rval = mb->list_entities( adj_sets );
    return MB_FAILURE;
  }
  // CHECK UNORDERED->ORDERED
  // add a duplicate vert
  rval = mb->add_entities( set, &verts.front(), 1);
  CHKERR(rval);
  // unordered sets cannot hold duplicates so size shouldn't change
  std::vector<EntityHandle> entities;
  rval = mb->get_entities_by_handle( set, entities );
  if(10 != entities.size()) {
    std::cerr << "set should not hold duplicate entities" << std::endl;
    return MB_FAILURE;
  }
  // change to an ordered set and check flags
  rval = mb->set_meshset_options( set, MESHSET_ORDERED);
  CHKERR(rval);
  rval = mb->get_meshset_options( set, flags );
  CHKERR(rval);
  if(MESHSET_SET&flags || MESHSET_TRACK_OWNER&flags || !MESHSET_ORDERED&flags){
    std::cerr << "set should be MESHSET_ORDERED only, flags=" << flags 
              << std::endl;
    return MB_FAILURE;
  }
  // swap the order with some entities to that the handles aren't ordered
  rval = mb->clear_meshset( &set, 1 );
  CHKERR(rval);
  entities.clear();
  entities[0] = verts[1];
  entities[1] = verts[0];
  rval = mb->add_entities( set, &entities[0], 2);
  CHKERR(rval);
  // check to ensure the entities keep their order
  entities.clear();
  rval = mb->get_entities_by_handle( set, entities );
  if(verts[0]!=entities[1] || verts[1]!=entities[0]) {
    std::cerr << "ordered set did not keep its order" << std::endl;
    return MB_FAILURE;
  }

  // CHECK ORDERED->UNORDERED
  // change to an unordered set and check flags
  rval = mb->set_meshset_options( set, MESHSET_SET);
  CHKERR(rval);
  rval = mb->get_meshset_options( set, flags );
  CHKERR(rval);
  if(!MESHSET_SET&flags || MESHSET_TRACK_OWNER&flags || MESHSET_ORDERED&flags){
    std::cerr << "set should be MESHSET_SET only, flags=" << flags 
              << std::endl;
    return MB_FAILURE;
  }
  // the entities in the set should now be ordered by handle
  entities.clear();
  rval = mb->get_entities_by_handle( set, entities );
  if(verts[0]!=entities[0] || verts[1]!=entities[1]) {
    std::cerr << "unordered set is still ordered" << std::endl;
    return MB_FAILURE;
  }
  return MB_SUCCESS;
}

  // number of entities of type MBVERTEX, MBEDGE, MBDTri, MBQUAD, MBTET, and MBHEX
  // in mbtest1.g  (all other values are 0.
static const EntityType types[] = { MBVERTEX, MBEDGE, MBTRI, MBQUAD, MBTET, MBHEX };
const int num_types = sizeof(types)/sizeof(types[0]);
static const unsigned int num_entities[num_types+1] = {47,12,18,8,22,8,0};

ErrorCode mb_delete_mesh_test()
{
  Core moab;
  Interface* gMB = &moab;
  ErrorCode error = load_file_one( gMB );
  if (MB_SUCCESS != error)
    return error;

    // Lets also test the global MB pointer (gMB) here.
  error = gMB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;

    // load the mesh again 
  error = load_file_one( gMB );
  if (error != MB_SUCCESS)
    return error;


    // step through each element type
  for (EntityType type = MBVERTEX; type != MBENTITYSET; type++)
  {
    // There should be entities
    Range entities;
    error = gMB->get_entities_by_type(0,  type, entities);
    if (error != MB_SUCCESS)
      return error;

    size_t idx = std::find( types, types+num_types, type ) - types;
    if ( entities.size() != num_entities[idx] )
      return MB_FAILURE;
  }

  return MB_SUCCESS;
}


ErrorCode mb_mesh_set_tracking_test()
{
  Core moab;
  Interface* MB = &moab;

    //read in a file so you have some data in the database
  
  ErrorCode error = load_file_one( MB );
  if (error != MB_SUCCESS)
    return error;

  EntityHandle ms1, ms2, ms3;

    //create meshsets 
  ErrorCode result = MB->create_meshset( MESHSET_TRACK_OWNER | MESHSET_ORDERED, ms1 ) ;
  if(result != MB_SUCCESS )
    return result;
  result = MB->create_meshset( MESHSET_SET | MESHSET_TRACK_OWNER, ms2 ) ;
  if(result != MB_SUCCESS )
    return result;
  result = MB->create_meshset( MESHSET_SET | MESHSET_TRACK_OWNER, ms3 ) ;
  if(result != MB_SUCCESS )
    return result;

    // get all hexes 
  Range hexes;
  result = MB->get_entities_by_type(0, MBHEX, hexes);
  if(result != MB_SUCCESS )
    return result;

    // get all tris 
  Range tris;
  result = MB->get_entities_by_type(0, MBTRI, tris );
  if(result != MB_SUCCESS )
    return result;

    // get all tets 
  Range temp_range;
  result = MB->get_entities_by_type(0, MBTET, temp_range);
  if(result != MB_SUCCESS )
    return result;

    //copy all the tets from the range into a vector 'tets'
  std::vector<EntityHandle> tets( temp_range.size() );
  std::copy(temp_range.begin(), temp_range.end(), tets.begin() );

    //Quick check on 'get_entities_by_dimension()'
  Range dim_3_range;
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

  Range::iterator iter;
  iter = tris.begin();

  std::vector< EntityHandle > temp_vec;

    //ask a tri which meshsets it is in
  result = MB->get_adjacencies( &(*iter), 1, 4, false, temp_vec ) ;
  if(result != MB_SUCCESS )
    return result;

    //cout<<"tris temp_vec.size() = "<<temp_vec.size()<<endl; 
  if( temp_vec.size() != 2 )
    return MB_FAILURE;

    //ask a tet which meshsets it is in
  temp_vec.clear();
  std::vector<EntityHandle>::iterator vec_iter = tets.begin();
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


// Compare internal representation of contents for a list (MESHSET_ORDERED)
// set to expected contents.  Assumes expected contents are correctly
// ordered.
static ErrorCode check_list_meshset_internal( Interface& mb,
                                              const EntityHandle* expected,
                                              int num_expected,
                                              const EntityHandle* contents,
                                              int length )
{
  bool okay = true;
  for (int i = 0; i < std::min(num_expected, length); ++i) {
    if (expected[i] != contents[i]) {
      std::cerr << "List set contents differ at index " << i 
                << ": expected " << expected[i] << " but got " 
                << contents[i] << std::endl;
      okay = false;
    }
  }
  if (num_expected > length) {
    std::cerr << "List set is missing " << num_expected - length << 
                 "handles" << std::endl;
    okay = false;
  }
  else if (length > num_expected) {
    std::cerr << "List set has " << num_expected - length << 
                 " extra handles" << std::endl;
    okay = false;
  }
  
  if (okay)
    return MB_SUCCESS;
  
  std::cerr << "Expected contents: ";
  if (!num_expected)
    std::cerr << "(empty)";
  else
    std::cerr << expected[0];
  for (int i = 1; i < num_expected; ++i)
    std::cerr << ", " << expected[i];
  std::cerr << std::endl;
  
  std::cerr << "Actual contents: ";
  if (!length)
    std::cerr << "(empty)";
  else
    std::cerr << contents[0];
  for (int i = 1; i < length; ++i)
    std::cerr << ", " << contents[i];
  std::cerr << std::endl;

  return MB_FAILURE;
}


// Compare internal representation of contents for a ranged (MESHSET_SET)
// set to expected contents.  Assumes expected contents are correctly
// ordered.
static ErrorCode check_ranged_meshset_internal( Interface& mb,
                                                const EntityHandle* expected,
                                                int num_expected,
                                                const EntityHandle* contents,
                                                int length )
{
  if (length % 2) {
    std::cerr << "Range set is corrupt.  Odd number of handles in content list" << std::endl;
    std::cerr << "Actual contents: " << contents[0];
    for (int i = 1; i < length; ++i)
      std::cerr << ", " << contents[i];
    std::cerr << std::endl;
    return MB_FAILURE;
  }
  bool okay = true;
    // check that all range pairs are valid (first handle less than or 
    // equal to second)
  for (int i = 0; i < length; i += 2) {
    if (contents[i] > contents[i+1]) {
      std::cerr << "Range set has invalid range pair at index " << i
                << ": [" << contents[i] << ',' << contents[i+1]
                << ']' << std::endl;
      okay = false;
    }
  }
    // check that all range pairs are sorted and non-overlapping
    // (start of a range must be greater than end of previous range)
  for (int i = 2; i < length; i += 2) {
    if (contents[i] < contents[i-1]) {
      std::cerr << "Range set has incorrectly ordered ranges at index " << i
                << ": [...," << contents[i-1] << "], [" << contents[i]
                << ",...]" << std::endl;
      okay = false;
    }
    if (contents[i] == contents[i-1]+1) {
      std::cerr << "Range set has pairs that should have been merged at index " << i
                << ": [...," << contents[i-1] << "], [" << contents[i]
                << ",...]" << std::endl;
      okay = false;
    }
  }
  if (!okay) {
    std::cerr << "Actual contents: ";
    if (!length)
      std::cerr << "(empty)";
    else
      std::cerr << '[' << contents[0] << ',' << contents[1] << ']';
    for (int i = 2; i < length; i += 2) 
      std::cerr << ", [" << contents[i] << ',' << contents[i+1] << ']';
    std::cerr << std::endl;
    return MB_FAILURE;
  }
  
  int j = 0;
  for (int i = 0; i < length; i += 2) {
    for (; j < num_expected && expected[j] < contents[i]; ++j) {
      std::cerr << "Range set missing expected handle: " << expected[j] << std::endl;
      okay = false;
    }
    int k = j;
    while (k < num_expected && expected[k] <= contents[i+1])
      ++k;
    if ((EntityHandle)(k-j) <= (contents[i+1]-contents[i])) {
      std::cerr << "Handle range [" << contents[i] << ',' << contents[i+1]
                << "] contains unexpected handles.  Expected handles: ";
      if (k == j)
        std::cerr << "(none)" << std::endl;
      else {
        std::cerr << expected[j];
        for (++j; j < k; ++j)
          std::cerr << ", " << expected[j];
        std::cerr << std::endl;
      }
      okay  = false;
    }
    j = k;
  }
  
  if (okay)
    return MB_SUCCESS;
    
  std::cerr << "Expected contents: ";
  if (!num_expected)
    std::cerr << "(empty)";
  else
    std::cerr << expected[0];
  for (int i = 1; i < num_expected; ++i)
    std::cerr << ", " << expected[i];
  std::cerr << std::endl;
  
  std::cerr << "Actual contents: ";
  if (!length)
    std::cerr << "(empty)";
  else
    std::cerr << '[' << contents[0] << ',' << contents[1] << ']';
  for (int i = 2; i < length; i += 2) 
    std::cerr << ", [" << contents[i] << ',' << contents[i+1] << ']';
  std::cerr << std::endl;
  
  return MB_FAILURE;
}

// Check the internal representation of a meshset
// to verify that it is correct.
static ErrorCode check_meshset_internal( Interface& mb,
                                         EntityHandle set,
                                         const EntityHandle* expected,
                                         int num_expected )
{
  ErrorCode rval;
  WriteUtilIface* tool = 0;
  rval = mb.query_interface( tool );
  CHKERR(rval);
  
  const EntityHandle* contents;
  int length;
  unsigned char flags;
  rval = tool->get_entity_list_pointers( &set, 1, &contents, 
                      WriteUtilIface::CONTENTS, &length, &flags );
  ErrorCode rval1 = mb.release_interface( tool );
  CHKERR(rval);
  CHKERR(rval1);
  
  if (flags & MESHSET_ORDERED)
    rval = check_list_meshset_internal( mb, expected, num_expected, contents, length );
  else
    rval = check_ranged_meshset_internal( mb, expected, num_expected, contents, length );
  CHKERR(rval);
  return MB_SUCCESS;
}

ErrorCode mb_mesh_set_set_add_remove_test()
{
  Core core;
  Interface& mb = core;
  EntityHandle set;
  ErrorCode rval = mb.create_meshset( MESHSET_SET, set );
  CHKERR(rval);
  
  EntityHandle list1[] = {10, 16, 18, 20, 24, 27};
  size_t len1 = sizeof(list1)/sizeof(list1[0]);
  EntityHandle list2[] = {10, 16, 17, 18, 19, 20, 24, 27};
  size_t len2 = sizeof(list2)/sizeof(list2[0]);
  rval = mb.add_entities( set, list1, len1 );
  CHKERR(rval);
  rval = check_meshset_internal( mb, set, list1, len2 );
  CHKERR(rval);
  rval = mb.add_entities( set, list2, len2 );
  CHKERR(rval);
  EntityHandle exp12[] = {10, 16, 17, 18, 19, 20, 24, 27};
  size_t len12 = sizeof(exp12)/sizeof(exp12[0]);
  rval = check_meshset_internal( mb, set, exp12, len12 );
  CHKERR(rval);
  
  EntityHandle list3[] = { 15, 16, 18, 20, 21, 24, 28 };
  size_t len3 = sizeof(list3)/sizeof(list3[0]);
  rval = mb.remove_entities( set, list3, len3 );
  CHKERR(rval);
  EntityHandle exp123[] = { 10, 17, 19, 27 };
  size_t len123 = sizeof(exp123)/sizeof(exp123[0]);
  rval = check_meshset_internal( mb, set, exp123, len123 );
  CHKERR(rval);
  
  EntityHandle list4[] = { 18, 10, 11, 12, 13, 14, 15, 16 };
  size_t len4 = sizeof(list4)/sizeof(list4[0]);
  rval = mb.add_entities( set, list4, len4 );
  CHKERR(rval);
  EntityHandle exp14[] = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 27};
  size_t len14 = sizeof(exp14)/sizeof(exp14[0]);
  rval = check_meshset_internal( mb, set, exp14, len14 );
  CHKERR(rval);
  
  EntityHandle list5[] = { 9, 10, 12, 13, 14, 15, 19, 20 };
  rval = mb.remove_entities( set, list5, sizeof(list5)/sizeof(list5[0]) );
  CHKERR(rval);
  EntityHandle exp5[] = { 11, 16, 17, 18, 27 };
  rval = check_meshset_internal( mb, set, exp5, sizeof(exp5)/sizeof(exp5[0]) );
  CHKERR(rval);
  
  EntityHandle list6[] = { 9, 10, 15, 16, 18, 19, 28 };
  rval = mb.add_entities( set, list6, sizeof(list6)/sizeof(list6[0]) );
  CHKERR(rval);
  EntityHandle exp6[] = { 9, 10, 11, 15, 16, 17, 18, 19, 27, 28 };
  rval = check_meshset_internal( mb, set, exp6, sizeof(exp6)/sizeof(exp6[0]) );
  CHKERR(rval);
  
  EntityHandle list7[] = { 13, 19, 27, 28 };
  rval = mb.add_entities( set, list7, sizeof(list7)/sizeof(list7[0]) );
  CHKERR(rval);
  EntityHandle exp7[] = { 9, 10, 11, 13, 15, 16, 17, 18, 19, 27, 28 };
  rval = check_meshset_internal( mb, set, exp7, sizeof(exp7)/sizeof(exp7[0]) );
  CHKERR(rval);
  
  EntityHandle list8[] = { 12, 14, 33 };
  rval = mb.add_entities( set, list8, sizeof(list8)/sizeof(list8[0]) );
  CHKERR(rval);
  EntityHandle exp8[] = { 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 27, 28, 33 };
  rval = check_meshset_internal( mb, set, exp8, sizeof(exp8)/sizeof(exp8[0]) );
  CHKERR(rval);
  
  EntityHandle list9[] = { 29, 30, 31, 32, 34 };
  rval = mb.remove_entities( set, list9, sizeof(list9)/sizeof(list9[0]) );
  CHKERR(rval);
  rval = check_meshset_internal( mb, set, exp8, sizeof(exp8)/sizeof(exp8[0]) );
  CHKERR(rval);
  
  EntityHandle list10[] = { 9, 11, 13, 17, 18, 19, 28, 33, 100 };
  rval = mb.remove_entities( set, list10, sizeof(list10)/sizeof(list10[0]) );
  CHKERR(rval);
  EntityHandle exp10[] = { 10, 12, 14, 15, 16, 27 };
  rval = check_meshset_internal( mb, set, exp10, sizeof(exp10)/sizeof(exp10[0]) );
  CHKERR(rval);
  
  EntityHandle list11[] = { 11, 12, 13, 14, 27, 28 };
  rval = mb.remove_entities( set, list11, sizeof(list11)/sizeof(list11[0]) );
  CHKERR(rval);
  EntityHandle exp11[] = { 10, 15, 16 };
  rval = check_meshset_internal( mb, set, exp11, sizeof(exp11)/sizeof(exp11[0]) );
  CHKERR(rval);
  
  EntityHandle list12[] = { 1, 10, 15, 16 };
  rval = mb.remove_entities( set, list12, sizeof(list12)/sizeof(list12[0]) );
  CHKERR(rval);
  rval = check_meshset_internal( mb, set, 0, 0 );
  CHKERR(rval);
  
  return MB_SUCCESS;
}


ErrorCode mb_higher_order_test()
{
  Core moab;
  Interface *MB = &moab;

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
  std::vector<EntityHandle> connectivity(8);
  EntityHandle node_handle;
  int i;
  for( i=0; i<7; i++)
  {
    if(MB->create_vertex( nodes_array[i], node_handle ) != MB_SUCCESS )
      return MB_FAILURE;
    connectivity[i] = node_handle;
  }

    //create the higher order tri
  EntityHandle tri_handle;
  ErrorCode result = MB->create_element(MBTRI, &connectivity[0], 6, tri_handle);
  if(result != MB_SUCCESS)
    return result;

    //create the higher order tri
  std::vector<EntityHandle> other_conn(3);

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

  EntityHandle other_tri_handle;
  result = MB->create_element(MBTRI, &other_conn[0], 3, other_tri_handle);
  if(result != MB_SUCCESS)
    return result;

    //get the connectivity now
  std::vector<EntityHandle> retrieved_conn; 

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
  std::vector<EntityHandle> topo_conn; 
  result = MB->get_connectivity(&other_tri_handle, 1, topo_conn, true) ;
  if(result != MB_SUCCESS)
    return result;

  if (topo_conn.size() != 3)
    return MB_FAILURE;

  for ( k=0; k<3; k++)
    if (topo_conn[k] != retrieved_conn[k] )
      return MB_FAILURE;

    // short check to make sure that Core::handle_from_id() works
  unsigned long handle_id = MB->id_from_handle( node_handle); 

  EntityHandle test_handle; 
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

ErrorCode mb_bit_tags_test()
{
  Core moab;
  Interface* MB = &moab;
  ErrorCode success = create_some_mesh( MB );
  if (MB_SUCCESS != success)
    return success;

  Tag bit_tag;
  Range entities;
  MB->get_entities_by_type(0, MBVERTEX, entities);

  if(MB->tag_get_handle("bit on vertex", 3, MB_TYPE_BIT, bit_tag, MB_TAG_CREAT) != MB_SUCCESS)
  {
    cout << "couldn't create bit tag" << endl;
    return MB_FAILURE;
  }

  Range::iterator iter;
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
  
  // test range-based query for all vertices
  std::vector<unsigned char> data(entities.size());
  success = MB->tag_get_data( bit_tag, entities, &data[0] );
  if (MB_SUCCESS != success) return success;
  std::vector<unsigned char>::iterator i = data.begin();
  for (iter = entities.begin(); iter != entities.end(); ++iter, ++i) 
    if (*i != ((*iter) & 0x7))
      return MB_FAILURE;

  // test vector-based query for all vertices
  std::vector<EntityHandle> verts(entities.begin(), entities.end());
  success = MB->tag_get_data( bit_tag, &verts[0], verts.size(), &data[0] );
  if (MB_SUCCESS != success) return success;
  i = data.begin();
  for (iter = entities.begin(); iter != entities.end(); ++iter, ++i) 
    if (*i != ((*iter) & 0x7))
      return MB_FAILURE;
  
  // test default value
  const unsigned char default_bits = '\005'; // 0000 0101
  Tag tag2;
  success = MB->tag_get_handle( "bit with default", 4, MB_TYPE_BIT, tag2, MB_TAG_CREAT, &default_bits );
  if (MB_SUCCESS != success) {
    cout << "Failed to create bit tag with default value" << std::endl;
    return success;
  }
  
  // set value to zero on a single vertex
  bits = 0;
  EntityHandle zh = verts[verts.size()/2];
  success = MB->tag_set_data( tag2, &zh, 1, &bits );
  if (MB_SUCCESS != success)
    return success;
  
  // get tag values for all entities
  data.clear();
  data.resize( verts.size(), 0x7A ); // initialize with 0111 1010
  success = MB->tag_get_data( tag2, entities, &data[0] );
  if (MB_SUCCESS != success)
    return success;
  
  // check values
  i = data.begin();
  for (iter = entities.begin(); iter != entities.end(); ++iter, ++i) 
    if (*iter == zh && *i) // the one we set to zero
      return MB_FAILURE;
    else if (*iter != zh && *i != default_bits)
      return MB_FAILURE;

  return MB_SUCCESS;
}

ErrorCode mb_tags_test()
{
  Core moab;
  Interface* MB = &moab;
  ErrorCode result = load_file_one( MB );
  if (MB_SUCCESS != result)
    return result;

  Tag stale_bits, stale_dense, stale_sparse;
  result = MB->tag_get_handle("stale data", 5, MB_TYPE_BIT, stale_bits, MB_TAG_CREAT);
  if (MB_SUCCESS != result)
    return result;
     
  int def_data = 9;
  result = MB->tag_get_handle("dense stale_data", 1, MB_TYPE_INTEGER, stale_dense, MB_TAG_DENSE|MB_TAG_EXCL, &def_data);
  if (MB_SUCCESS != result)
    return result;
  result = MB->tag_get_handle("sparse stale data", 1, MB_TYPE_INTEGER, stale_sparse, MB_TAG_SPARSE|MB_TAG_EXCL);
  if (MB_SUCCESS != result)
    return result;

  double coords[3] = { 0,0,0 };
  EntityHandle stale_handle1, stale_handle2;
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
    ErrorCode stale_result = MB->tag_get_data(stale_sparse, &stale_handle2, 1, &def_data);
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
  Range entities;
  int value = 1;
  const void *dum_ptr = &value;
  Tag material_tag;
  result = MB->tag_get_handle( MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER, material_tag);
  if (MB_SUCCESS != result)
    return result;
  if(MB->get_entities_by_type_and_tag( 0, MBENTITYSET, &material_tag, 
                                       &dum_ptr, 
                                       1, entities) != MB_SUCCESS)
    return MB_FAILURE;

  if( entities.size() != 1)
    return MB_FAILURE;
 
    //add a dense tag to hexes
  Tag junk_tag;
  if(MB->tag_get_handle( "junk_tag", 1, MB_TYPE_INTEGER, junk_tag, MB_TAG_DENSE|MB_TAG_EXCL) 
     != MB_SUCCESS)
    return MB_FAILURE;    

    //Set the dense tag on 5 hexes to 3489 
  Range test_range;
  result = MB->get_entities_by_type(0,  MBHEX, test_range ) ;
  if(result != MB_SUCCESS)
    return result;

  Range::iterator iter, end_iter;
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
  EntityHandle meshset;
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
    // an empty Range

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

ErrorCode mb_common_tag_test( TagType storage )
{
  Core moab;
  Interface* mb = &moab;
  ErrorCode result = create_some_mesh( mb );
  if (MB_SUCCESS != result)
    return result;

  char tagname[64];
  sprintf( tagname, "t%d", rand() );
  
  Tag tag;
  const EntityHandle def_val = ~(EntityHandle)0;
  ErrorCode rval = mb->tag_get_handle( tagname, 
                                       1, MB_TYPE_HANDLE,
                                       tag, storage|MB_TAG_EXCL,
                                       &def_val );
  if (MB_SUCCESS != rval)
    return rval;
  
  
  Range entities;
  mb->get_entities_by_handle( 0, entities );
  if (entities.empty())
    return MB_FAILURE;
  
    // set tag on every other entity to be the entities handle
  Range::const_iterator i;
  bool odd = true;
  for (i = entities.begin(); i != entities.end(); ++i, odd = !odd) {
    if (odd) {
      const EntityHandle h = *i;
      rval = mb->tag_set_data( tag, &h, 1, &h );
      if (MB_SUCCESS != rval)
        return rval;
    }
  }
  
    // check values on every entity -- expect default for every other entity
  odd = true;
  for (i = entities.begin(); i != entities.end(); ++i, odd = !odd) {
    EntityHandle val = 0;
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
  std::vector<EntityHandle> values( entities.size() );
  std::copy( entities.begin(), entities.end(), values.begin() );
  rval = mb->tag_set_data( tag, entities, &values[0] );
  if (MB_SUCCESS != rval)
    return rval;
  
    // check values on every entity -- expect default for every other entity
  for (i = entities.begin(); i != entities.end(); ++i) {
    EntityHandle val = 0;
    rval = mb->tag_get_data( tag, &*i, 1, &val );
    if (MB_SUCCESS != rval)
      return rval;
    if (val != *i)
      return MB_FAILURE;
  }
  
    // find each entity by tag value
  for (i = entities.begin(); i != entities.end(); ++i) {
    const EntityHandle h = *i;
    const EntityType type = mb->type_from_handle( h );
    const void* const tag_vals[] = { &h };
    Range result;
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


ErrorCode mb_dense_tag_test()
{
  return mb_common_tag_test( MB_TAG_DENSE );
}

ErrorCode mb_sparse_tag_test()
{
  return mb_common_tag_test( MB_TAG_SPARSE );
}
  
// class to offset hex center nodes
class OffsetHexCenterNodes : public Interface::HONodeAddedRemoved
{
public:
  OffsetHexCenterNodes(Interface* mb, double x, double y, double z)
      : gMB(mb)
    { 
      mOffset[0] = x; mOffset[1] = y; mOffset[2] = z; 
    }
    
  ~OffsetHexCenterNodes(){}

  void node_added(EntityHandle node, EntityHandle)
    {
      gMB->get_coords(&node, 1, mCoords);
      mCoords[0] += mOffset[0];
      mCoords[1] += mOffset[1];
      mCoords[2] += mOffset[2];
      gMB->set_coords(&node, 1, mCoords);
    }

    //do nothing
  void node_removed( EntityHandle /*node*/) {}

private:
  Interface* gMB;
  double mCoords[3];
  double mOffset[3];
};

ErrorCode mb_entity_conversion_test()
{
  ErrorCode error;
  Core moab;
  Interface* MB = &moab;

    //read in a file so you have some data in the database
  std::string file_name = TestDir + "/mbtest3.g";
  error = MB->load_mesh(file_name.c_str(), NULL, 0);
  if (error != MB_SUCCESS)
    return error;

  Range entities;
  EntityHandle meshset;
  MB->create_meshset(MESHSET_SET, meshset);
  
  MB->get_entities_by_type(0, MBHEX, entities);
  MB->add_entities(meshset, entities);


  OffsetHexCenterNodes function_object(MB,0.07, 0.15, 0);

  MB->convert_entities(meshset, false, false, true, &function_object);
  if (MB_SUCCESS != check_valid_connectivity( MB ))
    return MB_FAILURE;

  file_name = "hex_mid_volume_nodes.g";
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
  if (MB_SUCCESS != check_valid_connectivity( MB ))
    return MB_FAILURE;

  file_name = "hex_mid_edge_face_vol_nodes.g";
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
  if (MB_SUCCESS != check_valid_connectivity( MB ))
    return MB_FAILURE;

  file_name = "hex_mid_edge_nodes.g";
  error = MB->write_mesh(file_name.c_str());
  if (error != MB_SUCCESS)
    return error;
  
    // convert them back to hex8's
  MB->convert_entities(meshset, false, false, false);
  if (MB_SUCCESS != check_valid_connectivity( MB ))
    return MB_FAILURE;

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
  if (MB_SUCCESS != check_valid_connectivity( MB ))
    return MB_FAILURE;

  file_name = "tet_mid_edge_nodes.g";
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
  if (MB_SUCCESS != check_valid_connectivity( MB ))
    return MB_FAILURE;

  file_name = "tet_mid_face_nodes.g";
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
  if (MB_SUCCESS != check_valid_connectivity( MB ))
    return MB_FAILURE;

  file_name = "tet_mid_edge_face_nodes.g";
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
  for(Range::iterator tet_iter = entities.begin(); tet_iter != entities.end(); ++tet_iter)
  {
    std::vector<EntityHandle> adj;
    error = MB->get_adjacencies(&(*tet_iter), 1, 2, true, adj);
    if (MB_SUCCESS != error)
      return error;
    for(std::vector<EntityHandle>::iterator tri_iter = adj.begin();
        tri_iter != adj.end(); ++tri_iter)
    {
      std::vector<EntityHandle> up_adj;
      MB->get_adjacencies(&(*tri_iter), 1, 3, false, up_adj);
      if(up_adj.size() > 1) {
        error = MB->delete_entities(&(*tri_iter), 1);
        if (MB_SUCCESS != error)
          return error;
      }
    }
  }

    // create a meshset of the skin
  EntityHandle export_meshset;
  MB->create_meshset( MESHSET_SET, export_meshset);
  Tag material_tag;
  MB->tag_get_handle(MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER, material_tag);
  int block_id = 100;
  MB->tag_set_data(material_tag, &export_meshset, 1, &block_id);
  entities.clear();
  MB->get_entities_by_type(0, MBTRI, entities);
    // remove the first few tri's for fun
  Range tmp_ents;
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
  if (MB_SUCCESS != check_valid_connectivity( MB ))
    return MB_FAILURE;

    // make sure our first few tri's were untouched
  std::vector<EntityHandle> conn(3);
  for(Range::iterator kter=tmp_ents.begin(); kter != tmp_ents.end(); ++kter)
  {
    MB->get_connectivity(&(*kter), 1, conn);
    if(conn.size() != 3)
      return MB_FAILURE;
  }

    // output the skin
  file_name = "tri_mid_edge_face_nodes.g";
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

ErrorCode mb_forced_adjacencies_test()
{
    //! first clean up any existing mesh.
  ErrorCode error;
  Core moab;
  Interface* MB = &moab;


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

  EntityHandle node1, node2, node3, node4, node5, node6;
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

  std::vector<EntityHandle> conn(4);
    //! create the first quad
  EntityHandle              quad1;
  conn[0] = node1;
  conn[1] = node2;
  conn[2] = node5;
  conn[3] = node4;
  error = MB->create_element(MBQUAD, &conn[0], 4, quad1);
  if (error != MB_SUCCESS)
    return error;

    //! create the second quad
  EntityHandle              quad2;
  conn[0] = node2;
  conn[1] = node3;
  conn[2] = node6;
  conn[3] = node5;
  error = MB->create_element(MBQUAD, &conn[0], 4, quad2);
  if (error != MB_SUCCESS)
    return error;

    //! create the edges
  EntityHandle edge1;
  conn.resize(2);
  conn[0] = node1;
  conn[1] = node2;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge1);
  if (error != MB_SUCCESS)
    return error;

  EntityHandle edge2;
  conn[0] = node2;
  conn[1] = node3;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge2);
  if (error != MB_SUCCESS)
    return error;

  EntityHandle edge3;
  conn[0] = node1;
  conn[1] = node4;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge3);
  if (error != MB_SUCCESS)
    return error;

  EntityHandle edge4;
  conn[0] = node2;
  conn[1] = node5;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge4);
  if (error != MB_SUCCESS)
    return error;

  EntityHandle edge5;
  conn[0] = node2;
  conn[1] = node5;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge5);
  if (error != MB_SUCCESS)
    return error;

  EntityHandle edge6;
  conn[0] = node3;
  conn[1] = node6;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge6);
  if (error != MB_SUCCESS)
    return error;

  EntityHandle edge7;
  conn[0] = node4;
  conn[1] = node5;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge7);
  if (error != MB_SUCCESS)
    return error;

  EntityHandle edge8;
  conn[0] = node5;
  conn[1] = node6;
  error = MB->create_element(MBEDGE, &conn[0], 2, edge8);
  if (error != MB_SUCCESS)
    return error;


    //! Edge 4 and 5 share the same nodes, but should be different entities
  if (edge4 == edge5)
    return MB_FAILURE;

    //! Now that the geometry is created start adding the adjacency information
  std::vector<EntityHandle> edge_adjacencies1(4);
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

  std::vector<EntityHandle> edge_adjacencies2(4);
  edge_adjacencies2[0] = edge2;
  edge_adjacencies2[1] = edge6;
  edge_adjacencies2[2] = edge8;
  edge_adjacencies2[3] = edge5;
  error = MB->add_adjacencies(quad2, &edge_adjacencies2[0], edge_adjacencies2.size(), true);
  if (error != MB_SUCCESS)
    return error;

    //! now get the adjacencies of each quad.
  std::vector<EntityHandle> quad1_adjacencies;
  error = MB->get_adjacencies(&(quad1), 1, 1, false, quad1_adjacencies);
  if (error != MB_SUCCESS)
    return error;


  std::sort(quad1_adjacencies.begin(), quad1_adjacencies.end());
  std::sort(edge_adjacencies1.begin(), edge_adjacencies1.end());

  if (quad1_adjacencies != edge_adjacencies1)
    return MB_FAILURE;
  
  std::vector<EntityHandle> quad2_adjacencies;
  error = MB->get_adjacencies(&(quad2), 1, 1, false, quad2_adjacencies);
  if (error != MB_SUCCESS)
    return error;

  std::sort(quad2_adjacencies.begin(), quad2_adjacencies.end());
  std::sort(edge_adjacencies2.begin(), edge_adjacencies2.end());

  if (quad2_adjacencies != edge_adjacencies2)
    return MB_FAILURE;

    //! try getting the adjacency of edge1 (should be quad1)
  std::vector<EntityHandle> edge1_adjacencies;
  error = MB->get_adjacencies(&(edge1), 1, 2, false, edge1_adjacencies);
  if (error != MB_SUCCESS)
    return error;

    //! there should be only 1 entity adjacent to edge1
  if (edge1_adjacencies.size() != 1)
    return MB_FAILURE;

    //! and that entity should be quad1
  if (edge1_adjacencies[0] != quad1)
    return MB_FAILURE;

    //! try getting the adjacency of edge6 (should be one)
  std::vector<EntityHandle> edge6_adjacencies;
  error = MB->get_adjacencies(&(edge6), 1, 2, false, edge6_adjacencies);
  if (error != MB_SUCCESS)
    return error;

    //! there should be only 1 entity adjacent to edge6
  if (edge6_adjacencies.size() != 1)
    return MB_FAILURE;

    //! Now seal up the "gap" caused by edges 4 and 5.  Remove edge5
    //! from the adjacencies of quad2 and add edge 4 to quad2.

  std::vector<EntityHandle> edge5_adjacencies(1, edge5);
  error = MB->remove_adjacencies(quad2, &edge5_adjacencies[0], edge5_adjacencies.size());
  if (error != MB_SUCCESS)
    return error;


  std::vector<EntityHandle> edge4_adjacencies(1, edge4);
  error = MB->add_adjacencies(quad2, &edge4_adjacencies[0], edge4_adjacencies.size(), true);
  if (error != MB_SUCCESS)
    return error;
  
    //! get the adjacencies of edge4 and it should return both quads.
  std::vector<EntityHandle> quad_adjacencies;
  error = MB->get_adjacencies(&(edge4), 1, 2, false, quad_adjacencies);
  if (error != MB_SUCCESS)
    return error;

    //! there should be 2 entities adjacent to edge4
  if (quad_adjacencies.size() != 2)
    return MB_FAILURE;

    //! and they should be quad1 and quad2.  Note that we are not saying anything
    //! about order in the array.
  if ( (quad_adjacencies[0] != quad1 || quad_adjacencies[1] != quad2) &&
       (quad_adjacencies[0] != quad2 || quad_adjacencies[1] != quad1) ) 
    return MB_FAILURE;
  
    //! clean up on exit
  error = MB->delete_mesh();
  if (error != MB_SUCCESS)
    return error;


  return MB_SUCCESS;
}

/*bool lessnodesZ(const EntityHandle entity_handle1, const EntityHandle entity_handle2) 
  {
  double coords1[3], coords2[3];
  gMB->get_coords(entity_handle1, coords1);
  gMB->get_coords(entity_handle2, coords2);

  return coords2[2] < coords1[2];
  }*/

/*void sort_verts(Range vertices) 
  {
  std::vector<EntityHandle> vert_vec(vertices.size());
  Range::const_iterator iter;
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
ErrorCode find_coincident_nodes(Interface* gMB, Range vertices,
                                  std::vector< std::pair<EntityHandle,EntityHandle> > &coin_nodes)
{
  double first_coords[3], second_coords[3];
  Range::const_iterator iter, jter;
  std::pair<EntityHandle, EntityHandle> coincident_pair;
  ErrorCode result;

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

ErrorCode find_coincident_elements(Interface* gMB, Range entities, int num_nodes,
                                     std::vector< std::pair<EntityHandle,EntityHandle> > &coin)
{
  double coords1[8][3], coords2[8][3];
  Range::iterator iter, jter;
  std::vector<EntityHandle> conn(8);
  std::pair<EntityHandle, EntityHandle> coincident_pair;
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


ErrorCode mb_merge_test()
{ 
  Core moab;
  Interface* MB = &moab;
  
  time_t begin_time = clock();
  unsigned int i;
  ErrorCode result;
  Skinner Skinner_Obj(MB);

  std::string test_files[] = {std::string("cell1.gen"),
                              std::string("cell2.gen")};
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

    /*std::vector<Range> entities(sizeof(test_files));
      std::vector<Range> forward_lower(sizeof(test_files));
      std::vector<Range> reverse_lower(sizeof(test_files));
      std::vector<Range> nodes(sizeof(test_files));*/
  Range entities;
  Range forward_lower;
  Range reverse_lower;
  Range faces;
  Range nodes;

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
  Skinner_Obj.find_skin(entities,false,forward_lower,&reverse_lower);
  cout <<"num hexes = "<<entities.size()<<"\n";
  cout <<"fl = "<<forward_lower.size()<<" rl = "<<reverse_lower.size()<<"\n";
  
    //  Range::const_iterator iter;
  int dim = 0;
    //  int num_ents = 1;
  result = MB->get_adjacencies(forward_lower, dim, true, nodes, Interface::UNION);
  cout <<"nodes.size() = "<<nodes.size() <<"\n";
    
  if(result == MB_SUCCESS)
    cout << "---Success---";
  else
    cout << "---Failure---";
  cout << endl << endl;

    //  result = MB->get_entities_by_type(0, MBQUAD, faces);
    //  cout <<"num faces = "<<faces.size() <<"\n";

  std::vector<std::pair<EntityHandle, EntityHandle> > coin_nodes;
    //  cout <<"Begining sort...\n";
    //  std::sort(nodes.begin(),nodes.end(),lessnodesZ);
    //  cout <<"Ending sort...\n";
  result = find_coincident_nodes(MB,nodes, coin_nodes);
  cout <<"coin_nodes.size() = "<<coin_nodes.size() <<"\n";
  std::vector< std::pair<EntityHandle, EntityHandle> >::iterator n_iter;
  for (n_iter=coin_nodes.begin(); n_iter != coin_nodes.end(); n_iter++) {
    result = MB->merge_entities((*n_iter).first, (*n_iter).second, false, true);
    if (MB_SUCCESS != result)
      return result;
  }
    /*  std::vector<std::pair<EntityHandle, EntityHandle> > coin_faces;
        int nodes_per_elt = 4;
        result = find_coincident_elements(forward_lower, nodes_per_elt, coin_faces);
        if (result != MB_SUCCESS) cout <<"find_coincident_elements fail!\n";
        cout <<"coin_faces.size() = "<<coin_faces.size() <<"\n";
        std::vector< std::pair<EntityHandle, EntityHandle> >::iterator f_iter;
        for (f_iter=coin_faces.begin(); f_iter != coin_faces.end(); f_iter++)
        MB->merge_entities((*f_iter).first, (*f_iter).second, true, true);*/
    /*
      std::vector<std::pair<EntityHandle, EntityHandle> > coin_fl;
      nodes_per_elt = 4;
      result = find_coincident_elements(entities, nodes_per_elt, coin_fl);
      cout <<"coin_fl.size() = "<<coin_fl.size() <<"\n";
    */
  int num_ents;
  if ((MB_SUCCESS == MB->get_number_entities_by_dimension(0, 3, num_ents) &&
       0 != num_ents) ||
      (MB_SUCCESS == MB->get_number_entities_by_dimension(0, 2, num_ents) &&
       0 != num_ents))
    result = MB->write_mesh("merge_test.ncdf");
  ;
  

  double clocks_per_sec = (double) CLOCKS_PER_SEC;
  double real_time =  difftime(time(NULL), begin_time);
  cout <<"TIME: "<<(real_time/clocks_per_sec)<<" seconds.\n";
  return result;
}

ErrorCode mb_merge_update_test()
{
  Core moab;
  Interface* mb = &moab;
  ErrorCode rval;
  
    // create two quads with a coincident edge pair
  double coords[] = { 0, 0, 0,
                      1, 0, 0,
                      1, 1, 0,
                      0, 1, 0,
                      1, 1, 0,
                      1, 0, 0,
                      2, 0, 0,
                      2, 1, 0 };
  EntityHandle verts[8];
  for (int i = 0; i < 8; ++i) 
    mb->create_vertex( coords + 3*i, verts[i] );
  EntityHandle quad1, quad2, edge1, edge2;
  mb->create_element( MBQUAD, verts, 4, quad1 );
  mb->create_element( MBQUAD, verts+4, 4, quad2 );
  mb->create_element( MBEDGE, verts+1, 2, edge1 );
  mb->create_element( MBEDGE, verts+4, 2, edge2 );
  
    // create two tracking sets containing the vertices
    // and edge of each quad
  EntityHandle set1, set2;
  mb->create_meshset( MESHSET_TRACK_OWNER|MESHSET_SET, set1 );
  mb->create_meshset( MESHSET_TRACK_OWNER|MESHSET_ORDERED, set2 );
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
  Range r;
  mb->get_entities_by_type( 0, MBEDGE, r );
  if (r.size() != 1 || r.front() != edge1) {
    std::cerr << "Edge merge failed at " << __FILE__ << ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  std::vector<EntityHandle> exp(verts+1, verts+3), act;
  mb->get_connectivity( &edge1, 1, act );
  if (exp != act) {
    std::cerr << "Incorrect conn for edge at " << __FILE__ << ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  
    // check that quad connectivity is as expected
  exp = std::vector<EntityHandle>(verts, verts+4);
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
  exp = std::vector<EntityHandle>(verts, verts+4);
  exp.push_back( edge1 );
  act.clear();
  mb->get_entities_by_handle( set1, act );
  std::sort( exp.begin(), exp.end() );
  std::sort( act.begin(), act.end() );
  if (exp != act) {
    std::cerr << "Incorrect set contents at " << __FILE__ << ":" << __LINE__ << std::endl;
    std::cerr << "  Expected: ";
    std::copy( exp.begin(), exp.end(), std::ostream_iterator<EntityHandle>(std::cerr, " ") );
    std::cerr << std::endl << "  Actual  : ";
    std::copy( act.begin(), act.end(), std::ostream_iterator<EntityHandle>(std::cerr, " ") );
    std::cerr << std::endl;
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
  if (exp != act) {
    std::cerr << "Incorrect set contents at " << __FILE__ << ":" << __LINE__ << std::endl;
    std::cerr << "  Expected: ";
    std::copy( exp.begin(), exp.end(), std::ostream_iterator<EntityHandle>(std::cerr, " ") );
    std::cerr << std::endl << "  Actual  : ";
    std::copy( act.begin(), act.end(), std::ostream_iterator<EntityHandle>(std::cerr, " ") );
    std::cerr << std::endl;
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}

ErrorCode mb_stress_test()
{
  ErrorCode error;
  Core moab;
  Interface* MB = &moab;

  cout << "    Beginning Stress Test . . ." << endl;
  cout << "\n        Reading elements" << endl;
  clock_t start = clock();
  clock_t total_start = clock();

    //read in a file so you have some data in the database
  std::string file_name = "mb_big_test.g";
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

  Range hexes;
  error = MB->get_entities_by_type(0, MBHEX, hexes);
  if (error != MB_SUCCESS)
    return error;


  std::vector<EntityHandle> conn;
  Range::iterator iter;
  for (iter = hexes.begin(); iter != hexes.end(); iter++)
  {
    error = MB->get_connectivity(&(*iter), 1, conn);
    if (error != MB_SUCCESS)
      return error;

    double coords[3];
    int i = 0;
    std::vector<EntityHandle> vertex_handle(8);
    EntityHandle element_handle;
    std::vector<EntityHandle>::iterator jter;

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

  EntityHandle mesh_set;
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
  Tag tag_handle;
  ErrorCode result = MB->tag_get_handle( MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER, tag_handle ) ;
  if(result != MB_SUCCESS)
    return result;

  int id = 1;
  result = MB->tag_set_data( tag_handle, &mesh_set, 1, &id ) ;
  if(result != MB_SUCCESS)
    return result;

  std::vector<EntityHandle> output_list;
  output_list.push_back(mesh_set);

  file_name = "mb_stress_out.g";
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

ErrorCode mb_canon_number_test() 
{
  Core moab;
  Interface* MB = &moab;

    // various tests for canonical ordering

    // CN::AdjacentSubEntities
  std::vector<int> vec1, vec2;
  ErrorCode result;

  EntityType this_type;

  for (this_type = MBEDGE; this_type != MBKNIFE; this_type++) {
    
    for (int i = 0; i < CN::VerticesPerEntity(this_type); i++) {
        // test for edges and faces
      for (int dim = 1; dim <= CN::Dimension(this_type); dim++) {
          // get the sides adjacent to this vertex
        vec1.clear();
        int temp_result = CN::AdjacentSubEntities(this_type, &i, 1, 0, dim, vec1);
                                     
        if (0 != temp_result ||
            vec1.size() > (unsigned int) CN::NumSubEntities(this_type, dim)) {
          cout << "failed getting sides for type " << CN::EntityTypeName(this_type)
               << " dimension" << dim << endl;
          return MB_FAILURE;
        }
        

          // now get the vertices shared by these sides
        vec2.clear();
        temp_result = 
          CN::AdjacentSubEntities(this_type, &vec1[0], vec1.size(), dim, 0,
                                    vec2);

          // vertex side recovered should be i
        if (0 != temp_result || 
              // if dimension is same as DIMENSION(this_type), will get back all the
              // vertices in the entity
            (dim == CN::Dimension(this_type) && 
             vec2.size() != (unsigned int) CN::VerticesPerEntity(this_type)) ||
              // otherwise, we should get back only one vertex, and it should be the one
              // we started with
            (dim != CN::Dimension(this_type) &&
             (vec2.size() != 1 || vec2[0] != i))) {
          cout << "traversal from verts to sides to verts failed for " << endl
               << "vertex " << i << " type " << CN::EntityTypeName(this_type) 
               << " dimension " << dim << endl;
          return MB_FAILURE;
        }
      }
    }
  }
  
    // CN::side_number

    // create vertices to use later
  double xyz[3] = {0.0, 0.0, 0.0};
  EntityHandle vertex_handles[8];
  for (int i = 0; i < 8; i++) {
    result = MB->create_vertex(xyz, vertex_handles[i]);
    assert(result == MB_SUCCESS);
  }
  int side, sense, offset;
  
  EntityHandle this_entity;
  
  for (this_type = MBEDGE; this_type != MBKNIFE; this_type++) {
    
      // skip remainder of the loop for MBPOLYGONS and POLYHEDRA, which don't follow
      // the standard canonical ordering 
    if (this_type == MBPOLYGON || this_type == MBPOLYHEDRON)
      continue;
    
      // make an entity of this type
    result = MB->create_element(this_type, vertex_handles, 
                                CN::VerticesPerEntity(this_type),
                                this_entity);
    if (MB_SUCCESS != result || 0 == this_entity) {
      cout << "failed to create entity of type " 
           << CN::EntityTypeName(this_type) << endl;
      return MB_FAILURE;
    }

      // now get the connectivity vector *
    const EntityHandle *entity_vertices;
    int num_verts;
    result = MB->get_connectivity(this_entity, entity_vertices, num_verts);
    if (MB_SUCCESS != result || 
        num_verts != CN::VerticesPerEntity(this_type)) {
      cout << "failed to get connectivity for entity type " 
           << CN::EntityTypeName(this_type) << endl;
      return MB_FAILURE;
    }
    
      // for each dimension
    for (int dim = 1; dim <= CN::Dimension(this_type); dim++) {
        // for each side of this dimension
      const CN::ConnMap &cm = CN::mConnectivityMap[this_type][dim-1];
      int tmp_conn[moab::MAX_SUB_ENTITY_VERTICES];

      for (int side_no = 0; side_no < CN::NumSubEntities(this_type, dim); side_no++) {

        for (int j = 0; j < moab::MAX_SUB_ENTITY_VERTICES; j++) tmp_conn[j] = cm.conn[side_no][j];
        int temp_result = 
          CN::SideNumber(this_type,
                           tmp_conn, 
                           CN::VerticesPerEntity(CN::SubEntityType(this_type, dim, side_no)),
                           dim, side, sense, offset);
        if (0 != temp_result) {
          cout << "call to CN::side_number failed with non-success result"
               << " for type " 
               << CN::EntityTypeName(this_type) << " dimension " << dim 
               << " side no " << side_no << endl;
          return MB_FAILURE;
        }

          // side number should be the same as side_no, sense should be forward, offset should
          // be zero
        if (side != side_no || sense != 1 || offset != 0) {
          cout << "call to CN::side_number failed for type " 
               << CN::EntityTypeName(this_type) << " dimension " << dim << " side no " 
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

ErrorCode mb_poly_test() 
{
  Core moab;
  Interface* mb = &moab;

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

  EntityHandle verts[16];
  ErrorCode result;
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
  
  EntityHandle polygons[9], temp_connect[12];
  int idx = 0, nump = 0;
  ErrorCode tmp_result;
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
  const EntityHandle *connect;
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
  Range vert_range, poly_range;
  for (i = 0; i < 9; i++) poly_range.insert(polygons[i]);
  result = mb->get_adjacencies(poly_range, 0, false, vert_range, 
                               Interface::UNION);
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
  EntityHandle polyhedra[2];
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
  std::vector<EntityHandle> temp_verts;
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

ErrorCode mb_topo_util_test() 
{
  Core moab;
  Interface* gMB = &moab;
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
      
  EntityHandle grid_verts[18], grid_elems[4];
  ErrorCode result;
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
  EntityHandle connect[8];
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

  Range vert_range;
  std::copy(grid_verts, grid_verts+18, range_inserter(vert_range));
  
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
  Range all_hexes, middle_edge;
  std::copy(grid_elems, grid_elems+4, range_inserter(all_hexes));
    // get the shared edge
  result = gMB->get_adjacencies(all_hexes, 1, false, middle_edge);
  if (MB_SUCCESS != result || 1 != middle_edge.size()) {
    std::cout << "Bad result getting single shared edge." << std::endl;
    return MB_FAILURE;
  }
  
  std::vector<EntityHandle> star_faces, star_hexes;
  bool bdy_edge;
  result = mtu.star_entities(*middle_edge.begin(), star_faces, bdy_edge, 0, &star_hexes);
  if (MB_SUCCESS != result || bdy_edge || star_faces.size() != 4 || star_hexes.size() != 4) {
    std::cout << "Bad result from star_faces for non-bdy edge." << std::endl;
    return MB_FAILURE;
  }
  
    // now try for a different edge, which has to be on the bdy
  Range other_edges;
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

ErrorCode mb_split_test() 
{
  Core moab;
  Interface* gMB = &moab;
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
      
  EntityHandle grid_verts[27], grid_elems[8];
  ErrorCode result;
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
  EntityHandle connect[8];
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
  
  Range vert_range;
  std::copy(grid_verts, grid_verts+27, range_inserter(vert_range));
  
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
  Range split_faces, tmp_ents, tmp_faces;
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

  Range new_faces, new_regions;

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

ErrorCode mb_range_seq_intersect_test() 
{
  ErrorCode rval;
  SequenceManager sequences;
  RangeSeqIntersectIter iter( &sequences );
  Range range;

    // create some entity sequences
  EntitySequence *ts1, *ts2, *ts3, *qs1;
  EntityHandle th1, th2, th3, qh1;
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
  
    // construct an Range containing all valid handles;
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
    
  Range quads = range.subset_by_type( MBQUAD );
  EntityHandle removed = qs1->start_handle() + nq1/2;
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

  Range big;
  int junk;
  EntityHandle last = CREATE_HANDLE(MBQUAD+1, 0, junk);
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
  Error eh;
  EntityHandle ts1s  = ts1->start_handle();
  EntityHandle dead1 = ts1->start_handle() + 1;
  EntityHandle dead2 = ts1->start_handle() + 2;
  EntityHandle ts1e  = ts1->end_handle();
  EntityHandle dead3 = ts2->start_handle();
  EntityHandle dead4 = ts2->end_handle();
  EntityHandle qs1s  = qs1->start_handle();
  EntityHandle qs1e  = qs1->end_handle();
  EntityHandle dead5 = qs1->start_handle() + nq1/2;
  EntityHandle dead6 = dead5+1;
  rval = sequences.delete_entity( &eh, dead1 );
  if (MB_SUCCESS != rval) return rval;
  rval = sequences.delete_entity( &eh, dead2 );
  if (MB_SUCCESS != rval) return rval;
  rval = sequences.delete_entity( &eh, dead3 );
  if (MB_SUCCESS != rval) return rval;
  rval = sequences.delete_entity( &eh, dead4 );
  if (MB_SUCCESS != rval) return rval;
  rval = sequences.delete_entity( &eh, dead5 );
  if (MB_SUCCESS != rval) return rval;
  rval = sequences.delete_entity( &eh, dead6 );
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
    
std::ostream& operator<<( std::ostream& s, Range::const_iterator i ) {
  return s << *i;
}

ErrorCode mb_poly_adjacency_test()
{
  ErrorCode rval;
  Core moab;
  Interface *mbImpl = &moab;
  
    // make a couple polygons and a polyhedron
  double coords[3] = {0,1,2};
  EntityHandle verts[10], polygons[2], polyhedron;
  
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
  Range dum_range;
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

ErrorCode mb_memory_use_test() 
{
  Core mb;
  unsigned long init_total, total_with_elem, total_with_tag, total_with_tag_data;
  mb.estimated_memory_use(0,0,0,&init_total);
  
  double coords[12] = { 1, 2, 0, 3, 4, 0, 5, 6, 0, 7, 8, 0 };
  EntityHandle verts[4];
  for (int i = 0; i < 4; ++i)
    if (MB_SUCCESS != mb.create_vertex( coords + 3*i, verts[i] ))
      return MB_FAILURE;
  
  EntityHandle elem;
  mb.create_element( MBQUAD, verts, 4, elem );
  
  mb.estimated_memory_use(0,0,0,&total_with_elem);
  if (total_with_elem <= init_total)
    return MB_FAILURE;
  
  unsigned long min, am;
  Range r;
  r.insert( elem );
  mb.estimated_memory_use( r, &min, &am );
  if (min != 4*sizeof(EntityHandle))
    return MB_FAILURE;
  
  r.clear();
  r.insert( verts[0] );
  r.insert( verts[1] );
  mb.estimated_memory_use( r, &min, &am );
  if (min != 6*sizeof(double))
    return MB_FAILURE;
  
  Tag tag;
  if (MB_SUCCESS != mb.tag_get_handle( "TMP_TAG", 1, MB_TYPE_INTEGER, tag, MB_TAG_SPARSE|MB_TAG_EXCL ))
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

ErrorCode mb_skin_curve_test_common( bool use_adj );

ErrorCode mb_skin_curve_test()
  { return mb_skin_curve_test_common( false ); }

ErrorCode mb_skin_curve_adj_test()
  { return mb_skin_curve_test_common( true ); }

ErrorCode mb_skin_surface_test_common( bool use_adj );

ErrorCode mb_skin_surface_test()
  { return mb_skin_surface_test_common( false ); }

ErrorCode mb_skin_surface_adj_test()
  { return mb_skin_surface_test_common( true ); }

ErrorCode mb_skin_volume_test_common( bool use_adj );

ErrorCode mb_skin_volume_test()
  { return mb_skin_volume_test_common( false ); }

ErrorCode mb_skin_volume_adj_test()
  { return mb_skin_volume_test_common( true ); }

ErrorCode mb_skin_curve_test_common( bool use_adj )
{
  ErrorCode rval;
  Core moab;
  Interface *mb = &moab;
  
  std::vector<EntityHandle> verts;
  for (unsigned i = 0; i < 10; ++i) {
    double coords[] = { i, 0, 0 };
    EntityHandle h;
    mb->create_vertex( coords, h );
    verts.push_back(h);
  }
  Range edges;
  for (unsigned i = 1; i < verts.size(); ++i) {
    EntityHandle conn[] = { verts[i-1], verts[i] };
    EntityHandle h;
    mb->create_element( MBEDGE, conn, 2, h );
    edges.insert(h);
  }
  
  Range skin;
  Skinner tool(mb);
  rval = tool.find_skin( edges, 0, skin, use_adj );
  if (MB_SUCCESS != rval) {
    std::cerr << "Skinner failure at " __FILE__ ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  if (skin.size() != 2) {
    std::cerr << "Skinner bad result at " __FILE__ ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }

  if (verts.front() > verts.back())
    std::swap( verts.front(), verts.back() );
  if (skin.front() != verts.front() ||
      skin.back()  != verts.back()) {
    std::cerr << "Skinner bad result at " __FILE__ ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  
  
    // now test again with only one edge
  EntityHandle edge = edges.front();
  Range range(edge,edge);
  skin.clear();
  rval = tool.find_skin( range, 0, skin, use_adj );
  if (MB_SUCCESS != rval) {
    std::cerr << "Skinner failure at " __FILE__ ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  if (skin.size() != 2) {
    std::cerr << "Skinner bad result at " __FILE__ ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }

  Range verts2;
  mb->get_connectivity( &edge, 1, verts2 );
  if (skin.front() != verts2.front() ||
      skin.back()  != verts2.back()) {
    std::cerr << "Skinner bad result at " __FILE__ ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}

ErrorCode mb_skin_surface_test_common( bool use_adj )
{
  ErrorCode rval;
  Core moab;
  Interface *mb = &moab;
  
    /* Create 4 of 5 faces of a wedge: */
    /*
            4
           /|\
          / | \
         /  |  \
        /   2   \
       3.../.\...5
       |  /   \  |
       | /     \ |
       |/       \|
       0_________1
     */
     
  const double coords[][3] = { { 0, 0, 0 },
                               { 2, 0, 0 },
                               { 1, 0, 1 },
                               { 0, 2, 0 },
                               { 2, 2, 0 },
                               { 1, 2, 1 } };
  EntityHandle verts[6];
  for (unsigned i = 0; i < 6; ++i)
    mb->create_vertex( coords[i], verts[i] );

  EntityHandle faces[4];
  EntityHandle tri[] = { verts[0], verts[1], verts[2] };
  EntityHandle quad1[] = { verts[0], verts[1], verts[5], verts[3] };
  EntityHandle quad2[] = { verts[1], verts[5], verts[4], verts[2] };
  EntityHandle quad3[] = { verts[2], verts[4], verts[3], verts[0] };
  mb->create_element( MBTRI, tri, 3, faces[0] );
  mb->create_element( MBQUAD, quad1, 4, faces[1] );
  mb->create_element( MBQUAD, quad2, 4, faces[2] );
  mb->create_element( MBQUAD, quad3, 4, faces[3] );
  Range source;
  std::copy( faces, faces+4, range_inserter(source) );
  
    // Now skin the mesh.  The only missing face is the 
    // back triangle (verts 3, 4, & 5) so the skin should
    // be the edges bordering that face.
  
  Range skin;
  Skinner tool(mb);
  rval = tool.find_skin( source, 1, skin, use_adj );
  if (MB_SUCCESS != rval) {
    std::cerr << "Skinner failure at " __FILE__ ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  if (skin.size() != 3 || !skin.all_of_type(MBEDGE)) {
    std::cerr << "Skinner bad result at " __FILE__ ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }

    // Check each edge
  std::vector<EntityHandle> conn[3];
  rval = mb->get_connectivity( &skin.front(), 1, conn[0] );
  rval = mb->get_connectivity( &*++skin.begin(), 1, conn[1] );
  rval = mb->get_connectivity( &skin.back(), 1, conn[2] );
  for (int i = 0; i < 3; ++i)
    if (conn[i][0] > conn[i][1])
      std::swap(conn[i][0], conn[i][1]);

  for (int i = 0; i < 3; ++i) {
    EntityHandle s = verts[i+3], e = verts[(i+1)%3 + 3];
    if (s > e)
      std::swap(s,e);
    int j = 0; 
    for (j = 0; j < 3; ++j) 
      if (conn[j][0] == s && conn[j][1] == e)
        break;
    
    if (j == 3) {
      std::cerr << "Skin does not contain edge [" << s << "," << e << "] at " 
                << __FILE__ ":" << __LINE__ << std::endl;
      return MB_FAILURE;
    }
  }

  return MB_SUCCESS;
}

ErrorCode mb_skin_volume_test_common( bool use_adj )
{
  ErrorCode rval;
  Core moab;
  Interface *mb = &moab;
  
    /* A 2 adjacent hexes hexes */
    /*
          9-----10----11
         /     /     /|
        /     /     / |
       6-----7-----8..5
       | .   | .   | /
       |.    |.    |/
       0-----1-----2
     */
     
  const double coords[][3] = { { 0, 0, 0 },
                               { 1, 0, 0 },
                               { 2, 0, 0 },
                               { 0, 1, 0 },
                               { 1, 1, 0 },
                               { 2, 1, 0 },
                               { 0, 0, 1 },
                               { 1, 0, 1 },
                               { 2, 0, 1 },
                               { 0, 1, 1 },
                               { 1, 1, 1 },
                               { 2, 1, 1 } };
  EntityHandle verts[12];
  for (unsigned i = 0; i < 12; ++i)
    mb->create_vertex( coords[i], verts[i] );

  EntityHandle hex1c[] = { verts[0], verts[1], verts[4], verts[3],
                          verts[6], verts[7], verts[10], verts[9] };
  EntityHandle hex2c[] = { verts[1], verts[2], verts[5], verts[4],
                          verts[7], verts[8], verts[11], verts[10] };
  EntityHandle hex1, hex2;
  mb->create_element( MBHEX, hex1c, 8, hex1 );
  mb->create_element( MBHEX, hex2c, 8, hex2 );
  Range source;
  source.insert( hex1 );
  source.insert( hex2 );
    
    // get all quads and shared face
  Range tmp, all_faces;
  mb->get_adjacencies( source, 2, true, all_faces, Interface::UNION );
  mb->get_adjacencies( source, 2, true, tmp, Interface::INTERSECT );
  assert(tmp.size() == 1);
  Range non_shared = subtract( all_faces, tmp );
  
    // Now skin the mesh.  
  
  Range skin;
  Skinner tool(mb);
  rval = tool.find_skin( source, 2, skin, use_adj );
  if (MB_SUCCESS != rval) {
    std::cerr << "Skinner failure at " __FILE__ ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }
  if (skin != non_shared) {
    std::cerr << "Skinner bad result at " __FILE__ ":" << __LINE__ << std::endl;
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}

ErrorCode mb_skin_scd_test() 
{
    // make a 10x10x5 scd mesh
  Core moab;
  Interface *mb = &moab;
  ScdInterface *scdi;
  ErrorCode rval = mb->query_interface(scdi);
  if (MB_SUCCESS != rval) return rval;
  HomCoord low(0, 0, 0), high(10, 10, 5);
  ScdBox *this_box;
  rval = scdi->construct_box(low, high, NULL, 0, this_box);
  if (MB_SUCCESS != rval) return rval;
  
    // now skin it with the structured and original method, and compare results
  Skinner tool(mb);
  Range ents(this_box->start_element(), this_box->start_element()+this_box->num_elements()-1),
      scd_skin_ents, skin_ents;
  rval = tool.find_skin(ents, false, scd_skin_ents, NULL, true, true, true);
  if (MB_SUCCESS != rval) return rval;

  rval = tool.find_skin(ents, false, skin_ents, NULL, true, true, false);
  if (MB_SUCCESS != rval) return rval;

    // should be same number of entities
  if (scd_skin_ents.size() != skin_ents.size()) return MB_FAILURE;
  
  skin_ents.clear();
  scd_skin_ents.clear();
  
    // now test getting faces and vertices, also with existing faces now
  rval = tool.find_skin(ents, true, scd_skin_ents, NULL, true, true, true);
  if (MB_SUCCESS != rval) return rval;

  rval = tool.find_skin(ents, true, skin_ents, NULL, true, true, false);
  if (MB_SUCCESS != rval) return rval;

    // again, should have same numbers
  if (skin_ents.subset_by_type(MBVERTEX).size() != scd_skin_ents.subset_by_type(MBVERTEX).size()) 
    return MB_FAILURE;
  
  return MB_SUCCESS;
}

// It is a common problem for readers to incorrectly
// handle invalid/unknown file types and non-existant 
// files.  For either case, MOAB will try all readers
// (assuming it doesn't recongnize the file extension),
// so we can test all the readers w/out knowing which
// readers we have.
const char* argv0 = 0;
ErrorCode mb_read_fail_test()
{
  Core moab;
  Interface* mb = &moab;

  const char BAD_FILE_NAME[] = "non-existant-file.txt";
  ErrorCode rval;
  
  FILE* fptr = fopen(BAD_FILE_NAME,"r");
  if (fptr) {
    fclose(fptr);
    std::cout << "Test cannot proceed while file exists: " << BAD_FILE_NAME << std::endl;
    return MB_FAILURE;
  }
  
    // try reading a non-existant file
  
  rval = mb->load_file( BAD_FILE_NAME );
  if (MB_FILE_DOES_NOT_EXIST != rval)
    return MB_FAILURE;
  
    // try reading an invalid file
  if (!argv0)  // not set by main, oops!
    return MB_FAILURE;
  rval = mb->load_file( argv0 );
  if (MB_SUCCESS == rval)
    return MB_FAILURE;
  
  return MB_SUCCESS;
}

#define TEST_ERROR_CODE(E) \
  if (mb->get_error_string(E) != #E) { \
    std::cerr << "Invalid error string from get_error_string for " \
              << #E << ": " << mb->get_error_string(E) << std::endl;\
    return MB_FAILURE; \
  } 

ErrorCode mb_enum_string_test()
{
  Core moab;
  Interface* mb = &moab;

  TEST_ERROR_CODE( MB_SUCCESS );
  TEST_ERROR_CODE( MB_INDEX_OUT_OF_RANGE );
  TEST_ERROR_CODE( MB_TYPE_OUT_OF_RANGE );
  TEST_ERROR_CODE( MB_MEMORY_ALLOCATION_FAILED );
  TEST_ERROR_CODE( MB_ENTITY_NOT_FOUND );
  TEST_ERROR_CODE( MB_MULTIPLE_ENTITIES_FOUND );
  TEST_ERROR_CODE( MB_TAG_NOT_FOUND );
  TEST_ERROR_CODE( MB_FILE_DOES_NOT_EXIST );
  TEST_ERROR_CODE( MB_FILE_WRITE_ERROR );
  TEST_ERROR_CODE( MB_NOT_IMPLEMENTED );
  TEST_ERROR_CODE( MB_ALREADY_ALLOCATED );
  TEST_ERROR_CODE( MB_VARIABLE_DATA_LENGTH );
  TEST_ERROR_CODE( MB_INVALID_SIZE );
  TEST_ERROR_CODE( MB_UNSUPPORTED_OPERATION );
  TEST_ERROR_CODE( MB_UNHANDLED_OPTION );
  TEST_ERROR_CODE( MB_STRUCTURED_MESH );
  TEST_ERROR_CODE( MB_FAILURE );
  
  return MB_SUCCESS;
}

// Test basic skinning using vert-to-elem adjacencies
ErrorCode mb_skin_verts_common( unsigned dim, bool skin_elems );

ErrorCode mb_skin_surf_verts_test()
  { return mb_skin_verts_common( 2, false ); }

ErrorCode mb_skin_vol_verts_test()
  { return mb_skin_verts_common( 3, false ); }

ErrorCode mb_skin_surf_verts_elems_test()
  { return mb_skin_verts_common( 2, true ); }

ErrorCode mb_skin_vol_verts_elems_test()
  { return mb_skin_verts_common( 3, true ); }

ErrorCode mb_skin_verts_common( unsigned dim, bool skin_elems )
{
  const int INT = 10; // intervals+1
  const char* tmp_file = "structured.vtk";
  std::ofstream str(tmp_file);
  if (!str) {
    std::cerr << tmp_file << ": filed to create temp file" << std::endl;
    return MB_FAILURE;
  }
  str << "#vtk DataFile Version 2.0" << std::endl
      << "mb_skin_verts_common temp file" << std::endl
      << "ASCII" << std::endl
      << "DATASET STRUCTURED_POINTS" << std::endl
      << "DIMENSIONS " << INT << " " << (dim > 1 ? INT : 1) << " "  << (dim > 2 ? INT : 1) << std::endl
      << "ORIGIN 0 0 0" << std::endl
      << "SPACING 1 1 1" << std::endl;
  str.close();
  
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  
  rval = mb.load_file( tmp_file );
  remove( tmp_file );
  if (MB_SUCCESS != rval) 
    return rval;
  
  Range ents;
  rval = mb.get_entities_by_dimension( 0, dim, ents );
  if (MB_SUCCESS != rval)
    return rval;
  if (ents.empty())
    return MB_FAILURE;
  
  Skinner tool( &mb );
  
   // mesh is a structured quad/hex mesh, so we can
   // determine skin vertices from the number of
   // adjacent elements.
  unsigned interior_adj = 1;
  for (unsigned i = 0; i < dim; ++i)
    interior_adj *= 2;
  Range expected, verts;
  rval = mb.get_entities_by_dimension( 0, 0, verts );
  if (MB_SUCCESS != rval)
    return rval;
  Range::iterator h = expected.begin();
  std::vector<EntityHandle> adj;
  for (Range::iterator v = verts.begin(); v != verts.end(); ++v) {
    adj.clear();
    rval = mb.get_adjacencies( &*v, 1, dim, false, adj );
    if (MB_SUCCESS != rval)
      return rval;
    if (adj.size() < interior_adj)
      h = expected.insert( h, *v );
  }
  
    // Get skin vertices using skinner
  Range actual;
  rval = tool.find_skin( ents, !skin_elems, actual );
  if (MB_SUCCESS != rval)
    return rval;
 
  Range extra, missing;
  if (!skin_elems) {
      // Check that we got expected result
    extra = subtract( actual, expected );
    missing = subtract( expected, actual );
    if (!extra.empty() || !missing.empty()) {
      std::cout << "Extra vertices returned: " << extra << std::endl
                << "Missing vertices: " << missing << std::endl;
      return MB_FAILURE;
    }
    return MB_SUCCESS;
  }
  
    // test that no extra elements we're created
  extra.clear();
  rval = mb.get_entities_by_dimension( 0, dim-1, extra );
  if (MB_SUCCESS != rval)
    return rval;
  extra = subtract( extra, actual );
  if (!extra.empty()) {
    std::cout << "Extra/non-returned elements created: " << extra << std::endl;
    return MB_FAILURE;
  }
    
    // check that each skin vertex has the correct number of adjacent quads
  missing.clear(); extra.clear();
  for (Range::iterator i = expected.begin(); i != expected.end(); ++i) {
    std::vector<EntityHandle> elem, side;
    rval = mb.get_adjacencies( &*i, 1, dim, false, elem );
    if (MB_SUCCESS != rval) return rval;
    rval = mb.get_adjacencies( &*i, 1, dim-1, false, side );
    if (MB_SUCCESS != rval) return rval;
    if (elem.size() == 1) {
      if (side.size() < dim)
        missing.insert( *i );
      else if(side.size() > dim)
        extra.insert( *i );
    }
    else if (elem.size() == interior_adj) {
      if (!side.empty())
        extra.insert( *i );
    }
    else {
      if (side.size() < interior_adj/2)
        missing.insert( *i );
      else if (side.size() > interior_adj/2)
        extra.insert( *i );
    }
  }
  if (!missing.empty() || !extra.empty()) {
    std::cout << "Skin elements missing at vertices: " << missing << std::endl
              << "Extra skin elements at vertices: " << extra << std::endl;
    return MB_FAILURE;
  }
  
    // check that all returned elements are actually on the skin
  extra.clear();
  for (Range::iterator i = actual.begin(); i != actual.end(); ++i) {
    Range verts;
    rval = mb.get_adjacencies( &*i, 1, 0, false, verts );
    if (MB_SUCCESS != rval)
      return rval;
    verts = subtract( verts, expected );
    if (!verts.empty())
      extra.insert( *i );
  }
  if (!extra.empty()) {
    std::cout << "Skinner returned elements not on skin: " << extra << std::endl;
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;    
}

// Test that skinning of polyhedra works
ErrorCode mb_skin_poly_test()
{
  /* Create a mesh composed of 8 hexagonal prisms and 
     two hexahedra by extruding the following cross section
     two steps in the Z direction.
  
             0-----1
            /  (0)  \
       (11)/         \(1)
          /           \
        11             2
        / \     Y     / \
   (10)/   \(12)^(13)/   \(2)
      /     \   |   /     \
    10      12-----13      3
     |       |  |  |       |
  (9)|       |  +--|-->X   |(3)
     |       |     |       |
     9      15-----14      4
      \     /       \     /
    (8)\   /(15) (14)\   /(4)
        \ /           \ /
         8             5
          \           /
        (7)\         /(5)
            \  (6)  /
             7-----6
  */

  const double coords2D[][2] = { {-1, 5}, // 0
                                 { 1, 5},
                                 { 3, 3},
                                 { 5, 1},
                                 { 5,-1},
                                 { 3,-3}, // 5
                                 { 1,-5},
                                 {-1,-5},
                                 {-3,-3},
                                 {-5,-1},
                                 {-5, 1}, // 10
                                 {-3, 3},
                                 {-1, 1},
                                 { 1, 1},
                                 { 1,-1},
                                 {-1,-1}  // 15
                                 };
  const int polyconn[4][6] = { { 0,  1,  2, 13, 12, 11 },
                               { 2,  3,  4,  5, 14, 13 },
                               { 5,  6,  7,  8, 15, 14 },
                               { 8,  9, 10, 11, 12, 15 } };
  const int polyside[4][6] = { { 0,  1, 13, 16, 12, 11 },
                               { 2,  3,  4, 14, 17, 13 },
                               { 5,  6,  7, 15, 18, 14 },
                               { 8,  9, 10, 12, 19, 15 } };

  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  Range regions, faces, interior_faces;

    // create 48 vertices
  EntityHandle verts[3][16];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 16; ++j) {
      double coords[3] = { coords2D[j][0], coords2D[j][1], 2*i };
      rval = mb.create_vertex( coords, verts[i][j] );
      if (MB_SUCCESS != rval) return rval;
    }
  }
  
    // create two hexahedra
  EntityHandle hexes[2];
  for (int i = 0; i < 2; ++i) {
    EntityHandle conn[8] = { verts[i  ][15],
                               verts[i  ][14],
                               verts[i  ][13],
                               verts[i  ][12],
                               verts[i+1][15],
                               verts[i+1][14],
                               verts[i+1][13],
                               verts[i+1][12] };
    rval = mb.create_element( MBHEX, conn, 8, hexes[i] );
    if (MB_SUCCESS != rval) return rval;
    regions.insert(hexes[i]);
  }
  
    // create hexagonal faces
  EntityHandle hexagons[3][4];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 4; ++j) {
      EntityHandle conn[6];
      for (int k = 0; k < 6; ++k) 
        conn[k] = verts[i][polyconn[j][k]];
      rval = mb.create_element( MBPOLYGON, conn, 6, hexagons[i][j] );
      if (MB_SUCCESS != rval) return rval;
      faces.insert( hexagons[i][j] );
      if (i == 1)
        interior_faces.insert( hexagons[i][j] );
    }
  }
  
    // create quadrilateral faces
  EntityHandle quads[2][20];
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 20; ++j) {
      int c1, c2;
      if (j < 12) {
        c1 = j; c2 = (j+1)%12;
      }
      else if (j < 16) {
        c1 = j; c2 = 2 + 3*((j-9)%4);
      }
      else {
        c1 = j-4; c2 = 12 + (j-15)%4;
      }
      EntityHandle conn[4] = { verts[i  ][c1],
                                 verts[i  ][c2],
                                 verts[i+1][c2],
                                 verts[i+1][c1] };
      rval = mb.create_element( MBQUAD, conn, 4, quads[i][j] );
      if (MB_SUCCESS != rval) return rval;
      faces.insert( quads[i][j] );
      if (j > 11)
        interior_faces.insert( quads[i][j] );
    }
  }
  
    // create polyhedra
  EntityHandle poly[2][4];
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 4; ++j) {
      EntityHandle conn[8];
      for (int k = 0; k < 6; ++k) 
        conn[k] = quads[i][polyside[j][k]];
      conn[6] = hexagons[  i][j];
      conn[7] = hexagons[i+1][j];
      rval = mb.create_element( MBPOLYHEDRON, conn, 8, poly[i][j] );
      if (MB_SUCCESS != rval) return rval;
      regions.insert( poly[i][j] );
    }
  }
  
  Range interior_verts;
  interior_verts.insert( verts[1][12] );
  interior_verts.insert( verts[1][13] );
  interior_verts.insert( verts[1][14] );
  interior_verts.insert( verts[1][15] );
  
  Skinner tool(&mb);
  Range skin;
  rval = tool.find_skin( regions, true, skin, 0, true, false );
  if (MB_SUCCESS != rval) {
    std::cout << "Vertex skinning failed with: " << mb.get_error_string(rval) << std::endl;
    return rval;
  }
  
  Range all_verts, all_faces;
  rval = mb.get_entities_by_dimension( 0, 0, all_verts );
  if (MB_SUCCESS != rval) return rval;
  rval = mb.get_entities_by_dimension( 0, 2, all_faces );
  if (MB_SUCCESS != rval) return rval;
  
  Range expected = subtract( all_verts, interior_verts );
  if (expected != skin) {
    std::cout << "Skinner returned incorrect vertices." << std::endl;
    return MB_FAILURE;
  }
  if (all_faces != faces) {
    std::cout << "Skinner created/deleted faces for vertex-only skinning" << std::endl;
    return MB_FAILURE;
  }
  
  skin.clear();
  rval = tool.find_skin( regions, false, skin, 0, true, false );
  if (MB_SUCCESS != rval) {
    std::cout << "Non-create face skinning failed with: " << mb.get_error_string(rval) << std::endl;
    return rval;
  }
  expected = subtract( all_faces, interior_faces );
  if (expected != skin) {
    std::cout << "Skinner returned incorrect faces." << std::endl;
    return MB_FAILURE;
  }
  if (all_faces != faces) {
    std::cout << "Skinner created/deleted faces for no-create skinning" << std::endl;
    return MB_FAILURE;
  }
  
  skin.clear();
  rval = tool.find_skin( regions, false, skin, 0, true, true );
  if (MB_SUCCESS != rval) {
    std::cout << "Create face skinning failed with: " << mb.get_error_string(rval) << std::endl;
    return rval;
  }
  Range all_faces2;
  rval = mb.get_entities_by_dimension( 0, 2, all_faces2 );
  if (MB_SUCCESS != rval) return rval;
  Range difference = subtract( all_faces2, all_faces );
  if (difference.size() != 2) { // should have created two quads for hex top/bottom
    std::cout << "Skinner failed to create new quads or created to many." << std::endl;
    return MB_FAILURE;
  }
  expected.merge(difference);
  if (expected != skin) {
    std::cout << "Skinner returned incorrect faces." << std::endl;
    return MB_FAILURE;
  }
    // check that new faces are correct
  EntityHandle expected_conn[2][4] = {
    { verts[0][12],verts[0][13],verts[0][14],verts[0][15] },
    { verts[2][12],verts[2][13],verts[2][14],verts[2][15] } };
  EntityHandle nq[2] = { difference.front(), difference.back() };
  for (int i = 0; i < 2; ++i) {
    const EntityHandle* conn;
    int len;
    bool found = false;
    for (int j = 0; !found && j < 2; ++j) {
      rval = mb.get_connectivity( nq[j], conn, len );
      if (MB_SUCCESS != rval) return rval;
      int idx1 = std::find(conn,conn+len,expected_conn[i][0])-conn;
      if (idx1 == len) continue;
      found = true;
      for (int k = 1; k < 4; ++k)
        if (conn[(idx1+k)%4] != expected_conn[i][k])
          found = false;
      if (!found) {
        found = true;
        for (int k = 1; k < 4; ++k)
          if (conn[(idx1+4-k)%4] != expected_conn[i][k])
            found = false;
      }
    }
    if (!found) {
      std::cerr << "Skinner did not create & return expected quad " << i << std::endl;
      return MB_FAILURE;
    }
  }
  
  return MB_SUCCESS;
}

// Test that skinning of higher-order elements works
ErrorCode mb_skin_higher_order_faces_common( bool use_adj )
{
  /* Create mesh:
   
     0---1---2---3---4
     |       |      /
     |       |     /
     5   6   7  8 9
     |       |   /
     |       |  /
     10--11--12
  */

  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  
  double coords[13][3] = { 
   {0,4,0}, {2,4,0}, {4,4,0}, {6,4,0}, {8,4,0},
   {0,2,0}, {2,2,0}, {4,2,0}, {5,2,0}, {6,2,0},
   {0,0,0}, {2,0,0}, {4,0,0} };
  EntityHandle verts[13];
  for (int i = 0; i < 13; ++i) {
    rval = mb.create_vertex( coords[i], verts[i] );
    if (MB_SUCCESS != rval) return rval;
  }
  
  EntityHandle qconn[9] = {
    verts[0], verts[2], verts[12], verts[10],
    verts[1], verts[7], verts[11], verts[5],
    verts[6] };
  EntityHandle tconn[7] = {
    verts[2], verts[4], verts[12],
    verts[3], verts[9], verts[7],
    verts[8] };
  EntityHandle quad, tri;
  rval = mb.create_element( MBQUAD, qconn, 9, quad );
  if (MB_SUCCESS != rval) return rval;
  rval = mb.create_element( MBTRI, tconn, 7, tri );
  if (MB_SUCCESS != rval) return rval;
  
  Range faces;
  faces.insert(quad);
  faces.insert(tri);
  
  Range skin_verts;
  const int skin_vert_idx[] = { 0, 1, 2, 3, 4, 5, 9, 10, 11, 12 };
  for (size_t i = 0; i < sizeof(skin_vert_idx)/sizeof(skin_vert_idx[0]); ++i)
    skin_verts.insert( verts[skin_vert_idx[i]] );
  
  Skinner tool(&mb);
  Range skin;
  
  rval = tool.find_skin( faces, true, skin, 0, use_adj, false );
  if (MB_SUCCESS != rval) {
    std::cout << "Vertex skinning failed with: " << mb.get_error_string(rval) << std::endl;
    return rval;
  }
  if (skin != skin_verts) {
    std::cout << "Skinner returned incorrect vertices." << std::endl;
    return MB_FAILURE;
  }
  
  const int skin_edges[5][3] = {
    {0,1,2}, {2,3,4}, {4,9,12}, {12,11,10}, {10,5,0} };
  skin.clear();
  rval = tool.find_skin( faces, false, skin, 0, use_adj, true );
  if (MB_SUCCESS != rval) {
    std::cout << "Edge skinning failed with: " << mb.get_error_string(rval) << std::endl;
    return rval;
  }
  if (skin.size() != 5u) {
    std::cout << "Skinner returned " << skin.size() << " vertices.  Expected 5" << std::endl;
    return MB_FAILURE;
  }
  int num_quadratic = 0;
  const EntityHandle* conn;
  int len;
  for (Range::iterator i = skin.begin(); i != skin.end(); ++i) {
    rval = mb.get_connectivity( *i, conn, len, false );
    if (MB_SUCCESS != rval) return rval;
    if (len == 3)
      num_quadratic++;
    else if(len != 2) {
      std::cerr << "Skinner created edge with " << len << " vertices" << std::endl;
      return MB_FAILURE;
    }
  }
  if (num_quadratic != 5) {
    std::cerr << num_quadratic << " of 5 created edges were quadratic" << std::endl;
    return MB_FAILURE;
  }
  
  for (int i = 0; i < 5; ++i) {
    bool found = false;
    for (Range::iterator j = skin.begin(); j != skin.end(); ++j) {
      mb.get_connectivity( *j, conn, len, false );
      if (conn[2] == verts[skin_edges[i][1]]) {
        found = true;
        break;
      }
    }
    if (!found) {
      std::cerr << "One or more skin edges is incorrect" << std::endl;
      return MB_FAILURE;
    }
    if ((conn[0] != verts[skin_edges[i][0]] || conn[1] != verts[skin_edges[i][2]])
     && (conn[0] != verts[skin_edges[i][2]] || conn[1] != verts[skin_edges[i][0]])) {
      std::cerr << "Invalid skin edge connectivity" << std::endl;
      return MB_FAILURE;
    }
  }
   
  return MB_SUCCESS;
}
ErrorCode mb_skin_higher_order_faces_test()
  { return mb_skin_higher_order_faces_common( false ); }
ErrorCode mb_skin_adj_higher_order_faces_test()
  { return mb_skin_higher_order_faces_common( true ); }

// Test that skinning of higher-order elements works
ErrorCode mb_skin_higher_order_regions_common( bool use_adj )
{
  // create mesh containing two 27-node hexes
  /*
     0,2---1,2---2,2---3,2---4,2
      |           |           |
      |           |           |
     0,1   1,1   2,1   3,1   4,1
      |           |           |
      |           |           |
     0,0---1,0---2,0---3,0---4,0
  */

  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  Range hexes;
  
  
  EntityHandle verts[5][3][3];
  for (int i = 0; i < 5; ++i) 
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) {
        double coords[] = { i, j, k };
        rval = mb.create_vertex( coords, verts[i][j][k] );
        if (MB_SUCCESS != rval) return rval;
      }
  
  int hex_conn[][3] = {  // corners
                        {0,0,0}, {2,0,0}, {2,2,0}, {0,2,0},
                        {0,0,2}, {2,0,2}, {2,2,2}, {0,2,2},
                         // mid-edge
                        {1,0,0}, {2,1,0}, {1,2,0}, {0,1,0},
                        {0,0,1}, {2,0,1}, {2,2,1}, {0,2,1},
                        {1,0,2}, {2,1,2}, {1,2,2}, {0,1,2},
                        // mid-face
                        {1,0,1}, {2,1,1}, {1,2,1}, {0,1,1},
                        {1,1,0}, {1,1,2},
                        // mid-volume
                        {1,1,1} };
 
  EntityHandle hexverts[2][27];
  for (int i = 0; i < 2; ++i) {
    EntityHandle h;
    for (int j = 0; j < 27; ++j)
      hexverts[i][j] = verts[ 2*i+hex_conn[j][0] ][ hex_conn[j][1] ][ hex_conn[j][2] ];
    rval = mb.create_element( MBHEX, hexverts[i], 27, h );
    if (MB_SUCCESS != rval)
      return rval;
    hexes.insert( h );
  }
  
  Range interior_verts;
  interior_verts.insert( verts[1][1][1] ); // mid-node of hex 1
  interior_verts.insert( verts[3][1][1] ); // mid-node of hex 2
  interior_verts.insert( verts[2][1][1] ); // mid-node of shared face
  
  Skinner tool(&mb);
  Range skin;

  rval = tool.find_skin( hexes, true, skin, 0, use_adj, false );
  if (MB_SUCCESS != rval) {
    std::cout << "Vertex skinning failed with: " << mb.get_error_string(rval) 
              << std::endl;
    return rval;
  }
  Range extra = intersect( skin, interior_verts );
  if (!extra.empty()) {
    std::cout << "Skinner returned " << extra.size() << " interior vertices" 
              << std::endl;
    std::cout << extra << std::endl;
    return MB_FAILURE;
  }
  int num_vtx;
  mb.get_number_entities_by_dimension( 0, 0, num_vtx );
  size_t num_skin = num_vtx - interior_verts.size();
  if (skin.size() != num_skin) {
    std::cout << "Skinner returned " << skin.size() << " of " 
              << num_skin << " skin vertices" <<std::endl;
    return MB_FAILURE;
  }
  
  skin.clear();
  rval = tool.find_skin( hexes, false, skin, 0, use_adj, true );
  if (MB_SUCCESS != rval) {
    std::cout << "Element skinning failed with: " << mb.get_error_string(rval) << std::endl;
    return rval;
  }
  
  if (skin.size() > 10u) {
    std::cout << "Skinner created too many faces" << std::endl;
    return MB_FAILURE;
  }
  
  bool all_okay = true;
  bool faces[2][6] = { {0,0,0,0,0,0},{0,0,0,0,0,0} };
  const EntityHandle *conn;
  int len;
  for (Range::iterator it = skin.begin(); it != skin.end(); ++it) {
    rval = mb.get_connectivity( *it, conn, len );
    if (MB_SUCCESS != rval) return rval;
    if (len != 9) {
      std::cout << "Skinner created face with " << len << " nodes" << std::endl;
      all_okay = false;
      continue;
    }
    
    int valid, side, sense, offset;
    for (int h = 0; h < 2; ++h) {
      valid = CN::SideNumber( MBHEX, hexverts[h], conn, 4, 2, side, sense, offset );
      if (valid != 0)
        continue;
      if (sense != 1) {
        std::cout << "Skinner created reversed face for hex " 
                  << h << " side " << side << std::endl;
        all_okay = false;
        continue;
      }
      
      int idx[9], len2;
      EntityType sidetype;
      CN::SubEntityNodeIndices( MBHEX, 27, 2, side, sidetype, len2, idx );
      assert(sidetype == MBQUAD);
      assert(len2 == 9);
      if ( conn[   offset     ] != hexverts[h][idx[0]] ||
           conn[  (offset+1)%4] != hexverts[h][idx[1]] ||
           conn[  (offset+2)%4] != hexverts[h][idx[2]] ||
           conn[  (offset+3)%4] != hexverts[h][idx[3]] ||
           conn[4+ offset     ] != hexverts[h][idx[4]] ||
           conn[4+(offset+1)%4] != hexverts[h][idx[5]] ||
           conn[4+(offset+2)%4] != hexverts[h][idx[6]] ||
           conn[4+(offset+3)%4] != hexverts[h][idx[7]] ||
           conn[ 8            ] != hexverts[h][idx[8]]) {
        std::cout << "Side " << side << " of hex " << h 
                  << " has invalid connectivity" << std::endl;
        all_okay = false;
      }
           
      
      faces[h][side] = true;
    }
  }
  
  for (int h = 0; h < 2; ++h) {
    for (int i = 0; i < 6; ++i) {
      if ((h == 0 && i == 1) || (h == 1 && i == 3)) {
        if (faces[h][i]) {
          std::cout << "Skinner created interior face for side " 
                    << i << " of hex " << h << std::endl;
          all_okay = false;
        }
      }
      else if (!faces[h][i]) {
        std::cout << "Skinner failed to create side " 
                  << i << " of hex " << h << std::endl;
        all_okay = false;
      }
    }
  }

  return all_okay ? MB_SUCCESS : MB_FAILURE;
}

ErrorCode mb_skin_higher_order_regions_test()
  { return mb_skin_higher_order_regions_common(false); }
ErrorCode mb_skin_adj_higher_order_regions_test()
  { return mb_skin_higher_order_regions_common(true); }


ErrorCode mb_skin_reversed_common( int dim, bool use_adj )
{
  EntityType type, subtype;
  switch (dim) { 
    case 2: type = MBTRI; subtype = MBEDGE; break;
    case 3: type = MBTET; subtype = MBTRI;  break;
    default: assert(false); return MB_FAILURE;
  }
  
  /*      3
         /|\
        / | \
       /  |  \
      /   |   \
     0----1----2
  */

  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  Range hexes;
   
  double coords[][3] = { { 0, 0, 0 },
                         { 1, 0, 0 },
                         { 2, 0, 0 },
                         { 1, 2, 0 },
                         { 1, 2, 2 } };
  EntityHandle verts[5];
  const int nvtx = 2+dim;
  for (int i = 0; i < nvtx; ++i) {
    rval = mb.create_vertex( coords[i], verts[i] );
    if (MB_SUCCESS != rval) return rval;
  }
    // NOTE: order connectivity such that side 1 is shared!
  EntityHandle conn[2][4] = { 
    { verts[0], verts[1], verts[3], verts[4] },
    { verts[2], verts[3], verts[1], verts[4] } };
  const int conn_len = dim+1;
  Range elems;
  for (int i = 0; i < 2; ++i) {
    EntityHandle h;
    rval = mb.create_element( type, conn[i], conn_len, h );
    if (MB_SUCCESS != rval) return rval;
    elems.insert(h);
  }
  
    // create one reversed side
  EntityHandle side_conn[3];
  int side_indices[3] = {0,0,0};
  CN::SubEntityVertexIndices(type, dim-1, 0, side_indices );
  side_conn[0] = conn[0][side_indices[1]];
  side_conn[1] = conn[0][side_indices[0]];
  side_conn[2] = conn[0][side_indices[2]];
  EntityHandle side;
  rval = mb.create_element( subtype, side_conn, dim, side );
  if (MB_SUCCESS != rval) return rval;
  
  Range forward, reverse;
  Skinner tool(&mb);
  rval = tool.find_skin( elems, false, forward, &reverse, use_adj, true );
  if (MB_SUCCESS != rval) {
    std::cout << "Skinner failed." << std::endl;
    return rval;
  }
  
    // expect all faces created by the skinner to be forward,
    // so the only reversed side should be the one created above.
  if (reverse.size() != 1 || reverse.front() != side) {
    std::cout << "Expected 1 reversed side, got: " << reverse.size() << std::endl;
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}
ErrorCode mb_skin_faces_reversed_test()
  { return mb_skin_reversed_common( 2, false ); }
ErrorCode mb_skin_adj_faces_reversed_test()
  { return mb_skin_reversed_common( 2, true ); }
ErrorCode mb_skin_regions_reversed_test()
  { return mb_skin_reversed_common( 3, false ); }
ErrorCode mb_skin_adj_regions_reversed_test()
  { return mb_skin_reversed_common( 3, true ); }


ErrorCode mb_skin_subset_common( int dimension, bool use_adj )
{
  EntityType type;
  switch (dimension) { 
    case 2: type = MBTRI; break;
    case 3: type = MBPRISM; break;
    default: assert(false); return MB_FAILURE;
  }
  

  /*      0
         /|\
        / | \
       5  |  1
       |\ | /|
       | \|/ |
       |  6  |
       | /|\ |
       |/ | \|
       4  |  2
        \ | /
         \|/
          3
   */

  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  Range expected_verts;

  const double coords2D[7][2] = { {0,2}, {1,1}, {1,-1}, {0,-2}, {-1,-1}, {-1,1}, {0,0} };
  EntityHandle verts[2][7] = { {0,0,0,0,0,0,0}, {0,0,0,0,0,0,0} };
  for (int d = 1; d < dimension; ++d) {
    for (int i = 0; i < 7; ++i) {
      double coords[3] = { coords2D[i][0], coords2D[i][1], d-1 };
      rval = mb.create_vertex( coords, verts[d-1][i] );
      if (MB_SUCCESS != rval) return rval;
      if (i != 4 && i != 5)
        expected_verts.insert( verts[d-1][i] );
    }
  }
  
  EntityHandle elems[6];
  for (int i = 0; i < 6; ++i) {
    EntityHandle conn[6] = { verts[0][6], verts[0][(i+1)%6], verts[0][i],
                               verts[1][6], verts[1][(i+1)%6], verts[1][i] };
    rval = mb.create_element( type, conn, CN::VerticesPerEntity(type), elems[i] );
    if (MB_SUCCESS != rval) return rval;
  }
  
  Range input;
  input.insert( elems[0] );
  input.insert( elems[1] );
  input.insert( elems[2] );
  
  Range skin;
  Skinner tool(&mb);
  rval = tool.find_skin( input, true, skin, 0, use_adj, false );
  if (MB_SUCCESS != rval) {
    std::cout << "Skinner failed to find skin vertices" << std::endl;
    return MB_FAILURE;
  }
  if (skin != expected_verts) {
    std::cout << "Skinner returned incorrect skin vertices" << std::endl;
    return MB_FAILURE;
  }
  int n = 0;
  mb.get_number_entities_by_dimension( 0, dimension-1, n );
  if (n > 0) {
    std::cout << "Skinner created lower-dimension entities for vertex-only skinning" << std::endl;
    return MB_FAILURE;
  }
    
  std::vector<EntityHandle> sv( skin.begin(), skin.end() );
  std::vector<int> counts( sv.size(), 0 );
  skin.clear();
  rval = tool.find_skin( input, false, skin, 0, use_adj, true );
  if (MB_SUCCESS != rval) {
    std::cout << "Skinner failed to find skin elements" << std::endl;
    return MB_FAILURE;
  }
  for (Range::iterator i = skin.begin(); i != skin.end(); ++i) {
    const EntityHandle *conn;
    int len;
    rval = mb.get_connectivity( *i, conn, len );
    if (MB_SUCCESS != rval) return rval;
    for (int j = 0; j < len; ++j) {
      size_t idx = std::find(sv.begin(), sv.end(), conn[j]) - sv.begin();
      if (idx == sv.size()) {
        std::cout << "Skinner returned non-skin element" << std::endl;
        return MB_FAILURE;
      }
      counts[idx]++;
    }
  }
  for (size_t i = 0; i < counts.size(); ++i) {
    if (counts[i] < dimension) { // 2 for dim==2, {3,4,5} for dim==3
      std::cout << "Skinner did not return all skin elements" << std::endl;
      return MB_FAILURE;
    }
  }
  mb.get_number_entities_by_dimension( 0, dimension-1, n );
  if ((size_t)n != skin.size()) {
    std::cout << "Skinner created extra lower-dimension entities" << std::endl;
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

ErrorCode mb_skin_faces_subset_test()
  { return mb_skin_subset_common( 2, false ); }
ErrorCode mb_skin_adj_faces_subset_test()
  { return mb_skin_subset_common( 2, true ); }
ErrorCode mb_skin_regions_subset_test()
  { return mb_skin_subset_common( 3, false ); }
ErrorCode mb_skin_adj_regions_subset_test()
  { return mb_skin_subset_common( 3, true ); }
  
  
    
 
ErrorCode mb_skin_full_common( int dimension, bool use_adj )
{
  EntityType type;
  switch (dimension) { 
    case 2: type = MBQUAD; break;
    case 3: type = MBHEX; break;
    default: assert(false); return MB_FAILURE;
  }
  

  /*  
      3----4----5
      |    |    |
      |    |    |
      0----1----2
  */
  
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
 
    // create vertices
  const double coords2D[6][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1} };
  EntityHandle v[2][6] = { {0,0,0,0,0,0}, {0,0,0,0,0,0} };
  for (int d = 1; d < dimension; ++d) {
    for (int i = 0; i < 6; ++i) {
      double coords[3] = { coords2D[i][0], coords2D[i][1], d-1 };
      rval = mb.create_vertex( coords, v[d-1][i] );
      if (MB_SUCCESS != rval) return rval;
    }
  }
  
    // create elements
  Range input;
  EntityHandle elems[2], econn[2][8];;
  for (int i = 0; i < 2; ++i) {
    EntityHandle conn[8] = { v[0][i], v[0][i+1], v[0][i+4], v[0][i+3],
                               v[1][i], v[1][i+1], v[1][i+4], v[1][i+3] };
    memcpy( econn[i], conn, sizeof(conn) );
    rval = mb.create_element( type, conn, CN::VerticesPerEntity(type), elems[i] );
    if (MB_SUCCESS != rval) return rval;
    input.insert( elems[i] );
  }
  
    // create sides
    // NOTE: Shared side is element 0 side 1 and element 1 side 3
  Range expected;
  for (int i = 0; i < CN::NumSubEntities( type, dimension-1 ); ++i) {
    EntityType subtype;
    int len;
    const short* indices = CN::SubEntityVertexIndices( type, dimension-1,
                                                         i, subtype, len );
    EntityHandle conn[4];
    assert((size_t)len <= sizeof(conn)/sizeof(conn[0]));
    for (int j = 0; j < 2; ++j) {
      if (j == 1 && i == 3) // don't create shared face twice
        continue;
      for (int k = 0; k < len; ++k)
        conn[k] = econn[j][indices[k]];
      EntityHandle h;
      rval = mb.create_element( subtype, conn, len, h );
      if (MB_SUCCESS != rval) return rval;
      if (j != 0 || i != 1) // don't insert shared face
        expected.insert(h);
    }
  }
  
  Range skin;
  Skinner tool(&mb);
  rval = tool.find_skin( input, false, skin, 0, use_adj, true );
  if (MB_SUCCESS != rval) {
    std::cout << "Skinner failed to find skin elements" << std::endl;
    return MB_FAILURE;
  }
  if (skin != expected) {
    std::cout << "Skinner returned incorrect skin elements" << std::endl;
    return MB_FAILURE;
  }

  int n = 0;
  mb.get_number_entities_by_dimension( 0, dimension-1, n );
  if ((size_t)n != expected.size()+1) {
    std::cout << "Skinner created extra lower-dimension entities" << std::endl;
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

ErrorCode mb_skin_faces_full_test()
  { return mb_skin_full_common( 2, false ); }
ErrorCode mb_skin_adj_faces_full_test()
  { return mb_skin_full_common( 2, true ); }
ErrorCode mb_skin_regions_full_test()
  { return mb_skin_full_common( 3, false ); }
ErrorCode mb_skin_adj_regions_full_test()
  { return mb_skin_full_common( 3, true ); }
        
ErrorCode mb_skin_adjacent_surf_patches()
{
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  Skinner tool(&mb);
  
  /* Mesh with vertices and quads numbered. 
     Groups are indicated by letters {A,B,C,D}
  
     0----1----2----3----4----5
     | (0)| (1)| (2)| (3)| (4)|
     |  A |  A |  B |  B |  B |
     6----7----8----9---10---11
     | (5)| (6)| (7)| (8)| (9)|
     |  A |  A |  A |  B |  B |
    12---13---14---15---16---17
     |(10)|(11)|(12)|(13)|(14)|
     |  A |  C |  D |  D |  D |
    18---19---20---21---22---23          
     |(15)|(16)|(17)|(18)|(19)|
     |  C |  C |  C |  D |  D |
    24---25---26---27---28---29
  */
  
  const int num_vtx = 30;
  EntityHandle vtx[num_vtx];
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 5; ++j) {
      double coords[3] = {i, j, 0};
      rval = mb.create_vertex( coords, vtx[6*j+i] );
      if (MB_SUCCESS != rval) return rval;
    }
  }
  
  EntityHandle quads[20];
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 4; ++j) {
      int v = 6*j+i;
      EntityHandle conn[4] = { vtx[v+6], vtx[v+7], vtx[v+1], vtx[v+0] };
      rval = mb.create_element( MBQUAD, conn, 4, quads[5*j+i] );
      if (MB_SUCCESS != rval) return rval;
    }
  }
  
    // Define four groups of quads (as labeled above)
  const int Aquads[] = { 0, 1, 5, 6, 7, 10 };
  const int Aquads_size = sizeof(Aquads)/sizeof(Aquads[0]);
  const int Bquads[] = { 2, 3, 4, 8, 9 };
  const int Bquads_size = sizeof(Bquads)/sizeof(Bquads[0]);
  const int Cquads[] = { 11, 15, 16, 17 };
  const int Cquads_size = sizeof(Cquads)/sizeof(Cquads[0]);
  const int Dquads[] = { 12, 13, 14, 18, 19 };
  const int Dquads_size = sizeof(Dquads)/sizeof(Dquads[0]);
  
    // Define the results we expect for each group as a loop
    // of vertex indices
  const int Askin[] =   { 0, 1, 2, 8, 9, 15, 14, 13, 19, 18, 12, 6 };
  const int Ashared[] = {-1,-1, 1, 1, 1,  3,  2,  2,  2, -1, -1,-1 };
  const int Askin_size = sizeof(Askin)/sizeof(Askin[0]);
  const int Bskin[] =   { 2, 3, 4, 5, 11, 17, 16, 15, 9, 8 };
  const int Bshared[] = {-1,-1,-1,-1, -1,  3,  3,  0, 0, 0 };
  const int Bskin_size = sizeof(Bskin)/sizeof(Bskin[0]);
  const int Cskin[] =   { 18, 19, 13, 14, 20, 21, 27, 26, 25, 24 };
  const int Cshared[] = {  0,  0,  0,  3,  3,  3, -1, -1, -1, -1 };
  const int Cskin_size = sizeof(Cskin)/sizeof(Cskin[0]);
  const int Dskin[] =   { 14, 15, 16, 17, 23, 29, 28, 27, 21, 20 };
  const int Dshared[] = {  0,  1,  1, -1, -1, -1, -1,  2,  2,  2 };
  const int Dskin_size = sizeof(Dskin)/sizeof(Dskin[0]);
  
    // Make the above stuff indexable for easier looping
  const int* const gquads[4] = { Aquads, Bquads, Cquads, Dquads };
  const int gquads_size[4] = { Aquads_size, Bquads_size, Cquads_size, Dquads_size };
  const int* const skin[4] = { Askin, Bskin, Cskin, Dskin };
  const int* const shared[4] = { Ashared, Bshared, Cshared, Dshared };
  const int skin_size[4] = { Askin_size, Bskin_size, Cskin_size, Dskin_size };
  
    // Create an Range for each group of quads
  Range ranges[4];
  for (int grp = 0; grp < 4; ++grp) 
    for (int i = 0; i < gquads_size[grp]; ++i)
      ranges[grp].insert( quads[gquads[grp][i]] );
  
    // run test 4 times, one for each of:
    // o no adjacencies, no edges
    // o no adjacencies, edges
    // o adjacencies, edges
    // o adjacencies, no edges
    //
    // Be careful here: once we specify use_adj, we can't
    // unspecify it (the adjacencies will exist so the skinner
    // will use them, regardless of the passed flag.)
  int error_count = 0;
  for (int run = 0; run < 4; ++run) {
    const bool use_adj = run > 1;
    if (run == 3) {
      Range dead;
      mb.get_entities_by_type( 0, MBEDGE, dead );
      mb.delete_entities( dead );
    }
    
      // test each group
    Range edges[4];
    for (int grp = 0; grp < 4; ++grp) {
        // get the skin edges
      rval = tool.find_skin( ranges[grp], 1, edges[grp], use_adj );
      if (MB_SUCCESS != rval) {
        std::cout << "Skinner failed for run " << run << " group " << grp << std::endl;
        return rval;
      }
      
        // check that we have the expected result
      std::vector<bool> seen(skin_size[grp], false);
      for (Range::iterator e = edges[grp].begin(); e != edges[grp].end(); ++e) {
        const EntityHandle* conn;
        int len;
        rval = mb.get_connectivity( *e, conn, len );
        if (MB_SUCCESS != rval) return rval;
        if (len != 2) return MB_FAILURE;
        const int idx1 = std::find( vtx, vtx+num_vtx, conn[0] ) - vtx;
        const int idx2 = std::find( vtx, vtx+num_vtx, conn[1] ) - vtx;
        int pos = std::find( skin[grp], skin[grp]+skin_size[grp], idx1 ) - skin[grp];
        if (pos == skin_size[grp]) {
          std::cout << "Non-skin vertex in skin for run " << run << " group " << grp << std::endl;
          std::cout << "idx1 = " << idx1 << ", idx2 = " << idx2 << std::endl;
          ++error_count;
          continue;
        }
        
        
        if (skin[grp][(pos+skin_size[grp]-1)%skin_size[grp]] == idx2)
          pos = (pos + skin_size[grp] - 1)%skin_size[grp];
        else if (skin[grp][(pos+1)%skin_size[grp]] != idx2) {
          std::cout << "Non-skin edge in skin for run " << run << " group " << grp << std::endl;
          std::cout << "idx1 = " << idx1 << ", idx2 = " << idx2 << std::endl;
          ++error_count;
          continue;
        }
        
        if (seen[pos]) {
          std::cout << "Duplicate edge in skin for run " << run << " group " << grp << std::endl;
          std::cout << "idx1 = " << idx1 << ", idx2 = " << idx2 << std::endl;
          ++error_count;
        }
        seen[pos] = true;
        
        int shared_with = shared[grp][pos];
        if (shared_with < 0) // not shared with another group
          continue;
        if (shared_with > grp) // didn't skin other group yet
          continue;
        if (edges[shared_with].find(*e) == edges[shared_with].end()) {
          std::cout << "Skin edge duplicated for run " << run << " group " << grp << std::endl;
          std::cout << "idx1 = " << idx1 << ", idx2 = " << idx2 
                    << " not in skin for group " << shared_with << std::endl;
          ++error_count;
        }
      }
      
      int missing = std::count( seen.begin(), seen.end(), false );
      if (missing) {
        std::cout << "Missking " << missing << " skin edges for run " << run << " group " << grp << std::endl;
        error_count += missing;
      }
    }
  }
  
  return error_count ? MB_FAILURE : MB_SUCCESS;
}


static ErrorCode get_by_all_types_and_tag( Interface* mb,
                                           EntityHandle meshset,
                                           const Tag* tag_handles,
                                           const void* const* values,
                                           int num_tags,
                                           Range& result,
                                           int condition,
                                           bool recursive )
{
  ErrorCode rval;
  Range tmp;
  const EntityType LAST = recursive ? MBENTITYSET : MBMAXTYPE;
  for (EntityType t = MBVERTEX; t < LAST; ++t) {
    tmp.clear();
    rval = mb->get_entities_by_type_and_tag( meshset, t, tag_handles, values, num_tags, tmp, condition, recursive );
    if (MB_SUCCESS != rval)
      return rval;
    result.insert( tmp.begin(), tmp.end() );
  }
  return MB_SUCCESS;
}


/** Check that functions which accept a type return the
 *  result for all types when passed MBMAXTYPE
 */
ErrorCode mb_type_is_maxtype_test()
{
  Core moab;
  Interface* mb = &moab;
  ErrorCode rval = create_some_mesh( mb );
  if (MB_SUCCESS != rval)
    return rval;
  
  Range r1, r2;
  rval = mb->get_entities_by_type( 0, MBMAXTYPE, r1, false ); CHKERR(rval);
  rval = mb->get_entities_by_handle( 0, r2, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  
  std::vector<EntityHandle> v1, v2;
  rval = mb->get_entities_by_type( 0, MBMAXTYPE, v1, false ); CHKERR(rval);
  rval = mb->get_entities_by_handle( 0, v2, false ); CHKERR(rval);
  CHECK( v1 == v2 );
  
  int c1, c2;
  rval = mb->get_number_entities_by_type( 0, MBMAXTYPE, c1, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_handle( 0, c2, false ); CHKERR(rval);
  CHECK( c1 == c2 );
  
  Range h1, h2;
  Range::iterator it = r1.begin() + r1.size()/2;
  h1.insert( r1.begin(), it );
  if (it != r1.end())
    h2.insert( ++it, r1.end() );
  
  EntityHandle s1, s2;
  rval = mb->create_meshset( MESHSET_SET, s1 ); CHKERR(rval);
  rval = mb->create_meshset( MESHSET_ORDERED, s2 ); CHKERR(rval);
  rval = mb->add_entities( s1, r1 ); CHKERR(rval);
  rval = mb->add_entities( s2, r2 ); CHKERR(rval);
  rval = mb->add_entities( s2, &s1, 1 ); CHKERR(rval);
  
  r1.clear();
  r2.clear();
  rval = mb->get_entities_by_type( s1, MBMAXTYPE, r1, false ); CHKERR(rval);
  rval = mb->get_entities_by_handle( s1, r2, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  
  r1.clear();
  r2.clear();
  rval = mb->get_entities_by_type( s2, MBMAXTYPE, r1, false ); CHKERR(rval);
  rval = mb->get_entities_by_handle( s2, r2, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  
  r1.clear();
  r2.clear();
  rval = mb->get_entities_by_type( s2, MBMAXTYPE, r1, true ); CHKERR(rval);
  rval = mb->get_entities_by_handle( s2, r2, true ); CHKERR(rval);
  CHECK( r1 == r2 );
 
   
  v1.clear();
  v2.clear();
  rval = mb->get_entities_by_type( s1, MBMAXTYPE, v1, false ); CHKERR(rval);
  rval = mb->get_entities_by_handle( s1, v2, false ); CHKERR(rval);
  CHECK( v1 == v2 );
  
  v1.clear();
  v2.clear();
  rval = mb->get_entities_by_type( s2, MBMAXTYPE, v1, false ); CHKERR(rval);
  rval = mb->get_entities_by_handle( s2, v2, false ); CHKERR(rval);
  CHECK( v1 == v2 );
  
  v1.clear();
  v2.clear();
  rval = mb->get_entities_by_type( s2, MBMAXTYPE, v1, true ); CHKERR(rval);
  rval = mb->get_entities_by_handle( s2, v2, true ); CHKERR(rval);
  CHECK( v1 == v2 );
 
   
  rval = mb->get_number_entities_by_type( s1, MBMAXTYPE, c1, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_handle( s1, c2, false ); CHKERR(rval);
  CHECK( c1 == c2 );
  
  rval = mb->get_number_entities_by_type( s2, MBMAXTYPE, c1, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_handle( s2, c2, false ); CHKERR(rval);
  CHECK( c1 == c2 );
  
  rval = mb->get_number_entities_by_type( s2, MBMAXTYPE, c1, true ); CHKERR(rval);
  rval = mb->get_number_entities_by_handle( s2, c2, true ); CHKERR(rval);
  CHECK( c1 == c2 );
 
  r1.clear();
  rval = mb->get_entities_by_handle( s1, r1 ); CHKERR(rval);
  Tag t1;
  rval = mb->tag_get_handle( "maxtype1", 1, MB_TYPE_INTEGER, t1, MB_TAG_SPARSE|MB_TAG_EXCL );
  CHKERR(rval);
  std::vector<int> d1(r1.size());
  Range::iterator ri;
  std::vector<int>::iterator ii = d1.begin();
  for (ri = r1.begin(); ri != r1.end(); ++ri, ++ii)
    *ii = ((int)ID_FROM_HANDLE(*ri)) % 20;
  rval = mb->tag_set_data( t1, r1, &d1[0] ); CHKERR(rval);
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( 0, MBMAXTYPE, &t1, 0, 1, r1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( 0, MBMAXTYPE, &t1, 0, 1, c1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, 0, &t1, 0, 1, r2, Interface::INTERSECT, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s1, MBMAXTYPE, &t1, 0, 1, r1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s1, MBMAXTYPE, &t1, 0, 1, c1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s1, &t1, 0, 1, r2, Interface::INTERSECT, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s2, MBMAXTYPE, &t1, 0, 1, r1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s2, MBMAXTYPE, &t1, 0, 1, c1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s2, &t1, 0, 1, r2, Interface::INTERSECT, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s2, MBMAXTYPE, &t1, 0, 1, r1, Interface::INTERSECT, true ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s2, MBMAXTYPE, &t1, 0, 1, c1, Interface::INTERSECT, true ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s2, &t1, 0, 1, r2, Interface::INTERSECT, true ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  int value = 3;
  const void* vallist[2] = { &value, 0 };
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( 0, MBMAXTYPE, &t1, vallist, 1, r1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( 0, MBMAXTYPE, &t1, vallist, 1, c1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, 0, &t1, vallist, 1, r2, Interface::INTERSECT, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s1, MBMAXTYPE, &t1, vallist, 1, r1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s1, MBMAXTYPE, &t1, vallist, 1, c1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s1, &t1, vallist, 1, r2, Interface::INTERSECT, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s2, MBMAXTYPE, &t1, vallist, 1, r1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s2, MBMAXTYPE, &t1, vallist, 1, c1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s2, &t1, vallist, 1, r2, Interface::INTERSECT, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s2, MBMAXTYPE, &t1, vallist, 1, r1, Interface::INTERSECT, true ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s2, MBMAXTYPE, &t1, vallist, 1, c1, Interface::INTERSECT, true ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s2, &t1, vallist, 1, r2, Interface::INTERSECT, true ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_handle( s1, r1 ); CHKERR(rval);
  r2.insert( r1.back() );
  r1.clear();
  rval = mb->get_entities_by_handle( s2, r1 ); CHKERR(rval);
  r2.insert( r1.front() );
  
  Tag t2;
  rval = mb->tag_get_handle( "maxtype2", 1, MB_TYPE_INTEGER, t2, MB_TAG_DENSE|MB_TAG_EXCL );
  CHKERR(rval);
  d1.resize(r2.size());
  ii = d1.begin();;
  for (ri = r2.begin(); ri != r2.end(); ++ri, ++ii)
    *ii = ((int)ID_FROM_HANDLE(*ri)) % 2;
  rval = mb->tag_set_data( t2, r2, &d1[0] ); CHKERR(rval);
  
  Tag tags[] = { t1, t2 };
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( 0, MBMAXTYPE, tags, 0, 2, r1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( 0, MBMAXTYPE, tags, 0, 2, c1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, 0, tags, 0, 2, r2, Interface::INTERSECT, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s1, MBMAXTYPE, tags, 0, 2, r1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s1, MBMAXTYPE, tags, 0, 2, c1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s1, tags, 0, 2, r2, Interface::INTERSECT, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s2, MBMAXTYPE, tags, 0, 2, r1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s2, MBMAXTYPE, tags, 0, 2, c1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s2, tags, 0, 2, r2, Interface::INTERSECT, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s2, MBMAXTYPE, tags, 0, 2, r1, Interface::INTERSECT, true ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s2, MBMAXTYPE, tags, 0, 2, c1, Interface::INTERSECT, true ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s2, tags, 0, 2, r2, Interface::INTERSECT, true ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  int val2 = 1;
  vallist[1] = &val2;
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( 0, MBMAXTYPE, tags, vallist, 2, r1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( 0, MBMAXTYPE, tags, vallist, 2, c1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, 0, tags, vallist, 2, r2, Interface::INTERSECT, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s1, MBMAXTYPE, tags, vallist, 2, r1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s1, MBMAXTYPE, tags, vallist, 2, c1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s1, tags, vallist, 2, r2, Interface::INTERSECT, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s2, MBMAXTYPE, tags, vallist, 2, r1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s2, MBMAXTYPE, tags, vallist, 2, c1, Interface::INTERSECT, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s2, tags, vallist, 2, r2, Interface::INTERSECT, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s2, MBMAXTYPE, tags, vallist, 2, r1, Interface::INTERSECT, true ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s2, MBMAXTYPE, tags, vallist, 2, c1, Interface::INTERSECT, true ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s2, tags, vallist, 2, r2, Interface::INTERSECT, true ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
   
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( 0, MBMAXTYPE, tags, vallist, 2, r1, Interface::UNION, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( 0, MBMAXTYPE, tags, vallist, 2, c1, Interface::UNION, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, 0, tags, vallist, 2, r2, Interface::UNION, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s1, MBMAXTYPE, tags, vallist, 2, r1, Interface::UNION, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s1, MBMAXTYPE, tags, vallist, 2, c1, Interface::UNION, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s1, tags, vallist, 2, r2, Interface::UNION, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s2, MBMAXTYPE, tags, vallist, 2, r1, Interface::UNION, false ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s2, MBMAXTYPE, tags, vallist, 2, c1, Interface::UNION, false ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s2, tags, vallist, 2, r2, Interface::UNION, false ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  r1.clear(); r2.clear();
  rval = mb->get_entities_by_type_and_tag( s2, MBMAXTYPE, tags, vallist, 2, r1, Interface::UNION, true ); CHKERR(rval);
  rval = mb->get_number_entities_by_type_and_tag( s2, MBMAXTYPE, tags, vallist, 2, c1, Interface::UNION, true ); CHKERR(rval);
  rval = get_by_all_types_and_tag( mb, s2, tags, vallist, 2, r2, Interface::UNION, true ); CHKERR(rval);
  CHECK( r1 == r2 );
  CHECK( (unsigned)c1 == r2.size() );
  
  return MB_SUCCESS;
} 

/** Test behavior of various functions when passed the root set
 */
ErrorCode mb_root_set_test()
{
  ErrorCode rval;
  Core moab;
  Interface* mb = &moab;
  EntityHandle rs = mb->get_root_set();
  
    // expect root set to have zero handle
  CHECK(!rs);
  
    // check type
  EntityType type = mb->type_from_handle( rs );
  CHECK( MBENTITYSET == type );
  
    // Create a set to test with
  EntityHandle some_set;
  rval = mb->create_meshset( MESHSET_SET, some_set );
  CHKERR( rval );
  Range sets;
  sets.insert( some_set );
  
  int exp_dim = mb->dimension_from_handle( some_set );
  int dim = mb->dimension_from_handle( rs );
  CHECK( exp_dim == dim );
  
    // test some stuff that should fail
  rval = mb->clear_meshset( &rs, 1 );
  CHECK( MB_SUCCESS != rval );
  rval = mb->set_meshset_options( rs, MESHSET_ORDERED );
  CHECK( MB_SUCCESS != rval );
  rval = mb->subtract_meshset( rs, some_set );
  CHECK( MB_SUCCESS != rval );
  rval = mb->intersect_meshset( rs, some_set );
  CHECK( MB_SUCCESS != rval );
  rval = mb->unite_meshset( rs, some_set );
  CHECK( MB_SUCCESS != rval );
  
    // add/remove should fail for root set?
  rval = mb->add_entities( rs, &some_set, 1 );
  CHECK( MB_SUCCESS != rval );
  rval = mb->remove_entities( rs, &some_set, 1 );
  CHECK( MB_SUCCESS != rval );
  rval = mb->replace_entities( rs, &some_set, &some_set, 1 );
  CHECK( MB_SUCCESS != rval );
   
    // check flags
  unsigned flags;
  rval = mb->get_meshset_options( rs, flags );
  CHECK( flags & MESHSET_SET );
  CHECK( !(flags & MESHSET_ORDERED) );
  CHECK( flags & MESHSET_TRACK_OWNER );
  
    // contains tests
  bool c = mb->contains_entities( rs, &some_set, 1 );
  CHECK( c );
  
  Range sets2;
  rval = mb->get_contained_meshsets( rs, sets2 );
  CHKERR(rval);
  CHECK( sets == sets2 );
  
  int count;
  rval = mb->num_contained_meshsets( rs, &count );
  CHKERR(rval);
  CHECK( count == (int)sets.size() );
  
    // The expected behavior for parent/child queries on the root set
    // is to return an error.
    
  rval = mb->get_parent_meshsets( rs, sets2 );
  CHECK( MB_SUCCESS != rval );
  rval = mb->get_parent_meshsets( rs, sets2, 2 );
  CHECK( MB_SUCCESS != rval );
    
  rval = mb->get_child_meshsets( rs, sets2 );
  CHECK( MB_SUCCESS != rval );
  rval = mb->get_child_meshsets( rs, sets2, 2 );
  CHECK( MB_SUCCESS != rval );
    
  rval = mb->num_parent_meshsets( rs, &count );
  CHECK( MB_SUCCESS != rval );
  rval = mb->num_parent_meshsets( rs, &count, 2 );
  CHECK( MB_SUCCESS != rval );
    
  rval = mb->num_child_meshsets( rs, &count );
  CHECK( MB_SUCCESS != rval );
  rval = mb->num_child_meshsets( rs, &count, 2 );
  CHECK( MB_SUCCESS != rval );
  
  rval = mb->add_parent_meshset( rs, some_set );
  CHECK(MB_SUCCESS != rval);
  rval = mb->add_parent_meshsets( rs, &some_set, 1 );
  CHECK(MB_SUCCESS != rval);
  rval = mb->add_child_meshset( rs, some_set );
  CHECK(MB_SUCCESS != rval);
  rval = mb->add_child_meshsets( rs, &some_set, 1 );
  CHECK(MB_SUCCESS != rval);
  rval = mb->add_parent_child( rs, some_set );
  CHECK(MB_SUCCESS != rval);
  rval = mb->remove_parent_child( rs, some_set );
  CHECK(MB_SUCCESS != rval);
  rval = mb->remove_parent_meshset( rs, some_set );
  CHECK(MB_SUCCESS != rval);
  rval = mb->remove_child_meshset( rs, some_set );
  CHECK(MB_SUCCESS != rval);
  
  return MB_SUCCESS;
}


/* Create a regular 2x2x2 hex mesh */
ErrorCode create_some_mesh( Interface* iface )
{
  const double coords[] = { 0, 0, 0,
                            1, 0, 0,
                            2, 0, 0,
                            0, 1, 0,
                            1, 1, 0,
                            2, 1, 0,
                            0, 2, 0,
                            1, 2, 0,
                            2, 2, 0,
                            0, 0, 1,
                            1, 0, 1,
                            2, 0, 1,
                            0, 1, 1,
                            1, 1, 1,
                            2, 1, 1,
                            0, 2, 1,
                            1, 2, 1,
                            2, 2, 1,
                            0, 0, 2,
                            1, 0, 2,
                            2, 0, 2,
                            0, 1, 2,
                            1, 1, 2,
                            2, 1, 2,
                            0, 2, 2,
                            1, 2, 2,
                            2, 2, 2 };
  const size_t num_vtx = sizeof(coords)/sizeof(double)/3;
  assert(num_vtx == 27u);
  
  const int conn[] = {  0,  1,  4,  3,  9, 10, 13, 12,
                        1,  2,  5,  4, 10, 11, 14, 13,
                        3,  4,  7,  6, 12, 13, 16, 15,
                        4,  5,  8,  9, 13, 14, 17, 16,
                        9, 10, 13, 12, 18, 19, 22, 21,
                       10, 11, 14, 13, 19, 20, 23, 22,
                       12, 13, 16, 15, 21, 22, 25, 24,
                       13, 14, 17, 18, 22, 23, 26, 25 };
  const size_t num_elem = sizeof(conn)/sizeof(conn[0])/8;
  assert(num_elem == 8u);
                       
  EntityHandle verts[num_vtx], hexes[num_elem];
  for (size_t i = 0; i < num_vtx; ++i) {
    ErrorCode err = iface->create_vertex( coords + 3*i, verts[i] );
    if (MB_SUCCESS != err) return err;
  }
              
  for (size_t i = 0; i < num_elem; ++i) {
    EntityHandle c[8];
    for (int j = 0; j < 8; ++j) {
      assert(conn[8*i+j] < (int)num_vtx);
      c[j] = verts[conn[8*i+j]];
    }
    ErrorCode err = iface->create_element( MBHEX, c, 8, hexes[i] );
    if (MB_SUCCESS != err) return err;
  }
  
  return MB_SUCCESS;            
}

inline bool contained( const std::vector<EntityHandle>& list, EntityHandle h )
{
  std::vector<EntityHandle>::const_iterator i;
  i = std::lower_bound( list.begin(), list.end(), h );
  return i != list.end() && *i == h;
}

ErrorCode check_valid_connectivity( Interface* iface )
{
  ErrorCode rval;

  // get sorted array of vertex handles so that we can
  // check that all vertices in connectivity are valid
  // handles
  std::vector<EntityHandle> vertices, storage;
  rval = iface->get_entities_by_type( 0, MBVERTEX, vertices ); CHKERR(rval);
  std::sort( vertices.begin(), vertices.end() );
  
  // get all the elements
  Range elements, tmp;
  for (int d = 1; d < 4; ++d) {
    tmp.clear();
    rval = iface->get_entities_by_dimension( 0, d, tmp ); CHKERR(rval);
    elements.merge(tmp);
  }
  
  // check that all connectivity handles are valid
  Range::iterator it;
  ErrorCode result = MB_SUCCESS;
  for (it = elements.begin(); it != elements.end(); ++it) {
    const EntityHandle* conn;
    int len;
    rval = iface->get_connectivity( *it, conn, len, false, &storage ); CHKERR(rval);
    for (int i = 0; i < len; ++i) {
      if (!contained( vertices, conn[i] )) {
        printf("Invalid handle (%s %d) in connectivity of %s %d\n",
          CN::EntityTypeName(TYPE_FROM_HANDLE(conn[i])),
          (int)ID_FROM_HANDLE(conn[i]),
          CN::EntityTypeName(TYPE_FROM_HANDLE(*it)),
          (int)ID_FROM_HANDLE(*it));
        result = MB_FAILURE;
      }
    }
  }
  
  return result;
}


static void usage(const char* exe) {
  cerr << "Usage: " << exe << " [-nostress] [-d input_file_dir]\n";
  exit (1);
}


int number_tests = 0;
int number_tests_failed = 0;
#define RUN_TEST( A ) _run_test( (A), #A )


typedef ErrorCode (*TestFunc)();
static void _run_test( TestFunc func, const char* func_str ) 
{
  ++number_tests;
  cout << "   " << func_str << ": ";
  cout.flush();
  ErrorCode error = func( );
  
  if (MB_SUCCESS == error)
    std::cout << "Success" << std::endl;
  else if (MB_FAILURE == error)
    std::cout << "Failure" << std::endl;
  else {
    Core moab;
    std::cout << "Failed: " << moab.get_error_string( error ) << std::endl;
  }
  
  if (MB_SUCCESS != error) {
    ++number_tests_failed;
  }
}


/*!
  main routine for test harness 
*/


int main(int argc, char* argv[])
{
#ifdef USE_MPI
  int fail = MPI_Init(&argc, &argv);
  if (fail) return fail;
#endif

  argv0 = argv[0];

    // Check command line arg to see if we should avoid doing the stress test
  bool stress_test = true;

  std::cout << "Size of mConnMap = " << sizeof(CN::mConnectivityMap)
            << std::endl;
  std::cout << "Size of mUpConnMap = " << sizeof(CN::mUpConnMap)
            << std::endl;
  std::cout << "Size of CN = " << sizeof(CN)
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
  RUN_TEST( mb_adjacencies_create_delete_test );
  RUN_TEST( mb_upward_adjacencies_test );
  RUN_TEST( mb_adjacent_create_test );
  RUN_TEST( mb_vertex_coordinate_test );
  RUN_TEST( mb_vertex_tag_test );
  RUN_TEST( mb_temporary_test );
  RUN_TEST( mb_mesh_sets_set_test );
  RUN_TEST( mb_mesh_sets_list_test );
  RUN_TEST( mb_mesh_set_parent_child_test );
  RUN_TEST( mb_mesh_set_set_appends );
  RUN_TEST( mb_mesh_set_list_appends );
  RUN_TEST( mb_mesh_set_root_appends );
  RUN_TEST( mb_mesh_set_set_replace_test );
  RUN_TEST( mb_mesh_set_list_replace_test );
  RUN_TEST( mb_mesh_set_flag_test );
  RUN_TEST( mb_mesh_set_tracking_test );
  RUN_TEST( mb_mesh_set_set_add_remove_test );
  RUN_TEST( mb_tags_test );
  RUN_TEST( mb_dense_tag_test );
  RUN_TEST( mb_sparse_tag_test );
  RUN_TEST( mb_delete_mesh_test );
  RUN_TEST( mb_higher_order_test );
  RUN_TEST( mb_bit_tags_test );
  RUN_TEST( mb_entity_conversion_test );
  RUN_TEST( mb_forced_adjacencies_test );
  RUN_TEST( mb_canon_number_test );
  RUN_TEST( mb_poly_test );
  RUN_TEST( mb_topo_util_test );
  RUN_TEST( mb_split_test );
  RUN_TEST( mb_range_seq_intersect_test );
  RUN_TEST( mb_poly_adjacency_test );
  RUN_TEST( mb_memory_use_test );
  RUN_TEST( mb_skin_curve_test );
  RUN_TEST( mb_skin_curve_adj_test );
  RUN_TEST( mb_skin_surface_test );
  RUN_TEST( mb_skin_surface_adj_test );
  RUN_TEST( mb_skin_volume_test );
  RUN_TEST( mb_skin_volume_adj_test );
  RUN_TEST( mb_skin_surf_verts_test );
  RUN_TEST( mb_skin_vol_verts_test );
  RUN_TEST( mb_skin_surf_verts_elems_test );
  RUN_TEST( mb_skin_vol_verts_elems_test );
  RUN_TEST( mb_skin_poly_test );
  RUN_TEST( mb_skin_higher_order_faces_test );
  RUN_TEST( mb_skin_higher_order_regions_test );
  RUN_TEST( mb_skin_adj_higher_order_faces_test );
  RUN_TEST( mb_skin_adj_higher_order_regions_test );
  RUN_TEST( mb_skin_faces_reversed_test );
  RUN_TEST( mb_skin_adj_faces_reversed_test );
  RUN_TEST( mb_skin_regions_reversed_test );
  RUN_TEST( mb_skin_adj_regions_reversed_test );
  RUN_TEST( mb_skin_faces_subset_test );
  RUN_TEST( mb_skin_adj_faces_subset_test );
  RUN_TEST( mb_skin_regions_subset_test );
  RUN_TEST( mb_skin_adj_regions_subset_test );
  RUN_TEST( mb_skin_faces_full_test );
  RUN_TEST( mb_skin_adj_faces_full_test );
  RUN_TEST( mb_skin_regions_full_test );
  RUN_TEST( mb_skin_adj_regions_full_test );
  RUN_TEST( mb_skin_adjacent_surf_patches );
  RUN_TEST(mb_skin_scd_test);
  RUN_TEST( mb_read_fail_test );
  RUN_TEST( mb_enum_string_test );
  RUN_TEST( mb_merge_update_test );
  RUN_TEST( mb_type_is_maxtype_test );
  RUN_TEST( mb_root_set_test );
  RUN_TEST( mb_merge_test );
  if (stress_test) RUN_TEST( mb_stress_test );

    // summary

  cout << "\nMB TEST SUMMARY: \n"
       << "   Number Tests:           " << number_tests << "\n"
       << "   Number Successful:      " << number_tests - number_tests_failed << "\n"
       << "   Number Failed:          " << number_tests_failed 
       << "\n\n";

#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail) return fail;
#endif

  return number_tests_failed;
}
