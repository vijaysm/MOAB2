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

#include <stdio.h>
#include <math.h>

#include "MBCore.hpp"
#include "MBRange.hpp"

#define filename "h5test.h5m"

// Dense tag name
// Pick an ugly name to check special-char escaping in HDF5 file
#define tagname "foo/\\/\\"
#define bitname "bar\n"


MBInterface* iface;

void create();

bool compare();

void moab_error( const char* function );

int main()
{
  MBErrorCode rval;
  std::string msg;
  
  iface = new MBCore();
  
  // create a dodecahedron and inscribed hex
  fprintf( stderr, "creating... " );
  create();
  
  // write out the dodecahedron
  fprintf( stderr, "writing... " );
  rval = iface->write_mesh( filename );
  if (MB_SUCCESS != rval)
  {
    fprintf( stderr, "Failed to write \"%s\"\n", filename );
    if (MB_SUCCESS == iface->get_last_error( msg ))
      fprintf( stderr, "%s\n", msg.c_str() );
    return 1;
  }
  
  // Read back in as a copy of the original
  fprintf( stderr, "reading... " );
  rval = iface->load_mesh( filename );
  if (MB_SUCCESS != rval)
  {
    fprintf( stderr, "Failed to read \"%s\"\n", filename );
    if (MB_SUCCESS == iface->get_last_error( msg ))
      fprintf( stderr, "%s\n", msg.c_str() );
    return 1;
  }
  
  // Compare the two.
  fprintf( stderr, "comparing... " );
  if (!compare())
  {
    fprintf(stderr, "Comparison failed.\n");
    return 1;
  }
  fprintf( stderr, "success!\n" );
  
  // Write both the original and copy to a file
  fprintf( stderr, "writing... " );
  rval = iface->write_mesh( filename );
  if (MB_SUCCESS != rval)
  {
    fprintf( stderr, "Failed to write \"%s\"\n", filename );
    if (MB_SUCCESS == iface->get_last_error( msg ))
      fprintf( stderr, "%s\n", msg.c_str() );
    return 1;
  }
 
  // Delete the mesh
  fprintf( stderr, "clearing db... " );
  rval = iface->delete_mesh();
  if (MB_SUCCESS != rval)
    moab_error( "delete_mesh" );
  
  // Read the two dodecahedrons from the file
  fprintf( stderr, "reading... " );
  rval = iface->load_mesh( filename );
  if (MB_SUCCESS != rval)
  {
    fprintf( stderr, "Failed to read \"%s\"\n", filename );
    if (MB_SUCCESS == iface->get_last_error( msg ))
      fprintf( stderr, "%s\n", msg.c_str() );
    return 1;
  }

  // Compare them
  fprintf( stderr, "comparing... " );
  if (!compare())
  {
    fprintf(stderr, "Comparison failed.\n");
    return 1;
  }
  fprintf( stderr, "success!\n" );

  // Delete the mesh
  fprintf( stderr, "cleaning up... " );
  rval = iface->delete_mesh();
  if (MB_SUCCESS != rval)
    moab_error( "delete_mesh" );

  // Clean up the file.
  remove( filename );
  fprintf( stderr, "done.\n" );
  return 0;
}

MBEntityHandle vtx( double x, double y, double z )
{
  const double p[3] = { x, y, z };
  MBEntityHandle result;
  if (MB_SUCCESS != iface->create_vertex( p, result ))
      moab_error( "create_vertex" );
  return result;
}

MBEntityHandle pent( MBEntityHandle* vtx_list, int i1, int i2, int i3, int i4, int i5 )
{
  const MBEntityHandle conn[5] = { vtx_list[i1],
                                   vtx_list[i2],
                                   vtx_list[i3],
                                   vtx_list[i4],
                                   vtx_list[i5] };
  MBEntityHandle result;
  if (MB_SUCCESS != iface->create_element( MBPOLYGON, conn, 5, result ))
    moab_error( "create_element" );
  return result;
}

void create()
{
  // Create dodecahedron
  
  // radius
  const double r = 50.;  
  // center
  const double x = 0., y = 0., z = 0.;
  // length of edge of inscribed cube
  const double cube = r * (2.0 / sqrt( 3.0 ));  
  // length of dodecahedron edge
  const double edge = cube * (2.0 / (1.0  + sqrt( 5.0 )));
  // distance of projection of a dodecahedron vertex to
  // closest edge of inscribed cube
  const double p = (cube - edge) / 2.0;
  // distance of projection of a dodecahedron vertex to
  // closest face of inscribed cube
  const double d = sqrt( edge*edge - cube*cube / 4 - p );
  // coordinate values
  const double c = cube / 2.0;
  const double a = c + d;
  const double b = edge / 2.0;
  
  // list of vertex handles
  MBEntityHandle vertices[20];
  // list of pentagon handles
  MBEntityHandle faces[12];
  // Dodecahedron handle
  MBEntityHandle dodec;
  // Inscribed Hex handle
  MBEntityHandle hex;
  
  // Create vertices if inscribed cube
  vertices[ 0] = vtx( x-c, y+c, z+c );
  vertices[ 1] = vtx( x+c, y+c, z+c );
  vertices[ 2] = vtx( x+c, y-c, z+c );
  vertices[ 3] = vtx( x-c, y-c, z+c );
  vertices[ 4] = vtx( x-c, y+c, z-c );
  vertices[ 5] = vtx( x+c, y+c, z-c );
  vertices[ 6] = vtx( x+c, y-c, z-c );
  vertices[ 7] = vtx( x-c, y-c, z-c );
  
  // Create inscribed hex
  if (MB_SUCCESS != iface->create_element( MBHEX, vertices, 8, hex ))
    moab_error( "create_element" );
  
  // Create vertices, 2 "above" each face of inscribed cube
  // +z face
  vertices[ 8] = vtx( x-b, y  , z+a );
  vertices[ 9] = vtx( x+b, y  , z+a );
  // +x face
  vertices[10] = vtx( x+a, y+b, z   );
  vertices[11] = vtx( x+a, y-b, z   );
  // -z face
  vertices[12] = vtx( x-b, y  , z-a );
  vertices[13] = vtx( x+b, y  , z-a );
  // -x face
  vertices[14] = vtx( x-a, y+b, z   );
  vertices[15] = vtx( x-a, y-b, z   );
  // +y face
  vertices[16] = vtx( x  , y+a, z+b );
  vertices[17] = vtx( x  , y+a, z-b );
  // -y face
  vertices[18] = vtx( x  , y-a, z+b );
  vertices[19] = vtx( x  , y-a, z-b );
  
  // Create petagons
  faces[ 0] = pent( vertices,  0,  8,  9,  1, 16 );
  faces[ 1] = pent( vertices,  3, 18,  2,  9,  8 );
  faces[ 2] = pent( vertices,  2, 11, 10,  1,  9 );
  faces[ 3] = pent( vertices,  2, 18, 19,  6, 11 );
  faces[ 4] = pent( vertices,  1, 10,  5, 17, 16 );
  faces[ 5] = pent( vertices,  5, 10, 11,  6, 13 );
  faces[ 6] = pent( vertices,  4, 17,  5, 13, 12 );
  faces[ 7] = pent( vertices,  7, 12, 13,  6, 19 );
  faces[ 8] = pent( vertices,  4, 12,  7, 15, 14 );
  faces[ 9] = pent( vertices,  0, 16, 17,  4, 14 );
  faces[10] = pent( vertices,  3, 15,  7, 19, 18 );
  faces[11] = pent( vertices,  0, 14, 15,  3,  8 );
 
  // Create dodecahedron
  if (MB_SUCCESS != iface->create_element( MBPOLYHEDRON, faces, 12, dodec ))
    moab_error( "create_element" );
  
  // Create a dense tag 
  int zero = 0;
  MBTag tag;
  if (MB_SUCCESS != iface->tag_create( tagname, 4, MB_TAG_DENSE, tag, &zero ))
    moab_error( "tag_create" );
  
  // Put dense tag on all vertices.
  for (int i = 0; i < 20; ++i)
    if (MB_SUCCESS != iface->tag_set_data( tag, vertices + i, 1, &i ))
      moab_error( "tag_set_data" );
      
  // Create bit tag
  MBTag tag2;
  if (MB_SUCCESS != iface->tag_create( bitname, 2, MB_TAG_BIT, tag2, 0 ))
    moab_error( "tag_create" );
  
  // Set tag to 0 on Hex
  char two = '\002';
  if (MB_SUCCESS != iface->tag_set_data( tag2, &hex, 1, &two ))
    moab_error( "tag_set_data" );
  
  // set tag to 1 on dodecahedron
  char one  = '\001';
  if (MB_SUCCESS != iface->tag_set_data( tag2, &dodec, 1, &one ))
    moab_error( "tag_set_data" );
}

bool compare_conn( std::vector<MBEntityHandle>& conn1,
                   std::vector<MBEntityHandle>& conn2 )
{
  int i;
  
  if (conn1.size() != conn2.size() || conn1.size() == 0)
  {
    fprintf(stderr, "Error comparing connectivity: sizes %d and %d\n", 
                     conn1.size(), conn2.size());
    return false;
  }
  
  std::vector<double> coords[2];
  coords[0].resize( 3*conn1.size() ); coords[1].resize( 3*conn2.size() );
  if (MB_SUCCESS != iface->get_coords( &conn1[0], conn1.size(), &coords[0][0] ) ||
      MB_SUCCESS != iface->get_coords( &conn2[0], conn2.size(), &coords[1][0] ) ||
      coords[0].size() != coords[1].size())
    moab_error( "get_coords" );

  std::vector<double>::iterator citer1 = coords[0].begin(), citer2 = coords[1].begin();
  for (i = 0; i < conn1.size(); i++)
  {
    double x1 = *(citer1++), y1 = *(citer1++), z1 = *(citer1++);
    double x2 = *(citer2++), y2 = *(citer2++), z2 = *(citer2++);
    if (x1 != x2 || y1 != y2 || z1 != z2)
    {
      fprintf( stderr, "Vertex coords don't match: ( %f, %f, %f ) and ( %f, %f, %f )\n",
        x1, y1, z1, x2, y2, z2 );
      return false;
    }
  }
 
  
  std::vector<int> tags[2];
  tags[0].resize( conn1.size() ); tags[1].resize( conn2.size() );
  MBTag tag;
  if (MB_SUCCESS != iface->tag_get_handle( tagname, tag ))
    moab_error( "tag_get_handle" );
  if (MB_SUCCESS != iface->tag_get_data( tag, &conn1[0], conn1.size(), &tags[0][0] ) ||
      MB_SUCCESS != iface->tag_get_data( tag, &conn2[0], conn2.size(), &tags[1][0] ))
    moab_error( "tag_get_data" );
  
  std::vector<int>::iterator titer1 = tags[0].begin(), titer2 = tags[1].begin();
  for (i = 0; i < conn1.size(); i++)
  {
    int t1 = *(titer1++);
    int t2 = *(titer2++);
    if (t1 != t2)
    {
      fprintf( stderr, "Vertex tags don't match: %d != %d\n", t1, t2 );
      return false;
    }
  }
  
  return true;
}


bool compare()
{
  MBRange range;
  MBEntityHandle hex[2];
  MBEntityHandle dod[2];
  MBTag elemtag;
  MBRange::iterator iter;
  
  // get tag
  if (MB_SUCCESS != iface->tag_get_handle( bitname, elemtag ))
    moab_error( "tag_get_handle" );
  
  // get two hexes
  char two = '\002';
  const void* tarray[] = { &two };
  if (MB_SUCCESS != iface->
      get_entities_by_type_and_tag( 0, MBHEX, &elemtag, tarray, 1, range ))
    moab_error( "get_entities_by_type_and_tag" );
  if (range.size() != 2)
  {
    fprintf( stderr, "Expected 2 Hexes.  Got %d\n", range.size() );
    exit( 1 );
  }
  iter = range.begin();
  hex[0] = *iter;
  hex[1] = *++iter;
  
  // get two polyhedra
  range.clear();
  char one = '\001';
  const void* oarray[] = { &one };
  if (MB_SUCCESS != iface->
      get_entities_by_type_and_tag( 0, MBPOLYHEDRON, &elemtag, oarray, 1, range ))
    moab_error( "get_entities_by_type_and_tag" );
  if (range.size() != 2)
  {
    fprintf( stderr, "Expected 2 Polyhedra.  Got %d\n", range.size() );
    exit( 1 );
  }
  iter = range.begin();
  dod[0] = *iter;
  dod[1] = *++iter;
  
  // compare hexes
  std::vector<MBEntityHandle> conn[2];
  if (MB_SUCCESS != iface->get_connectivity( hex  , 1, conn[0], true ) ||
      MB_SUCCESS != iface->get_connectivity( hex+1, 1, conn[1], true ))
    moab_error( "get_connectivity" );
  if (!compare_conn( conn[0], conn[1] ))
    return false;
  
  // compare polyhedra
  
  std::vector<MBEntityHandle> face[2];
  conn[0].clear(); conn[1].clear();
  if (MB_SUCCESS != iface->get_connectivity( dod  , 1, conn[0], true ) ||
      MB_SUCCESS != iface->get_connectivity( dod+1, 1, conn[1], true ))
    moab_error( "get_connectivity" );
  if (conn[0].size() != 12 || conn[1].size() != 12)
  {
    fprintf(stderr, "Expected two dodecahedrons.  Got polyhedrons with "
                    "%d and %d faces respectively.\n", 
                    conn[0].size(), conn[1].size() );
    return false;
  }
  
  for (int i = 0; i < 12; ++i )
  {
    face[0].clear(); face[1].clear();
    if (MB_SUCCESS != iface->get_connectivity( &conn[0][i], 1, face[0], true) ||
        MB_SUCCESS != iface->get_connectivity( &conn[1][i], 1, face[1], true))
      moab_error( "get_connectivity" );
    if (!compare_conn( face[0], face[1] ))
      return false;
  }
  
  return true;
}

void moab_error( const char* str )
{
  std::string msg;
  fprintf( stderr, "%s() failed.\n", str );
  if (MB_SUCCESS == iface->get_last_error( msg ))
    fprintf( stderr, "%s\n", msg.c_str() );
  exit( 1 );
}
  
