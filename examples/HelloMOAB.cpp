/** @example HelloMOAB.cpp
 * Description: read a mesh, get the entities.\n
 * HelloMOAB is a simple test file which is used to read meshes from VTK file and test how many entities there are.\n
 *
 * To run: ./HelloMOAB [meshfile]\n
 * (default values can run if users don't specify a mesh file)
 */


#include "moab/Core.hpp"
#include <iostream>
#include <assert.h>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

string test_file_name = string(MESH_DIR) + string("/3k-tri-sphere.vtk");

int main( int argc, char** argv )
{
  Interface *iface = new Core;

    // need option handling here for input filename
  if (argc > 1){
    //user has input a mesh file
    test_file_name = argv[1];
  }  
    //load the mesh from vtk file
  ErrorCode rval = iface->load_mesh( test_file_name.c_str() );
  assert(rval == MB_SUCCESS);

    // get verts entities, by type
  Range verts;
  rval = iface->get_entities_by_type(0, MBVERTEX, verts);
  assert(rval == MB_SUCCESS);
    //get edge entities, by type
  Range edges;
  rval = iface->get_entities_by_type(0, MBEDGE, edges);
  assert(rval == MB_SUCCESS);

    // get faces, by dimension, so we stay generic to entity type
  Range faces;
  rval = iface->get_entities_by_dimension(0, 2, faces);
  assert(rval == MB_SUCCESS);

    //get regions, by dimension, so we stay generic to entity type
  Range elems;
  rval = iface->get_entities_by_dimension(0, 3, elems);
  assert(rval == MB_SUCCESS);

   //output the number of entities
  cout << "Number of vertices is " << verts.size() <<  endl;
  cout << "Number of edges is " << edges.size() <<  endl;
  cout << "Number of faces is " << faces.size() <<  endl;
  cout << "Number of elements is " << elems.size() <<  endl;
  
  return 0;
}
