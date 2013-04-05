/* \example HelloMOAB HelloMOAB.cpp
 * \Description: read a mesh, get the entities.
 * Prerequisite examples: none
 * better Doxygen-ized, standardized comment section
 *
 * general description: This is a simple file is used to read meshes from VTK file and test how many entities there are.
 *
 * To run: ./HelloMOAB <mesh file>
 * (default values can run if users don't specify a mesh file)
 */


#include "moab/Core.hpp"

using namespace moab;
using namespace std;

string test_file_name = string(MESH_DIR) + string("/3k-tri-sphere.vtk");

int main( int argc, char** argv[] )
{
  Interface *iface = new Core;

    // need option handling here for input filename
  if (argc > 1){
    //user has input a mesh file
    test_file_name = argv[1];
  }  
    //load the mesh from vtk file
  ErrorCode rval = iface->load_mesh( test_file_name );
  assert( rval == MB_SUCCESS);

    //get verts entities
  Range verts;
  rval = iface->get_entities_by_type(0, MBVERTEX, verts);
  assert( rval == MB_SUCCESS);

    //get triangular entities
  Range faces;
  rval = iface->get_entities_by_type(0, MBTRI, faces);
  assert( rval == MB_SUCCESS);

  cout << "Number of vertices is " << verts.size() << " and faces is " << faces.size() << endl;
  
  return 0;
}
