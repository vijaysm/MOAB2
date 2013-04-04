///Description: read a mesh, get the entities.
///Prerequisite examples: none

//general description: This is a simple file is used to read meshes from VTK file and test how many entities there are.
// Code

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/Interface.hpp"




#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <cassert>

using namespace moab;
using namespace std;

const char* test_file_name = "../MeshFiles/unittest/3k-tri-sphere.vtk";
const int fixed_num = 4457;
int main( int, char*  )
{

 Interface *iface = new Core;

 //load the mesh from vtk file
 ErrorCode rval = iface->load_mesh( test_file_name );
 assert( rval == MB_SUCCESS);

 //get the root set
 EntityHandle root_set = iface->get_root_set();

 //get node entities
 std::vector<EntityHandle> nodes;
 rval = iface->get_entities_by_type(root_set, MBVERTEX, nodes);
 assert( nodes.size() == 1487); 
 //get triangular entities
 std::vector<EntityHandle> faces;
 rval = iface->get_entities_by_type(root_set, MBTRI, faces);
 assert( nodes.size() == 2970); 

 return 0;
}


