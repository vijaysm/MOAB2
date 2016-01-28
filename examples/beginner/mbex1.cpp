/// \file moabuse1.cpp
/// 
/// \author Milad Fatenejad
/// 
/// \brief moabuse tutorial, example 1: Demonstrates
///        constructing/saving a simple 2x2x2 hex mesh
///
/// This example creates a 2x2x2 mesh (uniform grid) and writes it out
/// to a VTK file and an H5M file. Each cell is of size 1x1x1 and is
/// axis aligned. This example demonstrates how to manually create a
/// mesh using the unstructured mesh interface. The mesh is then
/// written out to file.

// The moab/Core.hpp header file is needed for all MOAB work...
#include "moab/Core.hpp"

#include <iostream>

// "mberr.hpp" contains the MBERR macro and is local to this tutorial
#include "mberr.hpp"

// ****************
// *              *
// *     main     *
// *              *
// ****************
int main()
{

  // MOAB functionality is accessed through an instance of the
  // moab::Interface class:
  moab::Core mbcore;

  // ***************************
  // *   Create the vertexes   *
  // ***************************
  
  // We are going to create 27 vertexes which will be used to define 8
  // hexahedron cells. First, we have to create an array to store the
  // coordinates of each vertex.

  const unsigned NUMVTX = 27; // The number of vertexes
  const unsigned NUMHEX = 8;  // The number of hexahedrons

  // This double array stores the x, y, z coordinate of each vertex.
  const double vertex_coords[3*NUMVTX] = { 0, 0, 0,  1, 0, 0,  2, 0, 0,
					   0, 1, 0,  1, 1, 0,  2, 1, 0,
					   0, 2, 0,  1, 2, 0,  2, 2, 0,
					   
					   0, 0, 1,  1, 0, 1,  2, 0, 1,
					   0, 1, 1,  1, 1, 1,  2, 1, 1,
					   0, 2, 1,  1, 2, 1,  2, 2, 1,
					   
					   0, 0, 2,  1, 0, 2,  2, 0, 2,
					   0, 1, 2,  1, 1, 2,  2, 1, 2,
					   0, 2, 2,  1, 2, 2,  2, 2, 2 };

  // Create the vertexes and store their entity handles in the
  // vertex_handles range. In MOAB, entities are defined using Entity
  // Handles (type moab::EntityHandle). An entity handle is a unique
  // integer that is used to identify/refer to specific entities (such
  // as vertexes, edges, hexahedrons, etc...) on the mesh. MOAB
  // guarentees that when multiple entities are created at the same
  // time, as in the create_vertices call below, that those entities
  // get adjacent handles. Thus, the second of the 27 vertexes we are
  // creating will have an entity handle that is 1 greater than the
  // handle of the first vertex. A range (type moab::Range) is a
  // container which stores entity handles. Of course, sets of handles
  // can also be stored in vectors or arrays, but ranges are much more
  // memory efficient, so use them when possible!
  moab::Range vertex_handles;
  moab::ErrorCode rval = mbcore.create_vertices( vertex_coords, NUMVTX, vertex_handles );
  MBERR("create_vertices", rval);

  // You can print out a range to see what elements it contains:
  std::cout << "Just created the following entities:" 
	    << vertex_handles << std::endl;
  
  // ******************************
  // *   Create the Hexahedrons   *
  // ******************************

  // The conn array stores the connectivity for each hex. It defines
  // which 8 vertexes connect to define a single hexahedron. The loop
  // below adds the handle of the first vertex (stored in
  // first_vertex_handle) to conn thereby ensuring that the entries of
  // conn are actual entity handles. This only works because MOAB
  // guarentees that entities (such as our vertexes) created at once
  // have adjacent entity handles.
  moab::EntityHandle conn[NUMHEX][8] = { {  0,  1,  4,  3,   9, 10, 13, 12 },
					 {  1,  2,  5,  4,  10, 11, 14, 13 }, 
					 {  3,  4,  7,  6,  12, 13, 16, 15 },
					 {  4,  5,  8,  7,  13, 14, 17, 16 },
					 {  9, 10, 13, 12,  18, 19, 22, 21 },
					 { 10, 11, 14, 13,  19, 20, 23, 22 },
					 { 12, 13, 16, 15,  21, 22, 25, 24 },
					 { 13, 14, 17, 16,  22, 23, 26, 25 } };

  // Lets get the handle for the first vertex. Note that we can use
  // the square brackets operator on ranges just like vectors or
  // arrays:
  moab::EntityHandle first_vertex_handle = vertex_handles[0];

  for (unsigned i = 0; i < NUMHEX; ++i) {
    for (unsigned j = 0; j < 8; ++j) {
      // Add first_vertex_handle to each elment of conn. This ensures
      // that the handles are specified properly (i.e. when
      // first_vertex_handle > 0)
      conn[i][j] = conn[i][j] + first_vertex_handle;
    }
  }


  // Now that the connectivity of each hex has been defined, we can
  // create each hex using a call to the create_element method which
  // gives back an entity handle for the hex that was created. We'll
  // then insert that entity into a range:
  moab::Range hexahedron_handles;
  moab::EntityHandle element;
  for (unsigned i = 0; i < NUMHEX; ++i) {
    rval = mbcore.create_element( moab::MBHEX, conn[i], 8, element );
    MBERR("create_element", rval);

    hexahedron_handles.insert(element);
  }

  // Let's see what entities we just created:
  std::cout << "Just created the following entities: "
	    << hexahedron_handles << std::endl;

  // ***************************
  // *   Write Mesh to Files   *
  // ***************************

  // Now that we've created this amazing mesh, we can write it out to
  // a file. Since MOAB can write out to a variety of standard file
  // formats, you can quickly visualize and manipulate your mesh using
  // standard tools, such as VisIt.

  // In these examples, I will stick to using the VTK file format
  // because it is text based and will work whether or not you've got
  // HDF5, NETCDF, etc... installed and is a fairly standard file
  // format so a lot of tools work with it out of the box. 
  rval = mbcore.write_file("moabuse1.vtk");
  MBERR("write_file(moabuse1.vtk)", rval);

  return 0;
}
