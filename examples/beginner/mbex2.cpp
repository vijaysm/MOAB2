/// \file moabuse2.cpp
/// 
/// \author Milad Fatenejad
/// 
/// \brief moabuse tutorial, example 2: Demonstrates loading a mesh
///        from a file, finding coordinate locations, connectivity
///        information...
///
/// In this example, we read in the VTK file (moabuse1.vtk) generated
/// in example 1, pick one of the hexahedrons, find out the vertexes
/// that define it (connectivity), and find the coordinates of those
/// vertexes. Finally, we will rotate the mesh a little and write it
/// out to a new file.

// The moab/Core.hpp header file is needed for all MOAB work...
#include "moab/Core.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>

// "mberr.hpp" contains the MBERR macro and is local to this tutorial
#include "mberr.hpp"

// ****************
// *              *
// *     main     *
// *              *
// ****************
int main()
{

  moab::Core mbcore;
  moab::Interface& mbint = mbcore;

  // First, lets load the file from the previous example. Note that
  // you must compile and run the previous example to get moabuse1.vtk!
  moab::ErrorCode rval = mbint.load_file("moabuse1.vtk");
  MBERR("load_file", rval);

  // *********************
  // *   Introspection   *
  // *********************

  // We can now access everything about the mesh through mbint. Lets
  // pick an element and find some information about it. In this case,
  // we know that all of the elements in the mesh are hexahedrons, so
  // we will ask MOAB for a range containing the handle for all hexes:
  moab::Range hex_range;
  rval = mbint.get_entities_by_type(0, moab::MBHEX, hex_range);
  MBERR("get_entities_by_type(HEX)", rval);

  std::cout << "Just loaded a mesh containing: " << hex_range << std::endl;

  // Let's analyze one of the hexes (in this case, I picked hex three):
  moab::EntityHandle handle = hex_range[3];

  // Find out the eight vertexes that define this particular hex:
  moab::Range connectivity;
  mbint.get_connectivity(&handle, 1, connectivity);

  // Connectivity should now contain the handles for the eight
  // vertexes of interest. Lets print the handle and coordinate of
  // each vertex. Note how we iterate through the range just like with
  // STL containers...
  double coord[3];
  moab::Range::iterator iter;

  std::cout << std::setw(6)  << "Handle" 
	    << std::setw(10) << "X" 
	    << std::setw(10) << "Y"
	    << std::setw(10) << "Z" << std::endl;

  for(iter = connectivity.begin(); iter != connectivity.end(); ++iter) {
    rval = mbint.get_coords(&(*iter), 1, coord);
    MBERR("get_coords", rval);
    
    // Print the entity handle followed by the x, y, z coordinate of
    // the vertex:
    std::cout << std::setw(6) << *iter
	      << std::setw(10) << coord[0]
	      << std::setw(10) << coord[1]
	      << std::setw(10) << coord[2] << std::endl;
  }

  // ***********************
  // *   Rotate the Mesh   *
  // ***********************

  // Now, let's rotate the mesh about the z-axis by pi/4 radians. Do
  // to this, we will iterate through all of the vertexes

  // Start by getting a range containing all of the vertex handles:
  moab::Range vertex_range;
  rval = mbint.get_entities_by_type(0, moab::MBVERTEX, vertex_range);
  MBERR("get_entities_by_type(VERTEX)", rval);

  // Get the coordinates of all of the vertexes:
  std::vector<double> vertex_coords(3*vertex_range.size());
  rval = mbint.get_coords(vertex_range, vertex_coords.data());
  MBERR("get_coords", rval);

  unsigned count = 0;
  const double PI = 3.14159265359;
  const double ANGLE = PI/4;
  for(iter = vertex_range.begin(); iter != vertex_range.end(); ++iter) {

    // Save the old coordinates:
    double x = vertex_coords[count+0];
    double y = vertex_coords[count+1];
    double z = vertex_coords[count+2];

    // Apply the rotation:
    vertex_coords[count+0] = x * std::cos(ANGLE) - y * std::sin(ANGLE);
    vertex_coords[count+1] = x * std::sin(ANGLE) + y * std::cos(ANGLE);

    count += 3;
  }

  // Now, the vertex_coords vector contains all of the updated
  // coordinates. Let's push these back into MOAB:
  mbint.set_coords(vertex_range, vertex_coords.data());
  MBERR("set_coords", rval);

  // **************************
  // *   Write Mesh to File   *
  // **************************

  rval = mbint.write_file("moabuse2.vtk");
  MBERR("write_file(moabuse2.vtk)", rval);

  // Now you can open your favorite visualization tool and look at the
  // original (moabuse1.vtk) and rotated (moabuse2.vtk) meshes.

  return 0;
}
