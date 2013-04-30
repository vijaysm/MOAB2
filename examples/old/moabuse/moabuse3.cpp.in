/// \file moabuse3.cpp
/// 
/// \author Milad Fatenejad
/// 
/// \brief moabuse tutorial, example 3: Demonstrates
///        constructing/saving a simple 2x2x2 hex mesh using the
///        structured mesh interface
///
/// In this example, we create a 2x2x2 mesh that is identical to the
/// previous example. Howver, in this case we will use the structured
/// mesh interface since the mesh we created is logically
/// structured. There are many advantages to using the structured mesh
/// interface...such as memory savings, speed, ease-of-use...
///
/// In the previous example, we had to create 27 vertexes manually,
/// define the connectivity, then manually create 8 hexahedrons. With
/// the structured mesh interface, we just have to create the 27
/// vertexes then tell MOAB that these define a 2x2x2 structured mesh
/// and everything else is taken care of for us!

// The moab/Core.hpp header file is needed for all MOAB work...
#include "moab/Core.hpp"

// The moab/ScdInterface.hpp contains code which defines the moab
// structured mesh interface
#include "moab/ScdInterface.hpp"

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
  moab::Core mbcore;
  moab::Interface& mbint = mbcore;

  // **********************************
  // *   Create the Structured Mesh   *
  // **********************************

  // As before, we have to create an array defining the coordinate of
  // each vertex:
  const unsigned NUMVTX = 27;
  const double vertex_coords[3*NUMVTX] = { 0, 0, 0,  1, 0, 0,  2, 0, 0,
					   0, 1, 0,  1, 1, 0,  2, 1, 0,
					   0, 2, 0,  1, 2, 0,  2, 2, 0,
					   
					   0, 0, 1,  1, 0, 1,  2, 0, 1,
					   0, 1, 1,  1, 1, 1,  2, 1, 1,
					   0, 2, 1,  1, 2, 1,  2, 2, 1,
					   
					   0, 0, 2,  1, 0, 2,  2, 0, 2,
					   0, 1, 2,  1, 1, 2,  2, 1, 2,
					   0, 2, 2,  1, 2, 2,  2, 2, 2 };

  // moab::ScdInterface is the structured mesh interface class for
  // MOAB.
  moab::ScdInterface *scdint;

  // The query_interface method will create a structured mesh instance
  // for mbcore and will point scdint to it. This is how you tell moab
  // that our moab::Core instance is going to represent a structured
  // mesh.
  moab::ErrorCode rval = mbint.query_interface(scdint);
  MBERR("mbint.query_interface", rval);

  // Structured meshes a divided into "boxes". Each box represents a
  // little structured mesh. A single mesh, for example a
  // block-structured mesh, can contain many individual boxes. In this
  // example, we want to create a single, 2x2x2 box. The construct_box
  // method will do this for us.
  moab::ScdBox *scdbox = NULL;
  rval = scdint->construct_box(moab::HomCoord(0,0,0), 
			       moab::HomCoord(2,2,2), 
			       vertex_coords, 
			       NUMVTX, 
			       scdbox);
  MBERR("scdint->construct_box", rval);

  // moab::HomCoord is a little class that is used to represent a
  // coordinate in logical space. Above, we told MOAB that we want
  // indexes which extend from 0,0,0 to 2,2,2 for vertexes. Since the
  // i,j,k start/end indexes are all different, MOAB knows that our
  // mesh consists of hexes. If we had gone from 0,0,0 to 1,1,0 then
  // MOAB would construct a 2x2x1 mesh of quadralaterals.


  // ***************
  // *   Inspect   *
  // ***************

  // Our mesh now exists and contains a single box! Lets loop over
  // that box and manually print out the handle of each vertex/hex:
  unsigned i,j,k;

  // Print out the entity handle associated with each hex:
  for(i = 0; i < 2; ++i)
    for(j = 0; j < 2; ++j)
      for(k = 0; k < 2; ++k) {
	std::cout << "Hex (" << i << "," << j << "," << k << ") "
		  << "has handle: " << scdbox->get_element(i,j,k) << std::endl;
      }

  // Print out the entity handle associated with each vertex:
  for(i = 0; i < 3; ++i)
    for(j = 0; j < 3; ++j)
      for(k = 0; k < 3; ++k) {
	std::cout << "Vertex (" << i << "," << j << "," << k << ") "
		  << "has handle: " << scdbox->get_vertex(i,j,k) << std::endl;
      }

  // Of course, you can still use all of the functionality defined in
  // mbint to get coordinates, connectivity, etc...
  
  // ***************************
  // *   Write Mesh to Files   *
  // ***************************

  rval = mbint.write_file("moabuse3.vtk");
  MBERR("write_file(moabuse3.vtk)", rval);

  return 0;
}
