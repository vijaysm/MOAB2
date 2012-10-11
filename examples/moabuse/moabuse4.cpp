/// \file moabuse4.cpp
/// 
/// \author Milad Fatenejad
/// 
/// \brief moabuse tutorial, example 3: Create a 2D structured mesh
///        and set some tag data on it
///
/// In this example, we create a 2D structured mesh (actually a 3D
/// mesh made of quads) and then we will actually set data on the
/// mesh. Tags represent data that is attached to entities. In this
/// example, we will create two tags:
///
/// -# A "temperature" tag that is a single double precision number
///    attached to each quad.
/// -# A "velocity" tag that is an array of 2 double precision numbers
///    attached to each vertex.
///
/// We will write the mesh out to a file, and you can visualize the
/// data using your favorite tool.
///
/// In this example, I am demonstrating these operations in the
/// clearest possible way - not using the most efficient method. MOAB
/// has been designed so that users do not need to sacrifice
/// performance - future examples will demonstrate how to
/// access/manipulate the mesh using the fastest methods possible.

// The moab/Core.hpp header file is needed for all MOAB work...
#include "moab/Core.hpp"

// The moab/ScdInterface.hpp contains code which defines the moab
// structured mesh interface
#include "moab/ScdInterface.hpp"

#include <iostream>
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

  // ***********************
  // *   Create the Mesh   *
  // ***********************

  // First, lets make the mesh. It will be a 100 by 100 uniform grid
  // (there will be 100x100 quads, 101x101 vertexes) with dx = dy =
  // 0.1. Unlike the previous example, we will first make the mesh,
  // then set the coordinates one at a time.

  const unsigned NI = 100;
  const unsigned NJ = 100;

  // moab::ScdInterface is the structured mesh interface class for
  // MOAB.
  moab::ScdInterface *scdint;

  // Tell MOAB that our mesh is structured:
  moab::ErrorCode rval = mbint.query_interface(scdint);
  MBERR("mbint.query_interface", rval);

  // Create the mesh:
  moab::ScdBox *scdbox = NULL;
  rval = scdint->construct_box(moab::HomCoord(0,0,0), 
			       moab::HomCoord(NI,NJ,0),
			       NULL, 
			       0, 
			       scdbox);
  MBERR("scdint->construct_box", rval);

  // MOAB knows to make quads instead of hexes because the last start
  // and end indexes are the same (0). Note that it is still a "3D"
  // mesh because each vertex coodinate is still defined using three
  // numbers althrough every element in the mesh is a quadralateral.

  // ******************************
  // *   Set Vertex Coordinates   *
  // ******************************

  // The "NULL" and "0" arguments in the call to construct_box are
  // where we could specify the vertex coordinates. Since we didn't
  // give any coordinates, every vertex is given a position of 0,0,0
  // by default. Now we will set the vertex coordinates...

  const double DX = 0.1;
  const double DY = 0.1;

  for(unsigned i = 0; i < NI+1; i++) 
    for(unsigned j = 0; j < NJ+1; j++) {
      // First, get the entity handle:
      moab::EntityHandle handle = scdbox->get_vertex(i,j);

      // Compute the coordinate:
      double coord[3] = {DX*i, DY*j, 0.0};

      // Change the coordinate of the vertex:
      mbint.set_coords(&handle, 1, coord);
    }


  // *******************
  // *   Attach Tags   *
  // *******************
  
  // The vertex coordinates have been defined, now let's attach some
  // data to the mesh. In MOAB this is done using "tags". Tags are
  // little bits of information that can be attached to any mesh
  // entity. In our example, we want to create two tags. The
  // "temperature" tag will be attached to each quad and will be 1
  // double. The "velocity" tag will be attached to each vertex and
  // will be an array of two doubles.

  // zero and twozeros represent the initial tag

  // Create the tags:
  moab::Tag temp_tag;
  double temp_default_value = 0.0;
  rval = mbint.tag_get_handle("temperature", 1, moab::MB_TYPE_DOUBLE, temp_tag, 
                              moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, 
			      &temp_default_value);
  MBERR("mbint.tag_get_handle(temperature)",rval);

  moab::Tag vel_tag;
  double vel_default_value[2] = {0.0,0.0};
  rval = mbint.tag_get_handle("velocity", 2, moab::MB_TYPE_DOUBLE, vel_tag, 
                              moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, 
			      vel_default_value);
  MBERR("mbint.tag_get_handle(velocity)",rval);

  // Note that when we created each tag, we specified two flags:
  //
  // The moab::MB_TAG_DENSE flag tells MOAB that this is a dense
  // tag. Dense tags will get automatically assigned to entities which
  // have continuous handles. For this example, this means that once
  // we set a tag on one vertex, memory will be allocated for
  // assigning the tag to all vertexes. The same is true of
  // quads. Dense tags are much more efficient when assigning tags to
  // lots of entities. If you only want to assign a tag to a few
  // entities, it is more efficient to use sparse tags
  // (moab::MB_TAG_SPARSE).
  // 
  // The moab::MB_TAG_CREAT flag tells MOAB to create the tag if it
  // doesn't already exist.

  // The tags have now been created, now we have to attach them to
  // entities and set their values. NOTE: I am going to do this in a
  // manner which emphasizes clarity - this is not the most efficient
  // approach - that will be saved for later tutorial examples.

  // Loop through each quad and set the temperature:
  for (unsigned i = 0; i < NI; i++)
    for (unsigned j = 0; j < NJ; j++) {
      // Get the handle for this quad:
      moab::EntityHandle handle = scdbox->get_element(i,j);
      
      // Compute the temperature...
      double xc = DX*(i+0.5);
      double yc = DY*(j+0.5);
      double r = std::sqrt(xc*xc + yc*yc);
      double temperature = std::exp(-0.5*r);

      // Set the temperature on a single quad:
      rval = mbint.tag_set_data(temp_tag, 
				&handle,
				1,
				&temperature);
      MBERR("mbint.tag_set_data(temp_tag)", rval);
    }


  // Loop through each vertex and set the velocity:
  for (unsigned i = 0; i < NI+1; i++)
    for (unsigned j = 0; j < NJ+1; j++) {
      // Get the handle for this vertex:
      moab::EntityHandle handle = scdbox->get_vertex(i,j);
      double velocity[2] = {i, j};

      // Set the velocity on a vertex:
      rval = mbint.tag_set_data(vel_tag, 
				&handle,
				1,
				velocity);
      MBERR("mbint.tag_set_data(vel_tag)", rval);
    }

  // ***************************
  // *   Write Mesh to Files   *
  // ***************************

  // NOTE: Some visualization software (such as VisIt) may not
  // interpret the velocity tag as a vector and you may not be able to
  // plot it. But you should be able to plot the temperature on top of
  // the mesh.

  rval = mbint.write_file("moabuse4.vtk");
  MBERR("write_file(moabuse4.vtk)", rval);


  return 0;
}
