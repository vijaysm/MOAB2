/* \example structuredmesh structuredmesh.cpp
 * \brief Show creation and query of structured mesh through MOAB's structured mesh interface.
 * This is a serial example showing creation and query of a 3D structured mesh.  A single
 * rectangular brick of hex elements is created, then referenced by its ijk parameterization.
 * 1D and 2D examples could be made simply by changing the dimension parameter and the number
 * of ijk parameters passed into the MOAB functions.
 *
 * This example:
 * 0. Instantiate MOAB and get the structured mesh interface
 * 1. Creates a IxJxK structured mesh, which includes I*J*K vertices and (I-1)*(J-1)*(K-1) hexes.
 * 2. Get the vertices and hexes from moab and check their numbers against I*J*K and (I-1)*(J-1)*(K-1), resp.
 * 3. Loop over elements in 3 nested loops over i, j, k; for each (i,j,k):
 * 3a. Get the element corresponding to (i,j,k)
 * 3b. Get the connectivity of the element
 * 3c. Get the coordinates of the vertices comprising that element
 * 4. Release the structured mesh interface and destroy the MOAB instance
 *
 * To run: ./structuredmesh <dim> <N>
 * (default values so can run w/ no user interaction)
 */

#include "moab/Core.hpp"
#include "moab/ScdInterface.hpp"
#include <iostream>
#include <vector>

using namespace moab;

int main(int argv, char **argv) 
{
  int N;
    // progoptions?
  std::cout << "Enter I, J, K... " << std::endl;
  std::cin >> I >> J >> K;
  
    // 0. Instantiate MOAB and get the structured mesh interface
  Interface *mb = new Core();
  ScdInterface *scdiface;
  ErrorCode rval = mb->query_interface(scdiface); // get a ScdInterface object through moab instance
  if (MB_SUCCESS != rval) return rval;
  
    // 1. Creates a IxJxK structured mesh, which includes I*J*K vertices and (I-1)*(J-1)*(K-1) hexes.
  ScdBox *box;
  rval = scdiface->construct_box(HomCoord(0, 0, 0), HomCoord(I-1, J-1, K-1), // low, high box corners in parametric space
                                 NULL, 0, // NULL coords vector and 0 coords (don't specify coords for now)
                                 box);    // box is the structured box object providing the parametric
                                          // structured mesh interface for this rectangle of elements
  if (MB_SUCCESS != rval) return rval;

    // 2. Get the vertices and hexes from moab and check their numbers against I*J*K and (I-1)*(J-1)*(K-1), resp.
  Range verts, hexes;
  rval = mb->get_entities_by_dimension(0, 0, verts); // first '0' specifies "root set", or entire MOAB instance, second the entity dimension being requested
  if (MB_SUCCESS != rval) return rval;
  rval = mb->get_entities_by_dimension(0, 3, hexes);
  if (MB_SUCCESS != rval) return rval;

  if ((I-1)*(J-1)*(K-1) == (int) hexes.size() && I*J*K == (int) verts.size())
    std::cout << "Created " << hexes.size() << " hexes and " << verts.size() << " vertices." << std::endl;
  else
    std::cout << "Created the wrong number of vertices or hexes!" << std::endl;
  
    // 3. Loop over elements in 3 nested loops over i, j, k; for each (i,j,k):
  std::vector<double> coords(8*3);
  std::vector<EntityHandle> connect;
  for (int k = 0; k < K-1; k++) {
    for (int j = 0; j < J-1; j++) {
      for (int i = 0; i < I-1; i++) {
          // 3a. Get the element corresponding to (i,j,k)
        EntityHandle ehandle = box->get_element(i, j, k);
        if (0 == ehandle) return MB_FAILURE;
          // 3b. Get the connectivity of the element
        rval = mb->get_connectivity(&ehandle, 1, connect); // get the connectivity, in canonical order
        if (MB_SUCCESS != rval) return rval;
          // 3c. Get the coordinates of the vertices comprising that element
        rval = mb->get_coords(connect.data(), connect.size(), coords.data()); // get the coordinates of those vertices
        if (MB_SUCCESS != rval) return rval;
      }
    }
  }

    // 4. Release the structured mesh interface and destroy the MOAB instance
  mb->release_interface(scdiface); // tell MOAB we're done with the ScdInterface
  delete mb;
  
  return 0;
}
