/** @example StructuredMeshSimple.cpp
 * \brief Show creation and query of structured mesh, serial or parallel, through MOAB's structured mesh interface.
 * This is an example showing creation and query of a 3D structured mesh.  In serial, a single N*N*N block of elements
 * is created; in parallel, each proc gets an N*N*N block, with blocks arranged in a 1d column, sharing vertices
 * and faces at their interfaces (proc 0 has no left neighbor and proc P-1 no right neighbor).
 * Each square block of hex elements is then referenced by its ijk parameterization.
 * 1D and 2D examples could be made simply by changing the dimension parameter passed into the MOAB functions. \n
 *
 * <b>This example </b>:
 *    -# Instantiate MOAB and get the structured mesh interface
 *    -# Decide what the local parameters of the mesh will be, based on parallel/serial and rank.
 *    -# Create a N^d structured mesh, which includes (N+1)^d vertices and N^d elements.
 *    -# Get the vertices and elements from moab and check their numbers against (N+1)^d and N^d, resp.
 *    -# Loop over elements in d nested loops over i, j, k; for each (i,j,k):
 *      -# Get the element corresponding to (i,j,k)
 *      -# Get the connectivity of the element
 *      -# Get the coordinates of the vertices comprising that element
 *    -# Release the structured mesh interface and destroy the MOAB instance
 *
 * <b> To run: </b> ./structuredmesh [d [N] ] \n
 * (default values so can run w/ no user interaction)
 */

#include "moab/Core.hpp"
#include "moab/ScdInterface.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CN.hpp"
#ifdef USE_MPI
#include "moab_mpi.h"
#endif
#include <iostream>
#include <vector>

using namespace moab;

int main(int argc, char **argv) 
{
  int N = 10, dim = 3;

#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif

  ProgOptions opts;
  opts.addOpt<int>(std::string("dim,d"), std::string("Dimension of mesh (default=3)"),
                   &dim);
  opts.addOpt<int>(std::string(",n"), std::string("Number of elements on a side (default=10)"),
                   &N);
  opts.parseCommandLine(argc, argv);
  
    // 0. Instantiate MOAB and get the structured mesh interface
  Interface *mb = new Core();
  ScdInterface *scdiface;
  ErrorCode rval = mb->query_interface(scdiface); // get a ScdInterface object through moab instance
  if (MB_SUCCESS != rval) return rval;

    // 1. Decide what the local parameters of the mesh will be, based on parallel/serial and rank.
  int ilow = 0, ihigh = N;
  int rank = 0, nprocs = 1;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs); MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ilow = rank*N; ihigh = ilow + N;
#endif  

    // 2. Create a N^d structured mesh, which includes (N+1)^d vertices and N^d elements.
  ScdBox *box;
  rval = scdiface->construct_box(HomCoord(ilow, (dim>1?0:-1), (dim>2?0:-1)), // use in-line logical tests to handle dimensionality
                                 HomCoord(ihigh, (dim>1?N:-1), (dim>2?N:-1)), 
                                 NULL, 0, // NULL coords vector and 0 coords (don't specify coords for now)
                                 box);    // box is the structured box object providing the parametric
                                          // structured mesh interface for this rectangle of elements
  if (MB_SUCCESS != rval) return rval;

    // 3. Get the vertices and elements from moab and check their numbers against (N+1)^d and N^d, resp.
  Range verts, elems;
  rval = mb->get_entities_by_dimension(0, 0, verts); // first '0' specifies "root set", or entire MOAB instance, second the entity dimension being requested
  if (MB_SUCCESS != rval) return rval;
  rval = mb->get_entities_by_dimension(0, dim, elems);
  if (MB_SUCCESS != rval) return rval;

#define MYSTREAM(a) if (!rank) std::cout << a << std::endl
  
  if (pow(N,dim) == (int) elems.size() && pow(N+1,dim) == (int) verts.size()) { // expected #e and #v are N^d and (N+1)^d, resp.
#ifdef USE_MPI
    MYSTREAM("Proc 0: ");
#endif
    MYSTREAM("Created " << elems.size() << " " << CN::EntityTypeName(mb->type_from_handle(*elems.begin())) 
             << " elements and " << verts.size() << " vertices." << std::endl);
  }
  else
    std::cout << "Created the wrong number of vertices or hexes!" << std::endl;
  
    // 4. Loop over elements in 3 nested loops over i, j, k; for each (i,j,k):
  std::vector<double> coords(3*pow(N+1,dim));
  std::vector<EntityHandle> connect;
  for (int k = 0; k < (dim>2?N:1); k++) {
    for (int j = 0; j < (dim>1?N:1); j++) {
      for (int i = 0; i < N-1; i++) {
          // 4a. Get the element corresponding to (i,j,k)
        EntityHandle ehandle = box->get_element(i, j, k);
        if (0 == ehandle) return MB_FAILURE;
          // 4b. Get the connectivity of the element
        rval = mb->get_connectivity(&ehandle, 1, connect); // get the connectivity, in canonical order
        if (MB_SUCCESS != rval) return rval;
          // 4c. Get the coordinates of the vertices comprising that element
        rval = mb->get_coords(connect.data(), connect.size(), coords.data()); // get the coordinates of those vertices
        if (MB_SUCCESS != rval) return rval;
      }
    }
  }

    // 5. Release the structured mesh interface and destroy the MOAB instance
  mb->release_interface(scdiface); // tell MOAB we're done with the ScdInterface
  delete mb;

#ifdef USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
