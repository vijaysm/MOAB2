/*This function tests the AHF datastructures on C*/
#include <iostream>
#include <vector>
#include <algorithm>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/HalfFacetRep.hpp"
#include "TestUtil.hpp"

#ifdef USE_MPI
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab/FileOptions.hpp"
#include "MBTagConventions.hpp"
#include "moab_mpi.h"
#endif

using namespace moab;

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

#ifdef USE_MPI
std::string read_options;
#endif

int number_tests_successful = 0;
int number_tests_failed = 0;

void handle_error_code(ErrorCode rv, int &number_failed, int &number_successful)
{
  if (rv == MB_SUCCESS) {
#ifdef USE_MPI
      int rank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank==0)
          std::cout << "Success";
#else
      std::cout << "Success";
#endif
    number_successful++;
  } else {
    std::cout << "Failure";
    number_failed++;
  }
}

ErrorCode testA()
{
  ErrorCode error;
  Core mb;
  Interface *mbImpl = &mb;
  MeshTopoUtil mtu(mbImpl);

  EntityHandle fset;
  error = mbImpl->create_meshset(MESHSET_SET, fset);CHECK_ERR(error);

  //Create a simple hex mesh
  error= create_simple_mesh(mbImpl, &fset, 1);CHECK_ERR(error);

  //Construct the AHF maps
  HalfFacetRep ahf(&mb, NULL, fset);

  error = ahf.initialize();CHECK_ERR(error);

  //Create another mesh
  error = create_simple_mesh(mbImpl, &fset, 2);CHECK_ERR(error);

  Range verts, hexes;
  error = mbImpl->get_entities_by_dimension(fset, 0, verts);CHECK_ERR(error);
  error = mbImpl->get_entities_by_dimension(fset, 3, hexes);CHECK_ERR(error);

  //Update the AHF maps for the new set and join the new mesh with the existing one
  Range setAhexes, setAverts, setBhexes, setBverts;
  setAhexes.insert(hexes.begin(), hexes.end()-3); setAhexes.print();
  setAverts.insert(verts.begin(), verts.end()-7); setAverts.print();
  setBhexes.insert(hexes.end()-2, hexes.end()); setBhexes.print();
  setBverts.insert(verts.end()-6, verts.end());setBverts.print();

  HalfFacetRep::ESet setA, setB;
  setA.verts = setAverts;
  setA.entities = setAhexes;
  setB.verts = setBverts;
  setB.entities = setBhexes;

  error = ahf.update_hf_maps(setB);CHECK_ERR(error);
  error = ahf.update_hf_maps(setA, setB);CHECK_ERR(error);

  //Check adjacencies
  std::vector<EntityHandle> adjents;
  Range mbents, ahfents;
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
      adjents.clear();
      error = ahf.get_up_adjacencies( *i, 3, adjents);
      CHECK_ERR(error);
      mbents.clear();
      error = mbImpl->get_adjacencies(&*i, 1, 3, false, mbents);
      CHECK_ERR(error);
      CHECK_EQUAL(adjents.size(), mbents.size());
      std::sort(adjents.begin(), adjents.end());
      std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
      mbents = subtract(mbents, ahfents);
      CHECK(!mbents.size());
    }

  for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
      adjents.clear();
      error = ahf.get_neighbor_adjacencies( *i, adjents);
      CHECK_ERR(error);
      mbents.clear();
      error = mtu.get_bridge_adjacencies( *i, 2, 3, mbents);
      CHECK_ERR(error);
      CHECK_EQUAL(adjents.size(), mbents.size());
      std::sort(adjents.begin(), adjents.end());
      std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
      mbents = subtract(mbents, ahfents);
      CHECK(!mbents.size());
    }

  error = ahf.deinitialize();CHECK_ERR(error);

   return MB_SUCCESS;

}

ErrorCode create_simple_mesh(Interface *mbImpl, EntityHandle *fileset, int set)
{
  ErrorCode error;
  if (set == 1){
      const double coords[] = {0,0,0,
                               1,0,0,
                               2,0,0,
                               2,1,0,
                               1,1,0,
                               0,1,0,
                               0,0,1,
                               1,0,1,
                               2,0,1,
                               2,1,1,
                               1,1,1,
                               0,1,1,
                               0,0,2,
                               1,0,2,
                               2,0,2,
                               2,1,2,
                               1,1,2,
                               0,1,2};

      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      const int conn[] = {0,1,4,5,6,7,10,11,
                          1,2,3,4,7,8,9,10,
                          6,7,10,11,12,13,16,17,
                          7,8,9,10,13,14,15,16};

      const size_t num_elems = sizeof(conn)/sizeof(int)/8;

      EntityHandle verts[num_vtx], cells[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[8];
          for (int j=0; j<8; j++)
            c[j] = verts[conn[8*i+j]];

          error = mbImpl->create_element(MBHEX, c, 8, cells[i]); CHECK_ERR(error);
        }

      Range vertices, hexes;
      error = mbImpl->get_entities_by_dimension(0, 0, vertices);CHECK_ERR(error);
      error = mbImpl->get_entities_by_dimension(0, 3, hexes);CHECK_ERR(error);

      error = mbImpl->add_entities(*fileset, vertices);CHECK_ERR(error);
      error = mbImpl->add_entities(*fileset, hexes);CHECK_ERR(error);
    }
  else if (set == 2)
    {
      Range everts;
      error = mbImpl->get_entities_by_dimension(*fileset, 0, everts);CHECK_ERR(error);

      const double coords[] = {3,0,0,
                               3,1,0,
                               3,0,1,
                               3,1,1,
                               3,0,2,
                               3,1,2};

      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      EntityHandle verts[num_vtx];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      EntityHandle c1[8] = {everts[2], verts[0], verts[1], everts[3], everts[8], verts[2], verts[3], everts[9]};
      EntityHandle c1[8] = {everts[8], verts[2], verts[3], everts[9], everts[14], verts[4], verts[5], everts[15]};

      EntityHandle hex1, hex2;

      error = mbImpl->create_element(MBHEX, c1, 8, hex1); CHECK_ERR(error);
      error = mbImpl->create_element(MBHEX, c2, 8, hex2); CHECK_ERR(error);

      error = mbImpl->add_entities(*fileset, verts, num_vtx);CHECK_ERR(error);
      error = mbImpl->add_entities(*fileset, &hex1, 1);CHECK_ERR(error);
      error = mbImpl->add_entities(*fileset, &hex2, 1);CHECK_ERR(error);
    }
  return MB_SUCCESS;
}


int main(int argc, char *argv[])
{
    ErrorCode result;

    std::cout<<"testA:";

    result = testA();
    handle_error_code(result, number_tests_failed, number_tests_successful);
    std::cout<<"\n";

    return number_tests_failed;
}

