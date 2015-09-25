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

      const size_t num_hexes = sizeof(conn)/sizeof(int)/8;

      EntityHandle verts[num_vtx], hexes[num_hexes];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      for (size_t i=0; i< num_hexes; ++i)
        {
          EntityHandle c[8];
          for (int j=0; j<8; j++)
            c[j] = verts[conn[8*i+j]];

          error = mbImpl->create_element(MBHEX, c, 8, hexes[i]); CHECK_ERR(error);
        }

      error = mbImpl->add_entities(*fileset, verts, num_vtx);CHECK_ERR(error);
      error = mbImpl->add_entities(*fileset, hexes, num_hexes);CHECK_ERR(error);
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
      EntityHandle c2[8] = {everts[8], verts[2], verts[3], everts[9], everts[14], verts[4], verts[5], everts[15]};

      EntityHandle hex1, hex2;

      error = mbImpl->create_element(MBHEX, c1, 8, hex1); CHECK_ERR(error);
      error = mbImpl->create_element(MBHEX, c2, 8, hex2); CHECK_ERR(error);

      error = mbImpl->add_entities(*fileset, verts, num_vtx);CHECK_ERR(error);
      error = mbImpl->add_entities(*fileset, &hex1, 1);CHECK_ERR(error);
      error = mbImpl->add_entities(*fileset, &hex2, 1);CHECK_ERR(error);
    }

  return MB_SUCCESS;
}

ErrorCode test_adjacencies(Interface *mbImpl , HalfFacetRep *ahf, Range &verts, Range &cells)
{
  ErrorCode error;
  MeshTopoUtil mtu(mbImpl);

  //Check adjacencies
  std::vector<EntityHandle> adjents;
  Range mbents, ahfents;
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
      adjents.clear();
      error = ahf->get_up_adjacencies( *i, 3, adjents);
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
      error = ahf->get_neighbor_adjacencies( *i, adjents);
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

  return MB_SUCCESS;
}

ErrorCode testA()
{
  ErrorCode error;
  Core mb;
  Interface *mbImpl = &mb;

  EntityHandle fset;
  error = mbImpl->create_meshset(MESHSET_SET, fset);CHECK_ERR(error);

  //Create a simple hex mesh
  error= create_simple_mesh(mbImpl, &fset, 1);CHECK_ERR(error);

  //Write out mesh for patchA
  std::stringstream fileA;
  fileA <<  "Mesh_patchA.vtk";
  std::string str = fileA.str();
  const char* output_file = str.c_str();
  error = mbImpl->write_file(output_file, 0, NULL, &fset, 1); CHECK_ERR(error);

  Range verts, hexes;
  error = mbImpl->get_entities_by_dimension(fset, 0, verts);CHECK_ERR(error);
  error = mbImpl->get_entities_by_dimension(fset, 3, hexes);CHECK_ERR(error);

  //Construct the AHF maps
  HalfFacetRep ahf(&mb, NULL, fset);

  error = ahf.initialize();CHECK_ERR(error);
  std::cout<<"HF maps: After initialize"<<std::endl;
  error = ahf.print_tags(3);CHECK_ERR(error);
  std::cout<<std::endl;

  //Test adjacencies of patchA
  error = test_adjacencies(mbImpl, &ahf, verts, hexes);CHECK_ERR(error);

  //Create another mesh
  error = create_simple_mesh(mbImpl, &fset, 2);CHECK_ERR(error);

  //Write out mesh for patchA+B
  std::stringstream fileB;
  fileB <<  "Mesh_patchAB.vtk";
  str = fileB.str();
  output_file = str.c_str();
  error = mbImpl->write_file(output_file, 0, NULL, &fset, 1); CHECK_ERR(error);

  //Update the AHF maps for the new set and join the new mesh with the existing one
  error = mbImpl->get_entities_by_dimension(fset, 0, verts);CHECK_ERR(error);
  error = mbImpl->get_entities_by_dimension(fset, 3, hexes);CHECK_ERR(error);

  std::vector<EntityHandle> patchA, patchB;
  patchA.insert(patchA.end(), hexes.begin(), hexes.end()-2);
  patchB.insert(patchB.end(), hexes.end()-2, hexes.end());

  for (int i=0; i<(int)patchA.size(); i++)
    std::cout<<"patchA["<<i<<"] = "<<patchA[i]<<std::endl;

  for (int i=0; i<(int)patchB.size(); i++)
    std::cout<<"patchB["<<i<<"] = "<<patchB[i]<<std::endl;

  error = ahf.update_hf_maps_patch(3, patchB);CHECK_ERR(error);

  std::cout<<"HF maps: After updating patch B"<<std::endl;
  error = ahf.print_tags(3);CHECK_ERR(error);
  std::cout<<std::endl;

  error = ahf.update_hf_maps_multimesh(3, patchA, patchB);CHECK_ERR(error);

   std::cout<<"HF maps: After joining patches A and B"<<std::endl;
  error = ahf.print_tags(3);CHECK_ERR(error);
  std::cout<<std::endl;

  //Test adjacencies of patchAB
  error = test_adjacencies(mbImpl, &ahf, verts, hexes);CHECK_ERR(error);

  error = ahf.deinitialize();CHECK_ERR(error);

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

