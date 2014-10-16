/*This function tests the AHF datastructures on CST meshes*/
#include <iostream>
#include <vector>
#include <algorithm>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/HalfFacetRep.hpp"
#include "../RefineMesh/moab/NestedRefine.hpp"
#include "TestUtil.hpp"

#ifdef USE_MPI
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "ReadParallel.hpp"
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

ErrorCode single_entity_single_level(Core *mb, EntityType type, int degree)
{
  ErrorCode error;

  NestedRefine uref(mb);
  EntityHandle set;
  error = uref.generate_mesh_hierarchy(&degree, 1, set);
  CHECK_ERR(error);

  Range verts, ents;

  error = mb->get_entities_by_type(set, MBVERTEX, verts);
  CHECK_ERR(error);

  error = mb->get_entities_by_type(set, type, ents);
  CHECK_ERR(error);
}

ErrorCode single_entity_multi_level(Core *mb, EntityType type, int *level_degrees, const int num_levels)
{
  ErrorCode error;

  NestedRefine uref(mb);
  EntityHandle set[num_levels];
  error = uref.generate_mesh_hierarchy(level_degrees, 1, set);
  CHECK_ERR(error);

  for (int l=0; l<num_levels; l++)
    {
      Range verts, ents;

      error = mb->get_entities_by_type(set[l], MBVERTEX, verts);
      CHECK_ERR(error);

      error = mb->get_entities_by_type(set[l], type, ents);
      CHECK_ERR(error);


    }

}

ErrorCode multiple_entities_single_level()
{

}

ErrorCode multiple_entities_multiple_levels()
{

}

ErrorCode create_mesh(Core *mb, EntityType type)
{
  ErrorCode error;
  Interface* mbImpl = mb;
  if (type == MBEDGE)
    {
      const double coords[] = {0,0,0,
                              1,0,0};
      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      const int conn[] = {0, 1};
      const size_t num_elems = sizeof(conn)/sizeof(int)/2;

      EntityHandle verts[num_vtx], edges[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]);
          if (error != MB_SUCCESS) return error;
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[2];
          c[0] = verts[conn[0]]; c[1] = verts[conn[1]];

          error = mbImpl->create_element(MBEDGE, c, 2, edges[i]);
          if (error != MB_SUCCESS) return error;
        }

    }
  else if (type == MBTRI)
    {
      const double coords[] = {0,0,0,
                              1,0,0,
                              0,1,0};
      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      const int conn[] = {0, 1, 2};
      const size_t num_elems = sizeof(conn)/sizeof(int)/3;

      EntityHandle verts[num_vtx], faces[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]);
          if (error != MB_SUCCESS) return error;
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[3];
          for (int j=0; j<3; j++)
            c[j] = verts[conn[3*i+j]];

          error = mbImpl->create_element(MBTRI, c, 3, faces[i]);
          if (error != MB_SUCCESS) return error;
        }
    }
  else if (type == MBQUAD)
    {
          const double coords[] = {0,0,0,
                                  1,0,0,
                                  1,1,0,
                                  0,1,0};
          const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

          const int conn[] = {0, 1, 2, 3};
          const size_t num_elems = sizeof(conn)/sizeof(int)/3;

          EntityHandle verts[num_vtx], faces[num_elems];
          for (size_t i=0; i< num_vtx; ++i)
            {
              error = mbImpl->create_vertex(coords+3*i, verts[i]);
              if (error != MB_SUCCESS) return error;
            }

          for (size_t i=0; i< num_elems; ++i)
            {
              EntityHandle c[4];
              for (int j=0; j<4; j++)
                c[j] = verts[conn[j]];

              error = mbImpl->create_element(MBQUAD, c, 4, faces[i]);
              if (error != MB_SUCCESS) return error;

            }
    }
  else if (type == MBTET)
    {
      const double coords[] = {0,0,0,
                              1,0,0,
                              0,1,0,
                              0,0,1};
      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      const int conn[] = {0, 1, 2, 3};
      const size_t num_elems = sizeof(conn)/sizeof(int)/3;

      EntityHandle verts[num_vtx], cells[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]);
          if (error != MB_SUCCESS) return error;
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[4];
          for (int j=0; j<4; j++)
            c[j] = verts[conn[j]];

          error = mbImpl->create_element(MBTET, c, 4, cells[i]);
          if (error != MB_SUCCESS) return error;
        }
    }
  else if (type == MBPRISM)
    {
      const double coords[] = {0,0,0,
                              1,0,0,
                              0,1,0,
                              0,0,1,
                              1,0,1,
                              0,1,1};
      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      const int conn[] = {0, 1, 2, 3, 4, 5};
      const size_t num_elems = sizeof(conn)/sizeof(int)/3;

      EntityHandle verts[num_vtx], cells[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]);
          if (error != MB_SUCCESS) return error;
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[6];
          for (int j=0; j<6; j++)
            c[j] = verts[conn[j]];

          error = mbImpl->create_element(MBPRISM, c, 6, cells[i]);
          if (error != MB_SUCCESS) return error;
        }
    }
  else if (type == MBHEX)
  {
    const double coords[] = {0,0,0,
                            1,0,0,
                            1,1,0,
                            0,1,0,
                            0,0,1,
                            1,0,1,
                            1,1,1,
                            0,1,1};
    const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

    const int conn[] = {0, 1, 2, 3, 4, 5, 6, 7};
    const size_t num_elems = sizeof(conn)/sizeof(int)/3;

    EntityHandle verts[num_vtx], cells[num_elems];
    for (size_t i=0; i< num_vtx; ++i)
      {
        error = mbImpl->create_vertex(coords+3*i, verts[i]);
        if (error != MB_SUCCESS) return error;
      }

    for (size_t i=0; i< num_elems; ++i)
      {
        EntityHandle c[8];
        for (int j=0; j<8; j++)
          c[j] = verts[conn[j]];

        error = mbImpl->create_element(MBHEX, c, 8, cells[i]);
        if (error != MB_SUCCESS) return error;
      }
    }
  return MB_SUCCESS;
}

ErrorCode test_single_entity(EntityType type, int *level_degrees, int num_levels)
{
  ErrorCode error;
  Core mb;
  error = create_mesh(&mb, type);
  if (error != MB_SUCCESS) return error;
  if (num_levels == 1)
    {
      error = single_entity_single_level(&mb, level_degrees[0]);
      if (error != MB_SUCCESS) return error;
    }
  else
    {
      error = single_entity_multi_level(&mb, level_degrees, num_levels);
      if (error != MB_SUCCESS) return error;
    }

}
ErrorCode test_all_entity_types()
{
  ErrorCode error;
  //Loop over all entity types and refine them to 1 level

  for (EntityType type = MBEDGE ; type != MBMAXTYPE; ++type)
    {
      int *degree;
      degree[0] = 2;
      error = test_single_entity(type, degree, 1);
      CHECK_ERR(error);

      degree[0] = 3;
      error = test_single_entity(type, degree, 1);
      CHECK_ERR(error);
    }

  return MB_SUCCESS;
}

ErrorCode ahf_test(const char* filename)
{

    Core moab;
    Interface* mbImpl = &moab;
    MeshTopoUtil mtu(mbImpl);
    ErrorCode error;

#ifdef USE_MPI
    int procs = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    if (procs > 1){
    read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=3.0.1.3";

    error = mbImpl->load_file(filename, 0, read_options.c_str());
    CHECK_ERR(error);
    }
    else if (procs == 1) {
#endif
    error = mbImpl->load_file(filename);
    CHECK_ERR(error);
#ifdef USE_MPI
    }
#endif

    /*Create ranges for handles of explicit elements of the mixed mesh*/
    Range verts, edges, faces, cells;
    error = mbImpl->get_entities_by_dimension( 0, 0, verts);
    error = mbImpl->get_entities_by_dimension( 0, 1, edges);
    error = mbImpl->get_entities_by_dimension( 0, 2, faces);
    error = mbImpl->get_entities_by_dimension( 0, 3, cells);

    // Create an ahf instance
    HalfFacetRep ahf(&moab);

    // Call the initialize function which creates the maps for each dimension
    ahf.initialize();

    //ahf.print_tags();

    //Perform queries
    std::vector<EntityHandle> adjents;
    Range mbents, ahfents;

    //1D Queries //
    //IQ1: For every vertex, obtain incident edges
    for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
        adjents.clear();
        error = ahf.get_up_adjacencies( *i, 1, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies( &*i, 1, 1, false, mbents );
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(),mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    //NQ1:  For every edge, obtain neighbor edges
    for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
        adjents.clear();
        error = ahf.get_neighbor_adjacencies( *i, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mtu.get_bridge_adjacencies( *i, 0, 1, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    // 2D Queries

    // IQ21: For every vertex, obtain incident faces
    for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
        adjents.clear();
        error = ahf.get_up_adjacencies( *i, 2, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies( &*i, 1, 2, false, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    //IQ22: For every edge, obtain incident faces
    for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
        adjents.clear();
        error = ahf.get_up_adjacencies( *i, 2, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies( &*i, 1, 2, false, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    //NQ2: For every face, obtain neighbor faces
    for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
        adjents.clear();
        error = ahf.get_neighbor_adjacencies( *i, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mtu.get_bridge_adjacencies( *i, 1, 2, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    //DQ 21: For every face, obtain its edges
    for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
        adjents.clear();
        error = ahf.get_down_adjacencies( *i, 1, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies( &*i, 1, 1, false, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    // 3D Queries
    //IQ 31: For every vertex, obtain incident cells
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

    // IQ 32: For every edge, obtain incident cells
    for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
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

    //IQ33: For every face, obtain incident cells
    for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
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

    //NQ3: For every cell, obtain neighbor cells
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


    //DQ 31: For every cell, obtain its edges
    for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
        adjents.clear();
        error = ahf.get_down_adjacencies( *i, 1, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies( &*i, 1, 1, false, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    //DQ 32: For every cell, obtain its faces
    for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
        adjents.clear();
        error = ahf.get_down_adjacencies( *i, 2, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies( &*i, 1, 2, false, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }



    ahf.deinitialize();

    return MB_SUCCESS;

}

int main(int argc, char *argv[])
{

#ifdef USE_MPI
    MPI_Init(&argc, &argv);

    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

 /*   const char* filename = 0;
#ifdef MESHDIR
 #ifdef HDF5_FILE
    filename = STRINGIFY(MESHDIR) "/32hex_ef.h5m";
 #else
    filename = STRINGIFY(MESHDIR) "/hexes_mixed.vtk";
 #endif
#else
 #ifdef HDF5_FILE
    filename = "32hex_ef.h5m";
 #else
    filename = "hexes_mixed.vtk";
 #endif
#endif */


    if (argc==1)
    {
#ifdef USE_MPI
        if (rank == 0)
            std::cout<<"Using default input file:"<<filename<<std::endl;
#else
        std::cout<<"Using default input file:"<<filename<<std::endl;
#endif
    }

    else if (argc==2)
        filename = argv[1];
    else {
            std::cerr << "Usage: " << argv[0] << " [filename]" << std::endl;
            return 1;
    }

    ErrorCode result;

#ifdef USE_MPI
    if (rank == 0)
        std::cout<<" para_ahf_test: ";
#else
    std::cout<<"ahf_test:";
#endif

    result = ahf_test(filename);
    handle_error_code(result, number_tests_failed, number_tests_successful);
    std::cout<<"\n";

#ifdef USE_MPI
    MPI_Finalize();
#endif

    return number_tests_failed;
}

