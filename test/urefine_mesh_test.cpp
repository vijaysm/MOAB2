/*This unit test is for the uniform refinement capability based on AHF datastructures*/
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "../RefineMesh/moab/NestedRefine.hpp"
#include "TestUtil.hpp"

#ifdef MOAB_HAVE_MPI
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

#ifdef MOAB_HAVE_MPI
std::string read_options;
#endif

int number_tests_successful = 0;
int number_tests_failed = 0;

void handle_error_code(ErrorCode rv, int &number_failed, int &number_successful)
{
  if (rv == MB_SUCCESS) {
#ifdef MOAB_HAVE_MPI
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

ErrorCode test_adjacencies(Interface *mbImpl, NestedRefine *nr, Range all_ents)
{
  MeshTopoUtil mtu(mbImpl);
  ErrorCode error;
  Range verts, edges, faces, cells;
  verts = all_ents.subset_by_dimension(0);
  edges = all_ents.subset_by_dimension(1);
  faces = all_ents.subset_by_dimension(2);
  cells = all_ents.subset_by_dimension(3);

  std::vector<EntityHandle> adjents;
  Range mbents, ahfents;

  if (!edges.empty())
    {
      //1D Queries //
      //IQ1: For every vertex, obtain incident edges
      for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
          adjents.clear(); mbents.clear(); ahfents.clear();
          error = nr->get_adjacencies( *i, 1, adjents);  CHECK_ERR(error);
          error = mbImpl->get_adjacencies( &*i, 1, 1, false, mbents ); CHECK_ERR(error);
          CHECK_EQUAL(adjents.size(),mbents.size());
          std::sort(adjents.begin(), adjents.end());
          std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
          mbents = subtract(mbents, ahfents);
          CHECK(!mbents.size());
      }

      //NQ1:  For every edge, obtain neighbor edges
      for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
          adjents.clear(); mbents.clear(); ahfents.clear();
          error = nr->get_adjacencies( *i, 1, adjents); CHECK_ERR(error);
          error = mtu.get_bridge_adjacencies( *i, 0, 1, mbents); CHECK_ERR(error);
          CHECK_EQUAL(adjents.size(), mbents.size());
          std::sort(adjents.begin(), adjents.end());
          std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
          mbents = subtract(mbents, ahfents);
          CHECK(!mbents.size());
      }
    }

  if (!faces.empty())
    {
      // IQ21: For every vertex, obtain incident faces
      for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
          adjents.clear(); mbents.clear(); ahfents.clear();
          error = nr->get_adjacencies( *i, 2, adjents); CHECK_ERR(error);
          error = mbImpl->get_adjacencies( &*i, 1, 2, false, mbents); CHECK_ERR(error);
          CHECK_EQUAL(adjents.size(), mbents.size());
          std::sort(adjents.begin(), adjents.end());
          std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
          mbents = subtract(mbents, ahfents);
          CHECK(!mbents.size());
        }

      //IQ22: For every edge, obtain incident faces
      if (!edges.empty()){
          for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
              adjents.clear(); mbents.clear(); ahfents.clear();
              error = nr->get_adjacencies( *i, 2, adjents);
              CHECK_ERR(error);
              error = mbImpl->get_adjacencies( &*i, 1, 2, false, mbents);
              CHECK_ERR(error);
              CHECK_EQUAL(adjents.size(), mbents.size());
              std::sort(adjents.begin(), adjents.end());
              std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
              mbents = subtract(mbents, ahfents);
              CHECK(!mbents.size());
            }
        }

      //NQ2: For every face, obtain neighbor faces
      for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
          adjents.clear(); mbents.clear(); ahfents.clear();
          error = nr->get_adjacencies( *i, 2, adjents); CHECK_ERR(error);
          error = mtu.get_bridge_adjacencies( *i, 1, 2, mbents); CHECK_ERR(error);
          CHECK_EQUAL(adjents.size(), mbents.size());
          std::sort(adjents.begin(), adjents.end());
          std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
          mbents = subtract(mbents, ahfents);
          CHECK(!mbents.size());
      }

      if (!edges.empty()){
          for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
              adjents.clear(); mbents.clear(); ahfents.clear();
              error = nr->get_adjacencies( *i, 1, adjents);
              CHECK_ERR(error);
              error = mbImpl->get_adjacencies( &*i, 1, 1, false, mbents);
              CHECK_ERR(error);
              CHECK_EQUAL(adjents.size(), mbents.size());
              std::sort(adjents.begin(), adjents.end());
              std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
              mbents = subtract(mbents, ahfents);
              CHECK(!mbents.size());
            }
        }
    }

  if (!cells.empty())
    {
      //IQ 31: For every vertex, obtain incident cells
      for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
          adjents.clear(); mbents.clear(); ahfents.clear();
          error = nr->get_adjacencies( *i, 3, adjents); CHECK_ERR(error);
          error = mbImpl->get_adjacencies(&*i, 1, 3, false, mbents); CHECK_ERR(error);
          CHECK_EQUAL(adjents.size(), mbents.size());
          std::sort(adjents.begin(), adjents.end());
          std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
          mbents = subtract(mbents, ahfents);
          CHECK(!mbents.size());
        }

      if (!edges.empty())
        {
          for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
              adjents.clear(); mbents.clear(); ahfents.clear();
              error = nr->get_adjacencies( *i, 3, adjents);
              CHECK_ERR(error);
              error = mbImpl->get_adjacencies(&*i, 1, 3, false, mbents);
              CHECK_ERR(error);
              CHECK_EQUAL(adjents.size(), mbents.size());
              std::sort(adjents.begin(), adjents.end());
              std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
              mbents = subtract(mbents, ahfents);
              CHECK(!mbents.size());
            }
        }

      if (!faces.empty())
        {
          for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
              adjents.clear(); mbents.clear(); ahfents.clear();
              error = nr->get_adjacencies( *i, 3, adjents);
              CHECK_ERR(error);
              error = mbImpl->get_adjacencies(&*i, 1, 3, false, mbents);
              CHECK_ERR(error);
              CHECK_EQUAL(adjents.size(), mbents.size());
              std::sort(adjents.begin(), adjents.end());
              std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
              mbents = subtract(mbents, ahfents);
              CHECK(!mbents.size());
            }
        }

      //NQ3: For every cell, obtain neighbor cells
      for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
          adjents.clear(); mbents.clear(); ahfents.clear();
          error = nr->get_adjacencies( *i, 3, adjents); CHECK_ERR(error);
          error = mtu.get_bridge_adjacencies( *i, 2, 3, mbents); CHECK_ERR(error);
          CHECK_EQUAL(adjents.size(), mbents.size());
          std::sort(adjents.begin(), adjents.end());
          std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
          mbents = subtract(mbents, ahfents);
          CHECK(!mbents.size());
      }

      if (!edges.empty())
        {
          for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
              adjents.clear(); mbents.clear(); ahfents.clear();
              error = nr->get_adjacencies( *i, 1, adjents);
              CHECK_ERR(error);
              error = mbImpl->get_adjacencies( &*i, 1, 1, false, mbents);
              CHECK_ERR(error);
              CHECK_EQUAL(adjents.size(), mbents.size());
              std::sort(adjents.begin(), adjents.end());
              std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
              mbents = subtract(mbents, ahfents);
              CHECK(!mbents.size());
            }
        }

      if (!faces.empty())
        {
          for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
              adjents.clear(); mbents.clear(); ahfents.clear();
              error = nr->get_adjacencies( *i, 2, adjents);
              CHECK_ERR(error);
              error = mbImpl->get_adjacencies( &*i, 1, 2, false, mbents);
              CHECK_ERR(error);
              CHECK_EQUAL(adjents.size(), mbents.size());
              std::sort(adjents.begin(), adjents.end());
              std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
              mbents = subtract(mbents, ahfents);
              CHECK(!mbents.size());
            }
        }
    }

  return MB_SUCCESS;
}


ErrorCode refine_entities(Interface *mb, int *level_degrees, const int num_levels, bool output)
{
  ErrorCode error;

  //Get the range of entities in the initial mesh
  Range init_ents[4];
  error = mb->get_entities_by_dimension(0, 0, init_ents[0]); CHECK_ERR(error);
  error = mb->get_entities_by_dimension(0, 1, init_ents[1]); CHECK_ERR(error);
  error = mb->get_entities_by_dimension(0, 2, init_ents[2]); CHECK_ERR(error);
  error = mb->get_entities_by_dimension(0, 3, init_ents[3]);  CHECK_ERR(error);

 int dim[3] = {1,2,3};

  if (output)
    {
    //  int inents = init_ents.size();
      std::stringstream file;
      file <<  "MESH_LEVEL_0.vtk";
      std::string str = file.str();
      const char* output_file = str.c_str();
      error = mb->write_file(output_file); CHECK_ERR(error);
    }

  //Create an hm object and generate the hierarchy
  std::cout<<"Creating a hm object"<<std::endl;
  NestedRefine uref(mb);
  EntityHandle *set = new EntityHandle[num_levels];

  std::cout<<"Starting hierarchy generation"<<std::endl;
  error = uref.generate_mesh_hierarchy( num_levels,level_degrees, set); CHECK_ERR(error);
  std::cout<<"Finished hierarchy generation"<<std::endl;

  std::cout<<std::endl;
  std::cout<<"Mesh size for level 0  :: inverts = "<<init_ents[0].size()<<", inedges = "<<init_ents[1].size()<<", infaces = "<<init_ents[2].size()<<", incells = "<<init_ents[3].size()<<std::endl;

  Range prev_ents[4];
  for (int i=0; i<4; i++)
    prev_ents[i] = init_ents[i];


  //Loop over each mesh level and check its topological properties
  for (int l=0; l<num_levels; l++)
    {
      Range all_ents;
      error = mb->get_entities_by_handle(set[l], all_ents); CHECK_ERR(error);

      Range ents[4];
      for (int k=0; k<4; k++)
        ents[k] = all_ents.subset_by_dimension(k);

      if (ents[0].empty() || all_ents.empty())
        std::cout<<"Something is not right"<<std::endl;

      std::cout<<std::endl;
      std::cout<<"Mesh size for level "<<l+1<<"  :: nverts = "<<ents[0].size()<<", nedges = "<<ents[1].size()<<", nfaces = "<<ents[2].size()<<", ncells = "<<ents[3].size()<<std::endl;

      //Check if the number of new entities created are correct.

      for (int type =0; type <3; type++)
        {
          int factor = 1;
          if (!ents[type+1].empty()){

              for (int p=0; p<=l; p++){
                  for (int d=0; d<dim[type]; d++)
                    factor *= level_degrees[p];
                }
              int  expected_nents = factor*init_ents[type+1].size();
              CHECK_EQUAL(expected_nents, (int)ents[type+1].size());
            }
        }

      //Check adjacencies
      error = test_adjacencies(mb, &uref, all_ents); CHECK_ERR(error);

      //Check interlevel child-parent query between previous and current level
      for (int type = 0; type < 3; type++)
        {
          if (!prev_ents[type+1].empty())
            {
              for (Range::iterator e = prev_ents[type+1].begin(); e != prev_ents[type+1].end(); e++)
                {
                  std::vector<EntityHandle> children;
                  error = uref.parent_to_child(*e, l, l+1, children); CHECK_ERR(error);
                  for (int i=0; i<(int)children.size(); i++)
                    {
                      EntityHandle parent;
                      error = uref.child_to_parent(children[i], l+1, l, &parent); CHECK_ERR(error);
                      assert(parent == *e);
                    }
                }
            }
        }

      for (int i=0; i<4; i++)
        prev_ents[i] = ents[i];

      //Print out the mesh
      if (output)
        {
          std::stringstream file;
          file <<  "MESH_LEVEL_" <<l+1<<".vtk";
          std::string str = file.str();
          const char* output_file = str.c_str();
          char * write_opts = NULL;
          error = mb->write_file(output_file, 0, write_opts, &set[l], 1); CHECK_ERR(error);
        }
    }

  //Check interlevel child-parent query between initial and most refined mesh
  for (int type = 0; type < 3; type++)
    {
      if (!init_ents[type+1].empty())
        {
          for (Range::iterator e= init_ents[type+1].begin(); e != init_ents[type+1].end(); e++)
            {
              std::vector<EntityHandle> children;
              error = uref.parent_to_child(*e, 0, num_levels, children); CHECK_ERR(error);
              for (int i=0; i<(int)children.size(); i++)
                {
                  EntityHandle parent;
                  error = uref.child_to_parent(children[i], num_levels, 0, &parent); CHECK_ERR(error);
                  assert(parent == *e);
                }
            }
        }
    }

  //Print out the whole hierarchy into a single file
  if (output)
    {
      std::stringstream file;
      file <<  "MESH_HIERARCHY.vtk";
      std::string str = file.str();
      const char* output_file = str.c_str();
      error = mb->write_file(output_file); CHECK_ERR(error);
    }

  delete [] set;
  return MB_SUCCESS;
}

ErrorCode create_single_entity(Interface *mbImpl, EntityType type)
{
  ErrorCode error;
  if (type == MBEDGE)
    {
      const double coords[] = {0.0,0.0,0.0,
                              1.0,0.0,0.0};
      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      const int conn[] = {0, 1};
      const size_t num_elems = sizeof(conn)/sizeof(conn[0])/2;

      std::cout<<"Specify verts and ents"<<std::endl;

      EntityHandle verts[num_vtx], edges[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      std::cout<<"Created vertices"<<std::endl;

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[2];
          c[0] = verts[conn[0]]; c[1] = verts[conn[1]];

          error = mbImpl->create_element(MBEDGE, c, 2, edges[i]); CHECK_ERR(error);
        }

      std::cout<<"Created ents"<<std::endl;

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

          error = mbImpl->create_element(MBTRI, c, 3, faces[i]); CHECK_ERR(error);
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
              error = mbImpl->create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
            }

          for (size_t i=0; i< num_elems; ++i)
            {
              EntityHandle c[4];
              for (int j=0; j<4; j++)
                c[j] = verts[conn[j]];

              error = mbImpl->create_element(MBQUAD, c, 4, faces[i]); CHECK_ERR(error);

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
      const size_t num_elems = sizeof(conn)/sizeof(int)/4;

      EntityHandle verts[num_vtx], cells[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[4];
          for (int j=0; j<4; j++)
            c[j] = verts[conn[j]];

          error = mbImpl->create_element(MBTET, c, 4, cells[i]); CHECK_ERR(error);
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
          c[j] = verts[conn[j]];

        error = mbImpl->create_element(MBHEX, c, 8, cells[i]); CHECK_ERR(error);
      }
    }
  return MB_SUCCESS;
}

ErrorCode create_mesh(Interface *mbImpl, EntityType type)
{
  ErrorCode error;
  if (type == MBEDGE)
    {
      const double coords[] = {0,0,0,
                              1,0,0,
                              0,1,0,
                              -1,0,0,
                               0,-1,0};
      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      const int conn[] = {1,0,
                         0,3,
                         2,0,
                         0,4};
      const size_t num_elems = sizeof(conn)/sizeof(int)/2;

      EntityHandle verts[num_vtx], edges[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[2];
          c[0] = verts[conn[2*i]]; c[1] = verts[conn[2*i+1]];

          error = mbImpl->create_element(MBEDGE, c, 2, edges[i]); CHECK_ERR(error);
        }

    }
  else if (type == MBTRI)
    {
     const double coords[] = {0,0,0,
                              1,-1,0,
                              1,1,0,
                              -1,1,0,
                              -1,-1,0,
                              0,0,1};

      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

     const int conn[] = {0, 1, 2,
                         0,2,3,
                         0,3,4,
                         0,4,1,
                         0,5,3,
                         0,2,5,
                         0,4,5,
                         0,5,1};

      const size_t num_elems = sizeof(conn)/sizeof(int)/3;

      EntityHandle verts[num_vtx], faces[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[3];
          for (int j=0; j<3; j++)
            c[j] = verts[conn[3*i+j]];

          error = mbImpl->create_element(MBTRI, c, 3, faces[i]); CHECK_ERR(error);
        }
    }
  else if (type == MBQUAD)
    {
        const double coords[] = {0,0,0,
                                  1,0,0,
                                  1,1,0,
                                  0,1,0,
                                  -1,1,0,
                                  -1,0,0,
                                  -1,-1,0,
                                  0,-1,0,
                                  1,-1,0,
                                  0,1,1,
                                  0,0,1,
                                  0,-1,1};

          const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

          const int conn[] = {0, 1, 2, 3,
                             0,3,4,5,
                             7,8,1,0,
                             6,7,0,5,
                             0,3,9,10,
                             0,10,11,7};

          const size_t num_elems = sizeof(conn)/sizeof(int)/4;

          EntityHandle verts[num_vtx], faces[num_elems];
          for (size_t i=0; i< num_vtx; ++i)
            {
              error = mbImpl->create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
            }

          for (size_t i=0; i< num_elems; ++i)
            {
              EntityHandle c[4];
              for (int j=0; j<4; j++)
                c[j] = verts[conn[4*i+j]];

              error = mbImpl->create_element(MBQUAD, c, 4, faces[i]); CHECK_ERR(error);

            }
    }
  else if (type == MBTET)
    {
      const double coords[] = {0,-1,0,
                               0,2,0,
                               1,0,0,
                               -0.5,0,0,
                               0,0,1,
                               0,0,-2};

      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      const int conn[] = {0,2,1,4,
                         3,0,1,4,
                         5,2,1,0,
                         5,0,1,3};

      const size_t num_elems = sizeof(conn)/sizeof(int)/4;

      EntityHandle verts[num_vtx], cells[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[4];
          for (int j=0; j<4; j++)
            c[j] = verts[conn[4*i+j]];

          error = mbImpl->create_element(MBTET, c, 4, cells[i]); CHECK_ERR(error);
        }
    }
  else if (type == MBHEX)
  {
    const double coords[] = {0,-1,0,
                            1,-1,0,
                            1, 1, 0,
                            0, 1, 0,
                           -1, 1, 0,
                            -1,-1,0,
                             0,-1,1,
                             1,-1,1,
                             1,1,1,
                             0,1,1,
                             -1,1,1,
                             -1,-1,1,
                             0,-1,-1,
                             1,-1,-1,
                             1,1,-1,
                             0,1,-1,
                             -1,1,-1,
                             -1,-1,-1};
    const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

    const int conn[] = {0, 1, 2, 3, 6, 7, 8, 9,
                       5,0,3,4,11,6,9,10,
                       12,13,14,15,0,1,2,3,
                       17,12,15,16,5,0,3,4};
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
    }
  return MB_SUCCESS;
}

ErrorCode create_simple_mesh(Interface *mbImpl, EntityType type)
{
  ErrorCode error;
  if (type == MBEDGE)
    {
      const double coords[] = {0,0,0,
                               1,0,0,
                               2,0,0,
                               3,0,0,
                               4,0,0,
                               4,1,0,
                               3,1,0,
                               2,1,0,
                               1,1,0,
                               0,1,0};
      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      const int conn[] = {1,0,
                         1,2,
                         2,3,
                         3,4,
                         4,5,
                         5,6,
                         6,7,
                         7,8,
                         8,9,
                         9,0};
      const size_t num_elems = sizeof(conn)/sizeof(int)/2;

      EntityHandle verts[num_vtx], edges[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[2];
          c[0] = verts[conn[2*i]]; c[1] = verts[conn[2*i+1]];

          error = mbImpl->create_element(MBEDGE, c, 2, edges[i]); CHECK_ERR(error);
        }

    }
  else if (type == MBTRI)
    {
     const double coords[] = {0,0,0,
                              1,0,0,
                              2,0,0,
                              2.5,1,0,
                              1.5,1,0,
                              0.5,1,0,
                              -0.5,1,0,
                              -0.5,-1,0,
                              0.5,-1,0,
                              1.5,-1,0,
                              2.5,-1,0};

      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

     const int conn[] = {0,5,6,
                        0,1,5,
                        1,4,5,
                        1,2,4,
                        2,3,4,
                        7,8,0,
                        8,1,0,
                        8,9,1,
                        9,2,1,
                        9,10,2};

      const size_t num_elems = sizeof(conn)/sizeof(int)/3;

      EntityHandle verts[num_vtx], faces[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[3];
          for (int j=0; j<3; j++)
            c[j] = verts[conn[3*i+j]];

          error = mbImpl->create_element(MBTRI, c, 3, faces[i]); CHECK_ERR(error);
        }
    }
  else if (type == MBQUAD)
    {
        const double coords[] = {0,0,0,
                                 1,0,0,
                                 2,0,0,
                                 3,0,0,
                                 4,0,0,
                                 5,0,0,
                                 0,1,0,
                                 1,1,0,
                                 2,1,0,
                                 3,1,0,
                                 4,1,0,
                                 5,1,0,
                                 0,2,0,
                                 1,2,0,
                                 2,2,0,
                                 3,2,0,
                                 4,2,0,
                                 5,2,0};

          const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

          const int conn[] = {0, 1,7,6,
                              1,2,8,7,
                              2,3,9,8,
                              3,4,10,9,
                              4,5,11,10,
                              6,7,13,12,
                              7,8,14,13,
                              8,9,15,14,
                              9,10,16,15,
                              10,11,17,16};

          const size_t num_elems = sizeof(conn)/sizeof(int)/4;

          EntityHandle verts[num_vtx], faces[num_elems];
          for (size_t i=0; i< num_vtx; ++i)
            {
              error = mbImpl->create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
            }

          for (size_t i=0; i< num_elems; ++i)
            {
              EntityHandle c[4];
              for (int j=0; j<4; j++)
                c[j] = verts[conn[4*i+j]];

              error = mbImpl->create_element(MBQUAD, c, 4, faces[i]); CHECK_ERR(error);

            }
    }
  else if (type == MBTET)
    {
      const double coords[] = {0,0,0,
                               1,0,0,
                               0,1,0,
                               -1,0,0,
                               0,-1,0,
                               0,0,1};

      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      const int conn[] = {0,1,2,5,
                         3,0,2,5,
                         4,1,0,5,
                         4,0,3,5};

      const size_t num_elems = sizeof(conn)/sizeof(int)/4;

      EntityHandle verts[num_vtx], cells[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mbImpl->create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[4];
          for (int j=0; j<4; j++)
            c[j] = verts[conn[4*i+j]];

          error = mbImpl->create_element(MBTET, c, 4, cells[i]); CHECK_ERR(error);
        }
    }
  else if (type == MBHEX)
  {
    const double coords[] = {0,0,0,
                             1,0,0,
                             2,0,0,
                             0,1,0,
                             1,1,0,
                             2,1,0,
                             0,2,0,
                             1,2,0,
                             2,2,0,
                             0,0,1,
                             1,0,1,
                             2,0,1,
                             0,1,1,
                             1,1,1,
                             2,1,1,
                             0,2,1,
                             1,2,1,
                             2,2,1};
    const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

    const int conn[] = {0,1,4,3,9,10,13,12,
                       1,2,5,4,10,11,14,13,
                       3,4,7,6,12,13,16,15,
                       4,5,8,7,13,14,17,16};
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
    }
  return MB_SUCCESS;
}


ErrorCode test_entities(int mesh_type, EntityType type, int *level_degrees, int num_levels, bool output)
{
  ErrorCode error;
  Core mb;
  Interface* mbimpl = &mb;

  //Create entities
  if (mesh_type == 1){
      error = create_single_entity(mbimpl, type);
      if (error != MB_SUCCESS) return error;
      std::cout<<"Entity created successfully"<<std::endl;
    }
  else if (mesh_type == 2)
    {
      error = create_mesh(mbimpl, type);
      if (error != MB_SUCCESS) return error;
      std::cout<<"Small mesh created successfully"<<std::endl;
    }
  else if (mesh_type == 3)
    {
      error = create_simple_mesh(mbimpl, type);
      if (error != MB_SUCCESS) return error;
      std::cout<<"Small simple mesh created successfully"<<std::endl;
    }

  //Generate hierarchy
  error = refine_entities(mbimpl, level_degrees, num_levels, output);
  if (error != MB_SUCCESS) return error;

  return MB_SUCCESS;

}
ErrorCode test_1D()
{
  ErrorCode error;

  std::cout<<"Testing EDGE"<<std::endl;
  EntityType type = MBEDGE;

  std::cout<<"Testing single entity"<<std::endl;
  int deg[3] = {2,3,5};
  int len = sizeof(deg) / sizeof(int);
  error = test_entities(1, type, deg, len, false); CHECK_ERR(error);

  std::cout<<std::endl;
  std::cout<<"Testing a small mesh"<<std::endl;
  error = test_entities(2, type, deg, len, false); CHECK_ERR(error);

  std::cout<<std::endl;
  std::cout<<"Testing a small simple mesh"<<std::endl;
  int degree[4] = {5,5,2,2};
  len = sizeof(degree) / sizeof(int);
  error = test_entities(3, type, degree, len, false); CHECK_ERR(error);

  return MB_SUCCESS;
}

ErrorCode test_2D()
{
  ErrorCode error;

 std::cout<<"Testing TRI"<<std::endl;
  EntityType type = MBTRI;

  std::cout<<"Testing single entity"<<std::endl;
  int deg[3] = {2,3,5};
  int len = sizeof(deg) / sizeof(int);
  error = test_entities(1, type, deg, len, false); CHECK_ERR(error);

  std::cout<<std::endl;
  std::cout<<"Testing a small mesh"<<std::endl;
  error = test_entities(2, type, deg, len, false); CHECK_ERR(error);

  std::cout<<std::endl;
  std::cout<<"Testing a small simple mesh"<<std::endl;
  int degree[2] = {5,2};
  int length = sizeof(degree) / sizeof(int);
  error = test_entities(3, type, degree, length, false); CHECK_ERR(error);

  std::cout<<std::endl;
  std::cout<<"Testing QUAD"<<std::endl;
  type = MBQUAD;

  std::cout<<"Testing single entity"<<std::endl;
  error = test_entities(1, type, deg, len, false); CHECK_ERR(error);

  std::cout<<std::endl;
  std::cout<<"Testing a small mesh"<<std::endl;
  error = test_entities(2, type, deg, len, false); CHECK_ERR(error);

  std::cout<<std::endl;
  std::cout<<"Testing a small simple mesh"<<std::endl;
  error = test_entities(3, type, degree, length, false); CHECK_ERR(error);

  return MB_SUCCESS;
}

ErrorCode test_3D()
{
  ErrorCode error;

  std::cout<<"Testing TET"<<std::endl;
  EntityType type = MBTET;
  int deg[2] = {2,3};
  int len = sizeof(deg) / sizeof(int);

  std::cout<<"Testing single entity"<<std::endl;
  error = test_entities(1, type, deg, len, false); CHECK_ERR(error);

  std::cout<<std::endl;
  std::cout<<"Testing a small mesh"<<std::endl;
  error = test_entities(2, type, deg, len, false); CHECK_ERR(error);

  std::cout<<std::endl;
  std::cout<<"Testing a small simple mesh"<<std::endl;
  int degree[4] = {2,2,2,2};
  int length = sizeof(degree) / sizeof(int);
  error = test_entities(3, type, degree, length, false); CHECK_ERR(error);

  std::cout<<std::endl;
  std::cout<<"Testing HEX"<<std::endl;
  type = MBHEX;

 std::cout<<"Testing single entity"<<std::endl;
  error = test_entities(1, type, deg, len, false); CHECK_ERR(error);

  std::cout<<std::endl;
  std::cout<<"Testing a small mesh"<<std::endl;
  error = test_entities(2, type, deg, len, false); CHECK_ERR(error);

  std::cout<<std::endl;
  std::cout<<"Testing a small simple mesh"<<std::endl;
  error = test_entities(3, type, degree, length, false); CHECK_ERR(error);

  return MB_SUCCESS;
}

ErrorCode test_mesh(const char* filename, int *level_degrees, int num_levels)
{
  Core moab;
  Interface* mbImpl = &moab;
  ErrorCode error;

#ifdef MOAB_HAVE_MPI
    int procs = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    if (procs > 1){
    read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;";

    error = mbImpl->load_file(filename, 0, read_options.c_str()); CHECK_ERR(error);
    }
    else if (procs == 1) {
#endif
    error = mbImpl->load_file(filename);  CHECK_ERR(error);
#ifdef MOAB_HAVE_MPI
    }
#endif

    Range verts, edges,faces,cells;
    error = mbImpl->get_entities_by_dimension(0, 0, verts); CHECK_ERR(error);
    error = mbImpl->get_entities_by_dimension(0, 1, edges); CHECK_ERR(error);
    error = mbImpl->get_entities_by_dimension(0, 2, faces); CHECK_ERR(error);
    error = mbImpl->get_entities_by_dimension(0, 3, cells); CHECK_ERR(error);

    //Generate hierarchy
    error = refine_entities(mbImpl, level_degrees, num_levels, false);  CHECK_ERR(error);

    return MB_SUCCESS;
}

int main(int argc, char *argv[])
{
#ifdef MOAB_HAVE_MPI
    MPI_Init(&argc, &argv);

    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    ErrorCode result;
    if (argc ==1){
        result = test_1D();
        handle_error_code(result, number_tests_failed, number_tests_successful);
        std::cout<<"\n";

        result = test_2D();
        handle_error_code(result, number_tests_failed, number_tests_successful);
        std::cout<<"\n";

        result = test_3D();
        handle_error_code(result, number_tests_failed, number_tests_successful);
        std::cout<<"\n";
  }
    else if (argc == 2)
      {
        const char* filename = argv[1];
        int deg[4] = {2,3,3,2};
        int len = sizeof(deg) / sizeof(int);
        result = test_mesh(filename, deg, len);
        handle_error_code(result, number_tests_failed, number_tests_successful);
        std::cout<<"\n";

     /*   deg[0] = 3; deg[1] = 3;
        result = test_mesh(filename, deg, len);
        handle_error_code(result, number_tests_failed, number_tests_successful);
        std::cout<<"\n";*/

    /*    deg = 5;
        result = test_mesh(filename, &deg, len);
        handle_error_code(result, number_tests_failed, number_tests_successful);
        std::cout<<"\n";*/
      }

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif

  return number_tests_failed;
}

