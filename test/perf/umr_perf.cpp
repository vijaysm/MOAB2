/*This function profiles the performance of the AHF datastructure */
#include <iostream>
#include <assert.h>
#include <time.h>
#include <vector>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/HalfFacetRep.hpp"
#include "RefineMesh/moab/NestedRefine.hpp"
#include "../TestUtil.hpp"
#include <sys/time.h>

#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab/FileOptions.hpp"
#include "MBTagConventions.hpp"
#include "moab_mpi.h"
#endif

using namespace moab;

#ifdef MOAB_HAVE_MPI
std::string read_options;
#endif

int number_tests_successful = 0;
int number_tests_failed = 0;

struct mesh_mem
{

 /* unsigned long long total_storage;
  unsigned long long amortized_total_storage;
  unsigned long long entity_storage;
  unsigned long long amortized_entity_storage;
  unsigned long long adjacency_storage;
  unsigned long long amortized_adjacency_storage;
  unsigned long long tag_storage;
  unsigned long long amortized_tag_storage;*/

  //Total storage
  unsigned long long total_storage;
  unsigned long long total_amortized_storage;

  //Tag storage
  unsigned long long tag_storage;
  unsigned long long amortized_tag_storage;

  //Vertex storage
   unsigned long long vertex_storage;
   unsigned long long amortized_vertex_storage;

   //Entity storage
   unsigned long long entity_storage;
   unsigned long long amortized_entity_storage;

};

enum OUTTYPE{
  TIME=0,
  MEM,
  BOTH
};

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

double wtime() {
  double y = -1;
  struct timeval cur_time;
  gettimeofday(&cur_time, NULL);
  y = (double)(cur_time.tv_sec) + (double)(cur_time.tv_usec)*1.e-6;
  return (y);
}

ErrorCode umr_perf_test(Core *mb, int *level_degrees, int num_levels, OUTTYPE output)
{
  ErrorCode error;
  Interface* mbImpl = mb;

  mesh_mem *umem = new mesh_mem[num_levels+2];
  double time_start, time_avg, time_total;

  //Create ranges of entities in the initial mesh
  Range inverts, inedges, infaces, incells;
  error = mbImpl->get_entities_by_dimension( 0, 0, inverts);
  error = mbImpl->get_entities_by_dimension( 0, 1, inedges);
  error = mbImpl->get_entities_by_dimension( 0, 2, infaces);
  error = mbImpl->get_entities_by_dimension( 0, 3, incells);

  Range init_ents;
  int dim=0;
  if (inedges.size())
    {
      dim = 1;
      init_ents = inedges;
    }
  else if (infaces.size())
    {
      dim = 2;
      init_ents = infaces;
    }
  else if (incells.size())
    {
      dim = 3;
      init_ents = incells;
    }

  std::cout<<std::endl;
  std::cout<<"Initial Mesh Size: "<<"NV = "<<inverts.size()<<", NE = "<<init_ents.size()<<std::endl;

  if (output == MEM || output == BOTH){
      //Storage Costs before mesh hierarchy generation
 /*     std::cout<<std::endl;
      unsigned long long sTotS, sTAS, sES, sAES, sAS, sAAS, sTS, sATS;
      sTotS = sTAS = sES = sAES = sAS = sAAS = sTS = sATS = 0;
      mbImpl->estimated_memory_use(NULL, 0, &sTotS, &sTAS, &sES, &sAES, &sAS, &sAAS, NULL, 0, &sTS, &sATS);

      umem[0].total_storage = sTotS;
      umem[0].amortized_total_storage = sTAS;
      umem[0].entity_storage = sES;
      umem[0].amortized_entity_storage = sAES;
      umem[0].adjacency_storage = sAS;
      umem[0].amortized_adjacency_storage = sAAS;
      umem[0].tag_storage = sTS;
      umem[0].amortized_tag_storage = sATS;*/

      unsigned long long vTotS, vTAS, vES, vAES, vAS, vAAS, vTS, vATS;
      vTotS = vTAS = vES = vAES = vAS = vAAS = vTS = vATS = 0;
      mbImpl->estimated_memory_use(inverts, &vTotS, &vTAS, &vES, &vAES, &vAS, &vAAS, NULL, 0, &vTS, &vATS);

      unsigned long long eTotS, eTAS, eES, eAES, eAS, eAAS, eTS, eATS;
      eTotS = eTAS = eES = eAES = eAS = eAAS = eTS = eATS = 0;
      mbImpl->estimated_memory_use(init_ents, &eTotS, &eTAS, &eES, &eAES, &eAS, &eAAS, NULL, 0, &eTS, &eATS);

      umem[0].total_storage = vTotS+eTotS;
      umem[0].vertex_storage = vES;
      umem[0].entity_storage = eES;
      umem[0].tag_storage = vTS+eTS;

      std::cout<<"MEMORY STORAGE:: Initial Mesh"<<std::endl;
      std::cout<<std::endl;
      std::cout<<"Total storage = "<<umem[0].total_storage<<std::endl;
      std::cout<<"Vertex storage = "<<umem[0].vertex_storage<<std::endl;
      std::cout<<"Entity storage = "<<umem[0].entity_storage<<std::endl;
      std::cout<<"Tag storage = "<< umem[0].tag_storage<<std::endl;


  /*    std::cout<<"Total storage = "<<umem[0].total_storage<<std::endl;
      std::cout<<"Total amortized storage = "<< umem[0].amortized_total_storage<<std::endl;
      std::cout<<"Entity storage = "<<umem[0].entity_storage<<std::endl;
      std::cout<<"Amortized entity storage = "<<umem[0].amortized_entity_storage<<std::endl;
      std::cout<<"Adjacency storage = "<< umem[0].adjacency_storage <<std::endl;
      std::cout<<"Amortized adjacency storage = "<<umem[0].amortized_adjacency_storage <<std::endl;
      std::cout<<"Tag storage = "<< umem[0].tag_storage<<std::endl;
      std::cout<<"Amortized tag storage = "<<umem[0].amortized_tag_storage <<std::endl;*/
      std::cout<<std::endl;
    }


  //Create an hm object and generate the hierarchy
  std::cout<<"Creating a hm object"<<std::endl;
  NestedRefine uref(mb);
  std::vector<EntityHandle> set;

  std::cout<<"Starting hierarchy generation"<<std::endl;
  time_start = wtime();
  error = uref.generate_mesh_hierarchy(num_levels, level_degrees, set);
  CHECK_ERR(error);

  time_total = wtime() - time_start;
  std::cout<<"Finished hierarchy generation"<<std::endl;

  if (output == TIME || output == BOTH){
      std::cout<<"Total time in secs:: generate mesh hierarchy:: L = "<<num_levels<<" :: "<<time_total<<std::endl;
      std::cout<<std::endl;
    }


  //Loop over each mesh level and check its topological properties
  for (int l=0; l<num_levels; l++)
    {
      //Get the current mesh level using its meshset
      Range verts, ents;
      error = mbImpl->get_entities_by_type(set[l+1], MBVERTEX, verts);
      CHECK_ERR(error);
      error = mbImpl->get_entities_by_dimension(set[l+1], dim, ents);
      CHECK_ERR(error);

      std::cout<<"Mesh size for level "<<l+1<<" :: deg = "<<level_degrees[l]<<" :: NV = "<<verts.size()<<", NE = "<<ents.size()<<std::endl;

      if (output == MEM || output == BOTH){
          //Storage Costs
          std::cout<<std::endl;
          unsigned long long vTotS, vTAS, vES, vAES, vAS, vAAS, vTS, vATS;
          vTotS = vTAS = vES = vAES = vAS = vAAS = vTS = vATS = 0;
          mbImpl->estimated_memory_use(verts, &vTotS, &vTAS, &vES, &vAES, &vAS, &vAAS, NULL, 0, &vTS, &vATS);
          unsigned long long eTotS, eTAS, eES, eAES, eAS, eAAS, eTS, eATS;
          eTotS = eTAS = eES = eAES = eAS = eAAS = eTS = eATS = 0;
          mbImpl->estimated_memory_use(ents, &eTotS, &eTAS, &eES, &eAES, &eAS, &eAAS, NULL, 0, &eTS, &eATS);

          umem[l+1].total_storage = vTotS+eTotS;
          umem[l+1].vertex_storage = vES;
          umem[l+1].entity_storage = eES;
          umem[l+1].tag_storage = vTS+eTS;

         /* umem[l+1].total_storage = vTotS+eTotS;
          umem[l+1].amortized_total_storage = vTAS+eTAS;
          umem[l+1].entity_storage = vES+eES;
          umem[l+1].amortized_entity_storage = vAES+eAES;
          umem[l+1].adjacency_storage = vAS+eAS;
          umem[l+1].amortized_adjacency_storage = vAAS+eAAS;
          umem[l+1].tag_storage = vTS+eTS;
          umem[l+1].amortized_tag_storage = vATS+eATS;*/

          std::cout<<"MEMORY STORAGE:: Mesh level "<<l+1<<std::endl;
          std::cout<<std::endl;
          std::cout<<"Total storage = "<<umem[l+1].total_storage<<std::endl;
          std::cout<<"Vertex storage = "<<umem[l+1].vertex_storage<<std::endl;
          std::cout<<"Entity storage = "<<umem[l+1].entity_storage<<std::endl;
          std::cout<<"Tag storage = "<< umem[l+1].tag_storage<<std::endl;

        /*  std::cout<<"Total storage = "<<umem[l+1].total_storage<<std::endl;
          std::cout<<"Total amortized storage = "<< umem[l+1].amortized_total_storage<<std::endl;
          std::cout<<"Entity storage = "<<umem[l+1].entity_storage<<std::endl;
          std::cout<<"Amortized entity storage = "<<umem[l+1].amortized_entity_storage<<std::endl;
          std::cout<<"Adjacency storage = "<< umem[l+1].adjacency_storage <<std::endl;
          std::cout<<"Amortized adjacency storage = "<<umem[l+1].amortized_adjacency_storage <<std::endl;
          std::cout<<"Tag storage = "<< umem[l+1].tag_storage<<std::endl;
          std::cout<<"Amortized tag storage = "<<umem[l+1].amortized_tag_storage <<std::endl;*/
          std::cout<<std::endl;
        }

      if (output == BOTH){
          //Loop over all vertices and get their coordinates
          time_start = wtime();
          for (Range::iterator i = verts.begin(); i != verts.end(); ++i)
            {
              double coords[3];
              EntityHandle vid = *i;
              error = uref.get_coordinates(&vid, 1, (int)l, &coords[0]);
              CHECK_ERR(error);
            }
          time_total = wtime() - time_start;
          time_avg = time_total/(double)verts.size();

          std::cout<<"Class NR :: OPERATION:: get_coordinates"<<" :: Time_total = "<<time_total<<" :: Time_avg = "<<time_avg<<std::endl;

          time_start = wtime();
          for (Range::iterator i = verts.begin(); i != verts.end(); ++i)
            {
              double coords[3];
              EntityHandle vid = *i;
              error = mbImpl->get_coords(&vid, 1, &coords[0]);
              CHECK_ERR(error);
            }
          time_total = wtime() - time_start;
          time_avg = time_total/(double)verts.size();

          std::cout<<"Class Core :: OPERATION:: get_coordinates"<<" :: Time_total = "<<time_total<<" :: Time_avg = "<<time_avg<<std::endl;


          //Loop over all entities and get their connectivity
          time_start = wtime();
          for (Range::iterator i = ents.begin(); i != ents.end(); ++i)
            {
              std::vector<EntityHandle> conn;
              error = uref.get_connectivity(*i, l, conn);
              CHECK_ERR(error);
            }
          time_total = wtime() - time_start;
          time_avg = time_total/(double)ents.size();
          std::cout<<std::endl;
          std::cout<<"Class NR :: OPERATION:: get_connectivity"<<" :: Time_total = "<<time_total<<" :: Time_avg = "<<time_avg<<std::endl;

          time_start = wtime();
          for (Range::iterator i = ents.begin(); i != ents.end(); ++i)
            {
              std::vector<EntityHandle> conn;
              error = mbImpl->get_connectivity(&*i, 1, conn);
              CHECK_ERR(error);
            }
          time_total = wtime() - time_start;
          time_avg = time_total/(double)ents.size();

          std::cout<<"Class Core :: OPERATION:: get_connectivity"<<" :: Time_total = "<<time_total<<" :: Time_avg = "<<time_avg<<std::endl;
          std::cout<<std::endl;
        }
    }


 /* if (output == MEM || output == BOTH){

      unsigned long long sTotS, sTAS, sES, sAES, sAS, sAAS, sTS, sATS;
      sTotS = sTAS = sES = sAES = sAS = sAAS = sTS = sATS = 0;
      mbImpl->estimated_memory_use(NULL, 0, &sTotS, &sTAS, &sES, &sAES, &sAS, &sAAS, NULL, 0, &sTS, &sATS);

      umem[num_levels+1].total_storage = sTotS;
      umem[num_levels+1].amortized_total_storage = sTAS;
      umem[num_levels+1].entity_storage = sES;
      umem[num_levels+1].amortized_entity_storage = sAES;
      umem[num_levels+1].adjacency_storage = sAS;
      umem[num_levels+1].amortized_adjacency_storage = sAAS;
      umem[num_levels+1].tag_storage = sTS;
      umem[num_levels+1].amortized_tag_storage = sATS;

      std::cout<<"MEMORY STORAGE:: WHOLE MESH HIERARCHY"<<std::endl;
      std::cout<<std::endl;
      std::cout<<"Total storage = "<<umem[num_levels+1].total_storage<<std::endl;
      std::cout<<"Total amortized storage = "<< umem[num_levels+1].amortized_total_storage<<std::endl;
      std::cout<<"Entity storage = "<<umem[num_levels+1].entity_storage<<std::endl;
      std::cout<<"Amortized entity storage = "<<umem[num_levels+1].amortized_entity_storage<<std::endl;
      std::cout<<"Adjacency storage = "<< umem[num_levels+1].adjacency_storage <<std::endl;
      std::cout<<"Amortized adjacency storage = "<<umem[num_levels+1].amortized_adjacency_storage <<std::endl;
      std::cout<<"Tag storage = "<< umem[num_levels+1].tag_storage<<std::endl;
      std::cout<<"Amortized tag storage = "<<umem[num_levels+1].amortized_tag_storage <<std::endl;
      std::cout<<std::endl;
    }*/

    return MB_SUCCESS;
}

ErrorCode create_simple_mesh(Core *mb, EntityType type)
{
  ErrorCode error;
  Interface* mbImpl = mb;
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
          error = mbImpl->create_vertex(coords+3*i, verts[i]);
          if (error != MB_SUCCESS) return error;
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[2];
          c[0] = verts[conn[2*i]]; c[1] = verts[conn[2*i+1]];

          error = mbImpl->create_element(MBEDGE, c, 2, edges[i]);
          if (error != MB_SUCCESS) return error;
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
              error = mbImpl->create_vertex(coords+3*i, verts[i]);
              if (error != MB_SUCCESS) return error;
            }

          for (size_t i=0; i< num_elems; ++i)
            {
              EntityHandle c[4];
              for (int j=0; j<4; j++)
                c[j] = verts[conn[4*i+j]];

              error = mbImpl->create_element(MBQUAD, c, 4, faces[i]);
              if (error != MB_SUCCESS) return error;

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
          error = mbImpl->create_vertex(coords+3*i, verts[i]);
          if (error != MB_SUCCESS) return error;
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[4];
          for (int j=0; j<4; j++)
            c[j] = verts[conn[4*i+j]];

          error = mbImpl->create_element(MBTET, c, 4, cells[i]);
          if (error != MB_SUCCESS) return error;
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
        error = mbImpl->create_vertex(coords+3*i, verts[i]);
        if (error != MB_SUCCESS) return error;
      }

    for (size_t i=0; i< num_elems; ++i)
      {
        EntityHandle c[8];
        for (int j=0; j<8; j++)
          c[j] = verts[conn[8*i+j]];

        error = mbImpl->create_element(MBHEX, c, 8, cells[i]);
        if (error != MB_SUCCESS) return error;
      }
    }
  return MB_SUCCESS;
}

ErrorCode test_mesh(EntityType type, int *level_degrees, int num_level)
{
  ErrorCode error;
  Core mb;

  error = create_simple_mesh(&mb, type);
  if (error != MB_SUCCESS) return error;

  OUTTYPE output = BOTH;
  error = umr_perf_test(&mb, level_degrees, num_level, output);
  if (error != MB_SUCCESS) return error;

  return MB_SUCCESS;
}


ErrorCode test_1D()
{
  ErrorCode error;
  std::cout<<std::endl;
  std::cout<<"EntityType = MBEDGE"<<std::endl;

  std::cout<<"Deg = 2"<<std::endl;
  int deg[7] = {5,5,2,2,2,2,2};
  int len = sizeof(deg) / sizeof(int);
  error = test_mesh(MBEDGE, deg, len);
  if (error != MB_SUCCESS) return error;

  std::cout<<std::endl;
  std::cout<<"Deg = 3"<<std::endl;
  deg[4] = 3; deg[5] = 3; deg[6] = 3;
  error = test_mesh(MBEDGE, deg, len);
  if (error != MB_SUCCESS) return error;

  std::cout<<std::endl;
  std::cout<<"Deg = 5"<<std::endl;
  deg[4] = 5; deg[5] = 5; deg[6] = 5;
  error = test_mesh(MBEDGE, deg, len);
  if (error != MB_SUCCESS) return error;

  return MB_SUCCESS;
}

ErrorCode test_2D()
{
  ErrorCode error;
  EntityType type = MBTRI;
  std::cout<<std::endl;
  std::cout<<"EntityType = MBTRI"<<std::endl;

  std::cout<<"Deg = 2"<<std::endl;
  int deg[5] = {5,2,2,2,2};
  int len = sizeof(deg) / sizeof(int);
  error = test_mesh(type, deg, len);
  if (error != MB_SUCCESS) return error;

  std::cout<<std::endl;
  std::cout<<"Deg = 3"<<std::endl;
  deg[2] = 3; deg[3] = 3; deg[4] = 3;
  error = test_mesh(type, deg, len);
  if (error != MB_SUCCESS) return error;

  std::cout<<std::endl;
  std::cout<<"Deg = 5"<<std::endl;
  deg[2] = 5; deg[3] = 5; deg[4] = 5;
  error = test_mesh(type, deg, len);
  if (error != MB_SUCCESS) return error;

  type = MBQUAD;
  std::cout<<std::endl;
  std::cout<<"EntityType = MBQUAD"<<std::endl;

  std::cout<<"Deg = 2"<<std::endl;
  deg[2] = 2; deg[3] = 2; deg[4] = 2;
  error = test_mesh(type, deg, len);
  if (error != MB_SUCCESS) return error;

  std::cout<<std::endl;
  std::cout<<"Deg = 3"<<std::endl;
  deg[2] = 3; deg[3] = 3; deg[4] = 3;
  error = test_mesh(type, deg, len);
  if (error != MB_SUCCESS) return error;

  std::cout<<std::endl;
  std::cout<<"Deg = 5"<<std::endl;
  deg[2] = 5; deg[3] = 5; deg[4] = 5;
  error = test_mesh(type, deg, len);
  if (error != MB_SUCCESS) return error;

  return MB_SUCCESS;
}

ErrorCode test_3D()
{
  ErrorCode error;
  EntityType type = MBTET;
  std::cout<<std::endl;
  std::cout<<"EntityType = MBTET"<<std::endl;

  std::cout<<"Deg = 2"<<std::endl;
  int deg[6] = {2,2,2,2,2,2};
  int len = sizeof(deg) / sizeof(int);
  error = test_mesh(type, deg, len);
  if (error != MB_SUCCESS) return error;

  std::cout<<std::endl;
  std::cout<<"Deg = 3"<<std::endl;
  deg[3] = 3; deg[4] = 3; deg[5] = 3;
  error = test_mesh(type, deg, len);
  if (error != MB_SUCCESS) return error;

  type = MBHEX;
  std::cout<<std::endl;
  std::cout<<"EntityType = MBHEX"<<std::endl;

  std::cout<<"Deg = 2"<<std::endl;
  deg[3] = 2; deg[4] = 2; deg[5] = 2;
  error = test_mesh(type, deg, len);
  if (error != MB_SUCCESS) return error;

  std::cout<<std::endl;
  std::cout<<"Deg = 3"<<std::endl;
  deg[3] = 3; deg[4] = 3; deg[5] = 3;
  error = test_mesh(type, deg, len);
  if (error != MB_SUCCESS) return error;

  return MB_SUCCESS;
}


ErrorCode perf_inmesh(const char* filename, int *level_degrees, int num_levels, OUTTYPE output)
{
  ErrorCode error;
  Core mb;
  Interface* mbImpl = &mb;

#ifdef MOAB_HAVE_MPI
    int procs = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    if (procs > 1){
    read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;";

    error = mbImpl->load_file(filename, 0, read_options.c_str());
    CHECK_ERR(error);
    }
    else if (procs == 1) {
#endif
    error = mbImpl->load_file(filename);
    CHECK_ERR(error);
#ifdef MOAB_HAVE_MPI
    }
#endif
  //  OUTTYPE output = MEM;
    error = umr_perf_test(&mb, level_degrees, num_levels, output);
    if (error != MB_SUCCESS) return error;

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

#ifdef MOAB_HAVE_MPI
    if (rank == 0)
        std::cout<<" para_umr_perf: ";
#else
    std::cout<<"umr_perf:";
#endif

    if (argc==1)
    {
        ErrorCode result;

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

    else if (argc ==2)
      {
        const char *filename = argv[1];
        ErrorCode result;

        OUTTYPE output = MEM;

        if (output == MEM){
            int deg[3] = {2,2,2};
            int len = sizeof(deg) / sizeof(int);
            result = perf_inmesh(filename, deg, len, output);
            handle_error_code(result, number_tests_failed, number_tests_successful);
            std::cout<<"\n";

            deg[0] = 3; deg[1] = 3; deg[2] = 3;
            result = perf_inmesh(filename, deg, len, output);
            handle_error_code(result, number_tests_failed, number_tests_successful);
            std::cout<<"\n";

          /*  deg[0] = 5; deg[1] = 5; deg[2] = 5;
            result = perf_inmesh(filename, deg, len, output);
            handle_error_code(result, number_tests_failed, number_tests_successful);
            std::cout<<"\n";*/
          }
        else if (output == TIME)
          {
            for (int L=0; L<3; L++)
              {
                int *level_degrees = new int[L+1];
                for (int i=0; i<L+1; i++)
                  {
                    level_degrees[i] = 3;
                  }

                result = perf_inmesh(filename, level_degrees, L+1, output);
                handle_error_code(result, number_tests_failed, number_tests_successful);
                std::cout<<"\n";

                delete [] level_degrees;
              }
          }

      }
    else {
            std::cerr << "Usage: " << argv[0] << " [filename]" << std::endl;
            return 1;
    }

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif

    return number_tests_failed;
}
