/*This unit test is for the uniform refinement capability based on AHF datastructures*/
#include <iostream>
#include <string>
#include <sstream>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "../RefineMesh/moab/RefineSlabs.hpp"
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

double wtime() {
  double y = -1;
  struct timeval cur_time;
  gettimeofday(&cur_time, NULL);
  y = (double)(cur_time.tv_sec) + (double)(cur_time.tv_usec)*1.e-6;
  return (y);
}


ErrorCode create_mesh_from_coords(Interface *mbImpl, const double *coords, const size_t num_vtx, const int *conn, const size_t num_elems, EntityHandle *verts, EntityHandle *cells )
{
  ErrorCode error = MB_SUCCESS;
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
  return error;
}

// copied from urefine_mesh_tst.cpp
ErrorCode create_simple_mesh(Interface *mbImpl, EntityType type)
{
  ErrorCode error;
  if (type == MBEDGE)
  {
    return MB_FAILURE;
  }
  else if (type == MBTRI)
  {
    return MB_FAILURE;
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
    return MB_FAILURE;
  }
  else if (type == MBHEX)
  {

    const int which_mesh = 0;
    // number of hexes in each direction. 
    // number of vertices is this plus one
    const int xx(5), yy(7), zz(9); 

    switch (which_mesh) 
    {
      case 0:
      {
        std::cout<<"Creating structured 2x2x1 hex mesh."<<std::endl;
        // 2 x 2 x 1 structured mesh
        const double coords0[]= {0,0,0,
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
        const size_t num_vtx0 = sizeof(coords0)/sizeof(double)/3;

        const int conn0[]= {0,1,4,3,9,10,13,12,
                           1,2,5,4,10,11,14,13,
                           3,4,7,6,12,13,16,15,
                           4,5,8,7,13,14,17,16};
        const size_t num_elems0 = sizeof(conn0)/sizeof(int)/8;

        EntityHandle verts0[num_vtx0], cells0[num_elems0];
      
        error = create_mesh_from_coords(mbImpl, coords0, num_vtx0, conn0, num_elems0, verts0, cells0 ); CHECK_ERR(error);
        // todo set pillow_set of hexes
      }
      break;
      case 1:
      {
        std::cout<<"Creating structured " << xx << "x" << yy << "x" << zz << " hex mesh."<<std::endl;
        // vertices
        const int num_vtx1 = (xx+1)*(yy+1)*(zz+1);
        double *coords1 = new double[num_vtx1];
        size_t i  = 0;
        for ( double z = 0; z <= zz; ++z)
        {
          for ( double y = 0; y <= yy; ++y)
          {
            for ( double x = 0; x <= xx; ++x)
            {
              coords1[i++] = x;
              coords1[i++] = y;
              coords1[i++] = z;
            }
          }
        }

        // hexes
        const int num_elems1 = xx*yy*zz;
        int conn1[ 8 * num_elems1];
        i = 0;
        for ( double z = 0; z < zz; ++z)
        {
          for ( double y = 0; y < yy; ++y)
          {
            for ( double x = 0; x < xx; ++x)
            {
              const int hex_id = z * yy * xx + y * xx + x;
              assert( hex_id * 8 == i );
              if ( hex_id * 8 != i )
                return MB_FAILURE;

              // lowest corner
              const double cx = x;
              const double cy = y;
              const double cz = z;

              // lower quad of the hex
              const double c0 = cz * (yy+1) * (xx+1) + cy * (xx+1) + cx;
              const double c1 = c0 + 1;
              const double c2 = c1 + (xx+1);
              const double c3 = c2 - 1;

              // upper quad of the hex
              const double c4 = c0 + (yy+1) * (xx+1);
              const double c5 = c4 + 1;
              const double c6 = c5 + (xx+1);
              const double c7 = c6 - 1;

              conn1[i++] = c0;
              conn1[i++] = c1;
              conn1[i++] = c2;
              conn1[i++] = c3;
              conn1[i++] = c4;
              conn1[i++] = c5;
              conn1[i++] = c6;
              conn1[i++] = c7;
            }
          }
        }

        EntityHandle verts1[num_vtx1], cells1[num_elems1];
        error = create_mesh_from_coords(mbImpl, coords1, num_vtx1, conn1, num_elems1, verts1, cells1 ); CHECK_ERR(error);

        // todo create quads on the outside

        // todo set pillow_set of hexes
      }
      break;
    }

    // todo
    // add a variant that creates surface quads, marked as being owned by a surface
    // and a test that pillows (some) of them as if they were in the interior of the set

  }
  return MB_SUCCESS;
}

ErrorCode test_B()
{
  ErrorCode error;
  Core mb;
  Interface* mbImpl = &mb;

  // create mesh
  error = create_simple_mesh(mbImpl, MBHEX);
  if (error != MB_SUCCESS) return error;
  std::cout<<"Small simple hex mesh created successfully"<<std::endl;

  // refine it 
  RefineSlabs refine_slabs(&mb);
  RefineSlabs::Entities coarse_hexes, coarse_quads, fine_hexes, fine_quads;
  // todo, put some hexes into the set
  error = refine_slabs.refine_mesh(coarse_hexes, coarse_quads, fine_hexes, fine_quads); CHECK_ERR(error);


  return MB_SUCCESS;
}


ErrorCode test_A(const char* filename)
{
  Core mb;
  Interface* mbImpl = &mb;
  EntityHandle fileset;

  // create up some mesh
  ErrorCode error = MB_SUCCESS;
  error = mbImpl->create_meshset(moab::MESHSET_SET, fileset); CHECK_ERR(error);

  error = mbImpl->load_file(filename, &fileset);  CHECK_ERR(error);

  // refine it 
  RefineSlabs refine_slabs(&mb);
  RefineSlabs::Entities coarse_hexes, coarse_quads, fine_hexes, fine_quads;
  // todo, put some hexes into the set
  error = refine_slabs.refine_mesh(coarse_hexes, coarse_quads, fine_hexes, fine_quads); CHECK_ERR(error);


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

    std::cout << "Hello there slab refine" << std::endl;

    ErrorCode result = MB_SUCCESS;
    if (argc == 1)
    {
      std::cout << "argc1 Hello there slab refine" << std::endl;
      result = test_B();
      handle_error_code(result, number_tests_failed, number_tests_successful);
      std::cout<<"\n";

    }
    else if (argc == 2)
    {
      std::cout << "argc2 Hello there slab refine" << std::endl;
      const char* filename = argv[1];
      result = test_A(filename);
      handle_error_code(result, number_tests_failed, number_tests_successful);
      std::cout<<"\n";

    }
    else if (argc == 3)
    {
      std::cout << "argc3 Hello there slab refine" << std::endl;
    }
    else
    {
      std::cout << "argc?? Hello there slab refine" << std::endl;      
    }


#ifdef USE_MPI
    MPI_Finalize();
#endif

  return number_tests_failed;
}

