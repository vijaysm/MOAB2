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

void report_elements( EntityHandle* hexes, int num_hexes, EntityHandle* vertices, int num_vertices )
{
  std::cout << "Reporting mesh" << std::endl;

  std::cout << "Hex handles: ";
  for (int h = 0; h < num_hexes; ++h )
  {
    EntityHandle handle = hexes[ h ];
    std::cout << handle << ", ";
  }
  std::cout << std::endl;

  std::cout << "Vertex handles: ";
  for (int v = 0; v < num_vertices; ++v )
  {
    EntityHandle handle = vertices[v];
    std::cout << handle << ", ";
  }
  std::cout << std::endl;
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

// works for a structured mesh of size xx by yy by zz
int get_hex_id( int x, int y, int z,  int xx, int yy, int zz) 
{ 
  if ( x >= xx || y >= yy || z >= zz || x < 0 || y < 0 || z < 0 )
  {
    std::cout << "Error, hex indices [" << x << " , " << y << " , " << z << "] out of bounds [" << xx << " , " << yy << " , " << zz << "]" << std::endl;
    return 0;
  }
  return z * yy * xx + y * xx + x; 
}
// to get the handle, access cells[ get_hex_id( ) ]

// copied from urefine_mesh_tst.cpp
ErrorCode create_simple_mesh(Interface *mbImpl, EntityType type, int &which_mesh, int &x, int &y, int &z, EntityHandle *&hexes, int &num_hexes, EntityHandle *&vertices, int &num_vertices)
{
  // number of hexes (quads) in each direction. 
  // number of vertices is this plus one
  // const int xx(5), yy(7), zz(9); 

  num_hexes = 0;
  num_vertices = 0;

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
    which_mesh = 0;
    x = 4; y = 1; z = 0;
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

    which_mesh = 1;

    switch (which_mesh) 
    {
      case 0:
      {
        x = 2; y = 2; z =1;
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

        EntityHandle *verts0 = new EntityHandle[num_vtx0];
        EntityHandle *cells0 = new EntityHandle[num_elems0];
      
        error = create_mesh_from_coords(mbImpl, coords0, num_vtx0, conn0, num_elems0, verts0, cells0 ); CHECK_ERR(error);
        hexes = cells0; num_hexes = num_elems0;
        vertices = verts0; num_vertices = num_vtx0;
        // todo set pillow_set of hexes
      }
      break;
      case 1:
      {
        const int xx(5), yy(7), zz(9); 
        x = xx; y = yy; z = zz;
        std::cout<<"Creating structured " << xx << "x" << yy << "x" << zz << " hex mesh."<<std::endl;
        // vertices
        const int num_vtx1 = (xx+1)*(yy+1)*(zz+1);
        double *coords1 = new double[3*num_vtx1];
        size_t i  = 0;
        for ( double za = 0; za <= zz; ++za)
        {
          for ( double ya = 0; ya <= yy; ++ya)
          {
            for ( double xa = 0; xa <= xx; ++xa)
            {
              coords1[i++] = xa;
              coords1[i++] = ya;
              coords1[i++] = za;
            }
          }
        }

        // hexes
        const int num_elems1 = xx*yy*zz;
        int conn1[ 8 * num_elems1];
        i = 0;
        for ( double za = 0; za < zz; ++za)
        {
          for ( double ya = 0; ya < yy; ++ya)
          {
            for ( double xa = 0; xa < xx; ++xa)
            {
              const int hex_id = get_hex_id( xa,ya,za, xx,yy,zz );
              assert( hex_id * 8 == (int) i );
              if ( hex_id * 8 != (int) i )
                return MB_FAILURE;

              // lowest corner
              const double cx = xa;
              const double cy = ya;
              const double cz = za;

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
              assert( i < 8 * num_elems1 );
              conn1[i++] = c7;
            }
          }
        }

        EntityHandle *verts1 = new EntityHandle[num_vtx1];
        EntityHandle *cells1 = new EntityHandle[num_elems1];

        error = create_mesh_from_coords(mbImpl, coords1, num_vtx1, conn1, num_elems1, verts1, cells1 ); CHECK_ERR(error);
        
        hexes = cells1; num_hexes = num_elems1;
        vertices = verts1; num_vertices = num_vtx1;

        // todo create explicit moab quads on the outside faces
      }
      break;
    }

    // todo
    // add a variant that creates surface quads, marked as being owned by a surface
    // and a test that pillows (some) of them as if they were in the interior of the set

  }

  report_elements( hexes, num_hexes, vertices, num_vertices );
  return MB_SUCCESS;
}


void list_block( Interface *mb, Range &range, int val )
{
  std::cout << "sets with tag value " << val << " are\n";
  if (range.empty())
    std::cout << " EMPTY." << std::endl;
  else
  {
    mb->list_entities(range);
  }
}
void block_range_to_vector( Interface *mb, Range &range, RefineSlabs::Entities &entities )
{
  if (range.empty())
    return;
  EntityHandle block_handle = *range.begin();
  Range ent_range;
  mb->get_entities_by_handle( block_handle, ent_range );
  entities.assign( ent_range.begin(), ent_range.end() );
}


ErrorCode load_mesh(std::string fname, Interface *mb,
                    RefineSlabs::Entities &hexes, RefineSlabs::Entities &vertices, 
                    RefineSlabs::Entities &coarse_hexes, RefineSlabs::Entities &coarse_quads,
                    RefineSlabs::Entities &surface_nodes, RefineSlabs::Entities &curve_nodes, RefineSlabs::Entities &vertex_nodes,
                    bool &is_sweep_edge, RefineSlabs::NodePair &sweep_edge )
{
  // Load the mesh from vtk file
  using namespace std;

  cout << "Reading mesh file " << fname << std::endl;
  ErrorCode rval = mb->load_mesh(fname.c_str());MB_CHK_ERR(rval);

    // Get regions, by dimension, so we stay generic to entity type
  Range elems;
  rval = mb->get_entities_by_dimension(0, 3, elems);MB_CHK_ERR(rval);
  cout << "Read number of elements is " << elems.size() << endl;
  hexes.assign( elems.begin(), elems.end() );

  Range verts;
  rval = mb->get_entities_by_dimension(0, 0, verts);MB_CHK_ERR(rval);
  cout << "Read number of vertices is " << verts.size() << endl;
  vertices.assign( verts.begin(), verts.end() );

  // parse material sets, which we use to define our refinement sets
  // and the geometry boundary

  Tag mtag;
  rval=mb->tag_get_handle("MATERIAL_SET", mtag);MB_CHK_ERR(rval);

  int val=1; // block 1 is the coarse_hexes
  void *vals[] = {&val};  
  Range refine_hex_block;
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &mtag, vals, 1, refine_hex_block, 1, false);
  std::cout << "Refinement hexes, ";
  list_block( mb, refine_hex_block, val );
  block_range_to_vector( mb, refine_hex_block, coarse_hexes );

  val = 2; // block 2 is the coarse_quads
  Range refine_quad_block;
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &mtag, vals, 1, refine_quad_block, 1, false);
  std::cout << "Refinement quads, ";
  list_block( mb, refine_quad_block, val );
  block_range_to_vector( mb, refine_quad_block, coarse_quads );

  // // get all Dirichlet sets
  // dsets.clear();  
  // rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &dtag, 0, 1, dsets, 1, false);
  // std::cout << "sets with dirichlet set tags are\n";
  // mb->list_entities(dsets);

  // get a Dirichlet set == cubit's node set  
  Tag dtag;
  rval=mb->tag_get_handle("DIRICHLET_SET", dtag);MB_CHK_ERR(rval);

  // example
  // Range dsets;
  // rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &dtag, vals, 1, dsets, 1, false);
  // std::cout << "sets with Dirichlet set (node set) tag value " << val << " are\n";
  // mb->list_entities(dsets);

  
  val = 10; // block 10 is the nodes on vertices 
  Range vertex_nodes_range;
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &dtag, vals, 1, vertex_nodes_range, 1, false);
  std::cout << "Vertex nodes, ";
  list_block( mb, vertex_nodes_range, val );
  block_range_to_vector( mb, vertex_nodes_range, vertex_nodes );

  val = 11; // block 11 is the nodes on vertices 
  Range curve_nodes_range;
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &dtag, vals, 1, curve_nodes_range, 1, false);
  std::cout << "Curve nodes, ";
  list_block( mb, curve_nodes_range, val );
  block_range_to_vector( mb, curve_nodes_range, curve_nodes );

  val = 12; // block 12 is the nodes on vertices 
  Range surface_nodes_range;
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &dtag, vals, 1, surface_nodes_range, 1, false);
  std::cout << "Surface nodes, ";
  list_block( mb, surface_nodes_range, val );
  block_range_to_vector( mb, surface_nodes_range, surface_nodes );


  val = 20; // block 20 is the head node of a sweep-direction edge 
  is_sweep_edge = false;
  Range head_nodes_range;
  RefineSlabs::Entities head_nodes;
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &dtag, vals, 1, head_nodes_range, 1, false);
  std::cout << "Head node, ";
  list_block( mb, head_nodes_range, val );
  block_range_to_vector( mb, head_nodes_range, head_nodes );
  if ( head_nodes.size() == 1 )
  {
    EntityHandle head_node = head_nodes.front();

    val = 21; // block 21 is the tail node of a sweep-direction edge 
      Range tail_nodes_range;
    RefineSlabs::Entities tail_nodes;
    rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &dtag, vals, 1, tail_nodes_range, 1, false);
    std::cout << "Tail node, ";
    list_block( mb, tail_nodes_range, val );
    block_range_to_vector( mb, tail_nodes_range, tail_nodes );
    if ( tail_nodes.size() == 1 )
    {
      EntityHandle tail_node = tail_nodes.front();
      sweep_edge.first = head_node;
      sweep_edge.second = tail_node;
      is_sweep_edge = true;
    }
  }



  return rval;
}

// maybe put something like this into RefineSlap
ErrorCode test_C()
{
  ErrorCode error;
  Core mb;
  Interface* mbImpl = &mb;

  // load a mesh
  RefineSlabs::Entities hexes, vertices, coarse_hexes, coarse_quads, surface_nodes, curve_nodes, vertex_nodes;
  bool is_sweep_edge;
  RefineSlabs::NodePair sweep_edge;

  error = load_mesh("amesh.exo", mbImpl, hexes, vertices, coarse_hexes, coarse_quads, surface_nodes, curve_nodes, vertex_nodes, is_sweep_edge, sweep_edge ); 
  if (error != MB_SUCCESS) return error;
  std::cout<<"Loaded mesh successfully"<<std::endl;

  // refine it 
  RefineSlabs refine_slabs(&mb);
  RefineSlabs::Entities fine_hexes, fine_quads;

  refine_slabs.register_entity_handles( &hexes[0], hexes.size(), &vertices[0], vertices.size() );
  refine_slabs.register_boundary_nodes( surface_nodes, curve_nodes, vertex_nodes );

  // write it out again to be sure we read it correctly
  refine_slabs.write_file_slab( hexes, "entire_mesh" );

  error = refine_slabs.refine_mesh(coarse_hexes, coarse_quads, fine_hexes, fine_quads, is_sweep_edge, sweep_edge); CHECK_ERR(error);

  return MB_SUCCESS;
}

ErrorCode test_B()
{
  ErrorCode error;
  Core mb;
  Interface* mbImpl = &mb;

  // create mesh
  int which_mesh;
  int xx, yy, zz;
  EntityHandle *hexes(0), *vertices(0);
  int num_hexes(0), num_vertices(0);
  error = create_simple_mesh(mbImpl, MBHEX, which_mesh, xx, yy, zz, hexes, num_hexes, vertices, num_vertices);
  if (error != MB_SUCCESS) return error;
  std::cout<<"Small simple hex mesh created successfully"<<std::endl;

  // refine it 
  RefineSlabs refine_slabs(&mb);
  RefineSlabs::Entities coarse_hexes, coarse_quads, fine_hexes, fine_quads;
  // todo, put some hexes into the set

  if ( xx < 2 || yy < 2 || zz < 2)
  {
    std::cout << " mesh is only one element thick, not enough to do slab refinement." << std::endl;
    return MB_FAILURE;
  }

  const int which_test = 0;
  std::cout << "refinement test set " << which_test << std::endl;
  switch (which_test)
  {
    case 0:
    {
      // 3 x 3 x 3 in the lower corner, plus one for buffer elements around it
      const int xs(3), ys(3), zs(3);
      const int xo(1), yo(1), zo(1);
      std::cout << "  " << xs << " by " << ys << " by " << zs << " block offset from lower corner " << xo << " by " << yo << " by " << zo << std::endl;
      for ( int x = xo; x < xs+xo; ++x)
        for ( int y = yo; y < ys+yo; ++y )
          for ( int z = zo; z < zs+zo; ++z)
          {
            int hex_id = get_hex_id( x,y,z, xx,yy,zz );
            EntityHandle hex_handle = hexes[ hex_id ];
            coarse_hexes.push_back( hex_handle );
          }
    }
    break;
    case 1:
    {
    }
    break;
    default:;
  }  
  refine_slabs.register_entity_handles( hexes, num_hexes, vertices, num_vertices );
  RefineSlabs::Entities hexvec(hexes, hexes + num_hexes);
  refine_slabs.write_file_slab( hexvec, "entire_mesh" );

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
      // result = test_A();
      // result = test_B();
      result = test_C();
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

