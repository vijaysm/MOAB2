
#include <iostream>
#include "moab/Interface.hpp"
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "InitCGMA.hpp"
#include "GeometryQueryTool.hpp"
#include "moab/MeshTopoUtil.hpp"
using namespace moab;

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)


#ifdef MESHDIR
static const char input_cube[] = STRINGIFY(MESHDIR) "/io/cube.sat";
#else
static const char input_cube[] = "/io/cube.sat";
#endif

// Function used to load the test file
void read_file( Interface* moab, const char* input_file );

// List of tests in this file
void cube_verts_connectivity_test();
void cube_tris_connectivity_test();

//Function used to match triangle connectivity and verts 
//void match_tri_connectivity( Range connectivity, 
//                             std::vector<EntityHandle> &reference_verts);


int main(int /* argc */, char** /* argv */)
{
  int result = 0;
 
  result += RUN_TEST(cube_verts_connectivity_test);
  result += RUN_TEST(cube_tris_connectivity_test);
 
  return result;
}



void read_file( Interface* moab, const char* input_file )
{
  InitCGMA::initialize_cgma();
  GeometryQueryTool::instance()->delete_geometry();

  ErrorCode rval = moab->load_file( input_file );
  CHECK_ERR(rval);
}

void cube_verts_connectivity_test()
{

  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube);

  //Get all vertex handles from the mesh
  Range verts;
  rval = mb->get_entities_by_type( 0, MBVERTEX, verts);
  CHECK_ERR(rval);

  //Check that each vertex connects to less than 4 triangles and no more than 6

  for(Range::const_iterator i = verts.begin(); i!=verts.end(); i++)
    {
      std::vector<EntityHandle> adj_tris;
      rval = mb->get_adjacencies( &(*i), 1, 2, false, adj_tris );
      CHECK_ERR(rval);

      int adj_size = adj_tris.size();
      CHECK( adj_size >= 4 && adj_size <= 6);
    }
    
}

void cube_tris_connectivity_test()
{
    ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube);

  //Get triangles from the mesh
  Range tris;
  rval = mb->get_entities_by_type( 0, MBTRI, tris);
  CHECK_ERR(rval);


  for(Range::const_iterator i = tris.begin()+1; i!=tris.end(); i++)
    {
      Range adj_tris;
      moab::MeshTopoUtil mu(mb);
      rval = mu.get_bridge_adjacencies( *i, 1, 2, adj_tris);
      CHECK_ERR(rval);
      int number_of_adj_tris=adj_tris.size();      
      CHECK_EQUAL( 3, number_of_adj_tris);
    }

}

