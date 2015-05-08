
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
#ifdef HAVE_OCC_STEP
static const char input_cube[] = STRINGIFY(MESHDIR) "/io/cube.stp";
#else
static const char input_cube[] = STRINGIFY(MESHDIR) "/io/cube.sat";
#endif
#else
#ifdef HAVE_OCC_STEP
static const char input_cube[] = "cube.stp";
#else
static const char input_cube[] = "cube.sat";
#endif
#endif


// Function used to load the test file
void read_file( Interface* moab, const char* input_file );

// List of tests in this file
void cube_verts_connectivity_test();
void cube_tris_connectivity_test();
void cube_tri_curve_coincidence_test();
void cube_edge_adjacencies_test();
void cube_tri_vertex_test();

//Other functions
void match_tri_edges_w_curve( Range tri_edges, Range curves );

int main(int /* argc */, char** /* argv */)
{
  int result = 0;
 
  result += RUN_TEST(cube_verts_connectivity_test);
  result += RUN_TEST(cube_tris_connectivity_test);
  result += RUN_TEST(cube_tri_curve_coincidence_test);
  result += RUN_TEST(cube_tri_vertex_test);
 
  return result;
}



void read_file( Interface* moab, const char* input_file )
{
  InitCGMA::initialize_cgma();
  GeometryQueryTool::instance()->delete_geometry();

  ErrorCode rval = moab->load_file( input_file );
  CHECK_ERR(rval);
}

// Checks the adjacency of each vertex entity in a simple cube file load
// to make sure it isn't adjacent to too many or too few triangles.
void cube_verts_connectivity_test()
{

  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );

  //Get all vertex handles from the mesh
  Range verts;
  rval = mb->get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);

  //Check that each vertex connects to at least 4 and no more than 6 triangles
  for(Range::const_iterator i = verts.begin(); i!=verts.end(); ++i)
    {
      std::vector<EntityHandle> adj_tris;
      rval = mb->get_adjacencies( &(*i), 1, 2, false, adj_tris );
      CHECK_ERR(rval);

      int adj_size = adj_tris.size();
      CHECK( adj_size >= 4 && adj_size <= 6 );
    }
    
}

// Check that each triangle in the mesh is adjacent to
// exactly three other triangles
void cube_tris_connectivity_test()
{
  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );

  //Get triangles from the mesh
  Range tris;
  rval = mb->get_entities_by_type( 0, MBTRI, tris );
  CHECK_ERR(rval);

  int expected_num_of_adj_tris = 3;

  for(Range::const_iterator i = tris.begin()+1; i!=tris.end(); ++i)
    {
      Range adj_tris;
      moab::MeshTopoUtil mu(mb);
      //Use Triangle edges to get all adjacent triangles
      rval = mu.get_bridge_adjacencies( *i, 1, 2, adj_tris );
      CHECK_ERR(rval);
      CHECK_EQUAL( expected_num_of_adj_tris, (int)adj_tris.size() );
      
      //Check that the entities we found from bridge_adjacencies
      //are triangles
      Range adj_tri_test = adj_tris.subset_by_type( MBTRI );
      CHECK_EQUAL( (int)adj_tris.size(), (int) adj_tri_test.size() );
    
    }

}


// Takes triangle edges and makes sure they match the EntityHandles of 
// curves in the case of a cube mesh
void cube_tri_curve_coincidence_test()
{

  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );

  //Get curves from the mesh
  Range curves;
  rval = mb->get_entities_by_type( 0, MBEDGE, curves );
  CHECK_ERR(rval);
  curves.print();

  //Get triangles from the mesh
  Range tris;
  rval = mb->get_entities_by_type( 0, MBTRI, tris );
  CHECK_ERR(rval);

  for(Range::const_iterator i=tris.begin(); i!=tris.end(); ++i)
    {
      //Get the any curve edges that are a part of the triangle
      Range tri_edges;
      rval = mb->get_adjacencies( &(*i), 1, 1, false, tri_edges );
      CHECK_ERR(rval);
      //Check that we've retrieved two edges from get_adjacencies
      //For a this file (cube), each triangle should have two curve
      //edges
      int num_of_tri_edges = tri_edges.size();
      CHECK_EQUAL( 2, num_of_tri_edges );
      match_tri_edges_w_curve( tri_edges, curves );
      CHECK_ERR(rval);
      }
}

void match_tri_edges_w_curve( Range tri_edges, Range curves )
{
  int match_counter=0;
  int num_of_tri_edges = tri_edges.size();
  CHECK(num_of_tri_edges);
  for(Range::const_iterator i=tri_edges.begin(); i!=tri_edges.end(); ++i)
    {
      for(Range::const_iterator j=curves.begin(); j!=curves.end(); ++j)
	{
          // If the edge handle matches a curve handle, increment the number
          // matches
          if( *i  == *j  ) match_counter++;
	}
    }
  //Make sure that each edge returned from triangle edges
  //has been matched to a curve
  CHECK_EQUAL( num_of_tri_edges, match_counter );
} 

// Ensures that each triangle edge is adjacent to no more than
// two triangles.
void cube_edge_adjacencies_test()
{
  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );

  //Get the curves 
  Range curves;
  rval = mb->get_entities_by_type( 0, MBEDGE, curves );
  CHECK_ERR(rval);

  for(Range::const_iterator i=curves.begin(); i!=curves.end(); ++i)
    {
      //Get triangle adjacent to each edge
      Range adj_tris;
      rval = mb->get_adjacencies( &(*i), 1, 2, false, adj_tris );
      CHECK_ERR(rval);
      
      int num_adj_tris = adj_tris.size();
      //Ensure that no edge is adjacent to more than two triangles
      CHECK( num_adj_tris <= 2 );
    }

}

// Checks, for each triangle, that none of the verices are the same
void cube_tri_vertex_test()
{
  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );
 
  //Get all triangles
  Range tris;
  rval = mb->get_entities_by_type( 0, MBTRI, tris );
  CHECK_ERR(rval);

  for(Range::const_iterator i=tris.begin(); i!=tris.end(); ++i)
    {
      //Get all triangle vertices
      Range verts;
      rval = mb->get_connectivity( &(*i), 1, verts );
      CHECK_ERR(rval);
      //Make sure that each vertex making up
      //the triangle is different
      int number_of_verts = verts.size();
      CHECK( 3 == number_of_verts );
      CHECK( verts[0]!=verts[1] );
      CHECK( verts[1]!=verts[2] );      
      CHECK( verts[2]!=verts[0] );
    } 
}
