
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
void read_cube_verts_test();
void read_cube_curves_test();
void read_cube_tris_test();
void read_cube_surfs_test();
void read_cube_vols_test();
void read_cube_vertex_pos_test();
void delete_mesh_test();


int main(int /* argc */, char** /* argv */)
{
  int result = 0;
 
  result += RUN_TEST( read_cube_verts_test );  
  result += RUN_TEST( read_cube_curves_test );  
  result += RUN_TEST( read_cube_tris_test );  
  result += RUN_TEST( read_cube_surfs_test );  
  result += RUN_TEST( read_cube_vols_test );  
  result += RUN_TEST( read_cube_vertex_pos_test );  
 
  return result;
}



void read_file( Interface* moab, const char* input_file )
{
  InitCGMA::initialize_cgma();
  GeometryQueryTool::instance()->delete_geometry();

  ErrorCode rval = moab->load_file( input_file );
  CHECK_ERR(rval);
}

void read_cube_verts_test()
{
  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );
   
  int number_of_vertices;
  rval = mb->get_number_entities_by_type(0, MBVERTEX, number_of_vertices);
  CHECK_ERR(rval);
 
  CHECK_EQUAL( 8, number_of_vertices);
}


void read_cube_tris_test()
{
  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );
   
  int number_of_tris;

  rval = mb->get_number_entities_by_type(0, MBTRI , number_of_tris);
  std::cout << "Number of Triangles = " << number_of_tris << std::endl;
  CHECK_ERR(rval);

  CHECK_EQUAL( 12, number_of_tris);  

}

void read_cube_curves_test()
{
  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );
   
  Tag geom_tag;

  rval = mb->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1,
				MB_TYPE_INTEGER, geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  CHECK_ERR(rval);
  
  Range curves;
  int dim = 1;
  void *val[] = {&dim};
  int number_of_curves;
  rval = mb->get_number_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
	  					    val, 1, number_of_curves );
  CHECK_ERR(rval);

  CHECK_EQUAL( 12, number_of_curves);  

} 

void read_cube_surfs_test()
{
  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );
  
  //Get geometry tag for pulling curve data from the mesh 
  Tag geom_tag;
  rval = mb->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, 
                             geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  CHECK_ERR(rval);
  
  Range curves;
  int dim = 2;
  void *val[] = {&dim};
  int number_of_surfs;
  rval = mb->get_number_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
	  					    val, 1, number_of_surfs );
  CHECK_ERR( rval );
  

  CHECK_EQUAL( 6, number_of_surfs );  

}

void read_cube_vols_test()
{
  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );
  
  //Get geometry tag for pulling curve data from the mesh 
  Tag geom_tag;
  rval = mb->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, 
                             geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  CHECK_ERR( rval );
  
  Range curves;
  int dim = 3;
  void *val[] = {&dim};
  int number_of_vols;
  rval = mb->get_number_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
	  					    val, 1, number_of_vols );
  CHECK_ERR( rval );
  

  CHECK_EQUAL( 1, number_of_vols );

}

void read_cube_vertex_pos_test()
{
  
  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );

  //First check that the correct number of vertices are present
  int number_of_verts;
  rval = mb->get_number_entities_by_type( 0, MBVERTEX, number_of_verts );
  CHECK_ERR( rval );

  CHECK_EQUAL( 8, number_of_verts );

  //Retrieve all vertex handles from the mesh
  Range verts;
  rval = mb->get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR( rval );

  //Get the vertex coordinates
  double x[verts.size()];
  double y[verts.size()];
  double z[verts.size()];
  rval = mb-> get_coords( verts, &x[0], &y[0], &z[0] );
  CHECK_ERR( rval );

  //Check against known locations of the vertices

  std::vector<double> x_ref;
  std::vector<double> y_ref;
  std::vector<double> z_ref;

  // Vertex 1

  x_ref.push_back( 5 );
  y_ref.push_back( -5 );
  z_ref.push_back( 5 );

  // Vertex 2
  x_ref.push_back( 5 );
  y_ref.push_back( 5 );
  z_ref.push_back( 5 );

  // Vertex 3
  x_ref.push_back( -5 );
  y_ref.push_back( 5 );
  z_ref.push_back( 5 );

  // Vertex 4
  x_ref.push_back( -5 );
  y_ref.push_back( -5 );
  z_ref.push_back( 5 );

  // Vertex 5
  x_ref.push_back( 5 );
  y_ref.push_back( 5 );
  z_ref.push_back( -5 );

  // Vertex 6
  x_ref.push_back( 5 );
  y_ref.push_back( -5 );
  z_ref.push_back( -5 );

  // Vertex 7
  x_ref.push_back( -5 );
  y_ref.push_back( -5 );
  z_ref.push_back( -5 );
 
  // Vertex 8
  x_ref.push_back( -5 );
  y_ref.push_back( 5 );
  z_ref.push_back( -5 );
 

  std::cout << verts.size() << std::endl;
  std::cout << x_ref.size() << std::endl;
  
  for (unsigned int i=0; i<verts.size(); i++)
    {
      for (unsigned int j=0; j<x_ref.size(); j++)
	{
	  if( x[i]==x_ref[j] && y[i]==y_ref[j] && z[i]==z_ref[j])
            {
              x_ref.erase( x_ref.begin()+j );
              y_ref.erase( y_ref.begin()+j );
              z_ref.erase( z_ref.begin()+j );
              
            }
	}
    }
  
  int leftovers = x_ref.size();
  CHECK_EQUAL(0, leftovers );

}

void delete_mesh_test()
{
 Core moab;
 Interface* mb = &moab;
 read_file( mb, input_cube );

 ErrorCode rval; 

 //Get geometry tag for pulling curve data from the mesh
 Tag geom_tag; 
 rval = mb->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, 
                            geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
 CHECK_ERR(rval);

 Range geom_sets[4];

 for(unsigned dim=0; dim<4; dim++) 
 {
	void *val[] = {&dim};
	rval = mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
	  					    val, 1, geom_sets[dim] );
        CHECK_ERR( rval ); 

        if( geom_sets[dim].size() == 0 ) std::cout << "Warning: No geom sets to begin with" << std::endl;

 }

 mb->delete_mesh();

 Range geom_sets_after[4];
 for(unsigned dim=0; dim<4; dim++) 
 {
	void *val_after[] = {&dim};
	rval = mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
	  					    val_after, 1, geom_sets_after[dim] );
        CHECK_ERR( rval ); 

        if( 0 != geom_sets_after[dim].size() ) rval = MB_FAILURE;

        CHECK_ERR( rval );
 }

}
 
