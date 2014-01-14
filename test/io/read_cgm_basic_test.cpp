
#include <iostream>
#include "moab/Interface.hpp"
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"

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

void read_file( Interface* moab, const char* input_file );
void read_cube_test();


void read_file( Interface* moab, const char* input_file )
{
  ErrorCode rval = moab->load_file( input_file );
  CHECK_ERR(rval);
}

void read_cube_verts_test()
{
  ErrorCode rval;
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
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );
   
  int number_of_tris;

  rval = mb->get_number_entities_by_type(0, MBTRI , number_of_tris);
  std::cout << "Number of Triangles = " << number_of_tris << std::endl;
  CHECK_ERR(rval);

  CHECK_EQUAL( 12, number_of_tris);  

}
 
int main(int /* argc */, char** /* argv */)
{
  int result = 0;

  result += RUN_TEST( read_cube_tris_test );  

  return result;
}


void delete_mesh_test()
{
 Core moab;
 Interface* mb = &moab;
 read_file( mb, input_cube );

 ErrorCode rval; 

 Tag geom_tag; 

 rval = mb->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1,
				MB_TYPE_INTEGER, geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
 CHECK_ERR(rval);

 Range geom_sets[4];

 for(unsigned dim=0; dim<4; dim++) 
 {
	void *val[] = {&dim};
	rval = mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
	  					    val, 1, geom_sets[dim] );
        CHECK_ERR(rval); 

        if( geom_sets[dim].size() == 0 ) std::cout << "Warning: No geom sets to begin with" << std::endl;

 }

 mb->delete_mesh();

 Range geom_sets_after[4];
 for(unsigned dim=0; dim<4; dim++) 
 {
	void *val_after[] = {&dim};
	rval = mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
	  					    val_after, 1, geom_sets_after[dim] );
        CHECK_ERR(rval); 

        if( 0 != geom_sets_after[dim].size() ) rval = MB_FAILURE;

        CHECK_ERR(rval);
 }

}

