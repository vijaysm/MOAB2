
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

void read_cube_test()
{
  ErrorCode rval;
  Core moab;
  Interface* mb = &moab;
  mb->delete_mesh();
  read_file( mb, input_cube );
   
  int number_of_tris;

  rval = mb->get_number_entities_by_type(0, MBTRI , number_of_tris);
  std::cout << "Number of Triangles = " << number_of_tris << std::endl;
  CHECK_ERR(rval);

  int number_of_vertices;
  rval = mb->get_number_entities_by_type(0, MBVERTEX, number_of_vertices);
  CHECK_ERR(rval);


  if( number_of_tris != 12) rval = MB_FAILURE; CHECK_ERR(rval);
   
  if( number_of_vertices !=8) rval = MB_FAILURE; CHECK_ERR(rval);


}
  
int main(int /* argc */, char** /* argv */)
{
  int result = 0;

  result += RUN_TEST( read_cube_test );
  result += RUN_TEST( read_cube_test );

  return result;
}
