#include <iostream>
#include "moab/Interface.hpp"
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"

#include "DagMC.hpp"

using namespace moab;

using moab::DagMC;

#define DAG DagMC::instance()

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)

#ifdef MESHDIR
static const char input_file[] = STRINGIFY(MESHDIR) "/dagmc/test_geom.h5m";
#else
static const char input_file[] = STRINGIFY(MESHDIR) "/dagmc/test_geom.h5m";
#endif


void dagmc_load_file() 
{
  ErrorCode rval = DAG->load_file(input_file); // open the Dag file
  CHECK_ERR(rval);
}

void dagmc_build_obb() 
{
  ErrorCode rval = DAG->init_OBBTree();
  CHECK_ERR(rval);
}
  
  
int main(int /* argc */, char** /* argv */)
{
  int result = 0;
  result += RUN_TEST( dagmc_load_file ); // test ray fire
  result += RUN_TEST( dagmc_build_obb ); // build the obb

  return result;
}
