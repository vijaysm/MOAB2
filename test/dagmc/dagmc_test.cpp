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

void dagmc_num_vols()
{
  int num_vols = DAG->num_entities(3); 
  CHECK_EQUAL(2,num_vols);
}

void dagmc_entity_handle()
{
  int num_vols = DAG->num_entities(3); 
  EntityHandle vol;
  for ( int i = 0 ; i < num_vols ; i++ )
    {
      vol = DAG->entity_by_index(3,i);
    }
  CHECK_EQUAL(12682136550675316765,vol);
}

void dagmc_point_in()
{
  int result = 0;
  double xyz[3]={0.0,0.0,0.0};
  ErrorCode = DAG->point_in_volume(12682136550675316765,pos,&result);
  CHECK_EQUAL(1,result);
}
  
int main(int /* argc */, char** /* argv */)
{
  int result = 0;
  result += RUN_TEST( dagmc_load_file ); // test ray fire
  result += RUN_TEST( dagmc_build_obb ); // build the obb
  result += RUN_TEST( dagmc_num_vols  ); // make sure the num of vols correct
  result += RUN_TEST( dagmc_entity_handle); // check the entity handle correct
  result += RUN_TEST( dagmc_point_in);
  return result;
}
