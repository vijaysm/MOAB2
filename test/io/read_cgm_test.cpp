
#include <iostream>
#include "moab/Interface.hpp"
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"

using namespace moab;

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)

#ifdef MESHDIR
static const char input_file[] = STRINGIFY(MESHDIR) "/io/dum.sat";
static const char input_file2[] = STRINGIFY(MESHDIR) "/io/dum.stp";
#else
static const char input_file[] = "dum.sat";
static const char input_file2[] = "dum.stp";
#endif

void read_multiple_test() 
{
  Core mb;

  ErrorCode rval = mb.load_file(input_file);
  if (rval!=MB_SUCCESS)
  {
    std::cout<<"try loading now an stp file, supported by occ\n";
    // try loading an stp file, maybe
    rval = mb.load_file(input_file2);
    CHECK_ERR(rval);
    // try to load it second time
    rval = mb.load_file(input_file2);
    CHECK_ERR(rval);
  }
  else
  {
    rval = mb.load_file(input_file);
    CHECK_ERR(rval);
  }
}
  
int main(int argc, char* argv[])
{
  int result = RUN_TEST( read_multiple_test );

  return result;
}
