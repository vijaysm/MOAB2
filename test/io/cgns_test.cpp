#include "TestUtil.hpp"
#include "moab/Core.hpp"
#define IS_BUILDING_MB
#include "moab/Range.hpp"
#include "moab/FileOptions.hpp"

using namespace moab;

#ifdef MESHDIR
static const char cgnsfile[] = STRINGIFY(MESHDIR) "/io/2d_naca0012.cgns";
static const char cgnsfilew[] = STRINGIFY(MESHDIR) "/io/test.cgns";
#else
static const char cgnsfile[] =  "2d_naca0012.cgns";
static const char cgnsfilew[] =  "test.cgns";
#endif

void test_read_write();

int main()
{
  int result = 0;

  result += RUN_TEST(test_read_write);

  return result;
}

void test_read_write()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  rval = mb.load_file(cgnsfile);
  if (MB_SUCCESS != rval) std::cerr << "Trouble reading file " << cgnsfile << std::endl;
  CHECK_ERR(rval);

  rval = mb.write_file(cgnsfilew);
  if (MB_SUCCESS != rval) std::cerr << "Trouble writing file " << cgnsfilew << std::endl;
  CHECK_ERR(rval);
}

