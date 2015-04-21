#include "TestUtil.hpp"
#include "moab/Core.hpp"
#define IS_BUILDING_MB
#include "moab/Range.hpp"
#include "moab/FileOptions.hpp"

using namespace moab;

#ifdef MESHDIR
static const char abaqusFile[] = STRINGIFY(MESHDIR) "/io/tet4.abq";
#else
static const char abaqusFile[] =  "tet4.abq";
#endif

void test_read();

int main()
{
  int result = 0;

  result += RUN_TEST(test_read);

  return result;
}

void test_read()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  rval = mb.load_file(abaqusFile);
  if (MB_SUCCESS != rval) std::cerr << "Trouble reading file " << abaqusFile << std::endl;
  CHECK_ERR(rval);
  // should elaborate more
}

