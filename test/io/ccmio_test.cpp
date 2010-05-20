#include "TestUtil.hpp"
#include "moab/Core.hpp"
#define IS_BUILDING_MB
#include "moab/Range.hpp"
#include "FileOptions.hpp"

using namespace moab;

#ifdef MESHDIR
static const char cubfile[] = STRINGIFY(MESHDIR) "/../singlecyl.cub";
static const char ccmgfile[] = STRINGIFY(MESHDIR) "/io/singlecyl.ccmg";
#else
static const char cubfile[] = STRINGIFY(MESHDIR) "singlecyl.cub";
static const char ccmgfile[] = STRINGIFY(MESHDIR) "singlecyl.ccmg";
#endif

void read_file( Interface& moab, const char* input_file );
void test_read();
void test_write();

int main()
{
  int result = 0;
  
  result += RUN_TEST(test_write);
  result += RUN_TEST(test_read);
  
  return result;
}

void test_write()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  rval = mb.load_file(cubfile);
  CHECK_ERR(rval);
  
  rval = mb.write_file(ccmgfile);
  CHECK_ERR(rval);
}  

void test_read()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  rval = mb.load_file(ccmgfile);
  CHECK_ERR(rval);
}  

