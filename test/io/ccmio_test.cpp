#include "TestUtil.hpp"
#include "moab/Core.hpp"
#define IS_BUILDING_MB
#include "moab/Range.hpp"
#include "FileOptions.hpp"

using namespace moab;

#ifdef MESHDIR
static const char cubfile[] = STRINGIFY(MESHDIR) "/io/singlecyl.cub";
static const char ccmgfiler[] = STRINGIFY(MESHDIR) "/io/singlecyl.ccmg";
static const char ccmgfilew[] = "singlecyl_tmp.ccmg";
#else
static const char cubfile[] =  "singlecyl.cub";
static const char ccmgfile[] = "singlecyl.ccmg";
static const char ccmgfilew[] = "singlecyl_tmp.ccmg";
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
  if (MB_SUCCESS != rval) std::cerr << "Trouble reading file " << cubfile << std::endl;
  CHECK_ERR(rval);
  
  rval = mb.write_file(ccmgfilew);
  if (MB_SUCCESS != rval) std::cerr << "Trouble writing file " << ccmgfilew << std::endl;
  CHECK_ERR(rval);
}  

void test_read()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  rval = mb.load_file(ccmgfiler);
  CHECK_ERR(rval);
}  

