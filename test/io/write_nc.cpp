#include "TestUtil.hpp"
#include "moab/Core.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example_eul[] = STRINGIFY(MESHDIR) "/io/camEul26x48x96.t3.nc";
#else
static const char example_eul[] = "/io/camEul26x48x96.t3.nc";
#endif

#ifdef USE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#endif

// CAM-EUL
void test_read_write_one_var();

ErrorCode get_options(std::string& opts);

int main(int argc, char* argv[])
{
  int result = 0;

#ifdef USE_MPI
  int fail = MPI_Init(&argc, &argv);
  if (fail)
    return 1;
#else
  argv[0] = argv[argc - argc]; // To remove the warnings in serial mode about unused variables
#endif

  result += RUN_TEST(test_read_write_one_var);

#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail)
    return 1;
#endif

  return result;
}

void test_read_write_one_var()
{
  Core moab;
  Interface& mb = moab;

  std::string orig;
  ErrorCode rval = get_options(orig);
  CHECK_ERR(rval);

  // we will have to read without mesh and without var first, to create some tags that we need
  // later when we write
  std::string opts;
  opts = orig + std::string(";DEBUG_IO=3;NOMESH;VARIABLE=");

  // Need a set for nomesh to work right
  // this set will have as tags a lot of header information, that will be used for writing back the file
  EntityHandle set;
  rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

  // load one variable, and the mesh
  opts= orig + std::string(";DEBUG_IO=3;VARIABLE=T");
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

  // first test will write information about one variable
  std::string writeopts;
  writeopts = std::string(";;VARIABLE=T;DEBUG_IO=2;");
  rval = mb.write_file( "testT.nc", 0, writeopts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

ErrorCode get_options(std::string& opts)
{
#ifdef USE_MPI
  // Use parallel options
  opts = std::string(";;PARALLEL=READ_PART;PARTITION_METHOD=SQIJ");
  return MB_SUCCESS;
#else
  opts = std::string(";;");
  return MB_SUCCESS;
#endif
}
