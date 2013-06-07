#include "TestUtil.hpp"
#include "moab/Core.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example_eul[] = STRINGIFY(MESHDIR) "/io/camEul26x48x96.t3.nc";
static const char example_fv[] = STRINGIFY(MESHDIR) "/io/fv26x46x72.t.3.nc";
#else
static const char example_eul[] = "/io/camEul26x48x96.t3.nc";
static const char example_fv[] = "/io/fv26x46x72.t.3.nc";
#endif

#ifdef USE_MPI
#include "moab_mpi.h"
#endif

// CAM-EUL
void test_read_eul_all();
void test_read_eul_onevar();
void test_read_eul_onetimestep();
void test_read_eul_nomesh();
void test_read_eul_novars();

// CAM-FV
void test_read_fv_all();
void test_read_fv_onevar();
void test_read_fv_onetimestep();
void test_read_fv_nomesh();
void test_read_fv_novars();

ErrorCode get_options(std::string& opts);

int main(int argc, char* argv[])
{
  int result = 0;

#ifdef USE_MPI
  int fail = MPI_Init(&argc, &argv);
  if (fail)
    return 1;
#else
  argv[0] = argv[argc - argc]; // to remove the warnings in serial mode about unused variables
#endif

  result += RUN_TEST(test_read_eul_all);
  result += RUN_TEST(test_read_eul_onevar);
  result += RUN_TEST(test_read_eul_onetimestep);
  result += RUN_TEST(test_read_eul_nomesh);
  result += RUN_TEST(test_read_eul_novars);

  // Exclude test_read_fv_all() since reading edge data is not implemented in MOAB yet
  //result += RUN_TEST(test_read_fv_all);
  result += RUN_TEST(test_read_fv_onevar);
  result += RUN_TEST(test_read_fv_onetimestep);
  result += RUN_TEST(test_read_fv_nomesh);
  result += RUN_TEST(test_read_fv_novars);

#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail) return 1;
#endif

  return result;
}

void test_read_eul_all()
{
  Core moab;
  Interface& mb = moab;

  std::string opts;
  ErrorCode rval = get_options(opts);
  CHECK_ERR(rval);

  rval = mb.load_file(example_eul, NULL, opts.c_str());
  CHECK_ERR(rval);

    // check for proper tags
  Tag Ttag0, Ttag1, coordTag;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("COORDS", 3, MB_TYPE_DOUBLE, coordTag);
  CHECK_ERR(rval);
}

void test_read_eul_onevar() 
{
  Core moab;
  Interface& mb = moab;
  std::string opts;
  ErrorCode rval = get_options(opts);
  CHECK_ERR(rval);

  opts += std::string(";VARIABLE=T");
  rval = mb.load_file(example_eul, NULL, opts.c_str());
  CHECK_ERR(rval);

    // check for proper tags
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_eul_onetimestep()
{
  Core moab;
  Interface& mb = moab;
  std::string opts;
  ErrorCode rval = get_options(opts);
  CHECK_ERR(rval);

  opts += std::string(";VARIABLE=T;TIMESTEP=1");
  rval = mb.load_file(example_eul, NULL, opts.c_str());
  CHECK_ERR(rval);

    // check for proper tags
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_eul_nomesh() 
{
  Core moab;
  Interface& mb = moab;

    // need a set for nomesh to work right
  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  std::string orig, opts;
  rval = get_options(orig);
  CHECK_ERR(rval);

  opts = orig + std::string(";VARIABLE=T;TIMESTEP=0");
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

    // check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

    // now read 2nd timestep with nomesh option
  opts = orig + std::string(";VARIABLE=T;TIMESTEP=1;NOMESH");
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

    // check for proper tag
  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_eul_novars() 
{
  Core moab;
  Interface& mb = moab;

    // need a set for nomesh to work right
  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  std::string orig, opts;
  rval = get_options(orig);
  CHECK_ERR(rval);

  opts = orig + std::string(";NOMESH;VARIABLE=");
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

  opts = orig + std::string(";VARIABLE=;TIMESTEP=0");
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

    // check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  opts = orig + std::string(";VARIABLE=T;TIMESTEP=0;NOMESH");
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

    // now read 2nd timestep with nomesh option
  opts = orig + std::string(";VARIABLE=T;TIMESTEP=1;NOMESH");
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

    // check for proper tag
  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_fv_all()
{
  Core moab;
  Interface& mb = moab;

  std::string opts;
  ErrorCode rval = get_options(opts);
  CHECK_ERR(rval);

  rval = mb.load_file(example_fv, NULL, opts.c_str());
  CHECK_ERR(rval);

    // check for proper tags
  Tag Ttag0, Ttag1, coordTag;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);

  rval=mb.tag_get_handle("COORDS", 3, MB_TYPE_DOUBLE, coordTag);
  CHECK_ERR(rval);
}

void test_read_fv_onevar() 
{
  Core moab;
  Interface& mb = moab;
  std::string opts;
  ErrorCode rval = get_options(opts);
  CHECK_ERR(rval);

  opts += std::string(";VARIABLE=T");
  rval = mb.load_file(example_fv, NULL, opts.c_str());
  CHECK_ERR(rval);

    // check for proper tags
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_fv_onetimestep()
{
  Core moab;
  Interface& mb = moab;
  std::string opts;
  ErrorCode rval = get_options(opts);
  CHECK_ERR(rval);

  opts += std::string(";VARIABLE=T;TIMESTEP=1");
  rval = mb.load_file(example_fv, NULL, opts.c_str());
  CHECK_ERR(rval);

    // check for proper tags
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_fv_nomesh() 
{
  Core moab;
  Interface& mb = moab;

    // need a set for nomesh to work right
  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  std::string orig, opts;
  rval = get_options(orig);
  CHECK_ERR(rval);

  opts = orig + std::string(";VARIABLE=T;TIMESTEP=0");
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);

    // check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

    // now read 2nd timestep with nomesh option
  opts = orig + std::string(";VARIABLE=T;TIMESTEP=1;NOMESH");
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);

    // check for proper tag
  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_fv_novars() 
{
  Core moab;
  Interface& mb = moab;

    // need a set for nomesh to work right
  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  std::string orig, opts;
  rval = get_options(orig);
  CHECK_ERR(rval);

  opts = orig + std::string(";NOMESH;VARIABLE=");
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);

  opts = orig + std::string(";VARIABLE=;TIMESTEP=0");
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);

    // check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  opts = orig + std::string(";VARIABLE=T;TIMESTEP=0;NOMESH");
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

    // now read 2nd timestep with nomesh option
  opts = orig + std::string(";VARIABLE=T;TIMESTEP=1;NOMESH");
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);

    // check for proper tag
  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

ErrorCode get_options(std::string& opts) 
{
#ifdef USE_MPI
    // use parallel options
  opts = std::string(";;PARALLEL=READ_PART;PARTITION_METHOD=SQIJ");
  return MB_SUCCESS;
#else
  opts = std::string(";;");
  return MB_SUCCESS;
#endif
}
