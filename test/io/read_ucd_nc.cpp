#include "TestUtil.hpp"
#include "moab/Core.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/homme26x3458.t.3.nc";
#else
static const char example[] = "/io/homme26x3458.t.3.nc";
#endif

#ifdef USE_MPI
#include "moab_mpi.h"
#endif

void read_file( Interface& moab, const char* input_file );
void test_read_all();
void test_read_onevar();
void test_read_onetimestep();
void test_read_nomesh();
void test_read_novars();
ErrorCode get_options(std::string &opts);

int main(int argc, char *argv[])
{
  int result = 0;

#ifdef USE_MPI
  int fail = MPI_Init(&argc, &argv);
  if (fail) return 1;
#endif
  
  result += RUN_TEST(test_read_all);
  result += RUN_TEST(test_read_onevar);
  result += RUN_TEST(test_read_onetimestep);
  result += RUN_TEST(test_read_nomesh);
  result += RUN_TEST(test_read_novars);
  
#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail) return 1;
#endif
  
  return result;
}


void test_read_all()
{
  Core moab;
  Interface& mb = moab;

  std::string opts;
  ErrorCode rval = get_options(opts);
  CHECK_ERR(rval);
  
  rval = mb.load_file( example, NULL, opts.c_str());
  CHECK_ERR(rval);
  
    // check for proper tags
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);
  
  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_onevar() 
{
  Core moab;
  Interface& mb = moab;
  std::string opts;
  ErrorCode rval = get_options(opts);
  CHECK_ERR(rval);

  opts += std::string(";VARIABLE=T");
  rval = mb.load_file( example, NULL, opts.c_str());
  CHECK_ERR(rval);
  
    // check for proper tags
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);
  
  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_onetimestep()
{
  Core moab;
  Interface& mb = moab;
  std::string opts;
  ErrorCode rval = get_options(opts);
  CHECK_ERR(rval);

  opts += std::string(";TIMESTEP=1");
  rval = mb.load_file( example, NULL, opts.c_str() );
  CHECK_ERR(rval);
  
    // check for proper tags
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);
  
  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_nomesh() 
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

  opts = orig + std::string(";TIMESTEP=0");
  rval = mb.load_file( example, &set, opts.c_str() );
  CHECK_ERR(rval);
  
    // check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);
  
  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

    // now read 2nd timestep with nomesh option
  opts = orig + std::string(";TIMESTEP=1;NOMESH");
  rval = mb.load_file( example, &set, opts.c_str() );
  CHECK_ERR(rval);
  
    // check for proper tag
  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_novars()
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
  rval = mb.load_file( example, &set, opts.c_str() );
  CHECK_ERR(rval);

  opts = orig + std::string(";VARIABLE=;TIMESTEP=0");
  rval = mb.load_file( example, &set, opts.c_str() );
  CHECK_ERR(rval);

    // check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  opts = orig + std::string(";VARIABLE=T;TIMESTEP=0;NOMESH");
  rval = mb.load_file( example, &set, opts.c_str() );
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

    // now read 2nd timestep with nomesh option
  opts = orig + std::string(";TIMESTEP=1;NOMESH");
  rval = mb.load_file( example, &set, opts.c_str() );
  CHECK_ERR(rval);

    // check for proper tag
  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

ErrorCode get_options(std::string &opts) 
{
#ifdef USE_MPI
    // use parallel options
  opts = std::string(";;TRIVIAL_PARTITION");
  return MB_SUCCESS;
#else
  opts = std::string(";;");
  return MB_SUCCESS;
#endif
}
