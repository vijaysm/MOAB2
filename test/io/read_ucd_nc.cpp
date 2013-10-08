#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/homme26x3458.t.3.nc";
#else
static const char example[] = "/io/homme26x3458.t.3.nc";
#endif

#ifdef USE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#endif

void test_read_all();
void test_read_onevar();
void test_read_onetimestep();
void test_read_nomesh();
void test_read_novars();

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

  result += RUN_TEST(test_read_all);
  result += RUN_TEST(test_read_onevar);
  result += RUN_TEST(test_read_onetimestep);
  result += RUN_TEST(test_read_nomesh);
  result += RUN_TEST(test_read_novars);

#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail)
    return 1;
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

  rval = mb.load_file(example, NULL, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tags
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
  // Create gather set
  opts += std::string(";GATHER_SET=");
  rval = mb.load_file(example, NULL, opts.c_str());
  CHECK_ERR(rval);

  // Check values of tag T0 at some strategically chosen places below
  int procs = 1;
#ifdef USE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  procs = pcomm->proc_config().proc_size();
#endif

  // Make check runs this test in one processor
  if (1 == procs) {
    // Check for proper tags
    Tag Ttag0, Ttag1;
    rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
    CHECK_ERR(rval);

    // Get vertices
    Range verts;
    rval = mb.get_entities_by_type(0, MBVERTEX, verts);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)6916, verts.size()); // Gather set vertices included

    // Get gather set
    EntityHandle gather_set;
    ReadUtilIface* readUtilIface;
    mb.query_interface(readUtilIface);
    rval = readUtilIface->get_gather_set(gather_set);
    CHECK_ERR(rval);

    // Get gather set entities
    Range gather_ents;
    rval = mb.get_entities_by_handle(gather_set, gather_ents);
    CHECK_ERR(rval);

    // Remove gather set vertices
    verts = subtract(verts, gather_ents);
    CHECK_EQUAL((size_t)3458, verts.size()); // Gather set vertices excluded

    // Get all values of tag T0
    int count;
    void* Tbuf;
    rval = mb.tag_iterate(Ttag0, verts.begin(), verts.end(), count, Tbuf);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)count, verts.size());

    const double eps = 0.0001;
    double* data = (double*) Tbuf;

    // Check first level values at some vertices
    CHECK_REAL_EQUAL(233.1136, data[0 * 26], eps); // First vert
    CHECK_REAL_EQUAL(236.1505, data[1728 * 26], eps); // Median vert
    CHECK_REAL_EQUAL(235.7722, data[1729 * 26], eps); // Median vert
    CHECK_REAL_EQUAL(234.0416, data[3457 * 26], eps); // Last vert
  }
}

void test_read_onetimestep()
{
  Core moab;
  Interface& mb = moab;
  std::string opts;
  ErrorCode rval = get_options(opts);
  CHECK_ERR(rval);

  opts += std::string(";TIMESTEP=1");
  rval = mb.load_file(example, NULL, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tags
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

  // Need a set for nomesh to work right
  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  std::string orig, opts;
  rval = get_options(orig);
  CHECK_ERR(rval);

  opts = orig + std::string(";TIMESTEP=0");
  rval = mb.load_file(example, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Now read 2nd timestep with nomesh option
  opts = orig + std::string(";TIMESTEP=1;NOMESH");
  rval = mb.load_file(example, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_novars()
{
  Core moab;
  Interface& mb = moab;

  // Need a set for nomesh to work right
  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  std::string orig, opts;
  rval = get_options(orig);
  CHECK_ERR(rval);

  opts = orig + std::string(";NOMESH;VARIABLE=");
  rval = mb.load_file(example, &set, opts.c_str());
  CHECK_ERR(rval);

  opts = orig + std::string(";VARIABLE=;TIMESTEP=0");
  rval = mb.load_file(example, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  opts = orig + std::string(";VARIABLE=T;TIMESTEP=0;NOMESH");
  rval = mb.load_file(example, &set, opts.c_str());
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Now read 2nd timestep with nomesh option
  opts = orig + std::string(";TIMESTEP=1;NOMESH");
  rval = mb.load_file(example, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

ErrorCode get_options(std::string& opts)
{
#ifdef USE_MPI
  // Use parallel options
  opts = std::string(";;PARALLEL=READ_PART;PARTITION_METHOD=TRIVIAL");
  return MB_SUCCESS;
#else
  opts = std::string(";;");
  return MB_SUCCESS;
#endif
}
