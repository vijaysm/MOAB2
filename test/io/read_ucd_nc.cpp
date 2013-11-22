#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include "TagInfo.hpp"

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
void test_read_dim_vars();

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
  result += RUN_TEST(test_read_dim_vars);

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
  EntityHandle file_set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, file_set);
  CHECK_ERR(rval);

  std::string orig, opts;
  rval = get_options(orig);
  CHECK_ERR(rval);

  opts = orig + std::string(";TIMESTEP=0");
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", 26, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Now read 2nd timestep with nomesh option
  opts = orig + std::string(";TIMESTEP=1;NOMESH");
  rval = mb.load_file(example, &file_set, opts.c_str());
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

void test_read_dim_vars()
{
  Core moab;
  Interface& mb = moab;

  EntityHandle file_set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, file_set);
  CHECK_ERR(rval);

  std::string orig, opts;
  rval = get_options(orig);
  CHECK_ERR(rval);

  opts = orig + std::string(";NOMESH;VARIABLE=");
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  std::string tag_name;
  int var_len;
  Tag var_tag;
  const void* var_data;

  // Check tag for regular dimension variable lev
  tag_name = "lev";
  var_len = 0;
  rval = mb.tag_get_handle(tag_name.c_str(), var_len, MB_TYPE_OPAQUE, var_tag, MB_TAG_SPARSE | MB_TAG_VARLEN);
  CHECK_ERR(rval);
  CHECK_EQUAL(true, var_tag->variable_length());
  CHECK_EQUAL(MB_TYPE_DOUBLE, var_tag->get_data_type());

  // Check lev tag size and values on file_set
  rval = mb.tag_get_by_ptr(var_tag, &file_set, 1, &var_data, &var_len);
  CHECK_ERR(rval);
  CHECK_EQUAL(26, var_len);
  double* lev_val = (double*)var_data;
  const double eps = 1e-10;
  CHECK_REAL_EQUAL(3.54463800000002, lev_val[0], eps);
  CHECK_REAL_EQUAL(992.556100000005, lev_val[25], eps);

  // Check tag for dummy dimension variable ncol
  tag_name = "ncol";
  var_len = 0;
  rval = mb.tag_get_handle(tag_name.c_str(), var_len, MB_TYPE_OPAQUE, var_tag, MB_TAG_SPARSE | MB_TAG_VARLEN);
  CHECK_ERR(rval);
  CHECK_EQUAL(true, var_tag->variable_length());
  CHECK_EQUAL(MB_TYPE_INTEGER, var_tag->get_data_type());

  // Check ncol tag size and values on file_set
  rval = mb.tag_get_by_ptr(var_tag, &file_set, 1, &var_data, &var_len);
  CHECK_ERR(rval);
  CHECK_EQUAL(1, var_len);
  int* ncol_val = (int*)var_data;
  CHECK_EQUAL(3458, ncol_val[0]);

  // Create another file set
  EntityHandle file_set2;
  rval = mb.create_meshset(MESHSET_SET, file_set2);
  CHECK_ERR(rval);

  // Read file again with file_set2
  rval = mb.load_file(example, &file_set2, opts.c_str());
  CHECK_ERR(rval);

  // Check tag for regular dimension variable lev
  tag_name = "lev";
  var_len = 0;
  rval = mb.tag_get_handle(tag_name.c_str(), var_len, MB_TYPE_OPAQUE, var_tag, MB_TAG_SPARSE | MB_TAG_VARLEN);
  CHECK_ERR(rval);
  CHECK_EQUAL(true, var_tag->variable_length());
  CHECK_EQUAL(MB_TYPE_DOUBLE, var_tag->get_data_type());

  // Check lev tag size and values on file_set2
  rval = mb.tag_get_by_ptr(var_tag, &file_set2, 1, &var_data, &var_len);
  CHECK_ERR(rval);
  CHECK_EQUAL(26, var_len);
  lev_val = (double*)var_data;
  CHECK_REAL_EQUAL(3.54463800000002, lev_val[0], eps);
  CHECK_REAL_EQUAL(992.556100000005, lev_val[25], eps);

  // Check tag for dummy dimension variable ncol
  tag_name = "ncol";
  var_len = 0;
  rval = mb.tag_get_handle(tag_name.c_str(), var_len, MB_TYPE_OPAQUE, var_tag, MB_TAG_SPARSE | MB_TAG_VARLEN);
  CHECK_ERR(rval);
  CHECK_EQUAL(true, var_tag->variable_length());
  CHECK_EQUAL(MB_TYPE_INTEGER, var_tag->get_data_type());

  // Check ncol tag size and values on file_set2
  rval = mb.tag_get_by_ptr(var_tag, &file_set2, 1, &var_data, &var_len);
  CHECK_ERR(rval);
  CHECK_EQUAL(1, var_len);
  ncol_val = (int*)var_data;
  CHECK_EQUAL(3458, ncol_val[0]);
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
