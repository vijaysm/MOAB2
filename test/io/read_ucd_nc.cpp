#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include "TagInfo.hpp"
#include "MBTagConventions.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/homme3x3458.t.3.nc";
static const char conn_fname[] = STRINGIFY(MESHDIR) "/io/HommeMapping.nc";
#else
static const char example[] = "/io/homme3x3458.t.3.nc";
static const char conn_fname[] = "io/HommeMapping.nc";
#endif

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#endif

void test_read_all();
void test_read_onevar();
void test_read_onetimestep();
void test_read_nomesh();
void test_read_novars();
void test_read_coord_vars(); // Test reading coordinate variables
void test_gather_onevar(); // Test gather set with one variable
void test_read_conn(); // Test reading connectivity file only

void get_options(std::string& opts);

const int levels = 3;

int main(int argc, char* argv[])
{
  int result = 0;

#ifdef MOAB_HAVE_MPI
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
  result += RUN_TEST(test_read_coord_vars);
  result += RUN_TEST(test_gather_onevar);
  result += RUN_TEST(test_read_conn);

#ifdef MOAB_HAVE_MPI
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
  get_options(opts);

  ErrorCode rval = mb.load_file(example, NULL, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tags
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_onevar() 
{
  Core moab;
  Interface& mb = moab;

  std::string opts;
  get_options(opts);

  // Read mesh and read vertex variable T at all timesteps
  opts += std::string(";VARIABLE=T");
  ErrorCode rval = mb.load_file(example, NULL, opts.c_str());
  CHECK_ERR(rval);

#ifdef MOAB_HAVE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  int procs = pcomm->proc_config().proc_size();
#else
  int procs = 1;
#endif

  // Make check runs this test on one processor
  if (1 == procs) {
    // Check for proper tags
    Tag Ttag0, Ttag1;
    rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
    CHECK_ERR(rval);

    // Get vertices
    Range verts;
    rval = mb.get_entities_by_type(0, MBVERTEX, verts);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)3458, verts.size());

    // Get all values of tag T0
    int count;
    void* Tbuf;
    rval = mb.tag_iterate(Ttag0, verts.begin(), verts.end(), count, Tbuf);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)count, verts.size());

    const double eps = 0.0001;
    double* data = (double*) Tbuf;

    // Check first level values on 4 strategically selected vertices
    CHECK_REAL_EQUAL(233.1136, data[0 * levels], eps); // First vert
    CHECK_REAL_EQUAL(236.1505, data[1728 * levels], eps); // Median vert
    CHECK_REAL_EQUAL(235.7722, data[1729 * levels], eps); // Median vert
    CHECK_REAL_EQUAL(234.0416, data[3457 * levels], eps); // Last vert
  }
}

void test_read_onetimestep()
{
  Core moab;
  Interface& mb = moab;

  std::string opts;
  get_options(opts);

  opts += std::string(";TIMESTEP=1");
  ErrorCode rval = mb.load_file(example, NULL, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tags
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
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
  get_options(orig);

  opts = orig + std::string(";TIMESTEP=0");
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Now read 2nd timestep with nomesh option
  opts = orig + std::string(";TIMESTEP=1;NOMESH");
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
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
  get_options(orig);

  opts = orig + std::string(";NOMESH;VARIABLE=");
  rval = mb.load_file(example, &set, opts.c_str());
  CHECK_ERR(rval);

  opts = orig + std::string(";VARIABLE=;TIMESTEP=0");
  rval = mb.load_file(example, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  opts = orig + std::string(";VARIABLE=T;TIMESTEP=0;NOMESH");
  rval = mb.load_file(example, &set, opts.c_str());
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Now read 2nd timestep with nomesh option
  opts = orig + std::string(";TIMESTEP=1;NOMESH");
  rval = mb.load_file(example, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_coord_vars()
{
  Core moab;
  Interface& mb = moab;

  EntityHandle file_set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, file_set);
  CHECK_ERR(rval);

  std::string orig, opts;
  get_options(orig);

  opts = orig + std::string(";NOMESH;VARIABLE=");
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  std::string tag_name;
  int var_len;
  Tag var_tag;
  const void* var_data;

  // Check tag for regular coordinate variable lev
  tag_name = "lev";
  var_len = 0;
  rval = mb.tag_get_handle(tag_name.c_str(), var_len, MB_TYPE_OPAQUE, var_tag, MB_TAG_SPARSE | MB_TAG_VARLEN);
  CHECK_ERR(rval);
  CHECK_EQUAL(true, var_tag->variable_length());
  CHECK_EQUAL(MB_TYPE_DOUBLE, var_tag->get_data_type());

  // Check lev tag size and values on file_set
  rval = mb.tag_get_by_ptr(var_tag, &file_set, 1, &var_data, &var_len);
  CHECK_ERR(rval);
  CHECK_EQUAL(levels, var_len);
  double* lev_val = (double*)var_data;
  const double eps = 1e-10;
  CHECK_REAL_EQUAL(3.54463800000002, lev_val[0], eps);
  CHECK_REAL_EQUAL(13.9672100000001, lev_val[levels - 1], eps);

  // Check tag for dummy coordinate variable ncol
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

  // Check tag for regular coordinate lev
  tag_name = "lev";
  var_len = 0;
  rval = mb.tag_get_handle(tag_name.c_str(), var_len, MB_TYPE_OPAQUE, var_tag, MB_TAG_SPARSE | MB_TAG_VARLEN);
  CHECK_ERR(rval);
  CHECK_EQUAL(true, var_tag->variable_length());
  CHECK_EQUAL(MB_TYPE_DOUBLE, var_tag->get_data_type());

  // Check lev tag size and values on file_set2
  rval = mb.tag_get_by_ptr(var_tag, &file_set2, 1, &var_data, &var_len);
  CHECK_ERR(rval);
  CHECK_EQUAL(levels, var_len);
  lev_val = (double*)var_data;
  CHECK_REAL_EQUAL(3.54463800000002, lev_val[0], eps);
  CHECK_REAL_EQUAL(13.9672100000001, lev_val[levels - 1], eps);

  // Check tag for dummy coordinate variable ncol
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

void test_gather_onevar()
{
  Core moab;
  Interface& mb = moab;

  EntityHandle file_set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, file_set);
  CHECK_ERR(rval);

  std::string opts;
  get_options(opts);

  // Read vertex variable T and create gather set on processor 0
  opts += ";VARIABLE=T;GATHER_SET=0";
#ifdef MOAB_HAVE_MPI
  opts += ";PARALLEL_RESOLVE_SHARED_ENTS";
#endif
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

#ifdef MOAB_HAVE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  int rank = pcomm->proc_config().proc_rank();

  Range verts, verts_owned;
  rval = mb.get_entities_by_type(file_set, MBVERTEX, verts);
  CHECK_ERR(rval);

  // Get local owned vertices
  rval = pcomm->filter_pstatus(verts, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &verts_owned);
  CHECK_ERR(rval);

  EntityHandle gather_set = 0;
  if (0 == rank) {
    // Get gather set
    ReadUtilIface* readUtilIface;
    mb.query_interface(readUtilIface);
    rval = readUtilIface->get_gather_set(gather_set);
    CHECK_ERR(rval);
    assert(gather_set != 0);
  }

  Tag Ttag0, gid_tag;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0, MB_TAG_DENSE);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid_tag, MB_TAG_DENSE);
  CHECK_ERR(rval);

  pcomm->gather_data(verts_owned, Ttag0, gid_tag, gather_set, 0);

  if (0 == rank) {
    // Get gather set vertices
    Range gather_set_verts;
    rval = mb.get_entities_by_type(gather_set, MBVERTEX, gather_set_verts);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)3458, gather_set_verts.size());

    // Get T0 tag values on 4 strategically selected gather set vertices
    double T0_val[4 * levels];
    EntityHandle vert_ents[] = {gather_set_verts[0], gather_set_verts[1728],
                                gather_set_verts[1729], gather_set_verts[3457]};
    rval = mb.tag_get_data(Ttag0, vert_ents, 4, T0_val);
    CHECK_ERR(rval);

    const double eps = 0.001;

    // Check first level values
    CHECK_REAL_EQUAL(233.1136, T0_val[0 * levels], eps); // First vert
    CHECK_REAL_EQUAL(236.1505, T0_val[1 * levels], eps); // Median vert
    CHECK_REAL_EQUAL(235.7722, T0_val[2 * levels], eps); // Median vert
    CHECK_REAL_EQUAL(234.0416, T0_val[3 * levels], eps); // Last vert
  }
#endif
}

void test_read_conn()
{
  Core moab;
  Interface& mb = moab;

  std::string opts;
  get_options(opts);

  ErrorCode rval = mb.load_file(conn_fname, NULL, opts.c_str());
  CHECK_ERR(rval);

#ifdef MOAB_HAVE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  int procs = pcomm->proc_config().proc_size();
#else
  int procs = 1;
#endif

  // Make check runs this test on one processor
  if (1 == procs) {
    // Get vertices
    Range verts;
    rval = mb.get_entities_by_type(0, MBVERTEX, verts);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)3458, verts.size());

    // Get cells
    Range cells;
    rval = mb.get_entities_by_type(0, MBQUAD, cells);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)3456, cells.size());
  }
}

void get_options(std::string& opts)
{
#ifdef MOAB_HAVE_MPI
  // Use parallel options
  opts = std::string(";;PARALLEL=READ_PART;PARTITION_METHOD=TRIVIAL");
#else
  opts = std::string(";;");
#endif
}
