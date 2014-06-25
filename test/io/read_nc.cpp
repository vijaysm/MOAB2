#include "TestUtil.hpp"
#include "moab/Core.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example_eul[] = STRINGIFY(MESHDIR) "/io/eul3x48x96.t.3.nc";
static const char example_fv[] = STRINGIFY(MESHDIR) "/io/fv3x46x72.t.3.nc";
#else
static const char example_eul[] = "/io/eul3x48x96.t.3.nc";
static const char example_fv[] = "/io/fv3x46x72.t.3.nc";
#endif

#ifdef USE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
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
#ifdef USE_MPI
void test_read_fv_ghosting();
#endif

ErrorCode get_options(std::string& opts);

const int levels = 3;

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
  // Before ghosting issues with ownership were fixed, this test failed on 4 processors
  result += RUN_TEST(test_read_fv_ghosting);
#endif

#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail)
    return 1;
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

  // Check for proper tags
  Tag Ttag0, Ttag1, coordTag;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("COORDS", 3, MB_TYPE_DOUBLE, coordTag);
  CHECK_ERR(rval);

  // Check for some tags with double underscore in the tag name
  Tag tempTag;
  rval = mb.tag_get_handle("__lon_LOC_MINMAX", 2, MB_TYPE_INTEGER, tempTag);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("__lon_LOC_VALS", 0, MB_TYPE_DOUBLE, tempTag, MB_TAG_VARLEN);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("__lon_GLOBAL_MINMAX", 2, MB_TYPE_INTEGER, tempTag);
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

  // Check for proper tags
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);

  // Check values of tag T0 (first level) at some strategically chosen places below
  int rank = 0;
  int procs = 1;
#ifdef USE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  rank = pcomm->proc_config().proc_rank();
  procs = pcomm->proc_config().proc_size();
#endif

  const double eps = 0.0001;
  double val[8 * levels];

  if (1 == procs) {
    Range global_quads;
    rval = mb.get_entities_by_type(0, MBQUAD, global_quads);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)4608, global_quads.size());

    EntityHandle gloabl_quad_ents[] = {global_quads[0], global_quads[2255], global_quads[2304], global_quads[4559],
                                       global_quads[48], global_quads[2303], global_quads[2352], global_quads[4607]};
    rval = mb.tag_get_data(Ttag0, &gloabl_quad_ents[0], 8, val);

    CHECK_REAL_EQUAL(252.8529, val[0 * levels], eps); // First global quad
    CHECK_REAL_EQUAL(234.8390, val[1 * levels], eps); // 2256th global quad
    CHECK_REAL_EQUAL(232.6458, val[2 * levels], eps); // 2305th global quad
    CHECK_REAL_EQUAL(205.3905, val[3 * levels], eps); // 4560th global quad
    CHECK_REAL_EQUAL(252.7116, val[4 * levels], eps); // 49th global quad
    CHECK_REAL_EQUAL(232.6670, val[5 * levels], eps); // 2304th global quad
    CHECK_REAL_EQUAL(234.6922, val[6 * levels], eps); // 2353th global quad
    CHECK_REAL_EQUAL(200.6828, val[7 * levels], eps); // Last global quad
  }
  else if (2 == procs) {
    Range local_quads;
    rval = mb.get_entities_by_type(0, MBQUAD, local_quads);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)2304, local_quads.size());

    EntityHandle local_quad_ents[] = {local_quads[0], local_quads[1151], local_quads[1152], local_quads[2303]};
    rval = mb.tag_get_data(Ttag0, &local_quad_ents[0], 4, val);

    if (0 == rank) {
      CHECK_REAL_EQUAL(252.8529, val[0 * levels], eps); // First local quad, first global quad
      CHECK_REAL_EQUAL(234.8390, val[1 * levels], eps); // Median local quad, 2256th global quad
      CHECK_REAL_EQUAL(232.6458, val[2 * levels], eps); // Median local quad, 2305th global quad
      CHECK_REAL_EQUAL(205.3905, val[3 * levels], eps); // Last local quad, 4560th global quad
    }
    else if (1 == rank) {
      CHECK_REAL_EQUAL(252.7116, val[0 * levels], eps); // First local quad, 49th global quad
      CHECK_REAL_EQUAL(232.6670, val[1 * levels], eps); // Median local quad, 2304th global quad
      CHECK_REAL_EQUAL(234.6922, val[2 * levels], eps); // Median local quad, 2353th global quad
      CHECK_REAL_EQUAL(200.6828, val[3 * levels], eps); // Last local quad, last global quad
    }
  }
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

  // Check for proper tags
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_eul_nomesh() 
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

  opts = orig + std::string(";VARIABLE=T;TIMESTEP=0");
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Now read 2nd timestep with nomesh option
  opts = orig + std::string(";VARIABLE=T;TIMESTEP=1;NOMESH");
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_eul_novars() 
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
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

  opts = orig + std::string(";VARIABLE=;TIMESTEP=0");
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  opts = orig + std::string(";VARIABLE=T;TIMESTEP=0;NOMESH");
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Now read 2nd timestep with nomesh option
  opts = orig + std::string(";VARIABLE=T;TIMESTEP=1;NOMESH");
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
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

  // Check for proper tags
  Tag Ttag0, Ttag1, coordTag;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("COORDS", 3, MB_TYPE_DOUBLE, coordTag);
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

  // Check for proper tags
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);

  // Check values of tag T0 (first level) at some strategically chosen places below
  int rank = 0;
  int procs = 1;
#ifdef USE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  rank = pcomm->proc_config().proc_rank();
  procs = pcomm->proc_config().proc_size();
#endif

  const double eps = 0.0001;
  double val[8 * levels];

  if (1 == procs) {
    Range global_quads;
    rval = mb.get_entities_by_type(0, MBQUAD, global_quads);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)3312, global_quads.size());

    EntityHandle gloabl_quad_ents[] = {global_quads[0], global_quads[1619], global_quads[1656], global_quads[3275],
                                       global_quads[36], global_quads[1655], global_quads[1692], global_quads[3311]};
    rval = mb.tag_get_data(Ttag0, &gloabl_quad_ents[0], 8, val);

    CHECK_REAL_EQUAL(253.6048, val[0 * levels], eps); // First global quad
    CHECK_REAL_EQUAL(232.2170, val[1 * levels], eps); // 1620th global quad
    CHECK_REAL_EQUAL(232.7454, val[2 * levels], eps); // 1657th global quad
    CHECK_REAL_EQUAL(210.2581, val[3 * levels], eps); // 3276th global quad
    CHECK_REAL_EQUAL(253.6048, val[4 * levels], eps); // 37th global quad
    CHECK_REAL_EQUAL(232.9553, val[5 * levels], eps); // 1656th global quad
    CHECK_REAL_EQUAL(232.1704, val[6 * levels], eps); // 1693th global quad
    CHECK_REAL_EQUAL(210.2581, val[7 * levels], eps); // Last global quad
  }
  else if (2 == procs) {
    Range local_quads;
    rval = mb.get_entities_by_type(0, MBQUAD, local_quads);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)1656, local_quads.size());

    EntityHandle local_quad_ents[] = {local_quads[0], local_quads[827], local_quads[828], local_quads[1655]};
    rval = mb.tag_get_data(Ttag0, &local_quad_ents[0], 4, val);

    if (0 == rank) {
      CHECK_REAL_EQUAL(253.6048, val[0 * levels], eps); // First local quad, first global quad
      CHECK_REAL_EQUAL(232.2170, val[1 * levels], eps); // Median local quad, 1620th global quad
      CHECK_REAL_EQUAL(232.7454, val[2 * levels], eps); // Median local quad, 1657th global quad
      CHECK_REAL_EQUAL(210.2581, val[3 * levels], eps); // Last local quad, 3276th global quad
    }
    else if (1 == rank) {
      CHECK_REAL_EQUAL(253.6048, val[0 * levels], eps); // First local quad, 37th global quad
      CHECK_REAL_EQUAL(232.9553, val[1 * levels], eps); // Median local quad, 1656th global quad
      CHECK_REAL_EQUAL(232.1704, val[2 * levels], eps); // Median local quad, 1693th global quad
      CHECK_REAL_EQUAL(210.2581, val[3 * levels], eps); // Last local quad, last global quad
    }
  }
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

  // Check for proper tags
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);

  // Check for some tags with double underscore in the tag name
  Tag tempTag;
  rval = mb.tag_get_handle("__lon_LOC_MINMAX", 2, MB_TYPE_INTEGER, tempTag);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("__lon_LOC_VALS", 0, MB_TYPE_DOUBLE, tempTag, MB_TAG_VARLEN);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("__lon_GLOBAL_MINMAX", 2, MB_TYPE_INTEGER, tempTag);
  CHECK_ERR(rval);
}

void test_read_fv_nomesh() 
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

  opts = orig + std::string(";VARIABLE=T;TIMESTEP=0");
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Now read 2nd timestep with nomesh option
  opts = orig + std::string(";VARIABLE=T;TIMESTEP=1;NOMESH");
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

void test_read_fv_novars() 
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
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);

  opts = orig + std::string(";VARIABLE=;TIMESTEP=0");
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  Tag Ttag0, Ttag1;
  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  opts = orig + std::string(";VARIABLE=T;TIMESTEP=0;NOMESH");
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T0", levels, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Now read 2nd timestep with nomesh option
  opts = orig + std::string(";VARIABLE=T;TIMESTEP=1;NOMESH");
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  rval = mb.tag_get_handle("T1", levels, MB_TYPE_DOUBLE, Ttag1);
  CHECK_ERR(rval);
}

#ifdef USE_MPI
void test_read_fv_ghosting()
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

  opts = std::string("PARALLEL=READ_PART;PARTITION;PARALLEL_GHOSTS=2.0.1;NOMESH;VARIABLE=;PARTITION_METHOD=SQIJ");
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);

  opts = std::string("PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=2.0.1;PARTITION_METHOD=SQIJ;VARIABLE=");
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);

  opts = std::string("PARALLEL=READ_PART;PARTITION;PARTITION_METHOD=SQIJ;VARIABLE=TOT_CLD_VISTAU;NOMESH;TIMESTEP=0;");
  rval = mb.load_file(example_fv, &set, opts.c_str());
  CHECK_ERR(rval);
}
#endif

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
