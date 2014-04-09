#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "TagInfo.hpp"

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
void test_read_write_T_gw();
void test_check_T_gw_values();

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

  result += RUN_TEST(test_read_write_T_gw);
  result += RUN_TEST(test_check_T_gw_values);

#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail)
    return 1;
#endif

  return result;
}

void test_read_write_T_gw()
{
  Core moab;
  Interface& mb = moab;

  std::string orig;
  ErrorCode rval = get_options(orig);
  CHECK_ERR(rval);

  EntityHandle set;
  rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Load non-set variable T, set variable gw, and the mesh
  std::string opts = orig + std::string(";DEBUG_IO=3;VARIABLE=T,gw");
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

  // This test will write information about variable T and gw
  // To load the output file with mesh, variable gw is required
  std::string writeopts;
  writeopts = std::string(";;VARIABLE=T,gw;DEBUG_IO=2;");
  rval = mb.write_file("testTgw.nc", 0, writeopts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

void test_check_T_gw_values()
{
  Core moab;
  Interface& mb = moab;

  std::string opts;
  ErrorCode rval = get_options(opts);
  CHECK_ERR(rval);

  EntityHandle set;
  rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Load non-set variable T, set variable gw, and the mesh
  opts += std::string(";VARIABLE=T,gw");
  rval = mb.load_file("testTgw.nc", &set, opts.c_str());
  CHECK_ERR(rval);

  int procs = 1;
#ifdef USE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  procs = pcomm->proc_config().proc_size();
#endif

  // Only test serial case for the time being
  if (1 == procs) {
    // Get tag gw
    Tag gw_tag;
    rval = mb.tag_get_handle("gw", 0, MB_TYPE_OPAQUE, gw_tag, MB_TAG_SPARSE | MB_TAG_VARLEN);
    CHECK_ERR(rval);

    // Check some values of tag gw
    const void* var_data;
    int var_len = 0;
    rval = mb.tag_get_by_ptr(gw_tag, &set, 1, &var_data, &var_len);
    CHECK_ERR(rval);
    CHECK_EQUAL(48, var_len);
    double* gw_val = (double*)var_data;
    double eps = 1e-10;
    CHECK_REAL_EQUAL(0.00315334605230584, gw_val[0], eps);
    CHECK_REAL_EQUAL(0.0647376968126839, gw_val[23], eps);
    CHECK_REAL_EQUAL(0.0647376968126839, gw_val[24], eps);
    CHECK_REAL_EQUAL(0.00315334605230584, gw_val[47], eps);

    // Get tag T0
    Tag Ttag0;
    rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
    CHECK_ERR(rval);

    // Check some values of tag T0
    double val[4 * 26];
    Range global_quads;
    rval = mb.get_entities_by_type(set, MBQUAD, global_quads);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)4608, global_quads.size());
    EntityHandle gloabl_quad_ents[] = {global_quads[0], global_quads[2303], global_quads[2304], global_quads[4607]};
    rval = mb.tag_get_data(Ttag0, &gloabl_quad_ents[0], 4, val);
    CHECK_ERR(rval);
    eps = 0.0001;
    CHECK_REAL_EQUAL(252.8529, val[0 * 26], eps); // First global quad
    CHECK_REAL_EQUAL(232.6670, val[1 * 26], eps); // 2304th global quad
    CHECK_REAL_EQUAL(232.6458, val[2 * 26], eps); // 2305th global quad
    CHECK_REAL_EQUAL(200.6828, val[3 * 26], eps); // Last global quad
  }
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
