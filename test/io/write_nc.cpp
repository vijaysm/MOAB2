#include "TestUtil.hpp"
#include "moab/Core.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example_eul[] = STRINGIFY(MESHDIR) "/io/camEul26x48x96.t3.nc";
static const char example_fv[] = STRINGIFY(MESHDIR) "/io/fv26x46x72.t.3.nc";
static const char example_homme[] = STRINGIFY(MESHDIR) "/io/homme26x3458.t.3.nc";
static const char example_mpas[] = STRINGIFY(MESHDIR) "/io/mpasx1.642.t.2.nc";
#else
static const char example_eul[] = "/io/camEul26x48x96.t3.nc";
static const char example_fv[] = "/io/fv26x46x72.t.3.nc";
static const char example_homme[] = "/io/homme26x3458.t.3.nc";
static const char example_mpas[] = "/io/mpasx1.642.t.2.nc";
#endif

#ifdef USE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#endif

#ifdef PNETCDF_FILE
#include "pnetcdf.h"
#define NCFUNC(func) ncmpi_ ## func
#define NCDF_SIZE MPI_Offset
#else
#include "netcdf.h"
#define NCFUNC(func) nc_ ## func
#define NCDF_SIZE size_t
#endif

// CAM-EUL
void test_eul_read_write_T();
void test_eul_check_T();

// CAM-FV
void test_fv_read_write_T();
void test_fv_check_T();

// CAM-SE (HOMME)
void test_homme_read_write_T();
void test_homme_check_T();

// MPAS
void test_mpas_read_write_vars();
void test_mpas_check_vars();

// Test timestep option
void test_eul_read_write_timestep();
void test_eul_check_timestep();

// Test append option
void test_eul_read_write_append();
void test_eul_check_append();

#ifdef USE_MPI
// Test mesh with ghosted entities
void test_eul_read_write_ghosting();
void test_eul_check_ghosting();
#endif

void get_eul_read_options(std::string& opts);
void get_fv_read_options(std::string& opts);
void get_homme_read_options(std::string& opts);
void get_mpas_read_options(std::string& opts);

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

  result += RUN_TEST(test_eul_read_write_T);
  result += RUN_TEST(test_eul_check_T);

  result += RUN_TEST(test_fv_read_write_T);
  result += RUN_TEST(test_fv_check_T);

  result += RUN_TEST(test_homme_read_write_T);
  result += RUN_TEST(test_homme_check_T);

  result += RUN_TEST(test_mpas_read_write_vars);
  result += RUN_TEST(test_mpas_check_vars);

  result += RUN_TEST(test_eul_read_write_timestep);
  result += RUN_TEST(test_eul_check_timestep);

  result += RUN_TEST(test_eul_read_write_append);
  result += RUN_TEST(test_eul_check_append);

#ifdef USE_MPI
  result += RUN_TEST(test_eul_read_write_ghosting);
  result += RUN_TEST(test_eul_check_ghosting);
#endif

#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail)
    return 1;
#endif

  return result;
}

// We also read and write gw (test writing a set variable without timesteps)
void test_eul_read_write_T()
{
  int procs = 1;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
#endif

// We will not test NC writer in parallel without pnetcdf support
#ifndef PNETCDF_FILE
  if (procs > 1)
    return;
#endif

  Core moab;
  Interface& mb = moab;

  std::string read_opts;
  get_eul_read_options(read_opts);

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Read non-set variable T and set variable gw
  read_opts += ";VARIABLE=T,gw;DEBUG_IO=0";
  rval = mb.load_file(example_eul, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variables T and gw
  std::string write_opts = ";;VARIABLE=T,gw;DEBUG_IO=0";
#ifdef USE_MPI
  // Use parallel options
  write_opts += ";PARALLEL=WRITE_PART";
#endif
  if (procs > 1)
    rval = mb.write_file("test_par_eul_T.nc", 0, write_opts.c_str(), &set, 1);
  else
    rval = mb.write_file("test_eul_T.nc", 0, write_opts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

// Check non-set variable T on some quads
// Also check set variable gw
void test_eul_check_T()
{
  int rank = 0;
  int procs = 1;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
#endif

// We will not test NC writer in parallel without pnetcdf support
#ifndef PNETCDF_FILE
  if (procs > 1)
    return;
#endif

  if (0 == rank) {
    int ncid;
    int success;

    std::string filename;
    if (procs > 1)
      filename = "test_par_eul_T.nc";
    else
      filename = "test_eul_T.nc";

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, filename.c_str(), NC_NOWRITE, MPI_INFO_NULL, &ncid);
#else
    success = NCFUNC(open)(filename.c_str(), NC_NOWRITE, &ncid);
#endif
    CHECK_EQUAL(0, success);

    int T_id;
    success = NCFUNC(inq_varid)(ncid, "T", &T_id);
    CHECK_EQUAL(0, success);

    int gw_id;
    success = NCFUNC(inq_varid)(ncid, "gw", &gw_id);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0, 0};
    NCDF_SIZE count[] = {2, 1, 48, 96}; // Read two timesteps and one level

    // Read variable T on 48 * 96 quads (first level)
    double T_vals_lev1[2 * 48 * 96];
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals_lev1);
    CHECK_EQUAL(0, success);

    // Read variable T on 48 * 96 quads (last level)
    double T_vals_lev26[2 * 48 * 96];
    start[1] = 25;
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals_lev26);
    CHECK_EQUAL(0, success);

    // Read variable gw on lat
    double gw_vals[48];
    count[0] = 48;
    success = NCFUNC(get_vara_double)(ncid, gw_id, start, count, gw_vals);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    double eps = 0.0001;

    // Check T values at some strategically chosen places (first level)
    // Timestep 0
    CHECK_REAL_EQUAL(252.8529, T_vals_lev1[0], eps); // First quad
    CHECK_REAL_EQUAL(232.6670, T_vals_lev1[2303], eps); // Median quad
    CHECK_REAL_EQUAL(232.6458, T_vals_lev1[2304], eps); // Median quad
    CHECK_REAL_EQUAL(200.6828, T_vals_lev1[4607], eps); // Last quad
    // Timestep 1
    CHECK_REAL_EQUAL(241.7352, T_vals_lev1[0 + 4608], eps); // First quad
    CHECK_REAL_EQUAL(234.7536, T_vals_lev1[2303 + 4608], eps); // Median quad
    CHECK_REAL_EQUAL(234.4739, T_vals_lev1[2304 + 4608], eps); // Median quad
    CHECK_REAL_EQUAL(198.2482, T_vals_lev1[4607 + 4608], eps); // Last quad

    // Check T values at some strategically chosen places (last level)
    // Timestep 0
    CHECK_REAL_EQUAL(253.1395, T_vals_lev26[0], eps); // First quad
    CHECK_REAL_EQUAL(299.0477, T_vals_lev26[2303], eps); // Median quad
    CHECK_REAL_EQUAL(300.0627, T_vals_lev26[2304], eps); // Median quad
    CHECK_REAL_EQUAL(241.1817, T_vals_lev26[4607], eps); // Last quad
    // Timestep 1
    CHECK_REAL_EQUAL(242.9252, T_vals_lev26[0 + 4608], eps); // First quad
    CHECK_REAL_EQUAL(299.9290, T_vals_lev26[2303 + 4608], eps); // Median quad
    CHECK_REAL_EQUAL(299.7614, T_vals_lev26[2304 + 4608], eps); // Median quad
    CHECK_REAL_EQUAL(241.1057, T_vals_lev26[4607 + 4608], eps); // Last quad

    eps = 1e-10;

    // Check gw values at some strategically chosen places
    CHECK_REAL_EQUAL(0.00315334605230584, gw_vals[0], eps);
    CHECK_REAL_EQUAL(0.0647376968126839, gw_vals[23], eps);
    CHECK_REAL_EQUAL(0.0647376968126839, gw_vals[24], eps);
    CHECK_REAL_EQUAL(0.00315334605230584, gw_vals[47], eps);

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);
  }
}

void test_fv_read_write_T()
{
  int procs = 1;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
#endif

// We will not test NC writer in parallel without pnetcdf support
#ifndef PNETCDF_FILE
  if (procs > 1)
    return;
#endif

  Core moab;
  Interface& mb = moab;

  std::string read_opts;
  get_fv_read_options(read_opts);

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Read non-set variable T
  read_opts += ";VARIABLE=T;DEBUG_IO=0";
  rval = mb.load_file(example_fv, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variable T
  std::string write_opts = ";;VARIABLE=T;DEBUG_IO=0";
#ifdef USE_MPI
  // Use parallel options
  write_opts += ";PARALLEL=WRITE_PART";
#endif
  if (procs > 1)
    rval = mb.write_file("test_par_fv_T.nc", 0, write_opts.c_str(), &set, 1);
  else
    rval = mb.write_file("test_fv_T.nc", 0, write_opts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

// Check non-set variable T on some quads
void test_fv_check_T()
{
  int rank = 0;
  int procs = 1;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
#endif

// We will not test NC writer in parallel without pnetcdf support
#ifndef PNETCDF_FILE
  if (procs > 1)
    return;
#endif

  if (0 == rank) {
    int ncid;
    int success;

    std::string filename;
    if (procs > 1)
      filename = "test_par_fv_T.nc";
    else
      filename = "test_fv_T.nc";

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, filename.c_str(), NC_NOWRITE, MPI_INFO_NULL, &ncid);
#else
    success = NCFUNC(open)(filename.c_str(), NC_NOWRITE, &ncid);
#endif
    CHECK_EQUAL(0, success);

    int T_id;
    success = NCFUNC(inq_varid)(ncid, "T", &T_id);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0, 0};
    NCDF_SIZE count[] = {2, 1, 46, 72}; // Read two timesteps and one level

    // Read variable T on 46 * 72 quads (first level)
    double T_vals_lev1[2 * 46 * 72];
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals_lev1);
    CHECK_EQUAL(0, success);

    // Read variable T on 46 * 72 quads (last level)
    double T_vals_lev26[2 * 46 * 72];
    start[1] = 25;
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals_lev26);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    const double eps = 0.0001;

    // Check T values at some strategically chosen places (first level)
    // Timestep 0
    CHECK_REAL_EQUAL(253.6048, T_vals_lev1[0], eps); // First quad
    CHECK_REAL_EQUAL(232.9553, T_vals_lev1[1655], eps); // Median quad
    CHECK_REAL_EQUAL(232.7454, T_vals_lev1[1656], eps); // Median quad
    CHECK_REAL_EQUAL(210.2581, T_vals_lev1[3311], eps); // Last quad
    // Timestep 1
    CHECK_REAL_EQUAL(242.4844, T_vals_lev1[0 + 3312], eps); // First quad
    CHECK_REAL_EQUAL(234.0176, T_vals_lev1[1655 + 3312], eps); // Median quad
    CHECK_REAL_EQUAL(233.8797, T_vals_lev1[1656 + 3312], eps); // Median quad
    CHECK_REAL_EQUAL(207.3904, T_vals_lev1[3311 + 3312], eps); // Last quad

    // Check T values at some strategically chosen places (last level)
    // Timestep 0
    CHECK_REAL_EQUAL(244.5516, T_vals_lev26[0], eps); // First quad
    CHECK_REAL_EQUAL(297.2558, T_vals_lev26[1655], eps); // Median quad
    CHECK_REAL_EQUAL(295.2663, T_vals_lev26[1656], eps); // Median quad
    CHECK_REAL_EQUAL(244.5003, T_vals_lev26[3311], eps); // Last quad
    // Timestep 1
    CHECK_REAL_EQUAL(238.8134, T_vals_lev26[0 + 3312], eps); // First quad
    CHECK_REAL_EQUAL(297.9755, T_vals_lev26[1655 + 3312], eps); // Median quad
    CHECK_REAL_EQUAL(296.1439, T_vals_lev26[1656 + 3312], eps); // Median quad
    CHECK_REAL_EQUAL(242.1957, T_vals_lev26[3311 + 3312], eps); // Last quad

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);
  }
}

// We also read and write lat (test writing a set variable without timesteps)
void test_homme_read_write_T()
{
  int procs = 1;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
#endif

// We will not test NC writer in parallel without pnetcdf support
#ifndef PNETCDF_FILE
  if (procs > 1)
    return;
#endif

  Core moab;
  Interface& mb = moab;

  std::string read_opts;
  get_homme_read_options(read_opts);

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Read non-set variable T and set variable lat
  read_opts += ";VARIABLE=T,lat;DEBUG_IO=0";
  if (procs > 1)
    read_opts += ";PARALLEL_RESOLVE_SHARED_ENTS";
  rval = mb.load_file(example_homme, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variables T and lat
  std::string write_opts = ";;VARIABLE=T,lat;DEBUG_IO=0";
#ifdef USE_MPI
  // Use parallel options
  write_opts += ";PARALLEL=WRITE_PART";
#endif
  if (procs > 1)
    rval = mb.write_file("test_par_homme_T.nc", 0, write_opts.c_str(), &set, 1);
  else
    rval = mb.write_file("test_homme_T.nc", 0, write_opts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

// Check non-set variable T on some vertices
// Also check set variable lat
void test_homme_check_T()
{
  int rank = 0;
  int procs = 1;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
#endif

// We will not test NC writer in parallel without pnetcdf support
#ifndef PNETCDF_FILE
  if (procs > 1)
    return;
#endif

  if (0 == rank) {
    int ncid;
    int success;

    std::string filename;
    if (procs > 1)
      filename = "test_par_homme_T.nc";
    else
      filename = "test_homme_T.nc";

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, filename.c_str(), NC_NOWRITE, MPI_INFO_NULL, &ncid);
#else
    success = NCFUNC(open)(filename.c_str(), NC_NOWRITE, &ncid);
#endif
    CHECK_EQUAL(0, success);

    int T_id;
    success = NCFUNC(inq_varid)(ncid, "T", &T_id);
    CHECK_EQUAL(0, success);

    int lat_id;
    success = NCFUNC(inq_varid)(ncid, "lat", &lat_id);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0};
    NCDF_SIZE count[] = {2, 1, 3458}; // Read two timesteps and one level

    // Read variable T on 3458 vertices (first level)
    double T_vals_lev1[2 * 3458];
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals_lev1);
    CHECK_EQUAL(0, success);

    // Read variable T on 3458 vertices (last level)
    double T_vals_lev26[2 * 3458];
    start[1] = 25;
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals_lev26);
    CHECK_EQUAL(0, success);

    // Read variable lat
    double lat_vals[3458];
    count[0] = 3458;
    success = NCFUNC(get_vara_double)(ncid, lat_id, start, count, lat_vals);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    double eps = 0.0001;

    // Check T values at some strategically chosen places (first level)
    // Timestep 0
    CHECK_REAL_EQUAL(233.1136, T_vals_lev1[0], eps); // First vertex
    CHECK_REAL_EQUAL(236.1505, T_vals_lev1[1728], eps); // Median vertex
    CHECK_REAL_EQUAL(235.7722, T_vals_lev1[1729], eps); // Median vertex
    CHECK_REAL_EQUAL(234.0416, T_vals_lev1[3457], eps); // Last vertex
    // Timestep 1
    CHECK_REAL_EQUAL(234.6015, T_vals_lev1[0 + 3458], eps); // First vertex
    CHECK_REAL_EQUAL(236.3139, T_vals_lev1[1728 + 3458], eps); // Median vertex
    CHECK_REAL_EQUAL(235.3373, T_vals_lev1[1729 + 3458], eps); // Median vertex
    CHECK_REAL_EQUAL(233.6020, T_vals_lev1[3457 + 3458], eps); // Last vertex

    // Check T values at some strategically chosen places (last level)
    // Timestep 0
    CHECK_REAL_EQUAL(281.5824, T_vals_lev26[0], eps); // First vertex
    CHECK_REAL_EQUAL(290.0607, T_vals_lev26[1728], eps); // Median vertex
    CHECK_REAL_EQUAL(285.7311, T_vals_lev26[1729], eps); // Median vertex
    CHECK_REAL_EQUAL(281.3145, T_vals_lev26[3457], eps); // Last vertex
    // Timestep 1
    CHECK_REAL_EQUAL(281.4600, T_vals_lev26[0 + 3458], eps); // First vertex
    CHECK_REAL_EQUAL(289.1411, T_vals_lev26[1728 + 3458], eps); // Median vertex
    CHECK_REAL_EQUAL(285.6183, T_vals_lev26[1729 + 3458], eps); // Median vertex
    CHECK_REAL_EQUAL(279.7856, T_vals_lev26[3457 + 3458], eps); // Last vertex

    eps = 1e-10;

    // Check lat values at some strategically chosen places
    CHECK_REAL_EQUAL(-35.2643896827547, lat_vals[0], eps); // First vertex
    CHECK_REAL_EQUAL(23.8854752772335, lat_vals[1728], eps); // Median vertex
    CHECK_REAL_EQUAL(29.8493120043874, lat_vals[1729], eps); // Median vertex
    CHECK_REAL_EQUAL(38.250274171077, lat_vals[3457], eps); // Last vertex

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);
  }
}

// Write vertex variable vorticity, edge variable u and cell veriable ke
void test_mpas_read_write_vars()
{
  int procs = 1;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
#endif

// We will not test NC writer in parallel without pnetcdf support
#ifndef PNETCDF_FILE
  if (procs > 1)
    return;
#endif

  Core moab;
  Interface& mb = moab;

  std::string read_opts;
  get_mpas_read_options(read_opts);

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Read non-set variables vorticity, u and ke
  read_opts += ";VARIABLE=vorticity,u,ke;DEBUG_IO=0";
  if (procs > 1)
    read_opts += ";PARALLEL_RESOLVE_SHARED_ENTS";
  rval = mb.load_file(example_mpas, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variables vorticity, u and ke
  std::string write_opts = ";;VARIABLE=vorticity,u,ke;DEBUG_IO=0";
#ifdef USE_MPI
  // Use parallel options
  write_opts += ";PARALLEL=WRITE_PART";
#endif
  if (procs > 1)
    rval = mb.write_file("test_par_mpas_vars.nc", 0, write_opts.c_str(), &set, 1);
  else
    rval = mb.write_file("test_mpas_vars.nc", 0, write_opts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

// Check vertex variable vorticity, edge variable u and cell veriable ke
void test_mpas_check_vars()
{
  int rank = 0;
  int procs = 1;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
#endif

// We will not test NC writer in parallel without pnetcdf support
#ifndef PNETCDF_FILE
  if (procs > 1)
    return;
#endif

  if (0 == rank) {
    int ncid;
    int success;

    std::string filename;
    if (procs > 1)
      filename = "test_par_mpas_vars.nc";
    else
      filename = "test_mpas_vars.nc";

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, filename.c_str(), NC_NOWRITE, MPI_INFO_NULL, &ncid);
#else
    success = NCFUNC(open)(filename.c_str(), NC_NOWRITE, &ncid);
#endif
    CHECK_EQUAL(0, success);

    int vorticity_id;
    success = NCFUNC(inq_varid)(ncid, "vorticity", &vorticity_id);
    CHECK_EQUAL(0, success);

    int u_id;
    success = NCFUNC(inq_varid)(ncid, "u", &u_id);
    CHECK_EQUAL(0, success);

    int ke_id;
    success = NCFUNC(inq_varid)(ncid, "ke", &ke_id);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0};
    NCDF_SIZE count[] = {2, 1, 1}; // Read two timesteps and one level

    // Read variable vorticity on 1st and 2nd vertices
    count[1] = 2;
    double vorticity_vals[4];
    success = NCFUNC(get_vara_double)(ncid, vorticity_id, start, count, vorticity_vals);
    CHECK_EQUAL(0, success);

    // Read variable u on 6th and 7th edges
    start[1] = 5;
    count[1] = 2;
    double u_vals[4];
    success = NCFUNC(get_vara_double)(ncid, u_id, start, count, u_vals);
    CHECK_EQUAL(0, success);

    // Read variable ke on all 642 cells
    start[1] = 0;
    count[1] = 642;
    double ke_vals[642 * 2];
    success = NCFUNC(get_vara_double)(ncid, ke_id, start, count, ke_vals);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    const double eps = 1e-10;

    // Check vorticity values on 1st and 2nd vertices
    // Timestep 0
    CHECK_REAL_EQUAL(1.1, vorticity_vals[0], eps);
    CHECK_REAL_EQUAL(1.2, vorticity_vals[1], eps);
    // Timestep 1
    CHECK_REAL_EQUAL(2.1, vorticity_vals[2], eps);
    CHECK_REAL_EQUAL(2.2, vorticity_vals[3], eps);

    // Check u values on 6th and 7th edges
    // Timestep 0
    CHECK_REAL_EQUAL(1.11313872154478, u_vals[0], eps);
    CHECK_REAL_EQUAL(-1.113138721930009, u_vals[1], eps);
    // Timestep 1
    CHECK_REAL_EQUAL(2.113138721544778, u_vals[2], eps);
    CHECK_REAL_EQUAL(-2.113138721930009, u_vals[3], eps);

    // Check ke values on first pentagon, last pentagon, first hexagon, and last hexagon
    // Timestep 0
    CHECK_REAL_EQUAL(15.001, ke_vals[0], eps);
    CHECK_REAL_EQUAL(15.012, ke_vals[11], eps);
    CHECK_REAL_EQUAL(16.013, ke_vals[12], eps);
    CHECK_REAL_EQUAL(16.642, ke_vals[641], eps);
    // Timestep 1
    CHECK_REAL_EQUAL(25.001, ke_vals[0 + 642], eps);
    CHECK_REAL_EQUAL(25.012, ke_vals[11 + 642], eps);
    CHECK_REAL_EQUAL(26.013, ke_vals[12 + 642], eps);
    CHECK_REAL_EQUAL(26.642, ke_vals[641 + 642], eps);

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);
  }
}

// Read non-set variable T on all 3 timesteps, and write only timestep 2
void test_eul_read_write_timestep()
{
  int procs = 1;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
#endif

// We will not test NC writer in parallel without pnetcdf support
#ifndef PNETCDF_FILE
  if (procs > 1)
    return;
#endif

  Core moab;
  Interface& mb = moab;

  std::string read_opts;
  get_eul_read_options(read_opts);

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Read non-set variable T
  read_opts += ";VARIABLE=T;DEBUG_IO=0";
  rval = mb.load_file(example_eul, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variable T on timestep 2
  std::string write_opts = ";;VARIABLE=T;TIMESTEP=2;DEBUG_IO=0";
#ifdef USE_MPI
  // Use parallel options
  write_opts += ";PARALLEL=WRITE_PART";
#endif
  if (procs > 1)
    rval = mb.write_file("test_par_eul_T2.nc", 0, write_opts.c_str(), &set, 1);
  else
    rval = mb.write_file("test_eul_T2.nc", 0, write_opts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

void test_eul_check_timestep()
{
  int rank = 0;
  int procs = 1;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
#endif

// We will not test NC writer in parallel without pnetcdf support
#ifndef PNETCDF_FILE
  if (procs > 1)
    return;
#endif

  if (0 == rank) {
    int ncid;
    int success;

    std::string filename;
    if (procs > 1)
      filename = "test_par_eul_T2.nc";
    else
      filename = "test_eul_T2.nc";

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, filename.c_str(), NC_NOWRITE, MPI_INFO_NULL, &ncid);
#else
    success = NCFUNC(open)(filename.c_str(), NC_NOWRITE, &ncid);
#endif
    CHECK_EQUAL(0, success);

    int T_id;
    success = NCFUNC(inq_varid)(ncid, "T", &T_id);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0, 0};
    NCDF_SIZE count[] = {1, 1, 48, 96}; // Read one timestep and one level

    // Read variable T on 48 * 96 quads (first level)
    double T_vals_lev1[48 * 96];
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals_lev1);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    double eps = 0.0001;

    // Check T values at some strategically chosen places (first level)
    // Timestep 2
    CHECK_REAL_EQUAL(224.1966, T_vals_lev1[0], eps); // First quad
    CHECK_REAL_EQUAL(236.1357, T_vals_lev1[2303], eps); // Median quad
    CHECK_REAL_EQUAL(235.9430, T_vals_lev1[2304], eps); // Median quad
    CHECK_REAL_EQUAL(218.7719, T_vals_lev1[4607], eps); // Last quad

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);
  }
}

// Read non-set variables T, U and V
// Write variable T, append U, and then append V (with a new name)
void test_eul_read_write_append()
{
  int procs = 1;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
#endif

// We will not test NC writer in parallel without pnetcdf support
#ifndef PNETCDF_FILE
  if (procs > 1)
    return;
#endif

  Core moab;
  Interface& mb = moab;

  std::string read_opts;
  get_eul_read_options(read_opts);

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Load non-set variables T, U, V, and the mesh
  read_opts += ";VARIABLE=T,U,V;DEBUG_IO=0";
  rval = mb.load_file(example_eul, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variable T
  std::string write_opts = ";;VARIABLE=T;DEBUG_IO=0";
#ifdef USE_MPI
  // Use parallel options
  write_opts += ";PARALLEL=WRITE_PART";
#endif
  if (procs > 1)
    rval = mb.write_file("test_par_eul_append.nc", 0, write_opts.c_str(), &set, 1);
  else
    rval = mb.write_file("test_eul_append.nc", 0, write_opts.c_str(), &set, 1);
  CHECK_ERR(rval);

  // Append to the file variable U
  write_opts = ";;VARIABLE=U;APPEND;DEBUG_IO=0";
#ifdef USE_MPI
  // Use parallel options
  write_opts += ";PARALLEL=WRITE_PART";
#endif
  if (procs > 1)
    rval = mb.write_file("test_par_eul_append.nc", 0, write_opts.c_str(), &set, 1);
  else
    rval = mb.write_file("test_eul_append.nc", 0, write_opts.c_str(), &set, 1);
  CHECK_ERR(rval);

  // Append to the file variable V, renamed to VNEWNAME
  write_opts = ";;VARIABLE=V;RENAME=VNEWNAME;APPEND;DEBUG_IO=0";
#ifdef USE_MPI
  // Use parallel options
  write_opts += ";PARALLEL=WRITE_PART";
#endif
  if (procs > 1)
    rval = mb.write_file("test_par_eul_append.nc", 0, write_opts.c_str(), &set, 1);
  else
    rval = mb.write_file("test_eul_append.nc", 0, write_opts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

void test_eul_check_append()
{
  int rank = 0;
  int procs = 1;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
#endif

// We will not test NC writer in parallel without pnetcdf support
#ifndef PNETCDF_FILE
  if (procs > 1)
    return;
#endif

  if (0 == rank) {
    int ncid;
    int success;

    std::string filename;
    if (procs > 1)
      filename = "test_par_eul_append.nc";
    else
      filename = "test_eul_append.nc";

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, filename.c_str(), NC_NOWRITE, MPI_INFO_NULL, &ncid);
#else
    success = NCFUNC(open)(filename.c_str(), NC_NOWRITE, &ncid);
#endif
    CHECK_EQUAL(0, success);

    int T_id;
    success = NCFUNC(inq_varid)(ncid, "T", &T_id);
    CHECK_EQUAL(0, success);

    int U_id;
    success = NCFUNC(inq_varid)(ncid, "U", &U_id);
    CHECK_EQUAL(0, success);

    int V_id;
    success = NCFUNC(inq_varid)(ncid, "VNEWNAME", &V_id);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0, 0};
    NCDF_SIZE count[] = {1, 1, 48, 96}; // Read one timestep and one level

    // Read variable T on 48 * 96 quads (first level)
    double T_vals_lev1[48 * 96];
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals_lev1);
    CHECK_EQUAL(0, success);

    // Read variable U on 48 * 96 quads (first level)
    double U_vals_lev1[48 * 96];
    success = NCFUNC(get_vara_double)(ncid, U_id, start, count, U_vals_lev1);
    CHECK_EQUAL(0, success);

    // Read variable V on 48 * 96 quads (first level)
    double V_vals_lev1[48 * 96];
    success = NCFUNC(get_vara_double)(ncid, V_id, start, count, V_vals_lev1);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    double eps = 0.0001;

    // Check T values at some strategically chosen places (first level)
    // Timestep 0
    CHECK_REAL_EQUAL(252.8529, T_vals_lev1[0], eps); // First quad
    CHECK_REAL_EQUAL(232.6670, T_vals_lev1[2303], eps); // Median quad
    CHECK_REAL_EQUAL(232.6458, T_vals_lev1[2304], eps); // Median quad
    CHECK_REAL_EQUAL(200.6828, T_vals_lev1[4607], eps); // Last quad

    eps = 1e-6;

    // Check U values at some strategically chosen places (first level)
    // Timestep 0
    CHECK_REAL_EQUAL(0.6060912, U_vals_lev1[0], eps); // First quad
    CHECK_REAL_EQUAL(-36.789986, U_vals_lev1[2303], eps); // Median quad
    CHECK_REAL_EQUAL(-31.429073, U_vals_lev1[2304], eps); // Median quad
    CHECK_REAL_EQUAL(-48.085426, U_vals_lev1[4607], eps); // Last quad

    eps = 1e-5;

    // Check V values at some strategically chosen places (first level)
    // Timestep 0
    CHECK_REAL_EQUAL(-0.44605, V_vals_lev1[0], eps); // First quad
    CHECK_REAL_EQUAL(0.89077, V_vals_lev1[2303], eps); // Median quad
    CHECK_REAL_EQUAL(1.141688, V_vals_lev1[2304], eps); // Median quad
    CHECK_REAL_EQUAL(-38.21262, V_vals_lev1[4607], eps); // Last quad

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);
  }
}

#ifdef USE_MPI
// NC writer should filter entities that are not owned, e.g. ghosted elements
void test_eul_read_write_ghosting()
{
  int procs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

// We will not test NC writer in parallel without pnetcdf support
#ifndef PNETCDF_FILE
  if (procs > 1)
    return;
#endif

  Core moab;
  Interface& mb = moab;

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  std::string read_opts = "PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=2.0.1;PARTITION_METHOD=SQIJ;VARIABLE=";
  rval = mb.load_file(example_eul, &set, read_opts.c_str());
  CHECK_ERR(rval);

  read_opts = "PARALLEL=READ_PART;PARTITION;PARTITION_METHOD=SQIJ;VARIABLE=T;NOMESH";
  rval = mb.load_file(example_eul, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variable T
  std::string write_opts = ";;PARALLEL=WRITE_PART;VARIABLE=T;DEBUG_IO=0";
  if (procs > 1)
    rval = mb.write_file("test_par_eul_ghosting.nc", 0, write_opts.c_str(), &set, 1);
  else
    rval = mb.write_file("test_eul_ghosting.nc", 0, write_opts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

void test_eul_check_ghosting()
{
  int rank = 0;
  int procs = 1;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
#endif

// We will not test NC writer in parallel without pnetcdf support
#ifndef PNETCDF_FILE
  if (procs > 1)
    return;
#endif

  if (0 == rank) {
    int ncid;
    int success;

    std::string filename;
    if (procs > 1)
      filename = "test_par_eul_ghosting.nc";
    else
      filename = "test_eul_ghosting.nc";

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, filename.c_str(), NC_NOWRITE, MPI_INFO_NULL, &ncid);
#else
    success = NCFUNC(open)(filename.c_str(), NC_NOWRITE, &ncid);
#endif
    CHECK_EQUAL(0, success);

    int T_id;
    success = NCFUNC(inq_varid)(ncid, "T", &T_id);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0, 0};
    NCDF_SIZE count[] = {1, 1, 48, 96}; // Read one timesteps and one level

    // Read variable T on 48 * 96 quads (first level)
    double T_vals_lev1[48 * 96];
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals_lev1);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    const double eps = 0.0001;

    // Check T values at some strategically chosen places (first level)
    // Timestep 0
    CHECK_REAL_EQUAL(252.8529, T_vals_lev1[0], eps); // First quad
    CHECK_REAL_EQUAL(232.6670, T_vals_lev1[2303], eps); // Median quad
    CHECK_REAL_EQUAL(232.6458, T_vals_lev1[2304], eps); // Median quad
    CHECK_REAL_EQUAL(200.6828, T_vals_lev1[4607], eps); // Last quad

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);
  }
}
#endif

void get_eul_read_options(std::string& opts)
{
#ifdef USE_MPI
  // Use parallel options
  opts = ";;PARALLEL=READ_PART;PARTITION_METHOD=SQIJ";
#else
  opts = ";;";
#endif
}

void get_fv_read_options(std::string& opts)
{
#ifdef USE_MPI
  // Use parallel options
  opts = ";;PARALLEL=READ_PART;PARTITION_METHOD=SQIJ";
#else
  opts = ";;";
#endif
}

void get_homme_read_options(std::string& opts)
{
#ifdef USE_MPI
  // Use parallel options
  opts = ";;PARALLEL=READ_PART;PARTITION_METHOD=TRIVIAL";
#else
  opts = ";;";
#endif
}

void get_mpas_read_options(std::string& opts)
{
#ifdef USE_MPI
  // Use parallel options
#ifdef HAVE_ZOLTAN
  opts = ";;PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN";
#else
  opts = ";;PARALLEL=READ_PART;PARTITION_METHOD=TRIVIAL";
#endif
#else
  opts = ";;";
#endif
}
