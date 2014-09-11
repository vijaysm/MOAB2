#include "TestUtil.hpp"
#include "moab/Core.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example_eul[] = STRINGIFY(MESHDIR) "/io/eul3x48x96.t.3.nc";
static const char example_eul_t0[] = STRINGIFY(MESHDIR) "/io/eul3x48x96.t0.nc";
static const char example_eul_t1[] = STRINGIFY(MESHDIR) "/io/eul3x48x96.t1.nc";
static const char example_eul_t2[] = STRINGIFY(MESHDIR) "/io/eul3x48x96.t2.nc";
static const char example_fv[] = STRINGIFY(MESHDIR) "/io/fv3x46x72.t.3.nc";
static const char example_homme[] = STRINGIFY(MESHDIR) "/io/homme3x3458.t.3.nc";
static const char example_mpas[] = STRINGIFY(MESHDIR) "/io/mpasx1.642.t.2.nc";
static const char example_gcrm[] = STRINGIFY(MESHDIR) "/io/gcrm_r3.nc";
#else
static const char example_eul[] = "/io/eul3x48x96.t.3.nc";
static const char example_eul_t0[] = "/io/eul3x48x96.t0.nc";
static const char example_eul_t1[] = "/io/eul3x48x96.t1.nc";
static const char example_eul_t2[] = "/io/eul3x48x96.t2.nc";
static const char example_fv[] = "/io/fv3x46x72.t.3.nc";
static const char example_homme[] = "/io/homme3x3458.t.3.nc";
static const char example_mpas[] = "/io/mpasx1.642.t.2.nc";
static const char example_gcrm[] = "/io/gcrm_r3.nc";
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

// GCRM
void test_gcrm_read_write_vars();
void test_gcrm_check_vars();

// Test timestep option
void test_eul_read_write_timestep();
void test_eul_check_timestep();

// Test append option
void test_eul_read_write_append();
void test_eul_check_append();

// Test writing variables with timesteps spread across files
void test_eul_read_write_across_files();
void test_eul_check_across_files();

#ifdef USE_MPI
// Test mesh with ghosted entities
void test_eul_read_write_ghosting();
void test_eul_check_ghosting();
#endif

void get_eul_read_options(std::string& opts);
void get_fv_read_options(std::string& opts);
void get_homme_read_options(std::string& opts);
void get_mpas_read_options(std::string& opts);

const double eps = 1e-10;
const int levels = 3;
const int mpas_levels = 1;

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

  result += RUN_TEST(test_gcrm_read_write_vars);
  result += RUN_TEST(test_gcrm_check_vars);

  result += RUN_TEST(test_eul_read_write_timestep);
  result += RUN_TEST(test_eul_check_timestep);

  result += RUN_TEST(test_eul_read_write_append);
  result += RUN_TEST(test_eul_check_append);

  result += RUN_TEST(test_eul_read_write_across_files);
  result += RUN_TEST(test_eul_check_across_files);

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
  read_opts += ";VARIABLE=T,gw;DEBUG_IO=3";
  rval = mb.load_file(example_eul, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variables T and gw
  std::string write_opts = ";;VARIABLE=T,gw;DEBUG_IO=3";
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

// Check non-set variable T
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
    int ncid_ref;
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

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, example_eul, NC_NOWRITE, MPI_INFO_NULL, &ncid_ref);
#else
    success = NCFUNC(open)(example_eul, NC_NOWRITE, &ncid_ref);
#endif
    CHECK_EQUAL(0, success);

    int T_id;
    success = NCFUNC(inq_varid)(ncid, "T", &T_id);
    CHECK_EQUAL(0, success);

    int T_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "T", &T_id_ref);
    CHECK_EQUAL(0, success);

    int gw_id;
    success = NCFUNC(inq_varid)(ncid, "gw", &gw_id);
    CHECK_EQUAL(0, success);

    int gw_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "gw", &gw_id_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0, 0};
    NCDF_SIZE count[] = {3, levels, 48, 96};
    const int size = 3 * levels * 48 * 96;

    // Read variable T from output file
    double T_vals[size];
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals);
    CHECK_EQUAL(0, success);

    // Read variable T from reference file
    double T_vals_ref[size];
    success = NCFUNC(get_vara_double)(ncid_ref, T_id_ref, start, count, T_vals_ref);
    CHECK_EQUAL(0, success);

    // Read variable gw (on lat) from output file
    count[0] = 48;
    double gw_vals[48];
    success = NCFUNC(get_vara_double)(ncid, gw_id, start, count, gw_vals);
    CHECK_EQUAL(0, success);

    // Read variable gw (on lat) from reference file
    double gw_vals_ref[48];
    success = NCFUNC(get_vara_double)(ncid_ref, gw_id_ref, start, count, gw_vals_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);

    success = NCFUNC(close)(ncid_ref);
    CHECK_EQUAL(0, success);

    // Check T values
    for (int i = 0; i < size; i++)
      CHECK_REAL_EQUAL(T_vals_ref[i], T_vals[i], eps);

    // Check gw values
    for (int i = 0; i < 48; i++)
      CHECK_REAL_EQUAL(gw_vals_ref[i], gw_vals[i], eps);
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

// Check non-set variable T
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
    int ncid_ref;
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

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, example_fv, NC_NOWRITE, MPI_INFO_NULL, &ncid_ref);
#else
    success = NCFUNC(open)(example_fv, NC_NOWRITE, &ncid_ref);
#endif
    CHECK_EQUAL(0, success);

    int T_id;
    success = NCFUNC(inq_varid)(ncid, "T", &T_id);
    CHECK_EQUAL(0, success);

    int T_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "T", &T_id_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0, 0};
    NCDF_SIZE count[] = {3, levels, 46, 72};
    const int size = 3 * levels * 46 * 72;

    // Read variable T from output file
    double T_vals[size];
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals);
    CHECK_EQUAL(0, success);

    // Read variable T from reference file
    double T_vals_ref[size];
    success = NCFUNC(get_vara_double)(ncid_ref, T_id_ref, start, count, T_vals_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);

    success = NCFUNC(close)(ncid_ref);
    CHECK_EQUAL(0, success);

    // Check T values
    for (int i = 0; i < size; i++)
      CHECK_REAL_EQUAL(T_vals_ref[i], T_vals[i], eps);
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
  if (procs > 1) {
    // Rotate trivial partition, otherwise localGidVertsOwned.psize() is always 1
    read_opts += ";PARALLEL_RESOLVE_SHARED_ENTS;TRIVIAL_PARTITION_SHIFT=1";
  }
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

// Check non-set variable T
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
    int ncid_ref;
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

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, example_homme, NC_NOWRITE, MPI_INFO_NULL, &ncid_ref);
#else
    success = NCFUNC(open)(example_homme, NC_NOWRITE, &ncid_ref);
#endif
    CHECK_EQUAL(0, success);

    int T_id;
    success = NCFUNC(inq_varid)(ncid, "T", &T_id);
    CHECK_EQUAL(0, success);

    int T_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "T", &T_id_ref);
    CHECK_EQUAL(0, success);

    int lat_id;
    success = NCFUNC(inq_varid)(ncid, "lat", &lat_id);
    CHECK_EQUAL(0, success);

    int lat_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "lat", &lat_id_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0};
    NCDF_SIZE count[] = {3, levels, 3458};
    const int size = 3 * levels * 3458;

    // Read variable T from output file
    double T_vals[size];
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals);
    CHECK_EQUAL(0, success);

    // Read variable T from reference file
    double T_vals_ref[size];
    success = NCFUNC(get_vara_double)(ncid_ref, T_id_ref, start, count, T_vals_ref);
    CHECK_EQUAL(0, success);

    // Read variable lat from output file
    count[0] = 3458;
    double lat_vals[3458];
    success = NCFUNC(get_vara_double)(ncid, lat_id, start, count, lat_vals);
    CHECK_EQUAL(0, success);

    // Read variable lat from reference file
    double lat_vals_ref[3458];
    success = NCFUNC(get_vara_double)(ncid_ref, lat_id_ref, start, count, lat_vals_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);

    success = NCFUNC(close)(ncid_ref);
    CHECK_EQUAL(0, success);

    // Check T values
    for (int i = 0; i < size; i++)
      CHECK_REAL_EQUAL(T_vals_ref[i], T_vals[i], eps);

    // Check gw values
    for (int i = 0; i < 3458; i++)
      CHECK_REAL_EQUAL(lat_vals_ref[i], lat_vals[i], eps);
  }
}

// Write vertex variable vorticity, edge variable u and cell variable ke
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

// Check vertex variable vorticity, edge variable u and cell variable ke
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
    int ncid_ref;
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

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, example_mpas, NC_NOWRITE, MPI_INFO_NULL, &ncid_ref);
#else
    success = NCFUNC(open)(example_mpas, NC_NOWRITE, &ncid_ref);
#endif
    CHECK_EQUAL(0, success);

    int vorticity_id;
    success = NCFUNC(inq_varid)(ncid, "vorticity", &vorticity_id);
    CHECK_EQUAL(0, success);

    int vorticity_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "vorticity", &vorticity_id_ref);
    CHECK_EQUAL(0, success);

    int u_id;
    success = NCFUNC(inq_varid)(ncid, "u", &u_id);
    CHECK_EQUAL(0, success);

    int u_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "u", &u_id_ref);
    CHECK_EQUAL(0, success);

    int ke_id;
    success = NCFUNC(inq_varid)(ncid, "ke", &ke_id);
    CHECK_EQUAL(0, success);

    int ke_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "ke", &ke_id_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0};
    NCDF_SIZE count[] = {2, 1, mpas_levels};
    const int size1 = 2 * 1280 * mpas_levels;
    const int size2 = 2 * 1920 * mpas_levels;
    const int size3 = 2 * 642 * mpas_levels;

    // Read vertex variable vorticity from output file
    count[1] = 1280;
    double vorticity_vals[size1];
    success = NCFUNC(get_vara_double)(ncid, vorticity_id, start, count, vorticity_vals);
    CHECK_EQUAL(0, success);

    // Read vertex variable vorticity from reference file
    double vorticity_vals_ref[size1];
    success = NCFUNC(get_vara_double)(ncid_ref, vorticity_id_ref, start, count, vorticity_vals_ref);
    CHECK_EQUAL(0, success);

    // Read edge variable u from output file
    count[1] = 1920;
    double u_vals[size2];
    success = NCFUNC(get_vara_double)(ncid, u_id, start, count, u_vals);
    CHECK_EQUAL(0, success);

    // Read edge variable u from reference file
    double u_vals_ref[size2];
    success = NCFUNC(get_vara_double)(ncid_ref, u_id_ref, start, count, u_vals_ref);
    CHECK_EQUAL(0, success);

    // Read cell variable ke from output file
    count[1] = 642;
    double ke_vals[size3];
    success = NCFUNC(get_vara_double)(ncid, ke_id, start, count, ke_vals);
    CHECK_EQUAL(0, success);

    // Read cell variable ke from reference file
    double ke_vals_ref[size3];
    success = NCFUNC(get_vara_double)(ncid_ref, ke_id_ref, start, count, ke_vals_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);

    success = NCFUNC(close)(ncid_ref);
    CHECK_EQUAL(0, success);

    // Check vorticity values
    for (int i = 0; i < size1; i++)
      CHECK_REAL_EQUAL(vorticity_vals_ref[i], vorticity_vals[i], eps);

    // Check u values
    for (int i = 0; i < size2; i++)
      CHECK_REAL_EQUAL(u_vals_ref[i], u_vals[i], eps);

    // Check ke values
    for (int i = 0; i < size3; i++)
      CHECK_REAL_EQUAL(ke_vals_ref[i], ke_vals[i], eps);
  }
}

// Write vertex variable u, edge variable wind, cell variable vorticity (on layers),
// and cell variable pressure (on interfaces)
void test_gcrm_read_write_vars()
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
  // we can use the same base options as mpas, because the zoltan can apply too
  get_mpas_read_options(read_opts);

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Read non-set variables u, wind, vorticity and pressure
  read_opts += ";VARIABLE=u,wind,vorticity,pressure;DEBUG_IO=0";
  if (procs > 1)
    read_opts += ";PARALLEL_RESOLVE_SHARED_ENTS";
  rval = mb.load_file(example_gcrm, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variables u, wind, vorticity and pressure
  std::string write_opts = ";;VARIABLE=u,wind,vorticity,pressure;DEBUG_IO=0";
#ifdef USE_MPI
  // Use parallel options
  write_opts += ";PARALLEL=WRITE_PART";
#endif
  if (procs > 1)
    rval = mb.write_file("test_par_gcrm_vars.nc", 0, write_opts.c_str(), &set, 1);
  else
    rval = mb.write_file("test_gcrm_vars.nc", 0, write_opts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

// Check vertex variable u, edge variable wind, cell variable vorticity (on layers),
// and cell variable pressure (on interfaces)
void test_gcrm_check_vars()
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
    int ncid_ref;
    int success;

    std::string filename;
    if (procs > 1)
      filename = "test_par_gcrm_vars.nc";
    else
      filename = "test_gcrm_vars.nc";

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, filename.c_str(), NC_NOWRITE, MPI_INFO_NULL, &ncid);
#else
    success = NCFUNC(open)(filename.c_str(), NC_NOWRITE, &ncid);
#endif
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, example_gcrm, NC_NOWRITE, MPI_INFO_NULL, &ncid_ref);
#else
    success = NCFUNC(open)(example_gcrm, NC_NOWRITE, &ncid_ref);
#endif
    CHECK_EQUAL(0, success);

    int u_id;
    success = NCFUNC(inq_varid)(ncid, "u", &u_id);
    CHECK_EQUAL(0, success);

    int u_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "u", &u_id_ref);
    CHECK_EQUAL(0, success);

    int wind_id;
    success = NCFUNC(inq_varid)(ncid, "wind", &wind_id);
    CHECK_EQUAL(0, success);

    int wind_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "wind", &wind_id_ref);
    CHECK_EQUAL(0, success);

    int vorticity_id;
    success = NCFUNC(inq_varid)(ncid, "vorticity", &vorticity_id);
    CHECK_EQUAL(0, success);

    int vorticity_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "vorticity", &vorticity_id_ref);
    CHECK_EQUAL(0, success);

    int pressure_id;
    success = NCFUNC(inq_varid)(ncid, "pressure", &pressure_id);
    CHECK_EQUAL(0, success);

    int pressure_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "pressure", &pressure_id_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0};
    NCDF_SIZE count[] = {2, 1, levels};
    const int size1 = 2 * 1280 * levels;
    const int size2 = 2 * 1920 * levels;
    const int size3 = 2 * 642 * levels;

    // Read vertex variable u from output file
    count[1] = 1280;
    double u_vals[size1];
    success = NCFUNC(get_vara_double)(ncid, u_id, start, count, u_vals);
    CHECK_EQUAL(0, success);

    // Read vertex variable u from reference file
    double u_vals_ref[size1];
    success = NCFUNC(get_vara_double)(ncid_ref, u_id_ref, start, count, u_vals_ref);
    CHECK_EQUAL(0, success);

    // Read edge variable wind from output file
    count[1] = 1920;
    double wind_vals[size2];
    success = NCFUNC(get_vara_double)(ncid, wind_id, start, count, wind_vals);
    CHECK_EQUAL(0, success);

    // Read edge variable wind from reference file
    double wind_vals_ref[size2];
    success = NCFUNC(get_vara_double)(ncid_ref, wind_id_ref, start, count, wind_vals_ref);
    CHECK_EQUAL(0, success);

    // Read cell variable vorticity from output file
    count[1] = 642;
    double vorticity_vals[size3];
    success = NCFUNC(get_vara_double)(ncid, vorticity_id, start, count, vorticity_vals);
    CHECK_EQUAL(0, success);

    // Read cell variable vorticity from reference file
    double vorticity_vals_ref[size3];
    success = NCFUNC(get_vara_double)(ncid_ref, vorticity_id_ref, start, count, vorticity_vals_ref);
    CHECK_EQUAL(0, success);

    // Read cell variable pressure from output file
    double pressure_vals[size3];
    success = NCFUNC(get_vara_double)(ncid, pressure_id, start, count, pressure_vals);
    CHECK_EQUAL(0, success);

    // Read cell variable pressure from reference file
    double pressure_vals_ref[size3];
    success = NCFUNC(get_vara_double)(ncid_ref, pressure_id_ref, start, count, pressure_vals_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);

    success = NCFUNC(close)(ncid_ref);
    CHECK_EQUAL(0, success);

    // Check u values
    for (int i = 0; i < size1; i++)
      CHECK_REAL_EQUAL(u_vals_ref[i], u_vals[i], eps);

    // Check wind values
    for (int i = 0; i < size2; i++)
      CHECK_REAL_EQUAL(wind_vals_ref[i], wind_vals[i], eps);

    // Check vorticity and pressure values
    for (int i = 0; i < size3; i++) {
      CHECK_REAL_EQUAL(vorticity_vals_ref[i], vorticity_vals[i], eps);
      CHECK_REAL_EQUAL(pressure_vals_ref[i], pressure_vals[i], eps);
    }
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
    int ncid_ref;
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

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, example_eul, NC_NOWRITE, MPI_INFO_NULL, &ncid_ref);
#else
    success = NCFUNC(open)(example_eul, NC_NOWRITE, &ncid_ref);
#endif
    CHECK_EQUAL(0, success);

    int T_id;
    success = NCFUNC(inq_varid)(ncid, "T", &T_id);
    CHECK_EQUAL(0, success);

    int T_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "T", &T_id_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0, 0};
    NCDF_SIZE count[] = {1, levels, 48, 96};
    const int size = levels * 48 * 96;

    // Read variable T from output file (timestep 0)
    double T_vals[size];
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals);
    CHECK_EQUAL(0, success);

    // Read variable T from reference file (timestep 2)
    start[0] = 2;
    double T_vals_ref[size];
    success = NCFUNC(get_vara_double)(ncid_ref, T_id_ref, start, count, T_vals_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);

    success = NCFUNC(close)(ncid_ref);
    CHECK_EQUAL(0, success);

    // Check T values
    for (int i = 0; i < size; i++)
      CHECK_REAL_EQUAL(T_vals_ref[i], T_vals[i], eps);
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
    int ncid_ref;
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

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, example_eul, NC_NOWRITE, MPI_INFO_NULL, &ncid_ref);
#else
    success = NCFUNC(open)(example_eul, NC_NOWRITE, &ncid_ref);
#endif
    CHECK_EQUAL(0, success);

    int T_id;
    success = NCFUNC(inq_varid)(ncid, "T", &T_id);
    CHECK_EQUAL(0, success);

    int T_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "T", &T_id_ref);
    CHECK_EQUAL(0, success);

    int U_id;
    success = NCFUNC(inq_varid)(ncid, "U", &U_id);
    CHECK_EQUAL(0, success);

    int U_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "U", &U_id_ref);
    CHECK_EQUAL(0, success);

    int V_id;
    success = NCFUNC(inq_varid)(ncid, "VNEWNAME", &V_id);
    CHECK_EQUAL(0, success);

    int V_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "V", &V_id_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0, 0};
    NCDF_SIZE count[] = {3, levels, 48, 96};
    const int size = 3 * levels * 48 * 96;

    // Read variable T from output file
    double T_vals[size];
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals);
    CHECK_EQUAL(0, success);

    // Read variable T from reference file
    double T_vals_ref[size];
    success = NCFUNC(get_vara_double)(ncid_ref, T_id_ref, start, count, T_vals_ref);
    CHECK_EQUAL(0, success);

    // Read variable U from output file
    double U_vals[size];
    success = NCFUNC(get_vara_double)(ncid, U_id, start, count, U_vals);
    CHECK_EQUAL(0, success);

    // Read variable U from reference file
    double U_vals_ref[size];
    success = NCFUNC(get_vara_double)(ncid_ref, U_id_ref, start, count, U_vals_ref);
    CHECK_EQUAL(0, success);

    // Read variable VNEWNAME from output file
    double V_vals[size];
    success = NCFUNC(get_vara_double)(ncid, V_id, start, count, V_vals);
    CHECK_EQUAL(0, success);

    // Read variable V from reference file
    double V_vals_ref[size];
    success = NCFUNC(get_vara_double)(ncid_ref, V_id_ref, start, count, V_vals_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);

    success = NCFUNC(close)(ncid_ref);
    CHECK_EQUAL(0, success);

    // Check T, U, and V values
    for (int i = 0; i < size; i++) {
      CHECK_REAL_EQUAL(T_vals_ref[i], T_vals[i], eps);
      CHECK_REAL_EQUAL(U_vals_ref[i], U_vals[i], eps);
      CHECK_REAL_EQUAL(V_vals_ref[i], V_vals[i], eps);
    }
  }
}

void test_eul_read_write_across_files()
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

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // This file contains single timestep 2 (option TIMESTEP=0 will be implicitly used)
  // Read T as tag T2 with option TIMESTEPBASE=2
  get_eul_read_options(read_opts);
  read_opts += ";VARIABLE=T;TIMESTEPBASE=2;DEBUG_IO=0";
  rval = mb.load_file(example_eul_t2, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // This file contains single timestep 0 (option TIMESTEP=0 will be implicitly used)
  // Read T as tag T0 with option TIMESTEPBASE=0
  get_eul_read_options(read_opts);
  read_opts += ";VARIABLE=T;TIMESTEPBASE=0;NOMESH;DEBUG_IO=0";
  rval = mb.load_file(example_eul_t0, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // This file contains single timestep 1 (option TIMESTEP=0 will be implicitly used)
  // Read T as tag T1 with option TIMESTEPBASE=1
  get_eul_read_options(read_opts);
  read_opts += ";VARIABLE=T;TIMESTEPBASE=1;NOMESH;DEBUG_IO=0";
  rval = mb.load_file(example_eul_t1, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variable T with 3 timesteps
  std::string write_opts = ";;VARIABLE=T;TIMESTEP=0,1,2;DEBUG_IO=0";
#ifdef USE_MPI
  // Use parallel options
  write_opts += ";PARALLEL=WRITE_PART";
#endif
  if (procs > 1)
    rval = mb.write_file("test_par_eul_across_files.nc", 0, write_opts.c_str(), &set, 1);
  else
    rval = mb.write_file("test_eul_across_files.nc", 0, write_opts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

void test_eul_check_across_files()
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
    int ncid_ref;
    int success;

    std::string filename;
    if (procs > 1)
      filename = "test_par_eul_across_files.nc";
    else
      filename = "test_eul_across_files.nc";

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, filename.c_str(), NC_NOWRITE, MPI_INFO_NULL, &ncid);
#else
    success = NCFUNC(open)(filename.c_str(), NC_NOWRITE, &ncid);
#endif
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, example_eul, NC_NOWRITE, MPI_INFO_NULL, &ncid_ref);
#else
    success = NCFUNC(open)(example_eul, NC_NOWRITE, &ncid_ref);
#endif
    CHECK_EQUAL(0, success);

    int T_id;
    success = NCFUNC(inq_varid)(ncid, "T", &T_id);
    CHECK_EQUAL(0, success);

    int T_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "T", &T_id_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0, 0};
    NCDF_SIZE count[] = {3, levels, 48, 96};
    const int size = 3 * levels * 48 * 96;

    // Read variable T from output file (with 3 timesteps)
    double T_vals[size];
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals);
    CHECK_EQUAL(0, success);

    // Read variable T from reference file (with 3 timesteps)
    double T_vals_ref[size];
    success = NCFUNC(get_vara_double)(ncid_ref, T_id_ref, start, count, T_vals_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);

    success = NCFUNC(close)(ncid_ref);
    CHECK_EQUAL(0, success);

    // Check T values
    for (int i = 0; i < size; i++)
      CHECK_REAL_EQUAL(T_vals_ref[i], T_vals[i], eps);
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
    int ncid_ref;
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

#ifdef PNETCDF_FILE
    success = NCFUNC(open)(MPI_COMM_SELF, example_eul, NC_NOWRITE, MPI_INFO_NULL, &ncid_ref);
#else
    success = NCFUNC(open)(example_eul, NC_NOWRITE, &ncid_ref);
#endif
    CHECK_EQUAL(0, success);

    int T_id;
    success = NCFUNC(inq_varid)(ncid, "T", &T_id);
    CHECK_EQUAL(0, success);

    int T_id_ref;
    success = NCFUNC(inq_varid)(ncid_ref, "T", &T_id_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0, 0};
    NCDF_SIZE count[] = {3, levels, 48, 96};
    const int size = 3 * levels * 48 * 96;

    // Read variable T from output file
    double T_vals[size];
    success = NCFUNC(get_vara_double)(ncid, T_id, start, count, T_vals);
    CHECK_EQUAL(0, success);

    // Read variable T from reference file
    double T_vals_ref[size];
    success = NCFUNC(get_vara_double)(ncid_ref, T_id_ref, start, count, T_vals_ref);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid_ref);
    CHECK_EQUAL(0, success);
#endif

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);

    success = NCFUNC(close)(ncid_ref);
    CHECK_EQUAL(0, success);

    // Check T values
    for (int i = 0; i < size; i++)
      CHECK_REAL_EQUAL(T_vals_ref[i], T_vals[i], eps);
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
