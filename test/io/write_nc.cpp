#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "TagInfo.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example_eul[] = STRINGIFY(MESHDIR) "/io/camEul26x48x96.t3.nc";
static const char example_fv[] = STRINGIFY(MESHDIR) "/io/fv26x46x72.t.3.nc";
static const char example_homme[] = STRINGIFY(MESHDIR) "/io/homme26x3458.t.3.nc";
static const char example_homme_mapping[] = STRINGIFY(MESHDIR) "/io/HommeMapping.nc";
static const char example_mpas[] = STRINGIFY(MESHDIR) "/io/mpasx1.642.t.2.nc";
#else
static const char example_eul[] = "/io/camEul26x48x96.t3.nc";
static const char example_fv[] = "/io/fv26x46x72.t.3.nc";
static const char example_homme[] = "/io/homme26x3458.t.3.nc";
static const char example_homme_mapping[] = "/io/HommeMapping.nc";
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
void test_eul_read_write_append();

// CAM-FV
void test_fv_read_write_T();
void test_fv_check_T();

// CAM-SE (HOMME)
void test_homme_read_write_T();
void test_homme_check_T();

// MPAS
void test_mpas_read_write_vars();
void test_mpas_check_vars();

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
  result += RUN_TEST(test_eul_read_write_append);
  result += RUN_TEST(test_fv_read_write_T);
  result += RUN_TEST(test_fv_check_T);
  result += RUN_TEST(test_homme_read_write_T);
  result += RUN_TEST(test_homme_check_T);
  result += RUN_TEST(test_mpas_read_write_vars);
  result += RUN_TEST(test_mpas_check_vars);

#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail)
    return 1;
#endif

  return result;
}

// We also read and write set variable gw, which is required to create the mesh
// In test_eul_check_T(), we need to load the output file with mesh
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

  // Load non-set variable T, set variable gw, and the mesh
  read_opts += ";DEBUG_IO=0;VARIABLE=T,gw";
  rval = mb.load_file(example_eul, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variables T and gw
  std::string write_opts;
  write_opts = std::string(";;VARIABLE=T,gw;DEBUG_IO=0");
#ifdef USE_MPI
  // Use parallel options
  write_opts += std::string(";PARALLEL=WRITE_PART");
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

  Core moab;
  Interface& mb = moab;

  std::string read_opts;
  get_eul_read_options(read_opts);

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Load non-set variable T, set variable gw, and the mesh
  read_opts += ";VARIABLE=T,gw";
  if (procs > 1)
    rval = mb.load_file("test_par_eul_T.nc", &set, read_opts.c_str());
  else
    rval = mb.load_file("test_eul_T.nc", &set, read_opts.c_str());
  CHECK_ERR(rval);

  double eps = 1e-10;

  // Only check tag gw values on the root processor
  if (0 == rank) {
    // Get tag gw
    Tag gw_tag;
    rval = mb.tag_get_handle("gw", 0, MB_TYPE_OPAQUE, gw_tag, MB_TAG_SPARSE | MB_TAG_VARLEN);
    CHECK_ERR(rval);

    const void* var_data;
    int var_len = 0;
    rval = mb.tag_get_by_ptr(gw_tag, &set, 1, &var_data, &var_len);
    CHECK_ERR(rval);
    CHECK_EQUAL(48, var_len);
    double* gw_val = (double*)var_data;
    CHECK_REAL_EQUAL(0.00315334605230584, gw_val[0], eps);
    CHECK_REAL_EQUAL(0.0647376968126839, gw_val[23], eps);
    CHECK_REAL_EQUAL(0.0647376968126839, gw_val[24], eps);
    CHECK_REAL_EQUAL(0.00315334605230584, gw_val[47], eps);
  }

  // Get tag T0
  Tag Ttag0;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  eps = 0.0001;
  double val[8 * 26];

  if (1 == procs) {
    Range global_quads;
    rval = mb.get_entities_by_type(0, MBQUAD, global_quads);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)4608, global_quads.size());

    EntityHandle gloabl_quad_ents[] = {global_quads[0], global_quads[2255], global_quads[2304], global_quads[4559],
                                       global_quads[48], global_quads[2303], global_quads[2352], global_quads[4607]};
    rval = mb.tag_get_data(Ttag0, &gloabl_quad_ents[0], 8, val);

    CHECK_REAL_EQUAL(252.8529, val[0 * 26], eps); // First global quad
    CHECK_REAL_EQUAL(234.8390, val[1 * 26], eps); // 2256th global quad
    CHECK_REAL_EQUAL(232.6458, val[2 * 26], eps); // 2305th global quad
    CHECK_REAL_EQUAL(205.3905, val[3 * 26], eps); // 4560th global quad
    CHECK_REAL_EQUAL(252.7116, val[4 * 26], eps); // 49th global quad
    CHECK_REAL_EQUAL(232.6670, val[5 * 26], eps); // 2304th global quad
    CHECK_REAL_EQUAL(234.6922, val[6 * 26], eps); // 2353th global quad
    CHECK_REAL_EQUAL(200.6828, val[7 * 26], eps); // Last global quad
  }
  else if (2 == procs) {
    Range local_quads;
    rval = mb.get_entities_by_type(0, MBQUAD, local_quads);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)2304, local_quads.size());

    EntityHandle local_quad_ents[] = {local_quads[0], local_quads[1151], local_quads[1152], local_quads[2303]};
    rval = mb.tag_get_data(Ttag0, &local_quad_ents[0], 4, val);

    if (0 == rank) {
      CHECK_REAL_EQUAL(252.8529, val[0 * 26], eps); // First local quad, first global quad
      CHECK_REAL_EQUAL(234.8390, val[1 * 26], eps); // Median local quad, 2256th global quad
      CHECK_REAL_EQUAL(232.6458, val[2 * 26], eps); // Median local quad, 2305th global quad
      CHECK_REAL_EQUAL(205.3905, val[3 * 26], eps); // Last local quad, 4560th global quad
    }
    else if (1 == rank) {
      CHECK_REAL_EQUAL(252.7116, val[0 * 26], eps); // First local quad, 49th global quad
      CHECK_REAL_EQUAL(232.6670, val[1 * 26], eps); // Median local quad, 2304th global quad
      CHECK_REAL_EQUAL(234.6922, val[2 * 26], eps); // Median local quad, 2353th global quad
      CHECK_REAL_EQUAL(200.6828, val[3 * 26], eps); // Last local quad, last global quad
    }
  }
}

// We read and write variables T and U; U after we write T, so we append
// we will also write gw, just to test the file after writing
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

  // Load non-set variable T, set variable gw, and the mesh
  read_opts += ";DEBUG_IO=0;VARIABLE=T,U,gw";
  rval = mb.load_file(example_eul, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variables T and gw
  std::string write_opts;
  write_opts = std::string(";;VARIABLE=T,gw;DEBUG_IO=0");
#ifdef USE_MPI
  // Use parallel options
  write_opts += std::string(";PARALLEL=WRITE_PART");
#endif
  if (procs > 1)
    rval = mb.write_file("test_par_eul_TU.nc", 0, write_opts.c_str(), &set, 1);
  else
    rval = mb.write_file("test_eul_TU.nc", 0, write_opts.c_str(), &set, 1);
  CHECK_ERR(rval);
  // append to the file variable U
  std::string write_opts2;
  write_opts2 = std::string(";;VARIABLE=U;DEBUG_IO=0;APPEND");
#ifdef USE_MPI
  // Use parallel options
  write_opts2 += std::string(";PARALLEL=WRITE_PART");
#endif
  if (procs > 1)
    rval = mb.write_file("test_par_eul_TU.nc", 0, write_opts2.c_str(), &set, 1);
  else
    rval = mb.write_file("test_eul_TU.nc", 0, write_opts2.c_str(), &set, 1);
  CHECK_ERR(rval);
}

// We also write coordinate variables slat and slon to the output file, so that
// it can be recognized by NC reader later in test_fv_check_T()
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

  // Load non-set variable T, and the mesh
  read_opts += ";DEBUG_IO=0;VARIABLE=T";
  rval = mb.load_file(example_fv, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variables T, slat and slon
  std::string write_opts;
  write_opts = std::string(";;VARIABLE=T,slat,slon;DEBUG_IO=0");
#ifdef USE_MPI
  // Use parallel options
  write_opts += std::string(";PARALLEL=WRITE_PART");
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

  Core moab;
  Interface& mb = moab;

  std::string read_opts;
  get_eul_read_options(read_opts);

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Load non-set variable T and the mesh
  read_opts += ";VARIABLE=T";
  if (procs > 1)
    rval = mb.load_file("test_par_fv_T.nc", &set, read_opts.c_str());
  else
    rval = mb.load_file("test_fv_T.nc", &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Get tag T0
  Tag Ttag0;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  double eps = 0.0001;
  double val[8 * 26];

  if (1 == procs) {
    Range global_quads;
    rval = mb.get_entities_by_type(0, MBQUAD, global_quads);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)3312, global_quads.size());

    EntityHandle gloabl_quad_ents[] = {global_quads[0], global_quads[1619], global_quads[1656], global_quads[3275],
                                       global_quads[36], global_quads[1655], global_quads[1692], global_quads[3311]};
    rval = mb.tag_get_data(Ttag0, &gloabl_quad_ents[0], 8, val);

    CHECK_REAL_EQUAL(253.6048, val[0 * 26], eps); // First global quad
    CHECK_REAL_EQUAL(232.2170, val[1 * 26], eps); // 1620th global quad
    CHECK_REAL_EQUAL(232.7454, val[2 * 26], eps); // 1657th global quad
    CHECK_REAL_EQUAL(210.2581, val[3 * 26], eps); // 3276th global quad
    CHECK_REAL_EQUAL(253.6048, val[4 * 26], eps); // 37th global quad
    CHECK_REAL_EQUAL(232.9553, val[5 * 26], eps); // 1656th global quad
    CHECK_REAL_EQUAL(232.1704, val[6 * 26], eps); // 1693th global quad
    CHECK_REAL_EQUAL(210.2581, val[7 * 26], eps); // Last global quad
  }
  else if (2 == procs) {
    Range local_quads;
    rval = mb.get_entities_by_type(0, MBQUAD, local_quads);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)1656, local_quads.size());

    EntityHandle local_quad_ents[] = {local_quads[0], local_quads[827], local_quads[828], local_quads[1655]};
    rval = mb.tag_get_data(Ttag0, &local_quad_ents[0], 4, val);

    if (0 == rank) {
      CHECK_REAL_EQUAL(253.6048, val[0 * 26], eps); // First local quad, first global quad
      CHECK_REAL_EQUAL(232.2170, val[1 * 26], eps); // Median local quad, 1620th global quad
      CHECK_REAL_EQUAL(232.7454, val[2 * 26], eps); // Median local quad, 1657th global quad
      CHECK_REAL_EQUAL(210.2581, val[3 * 26], eps); // Last local quad, 3276th global quad
    }
    else if (1 == rank) {
      CHECK_REAL_EQUAL(253.6048, val[0 * 26], eps); // First local quad, 37th global quad
      CHECK_REAL_EQUAL(232.9553, val[1 * 26], eps); // Median local quad, 1656th global quad
      CHECK_REAL_EQUAL(232.1704, val[2 * 26], eps); // Median local quad, 1693th global quad
      CHECK_REAL_EQUAL(210.2581, val[3 * 26], eps); // Last local quad, last global quad
    }
  }
}

// We also read and write set variables lat and lon, which are required to create the mesh
// In test_homme_check_T(), we need to load the output file with mesh
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

  // Load non-set variable T, set variable lat, set variable lon, and the mesh
  read_opts += ";DEBUG_IO=0;VARIABLE=T,lat,lon";
  if (procs > 1)
    read_opts += ";PARALLEL_RESOLVE_SHARED_ENTS";
  rval = mb.load_file(example_homme, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variables T, lat and lon
  std::string write_opts = ";;VARIABLE=T,lat,lon;DEBUG_IO=0";
#ifdef USE_MPI
  // Use parallel options
  write_opts += std::string(";PARALLEL=WRITE_PART");
#endif
  if (procs > 1)
    rval = mb.write_file("test_par_homme_T.nc", 0, write_opts.c_str(), &set, 1);
  else
    rval = mb.write_file("test_homme_T.nc", 0, write_opts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

// Check non-set variable T on some vertices
// Also check set variables lat and lon
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

  Core moab;
  Interface& mb = moab;

  std::string read_opts;
  get_homme_read_options(read_opts);

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Load non-set variable T, set variable lat, set variable lon, and the mesh
  read_opts += ";VARIABLE=T,lat,lon";
  read_opts += ";CONN=";
  read_opts += example_homme_mapping;
  if (procs > 1)
    rval = mb.load_file("test_par_homme_T.nc", &set, read_opts.c_str());
  else
    rval = mb.load_file("test_homme_T.nc", &set, read_opts.c_str());
  CHECK_ERR(rval);

  double eps = 1e-10;

  if (1 == procs) {
    // Get tag lat
    Tag lat_tag;
    rval = mb.tag_get_handle("lat", 0, MB_TYPE_OPAQUE, lat_tag, MB_TAG_SPARSE | MB_TAG_VARLEN);
    CHECK_ERR(rval);

    // Check some values of tag lat
    const void* var_data;
    int var_len = 0;
    rval = mb.tag_get_by_ptr(lat_tag, &set, 1, &var_data, &var_len);
    CHECK_ERR(rval);
    CHECK_EQUAL(3458, var_len);
    double* lat_val = (double*)var_data;
    CHECK_REAL_EQUAL(-35.2643896827547, lat_val[0], eps);
    CHECK_REAL_EQUAL(23.8854752772335, lat_val[1728], eps);
    CHECK_REAL_EQUAL(29.8493120043874, lat_val[1729], eps);
    CHECK_REAL_EQUAL(38.250274171077, lat_val[3457], eps);

    // Get tag lon
    Tag lon_tag;
    rval = mb.tag_get_handle("lon", 0, MB_TYPE_OPAQUE, lon_tag, MB_TAG_SPARSE | MB_TAG_VARLEN);
    CHECK_ERR(rval);

    // Check some values of tag lon
    var_len = 0;
    rval = mb.tag_get_by_ptr(lon_tag, &set, 1, &var_data, &var_len);
    CHECK_ERR(rval);
    CHECK_EQUAL(3458, var_len);
    double* lon_val = (double*)var_data;
    CHECK_REAL_EQUAL(315, lon_val[0], eps);
    CHECK_REAL_EQUAL(202.5, lon_val[1728], eps);
    CHECK_REAL_EQUAL(194.359423525313, lon_val[1729], eps);
    CHECK_REAL_EQUAL(135, lon_val[3457], eps);
  }

  // Get tag T0
  Tag Ttag0;
  rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
  CHECK_ERR(rval);

  // Get vertices
  Range verts;
  rval = mb.get_entities_by_type(0, MBVERTEX, verts);
  CHECK_ERR(rval);

  // Get all values of tag T0
  int count;
  void* Tbuf;
  rval = mb.tag_iterate(Ttag0, verts.begin(), verts.end(), count, Tbuf);
  CHECK_ERR(rval);
  CHECK_EQUAL((size_t)count, verts.size());

  double* data = (double*) Tbuf;
  eps = 0.0001;

  if (1 == procs) {
    CHECK_EQUAL((size_t)3458, verts.size());
    CHECK_REAL_EQUAL(233.1136, data[0 * 26], eps); // First vert
    CHECK_REAL_EQUAL(236.1505, data[1728 * 26], eps); // Median vert
    CHECK_REAL_EQUAL(235.7722, data[1729 * 26], eps); // Median vert
    CHECK_REAL_EQUAL(234.0416, data[3457 * 26], eps); // Last vert
  }
  else if (2 == procs) {
    if (0 == rank) {
      CHECK_EQUAL((size_t)1825, verts.size());
      CHECK_REAL_EQUAL(233.1136, data[0 * 26], eps); // First vert
      CHECK_REAL_EQUAL(237.1977, data[912 * 26], eps); // Median vert
      CHECK_REAL_EQUAL(234.9711, data[1824 * 26], eps); // Last vert
    }
    else if (1 == rank) {
      CHECK_EQUAL((size_t)1825, verts.size());
      CHECK_REAL_EQUAL(233.1136, data[0 * 26], eps); // First vert
      CHECK_REAL_EQUAL(231.0446, data[912 * 26], eps); // Median vert
      CHECK_REAL_EQUAL(234.0416, data[1824 * 26], eps); // Last vert
    }
  }
}

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

  // Load non-set variable u, non-set variable ke, non-set variable vorticity, and the mesh
  read_opts += ";DEBUG_IO=0;VARIABLE=u,ke,vorticity";
  if (procs > 1)
    read_opts += ";PARALLEL_RESOLVE_SHARED_ENTS";
  rval = mb.load_file(example_mpas, &set, read_opts.c_str());
  CHECK_ERR(rval);

  // Write variables u, ke and vorticity (no mesh information)
  std::string write_opts = ";;VARIABLE=u,ke,vorticity;DEBUG_IO=0";
#ifdef USE_MPI
  // Use parallel options
  write_opts += std::string(";PARALLEL=WRITE_PART");
#endif
  if (procs > 1)
    rval = mb.write_file("test_par_mpas_vars.nc", 0, write_opts.c_str(), &set, 1);
  else
    rval = mb.write_file("test_mpas_vars.nc", 0, write_opts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

// Check variables without using NC reader (the output file does not contain mesh information)
// For the time being, we only check cell variable ke
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
    int ncid, varid;
    int success;
    const double eps = 1e-10;

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

    success = NCFUNC(inq_varid)(ncid, "ke", &varid);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // Enter independent I/O mode
    success = NCFUNC(begin_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    NCDF_SIZE start[] = {0, 0, 0};
    NCDF_SIZE count[] = {1, 642, 1};
    double vals[642];
    success = NCFUNC(get_vara_double)(ncid, varid, start, count, vals);
    CHECK_EQUAL(0, success);

#ifdef PNETCDF_FILE
    // End independent I/O mode
    success = NCFUNC(end_indep_data)(ncid);
    CHECK_EQUAL(0, success);
#endif

    CHECK_REAL_EQUAL(15.001, vals[0], eps);
    CHECK_REAL_EQUAL(15.012, vals[11], eps);
    CHECK_REAL_EQUAL(16.013, vals[12], eps);
    CHECK_REAL_EQUAL(16.642, vals[641], eps);

    success = NCFUNC(close)(ncid);
    CHECK_EQUAL(0, success);
  }
}

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
