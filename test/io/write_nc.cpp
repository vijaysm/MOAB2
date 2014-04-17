#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "TagInfo.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example_eul[] = STRINGIFY(MESHDIR) "/io/camEul26x48x96.t3.nc";
static const char example_homme[] = STRINGIFY(MESHDIR) "/io/homme26x3458.t.3.nc";
static const char example_homme_mapping[] = STRINGIFY(MESHDIR) "/io/HommeMapping.nc";
#else
static const char example_eul[] = "/io/camEul26x48x96.t3.nc";
static const char example_homme[] = "/io/homme26x3458.t.3.nc";
static const char example_homme_mapping[] = "/io/HommeMapping.nc";
#endif

#ifdef USE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#endif

// CAM-EUL
void test_eul_read_write_T();
void test_eul_check_T();

// CAM-SE (HOMME)
void test_homme_read_write_T();
void test_homme_check_T();

void get_eul_options(std::string& opts);
void get_homme_options(std::string& opts);

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
  result += RUN_TEST(test_homme_read_write_T);
  result += RUN_TEST(test_homme_check_T);

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
  Core moab;
  Interface& mb = moab;

  std::string orig;
  get_eul_options(orig);

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Load non-set variable T, set variable gw, and the mesh
  std::string opts = orig + std::string(";DEBUG_IO=0;VARIABLE=T,gw");
  rval = mb.load_file(example_eul, &set, opts.c_str());
  CHECK_ERR(rval);

  // Write variables T and gw
  std::string writeopts;
  writeopts = std::string(";;VARIABLE=T,gw;DEBUG_IO=0;");
  rval = mb.write_file("test_eul_T.nc", 0, writeopts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

// Check non-set variable T on some quads
// Also check set variable gw
void test_eul_check_T()
{
  Core moab;
  Interface& mb = moab;

  std::string opts;
  get_eul_options(opts);

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Load non-set variable T, set variable gw, and the mesh
  opts += std::string(";VARIABLE=T,gw");
  rval = mb.load_file("test_eul_T.nc", &set, opts.c_str());
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

    // Check some values of tag T0 on first level
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

// We also read and write set variables lat and lon, which are are required to create the mesh
// In test_homme_check_T(), we need to load the output file with mesh
void test_homme_read_write_T()
{
  Core moab;
  Interface& mb = moab;

  std::string orig;
  get_homme_options(orig);

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Load non-set variable T, set variable lat, set variable lon, and the mesh
  std::string opts = orig + std::string(";DEBUG_IO=0;VARIABLE=T,lat,lon");
  rval = mb.load_file(example_homme, &set, opts.c_str());
  CHECK_ERR(rval);

  // Write variables T, lat and lon
  std::string writeopts;
  writeopts = std::string(";;VARIABLE=T,lat,lon;DEBUG_IO=0;");
  rval = mb.write_file("test_homme_T.nc", 0, writeopts.c_str(), &set, 1);
  CHECK_ERR(rval);
}

// Check non-set variable T on some vertices
// Also check set variables lat and lon
void test_homme_check_T()
{
  Core moab;
  Interface& mb = moab;

  std::string opts;
  get_homme_options(opts);

  EntityHandle set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, set);
  CHECK_ERR(rval);

  // Load non-set variable T, set variable lat, set variable lon, and the mesh
  opts += std::string(";VARIABLE=T,lat,lon");
  opts += ";CONN=";
  opts += example_homme_mapping;
  rval = mb.load_file("test_homme_T.nc", &set, opts.c_str());
  CHECK_ERR(rval);

  int procs = 1;
#ifdef USE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  procs = pcomm->proc_config().proc_size();
#endif

  // Only test serial case for the time being
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
    double eps = 1e-10;
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

    // Get tag T0
    Tag Ttag0;
    rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0);
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

    // Check some values of tag T0 on first level
    eps = 0.0001;
    double* data = (double*) Tbuf;
    CHECK_REAL_EQUAL(233.1136, data[0 * 26], eps); // First vert
    CHECK_REAL_EQUAL(236.1505, data[1728 * 26], eps); // Median vert
    CHECK_REAL_EQUAL(235.7722, data[1729 * 26], eps); // Median vert
    CHECK_REAL_EQUAL(234.0416, data[3457 * 26], eps); // Last vert
  }
}

void get_eul_options(std::string& opts)
{
#ifdef USE_MPI
  // Use parallel options
  opts = std::string(";;PARALLEL=READ_PART;PARTITION_METHOD=SQIJ");
#else
  opts = std::string(";;");
#endif
}

void get_homme_options(std::string& opts)
{
#ifdef USE_MPI
  // Use parallel options
  opts = std::string(";;PARALLEL=READ_PART;PARTITION_METHOD=TRIVIAL");
#else
  opts = std::string(";;");
#endif
}
