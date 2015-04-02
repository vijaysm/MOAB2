#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/gcrm_r3.nc";
#else
static const char example[] = "/io/gcrm_r3.nc";
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
void test_read_no_edges(); // Test read option NO_EDGES
void test_gather_onevar(); // Test gather set with one variable

void get_options(std::string& opts);

const double eps = 1e-6;
const int layers = 3;
const int interfaces = 3;

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
  result += RUN_TEST(test_read_no_edges);
  result += RUN_TEST(test_gather_onevar);

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

  // Read mesh and read all variables at all timesteps
  ErrorCode rval = mb.load_file(example, 0, opts.c_str());
  CHECK_ERR(rval);

  int procs = 1;
#ifdef MOAB_HAVE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  procs = pcomm->proc_config().proc_size();
#endif

  // Make check runs this test on one processor
  if (1 == procs) {
    // For u, wind and vorticity, check tag values on two entities
    double val[2 * layers];

    // Check tags for vertex variable u
    Tag u_tag0, u_tag1;
    rval = mb.tag_get_handle("u0", layers, MB_TYPE_DOUBLE, u_tag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("u1", layers, MB_TYPE_DOUBLE, u_tag1);
    CHECK_ERR(rval);

    // Get vertices (1280 vertices)
    Range verts;
    rval = mb.get_entities_by_type(0, MBVERTEX, verts);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)1280, verts.size());
    CHECK_EQUAL((size_t)1, verts.psize());

    // Check u tag values on first and last vertices
    EntityHandle vert_ents[] = {verts[0], verts[1279]};

    // Only check first two layers
    // Timestep 0
    rval = mb.tag_get_data(u_tag0, vert_ents, 2, val);
    CHECK_ERR(rval);
    // Layer 0
    CHECK_REAL_EQUAL(-4.839992, val[0 * layers], eps);
    CHECK_REAL_EQUAL(-3.699257, val[1 * layers], eps);
    // Layer 1
    CHECK_REAL_EQUAL(-4.839925, val[0 * layers + 1], eps);
    CHECK_REAL_EQUAL(-3.699206, val[1 * layers + 1], eps);

    // Timestep 1
    rval = mb.tag_get_data(u_tag1, vert_ents, 2, val);
    CHECK_ERR(rval);
    // Layer 0
    CHECK_REAL_EQUAL(-4.712473, val[0 * layers], eps);
    CHECK_REAL_EQUAL(-3.601793, val[1 * layers], eps);
    // Layer 1
    CHECK_REAL_EQUAL(-4.712409, val[0 * layers + 1], eps);
    CHECK_REAL_EQUAL(-3.601743, val[1 * layers + 1], eps);

    // Check tags for edge variable wind
    Tag wind_tag0, wind_tag1;
    rval = mb.tag_get_handle("wind0", layers, MB_TYPE_DOUBLE, wind_tag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("wind1", layers, MB_TYPE_DOUBLE, wind_tag1);
    CHECK_ERR(rval);

    // Get edges (1920 edges)
    Range edges;
    rval = mb.get_entities_by_type(0, MBEDGE, edges);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)1920, edges.size());
    CHECK_EQUAL((size_t)1, edges.psize());

    // Check wind tag values on first and last edges
    EntityHandle edge_ents[] = {edges[0], edges[1919]};

    // Only check first two layers
    // Timestep 0
    rval = mb.tag_get_data(wind_tag0, edge_ents, 2, val);
    CHECK_ERR(rval);
    // Layer 0
    CHECK_REAL_EQUAL(-5.081991, val[0 * layers], eps);
    CHECK_REAL_EQUAL(-6.420274, val[1 * layers], eps);
    // Layer 1
    CHECK_REAL_EQUAL(-5.081781, val[0 * layers + 1], eps);
    CHECK_REAL_EQUAL(-6.419831, val[1 * layers + 1], eps);

    // Timestep 1
    rval = mb.tag_get_data(wind_tag1, edge_ents, 2, val);
    CHECK_ERR(rval);
    // Layer 0
    CHECK_REAL_EQUAL(-4.948097, val[0 * layers], eps);
    CHECK_REAL_EQUAL(-6.251121, val[1 * layers], eps);
    // Layer 1
    CHECK_REAL_EQUAL(-4.947892, val[0 * layers + 1], eps);
    CHECK_REAL_EQUAL(-6.250690, val[1 * layers + 1], eps);

    // Check tags for cell variable vorticity
    Tag vorticity_tag0, vorticity_tag1;
    rval = mb.tag_get_handle("vorticity0", layers, MB_TYPE_DOUBLE, vorticity_tag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("vorticity1", layers, MB_TYPE_DOUBLE, vorticity_tag1);
    CHECK_ERR(rval);

    // Get cells (12 pentagons and 630 hexagons)
    Range cells;
    rval = mb.get_entities_by_type(0, MBPOLYGON, cells);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)642, cells.size());

    // GCRM pentagons are always padded to hexagons
    CHECK_EQUAL((size_t)1, cells.psize());

    // Check vorticity tag values on first and last cells
    EntityHandle cell_ents[] = {cells[0], cells[641]};

    // Only check first two layers
    // Timestep 0
    rval = mb.tag_get_data(vorticity_tag0, cell_ents, 2, val);
    CHECK_ERR(rval);
    // Layer 0
    CHECK_REAL_EQUAL(3.629994, val[0 * layers], eps);
    CHECK_REAL_EQUAL(-0.554888, val[1 * layers], eps);
    // Layer 1
    CHECK_REAL_EQUAL(3.629944, val[0 * layers + 1], eps);
    CHECK_REAL_EQUAL(-0.554881, val[1 * layers + 1], eps);

    // Timestep 1
    rval = mb.tag_get_data(vorticity_tag1, cell_ents, 2, val);
    CHECK_ERR(rval);
    // Layer 0
    CHECK_REAL_EQUAL(3.534355, val[0 * layers], eps);
    CHECK_REAL_EQUAL(-0.540269, val[1 * layers], eps);
    // Layer 1
    CHECK_REAL_EQUAL(3.534306, val[0 * layers + 1], eps);
    CHECK_REAL_EQUAL(-0.540262, val[1 * layers + 1], eps);

    // Check tags for cell variable pressure
    Tag pressure_tag0, pressure_tag1;
    rval = mb.tag_get_handle("pressure0", interfaces, MB_TYPE_DOUBLE, pressure_tag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("pressure1", interfaces, MB_TYPE_DOUBLE, pressure_tag1);
    CHECK_ERR(rval);

    // For pressure, check tag values on two cells
    double pressure_val[2 * interfaces];

    // Check pressure tag values on first and last cells
    // Only check first two interfaces
    // Timestep 0
    rval = mb.tag_get_data(pressure_tag0, cell_ents, 2, pressure_val);
    CHECK_ERR(rval);
    // Interface 0
    CHECK_REAL_EQUAL(4.44234e-06, pressure_val[0 * interfaces], 1e-11);
    CHECK_REAL_EQUAL(0.2486804, pressure_val[1 * interfaces], 1e-7);
    // Interface 1
    CHECK_REAL_EQUAL(4.44234e-06, pressure_val[0 * interfaces + 1], 1e-11);
    CHECK_REAL_EQUAL(0.2486804, pressure_val[1 * interfaces + 1], 1e-7);

    // Timestep 1
    rval = mb.tag_get_data(pressure_tag1, cell_ents, 2, pressure_val);
    CHECK_ERR(rval);
    // Interface 0
    CHECK_REAL_EQUAL(2.365176e-07, pressure_val[0 * interfaces], 1e-13);
    CHECK_REAL_EQUAL(0.02234409, pressure_val[1 * interfaces], 1e-8);
    // Interface 1
    CHECK_REAL_EQUAL(2.365176e-07, pressure_val[0 * interfaces + 1], 1e-13);
    CHECK_REAL_EQUAL(0.02234409, pressure_val[1 * interfaces + 1], 1e-8);
  }
}

void test_read_onevar()
{
  Core moab;
  Interface& mb = moab;

  std::string opts;
  get_options(opts);

  // Read mesh and read cell variable vorticity at all timesteps
  opts += ";VARIABLE=vorticity";
  ErrorCode rval = mb.load_file(example, NULL, opts.c_str());
  CHECK_ERR(rval);

  int procs = 1;
#ifdef MOAB_HAVE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  procs = pcomm->proc_config().proc_size();
#endif

  // Make check runs this test on one processor
  if (1 == procs) {
    // Check vorticity tags
    Tag vorticity_tag0, vorticity_tag1;
    rval = mb.tag_get_handle("vorticity0", layers, MB_TYPE_DOUBLE, vorticity_tag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("vorticity1", layers, MB_TYPE_DOUBLE, vorticity_tag1);
    CHECK_ERR(rval);

    // Get cells (12 pentagons and 630 hexagons)
    Range cells;
    rval = mb.get_entities_by_type(0, MBPOLYGON, cells);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)642, cells.size());

    // GCRM pentagons are always padded to hexagons
    CHECK_EQUAL((size_t)1, cells.psize());

    // Check vorticity tag values on 4 cells: first cell, two median cells, and last cell
    EntityHandle cell_ents[] = {cells[0], cells[320], cells[321], cells[641]};
    double vorticity_val[4 * layers];

    // Only check first two layers
    // Timestep 0
    rval = mb.tag_get_data(vorticity_tag0, cell_ents, 4, vorticity_val);
    CHECK_ERR(rval);
    // Layer 0
    CHECK_REAL_EQUAL(3.629994, vorticity_val[0 * layers], eps);
    CHECK_REAL_EQUAL(0.131688, vorticity_val[1 * layers], eps);
    CHECK_REAL_EQUAL(-0.554888, vorticity_val[2 * layers], eps);
    CHECK_REAL_EQUAL(-0.554888, vorticity_val[3 * layers], eps);
    // Layer 1
    CHECK_REAL_EQUAL(3.629944, vorticity_val[0 * layers + 1], eps);
    CHECK_REAL_EQUAL(0.131686, vorticity_val[1 * layers + 1], eps);
    CHECK_REAL_EQUAL(-0.554881, vorticity_val[2 * layers + 1], eps);
    CHECK_REAL_EQUAL(-0.554881, vorticity_val[3 * layers + 1], eps);

    // Timestep 1
    rval = mb.tag_get_data(vorticity_tag1, cell_ents, 4, vorticity_val);
    CHECK_ERR(rval);
    // Layer 0
    CHECK_REAL_EQUAL(3.534355, vorticity_val[0 * layers], eps);
    CHECK_REAL_EQUAL(0.128218, vorticity_val[1 * layers], eps);
    CHECK_REAL_EQUAL(-0.540269, vorticity_val[2 * layers], eps);
    CHECK_REAL_EQUAL(-0.540269, vorticity_val[3 * layers], eps);
    // Layer 1
    CHECK_REAL_EQUAL(3.534306, vorticity_val[0 * layers + 1], eps);
    CHECK_REAL_EQUAL(0.128216, vorticity_val[1 * layers + 1], eps);
    CHECK_REAL_EQUAL(-0.540262, vorticity_val[2 * layers + 1], eps);
    CHECK_REAL_EQUAL(-0.540262, vorticity_val[3 * layers + 1], eps);
  }
}

void test_read_onetimestep()
{
  Core moab;
  Interface& mb = moab;

  std::string opts;
  get_options(opts);

  // Read mesh and read all variables at 2nd timestep
  opts += ";TIMESTEP=1";
  ErrorCode rval = mb.load_file(example, NULL, opts.c_str());
  CHECK_ERR(rval);

  // Check vorticity tags
  Tag vorticity_tag0, vorticity_tag1;
  rval = mb.tag_get_handle("vorticity0", layers, MB_TYPE_DOUBLE, vorticity_tag0);
  // Tag vorticity0 should not exist
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);
  rval = mb.tag_get_handle("vorticity1", layers, MB_TYPE_DOUBLE, vorticity_tag1);
  CHECK_ERR(rval);
}

void test_read_nomesh()
{
  Core moab;
  Interface& mb = moab;

  // Need a file set for nomesh to work right
  EntityHandle file_set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, file_set);
  CHECK_ERR(rval);

  std::string orig, opts;
  get_options(orig);

  // Read mesh and read all variables at 1st timestep
  opts = orig + ";TIMESTEP=0";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Check u tags
  Tag u_tag0, u_tag1;
  rval = mb.tag_get_handle("u0", layers, MB_TYPE_DOUBLE, u_tag0);
  CHECK_ERR(rval);
  // Tag u1 should not exist
  rval = mb.tag_get_handle("u1", layers, MB_TYPE_DOUBLE, u_tag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Read all variables at 2nd timestep 0, no need to read mesh
  opts = orig + ";TIMESTEP=1;NOMESH";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Check tag u1 again
  rval = mb.tag_get_handle("u1", layers, MB_TYPE_DOUBLE, u_tag1);
  // Tag u1 should exist at this time
  CHECK_ERR(rval);
}

void test_read_novars()
{
  Core moab;
  Interface& mb = moab;

  // Need a file set for nomesh to work right
  EntityHandle file_set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, file_set);
  CHECK_ERR(rval);

  std::string orig, opts;
  get_options(orig);
  CHECK_ERR(rval);

  // Read header info only, no mesh, no variables
  opts = orig + ";NOMESH;VARIABLE=";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Read mesh, but still no variables
  opts = orig + ";VARIABLE=";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Check vorticity tags
  Tag vorticity_tag0, vorticity_tag1;
  rval = mb.tag_get_handle("vorticity0", layers, MB_TYPE_DOUBLE, vorticity_tag0);
  // Tag vorticity0 should not exist
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);
  rval = mb.tag_get_handle("vorticity1", layers, MB_TYPE_DOUBLE, vorticity_tag1);
  // Tag vorticity1 should not exist
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Read vorticity at 1st timestep, no need to read mesh
  opts = orig + ";VARIABLE=vorticity;TIMESTEP=0;NOMESH";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Check vorticity tags again
  rval = mb.tag_get_handle("vorticity0", layers, MB_TYPE_DOUBLE, vorticity_tag0);
  // Tag vorticity0 should exist at this time
  CHECK_ERR(rval);
  // Tag vorticity1 should still not exist
  rval = mb.tag_get_handle("vorticity1", layers, MB_TYPE_DOUBLE, vorticity_tag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Read vorticity at 2nd timestep, no need to read mesh
  opts = orig + ";VARIABLE=vorticity;TIMESTEP=1;NOMESH";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Check tag vorticity1 again
  rval = mb.tag_get_handle("vorticity1", layers, MB_TYPE_DOUBLE, vorticity_tag1);
  // Tag vorticity1 should exist at this time
  CHECK_ERR(rval);

  int procs = 1;
#ifdef MOAB_HAVE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  procs = pcomm->proc_config().proc_size();
#endif

  // Make check runs this test on one processor
  if (1 == procs) {
    // Get cells (12 pentagons and 630 hexagons)
    Range cells;
    rval = mb.get_entities_by_type(file_set, MBPOLYGON, cells);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)642, cells.size());

    // GCRM pentagons are always padded to hexagons
    CHECK_EQUAL((size_t)1, cells.psize());

    // Check vorticity tag values on 4 cells: first cell, two median cells, and last cell
    EntityHandle cell_ents[] = {cells[0], cells[320], cells[321], cells[641]};
    double vorticity_val[4 * layers];

    // Only check first two layers
    // Timestep 0
    rval = mb.tag_get_data(vorticity_tag0, cell_ents, 4, vorticity_val);
    CHECK_ERR(rval);
    // Layer 0
    CHECK_REAL_EQUAL(3.629994, vorticity_val[0 * layers], eps);
    CHECK_REAL_EQUAL(0.131688, vorticity_val[1 * layers], eps);
    CHECK_REAL_EQUAL(-0.554888, vorticity_val[2 * layers], eps);
    CHECK_REAL_EQUAL(-0.554888, vorticity_val[3 * layers], eps);
    // Layer 1
    CHECK_REAL_EQUAL(3.629944, vorticity_val[0 * layers + 1], eps);
    CHECK_REAL_EQUAL(0.131686, vorticity_val[1 * layers + 1], eps);
    CHECK_REAL_EQUAL(-0.554881, vorticity_val[2 * layers + 1], eps);
    CHECK_REAL_EQUAL(-0.554881, vorticity_val[3 * layers + 1], eps);

    // Timestep 1
    rval = mb.tag_get_data(vorticity_tag1, cell_ents, 4, vorticity_val);
    CHECK_ERR(rval);
    // Layer 0
    CHECK_REAL_EQUAL(3.534355, vorticity_val[0 * layers], eps);
    CHECK_REAL_EQUAL(0.128218, vorticity_val[1 * layers], eps);
    CHECK_REAL_EQUAL(-0.540269, vorticity_val[2 * layers], eps);
    CHECK_REAL_EQUAL(-0.540269, vorticity_val[3 * layers], eps);
    // Layer 1
    CHECK_REAL_EQUAL(3.534306, vorticity_val[0 * layers + 1], eps);
    CHECK_REAL_EQUAL(0.128216, vorticity_val[1 * layers + 1], eps);
    CHECK_REAL_EQUAL(-0.540262, vorticity_val[2 * layers + 1], eps);
    CHECK_REAL_EQUAL(-0.540262, vorticity_val[3 * layers + 1], eps);
  }
}

void test_read_no_edges()
{
  Core moab;
  Interface& mb = moab;

  std::string opts;
  get_options(opts);

  opts += ";NO_EDGES;VARIABLE=";
  ErrorCode rval = mb.load_file(example, NULL, opts.c_str());
  CHECK_ERR(rval);

  // Get edges
  Range edges;
  rval = mb.get_entities_by_type(0, MBEDGE, edges);
  CHECK_ERR(rval);
  CHECK_EQUAL((size_t)0, edges.size());
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

  // Read cell variable vorticity and create gather set on processor 0
  opts += ";VARIABLE=vorticity;GATHER_SET=0";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

#ifdef MOAB_HAVE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  int rank = pcomm->proc_config().proc_rank();

  Range cells, cells_owned;
  rval = mb.get_entities_by_type(file_set, MBPOLYGON, cells);
  CHECK_ERR(rval);

  // Get local owned cells
  rval = pcomm->filter_pstatus(cells, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &cells_owned);
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

  Tag vorticity_tag0, gid_tag;
  rval = mb.tag_get_handle("vorticity0", layers, MB_TYPE_DOUBLE, vorticity_tag0, MB_TAG_DENSE);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid_tag, MB_TAG_DENSE);
  CHECK_ERR(rval);

  pcomm->gather_data(cells_owned, vorticity_tag0, gid_tag, gather_set, 0);

  if (0 == rank) {
    // Get gather set cells
    Range gather_set_cells;
    rval = mb.get_entities_by_type(gather_set, MBPOLYGON, gather_set_cells);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)642, gather_set_cells.size());
    CHECK_EQUAL((size_t)1, gather_set_cells.psize());

    // Check vorticity0 tag values on 4 gather set cells: first cell, two median cells, and last cell
    EntityHandle cell_ents[] = {gather_set_cells[0], gather_set_cells[320],
                                gather_set_cells[321], gather_set_cells[641]};
    double vorticity0_val[4 * layers];
    rval = mb.tag_get_data(vorticity_tag0, cell_ents, 4, vorticity0_val);
    CHECK_ERR(rval);

    // Only check first two layers
    // Layer 0
    CHECK_REAL_EQUAL(3.629994, vorticity0_val[0 * layers], eps);
    CHECK_REAL_EQUAL(0.131688, vorticity0_val[1 * layers], eps);
    CHECK_REAL_EQUAL(-0.554888, vorticity0_val[2 * layers], eps);
    CHECK_REAL_EQUAL(-0.554888, vorticity0_val[3 * layers], eps);
    // Layer 1
    CHECK_REAL_EQUAL(3.629944, vorticity0_val[0 * layers + 1], eps);
    CHECK_REAL_EQUAL(0.131686, vorticity0_val[1 * layers + 1], eps);
    CHECK_REAL_EQUAL(-0.554881, vorticity0_val[2 * layers + 1], eps);
    CHECK_REAL_EQUAL(-0.554881, vorticity0_val[3 * layers + 1], eps);
  }
#endif
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
