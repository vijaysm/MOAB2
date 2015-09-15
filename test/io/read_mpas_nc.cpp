#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/mpasx1.642.t.2.nc";
#else
static const char example[] = "/io/mpasx1.642.t.2.nc";
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
void test_read_no_mixed_elements(); // Test read option NO_MIXED_ELEMENTS
void test_read_no_edges(); // Test read option NO_EDGES
void test_gather_onevar(); // Test gather set with one variable

void get_options(std::string& opts);

const double eps = 1e-20;

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
  result += RUN_TEST(test_read_no_mixed_elements);
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
  ErrorCode rval = mb.load_file(example, NULL, opts.c_str());
  CHECK_ERR(rval);

  int procs = 1;
#ifdef MOAB_HAVE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  procs = pcomm->proc_config().proc_size();
#endif

  // Make check runs this test on one processor
  if (1 == procs) {
    // For each tag, check two values
    double val[2];

    // Check tags for vertex variable vorticity
    Tag vorticity_tag0, vorticity_tag1;
    rval = mb.tag_get_handle("vorticity0", 1, MB_TYPE_DOUBLE, vorticity_tag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("vorticity1", 1, MB_TYPE_DOUBLE, vorticity_tag1);
    CHECK_ERR(rval);

    // Get vertices (1280 edges)
    Range verts;
    rval = mb.get_entities_by_type(0, MBVERTEX, verts);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)1280, verts.size());
    CHECK_EQUAL((size_t)1, verts.psize());

    // Check vorticity tag values on first two vertices
    EntityHandle vert_ents[] = {verts[0], verts[1]};
    rval = mb.tag_get_data(vorticity_tag0, vert_ents, 2, val);
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL(1.1, val[0], eps);
    CHECK_REAL_EQUAL(1.2, val[1], eps);
    rval = mb.tag_get_data(vorticity_tag1, vert_ents, 2, val);
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL(2.1, val[0], eps);
    CHECK_REAL_EQUAL(2.2, val[1], eps);

    // Check tags for edge variable u
    Tag u_tag0, u_tag1;
    rval = mb.tag_get_handle("u0", 1, MB_TYPE_DOUBLE, u_tag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("u1", 1, MB_TYPE_DOUBLE, u_tag1);
    CHECK_ERR(rval);

    // Get edges (1920 edges)
    Range edges;
    rval = mb.get_entities_by_type(0, MBEDGE, edges);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)1920, edges.size());
    CHECK_EQUAL((size_t)1, edges.psize());

    // Check u tag values on two specified edges
    EntityHandle edge_ents[] = {edges[5], edges[6]};
    rval = mb.tag_get_data(u_tag0, edge_ents, 2, val);
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL(1.113138721544778, val[0], eps);
    CHECK_REAL_EQUAL(-1.113138721930009, val[1], eps);
    rval = mb.tag_get_data(u_tag1, edge_ents, 2, val);
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL(2.113138721544778, val[0], eps);
    CHECK_REAL_EQUAL(-2.113138721930009, val[1], eps);

    // Check tags for cell variable ke
    Tag ke_tag0, ke_tag1;
    rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("ke1", 1, MB_TYPE_DOUBLE, ke_tag1);
    CHECK_ERR(rval);

    // Get cells (12 pentagons and 630 hexagons)
    Range cells;
    rval = mb.get_entities_by_type(0, MBPOLYGON, cells);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)642, cells.size());
    CHECK_EQUAL((size_t)1, cells.psize());

    // Check ke tag values on first pentagon and first hexagon
    EntityHandle cell_ents[] = {cells[0], cells[12]};
    rval = mb.tag_get_data(ke_tag0, cell_ents, 2, val);
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL(15.001, val[0], eps);
    CHECK_REAL_EQUAL(16.013, val[1], eps);
    rval = mb.tag_get_data(ke_tag1, cell_ents, 2, val);
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL(25.001, val[0], eps);
    CHECK_REAL_EQUAL(26.013, val[1], eps);
  }
}

void test_read_onevar() 
{
  Core moab;
  Interface& mb = moab;

  std::string opts;
  get_options(opts);

  // Read mesh and read cell variable ke at all timesteps
  opts += ";VARIABLE=ke";
  ErrorCode rval = mb.load_file(example, NULL, opts.c_str());
  CHECK_ERR(rval);

  int procs = 1;
#ifdef MOAB_HAVE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  procs = pcomm->proc_config().proc_size();
#endif

  // Make check runs this test on one processor
  if (1 == procs) {
    // Check ke tags
    Tag ke_tag0, ke_tag1;
    rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("ke1", 1, MB_TYPE_DOUBLE, ke_tag1);
    CHECK_ERR(rval);

    // Get cells (12 pentagons and 630 hexagons)
    Range cells;
    rval = mb.get_entities_by_type(0, MBPOLYGON, cells);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)642, cells.size());
    CHECK_EQUAL((size_t)1, cells.psize());

    // Check ke tag values on 4 cells: first pentagon, last pentagon,
    // first hexagon, and last hexagon
    EntityHandle cell_ents[] = {cells[0], cells[11], cells[12], cells[641]};
    double ke0_val[4];
    rval = mb.tag_get_data(ke_tag0, cell_ents, 4, ke0_val);
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL(15.001, ke0_val[0], eps);
    CHECK_REAL_EQUAL(15.012, ke0_val[1], eps);
    CHECK_REAL_EQUAL(16.013, ke0_val[2], eps);
    CHECK_REAL_EQUAL(16.642, ke0_val[3], eps);
    double ke1_val[4];
    rval = mb.tag_get_data(ke_tag1, cell_ents, 4, ke1_val);
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL(25.001, ke1_val[0], eps);
    CHECK_REAL_EQUAL(25.012, ke1_val[1], eps);
    CHECK_REAL_EQUAL(26.013, ke1_val[2], eps);
    CHECK_REAL_EQUAL(26.642, ke1_val[3], eps);
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
  rval = mb.tag_get_handle("vorticity0", 1, MB_TYPE_DOUBLE, vorticity_tag0);
  // Tag vorticity0 should not exist
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);
  rval = mb.tag_get_handle("vorticity1", 1, MB_TYPE_DOUBLE, vorticity_tag1);
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
  rval = mb.tag_get_handle("u0", 1, MB_TYPE_DOUBLE, u_tag0);
  CHECK_ERR(rval);
  // Tag u1 should not exist
  rval = mb.tag_get_handle("u1", 1, MB_TYPE_DOUBLE, u_tag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Read all variables at 2nd timestep 0, no need to read mesh
  opts = orig + ";TIMESTEP=1;NOMESH";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Check tag u1 again
  rval = mb.tag_get_handle("u1", 1, MB_TYPE_DOUBLE, u_tag1);
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
  opts = orig + ";VARIABLE=;TIMESTEP=0";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Check ke tags
  Tag ke_tag0, ke_tag1;
  rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
  // Tag ke0 should not exist
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);
  rval = mb.tag_get_handle("ke1", 1, MB_TYPE_DOUBLE, ke_tag1);
  // Tag ke1 should not exist
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Read ke at 1st timestep, no need to read mesh
  opts = orig + ";VARIABLE=ke;TIMESTEP=0;NOMESH";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Check ke tags again
  rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
  // Tag ke0 should exist at this time
  CHECK_ERR(rval);
  // Tag ke1 should still not exist
  rval = mb.tag_get_handle("ke1", 1, MB_TYPE_DOUBLE, ke_tag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Read ke at 2nd timestep, no need to read mesh
  opts = orig + ";VARIABLE=ke;TIMESTEP=1;NOMESH";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Check tag ke1 again
  rval = mb.tag_get_handle("ke1", 1, MB_TYPE_DOUBLE, ke_tag1);
  // Tag ke1 should exist at this time
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
    CHECK_EQUAL((size_t)1, cells.psize());

    // Check ke tag values on 4 cells: first pentagon, last pentagon,
    // first hexagon, and last hexagon
    EntityHandle cell_ents[] = {cells[0], cells[11], cells[12], cells[641]};
    double ke0_val[4];
    rval = mb.tag_get_data(ke_tag0, cell_ents, 4, ke0_val);
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL(15.001, ke0_val[0], eps);
    CHECK_REAL_EQUAL(15.012, ke0_val[1], eps);
    CHECK_REAL_EQUAL(16.013, ke0_val[2], eps);
    CHECK_REAL_EQUAL(16.642, ke0_val[3], eps);
    double ke1_val[4];
    rval = mb.tag_get_data(ke_tag1, cell_ents, 4, ke1_val);
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL(25.001, ke1_val[0], eps);
    CHECK_REAL_EQUAL(25.012, ke1_val[1], eps);
    CHECK_REAL_EQUAL(26.013, ke1_val[2], eps);
    CHECK_REAL_EQUAL(26.642, ke1_val[3], eps);
  }
}

void test_read_no_mixed_elements()
{
  Core moab;
  Interface& mb = moab;

  std::string opts;
  get_options(opts);

  // Read mesh with no mixed elements and read all variables at all timesteps
  opts += ";NO_MIXED_ELEMENTS";
  ErrorCode rval = mb.load_file(example, NULL, opts.c_str());
  CHECK_ERR(rval);

  int procs = 1;
#ifdef MOAB_HAVE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  procs = pcomm->proc_config().proc_size();
#endif

  // Make check runs this test on one processor
  if (1 == procs) {
    // Check ke tags
    Tag ke_tag0, ke_tag1;
    rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("ke1", 1, MB_TYPE_DOUBLE, ke_tag1);
    CHECK_ERR(rval);

    // Get cells (12 pentagons and 630 hexagons)
    Range cells;
    rval = mb.get_entities_by_type(0, MBPOLYGON, cells);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)642, cells.size());
    // Only one group of cells (pentagons are padded to hexagons,
    // e.g. connectivity [1 2 3 4 5] => [1 2 3 4 5 5])
    CHECK_EQUAL((size_t)1, cells.psize());

    // Check ke tag values on 4 cells: first pentagon, last pentagon,
    // first hexagon, and last hexagon
    EntityHandle cell_ents[] = {cells[0], cells[11], cells[12], cells[641]};
    double ke0_val[4];
    rval = mb.tag_get_data(ke_tag0, cell_ents, 4, ke0_val);
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL(15.001, ke0_val[0], eps);
    CHECK_REAL_EQUAL(15.012, ke0_val[1], eps);
    CHECK_REAL_EQUAL(16.013, ke0_val[2], eps);
    CHECK_REAL_EQUAL(16.642, ke0_val[3], eps);
    double ke1_val[4];
    rval = mb.tag_get_data(ke_tag1, cell_ents, 4, ke1_val);
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL(25.001, ke1_val[0], eps);
    CHECK_REAL_EQUAL(25.012, ke1_val[1], eps);
    CHECK_REAL_EQUAL(26.013, ke1_val[2], eps);
    CHECK_REAL_EQUAL(26.642, ke1_val[3], eps);
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

  // Read cell variable ke and create gather set on processor 0
  opts += ";VARIABLE=ke;GATHER_SET=0";
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

  Tag ke_tag0, gid_tag;
  rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0, MB_TAG_DENSE);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid_tag, MB_TAG_DENSE);
  CHECK_ERR(rval);

  pcomm->gather_data(cells_owned, ke_tag0, gid_tag, gather_set, 0);

  if (0 == rank) {
    // Get gather set cells
    Range gather_set_cells;
    rval = mb.get_entities_by_type(gather_set, MBPOLYGON, gather_set_cells);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)642, gather_set_cells.size());
    CHECK_EQUAL((size_t)1, gather_set_cells.psize());

    // Check ke0 tag values on 4 gather set cells: first pentagon, last pentagon,
    // first hexagon, and last hexagon
    double ke0_val[4];
    EntityHandle cell_ents[] = {gather_set_cells[0], gather_set_cells[11],
                                gather_set_cells[12], gather_set_cells[641]};
    rval = mb.tag_get_data(ke_tag0, cell_ents, 4, ke0_val);
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL(15.001, ke0_val[0], eps);
    CHECK_REAL_EQUAL(15.012, ke0_val[1], eps);
    CHECK_REAL_EQUAL(16.013, ke0_val[2], eps);
    CHECK_REAL_EQUAL(16.642, ke0_val[3], eps);
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
