#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/mpasx1.642.t.2.nc";
#else
static const char example[] = "/io/mpasx1.642.t.2.nc";
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
void test_read_no_mixed_elements();

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
  // Test read option NO_MIXED_ELEMENTS
  result += RUN_TEST(test_read_no_mixed_elements);

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

  int procs = 1;
#ifdef USE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  procs = pcomm->proc_config().proc_size();
#endif

  // Make check runs this test in one processor
  if (1 == procs) {
    // Check tags for vertex variable vorticity
    Tag vorticity_tag0, vorticity_tag1;
    rval = mb.tag_get_handle("vorticity0", 1, MB_TYPE_DOUBLE, vorticity_tag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("vorticity1", 1, MB_TYPE_DOUBLE, vorticity_tag1);
    CHECK_ERR(rval);

    // Get vertices (1280 edges)
    Range verts;
    rval = mb.get_entities_by_type(0, MBVERTEX, verts);
    assert(rval == MB_SUCCESS);
    CHECK_EQUAL((size_t)1280, verts.size());
    CHECK_EQUAL((size_t)1, verts.psize());

    const double eps = 1e-20;
    double val[2];

    // Check vorticity tag values on first two vertices
    EntityHandle vert_ents[] = {verts[0], verts[1]};
    rval = mb.tag_get_data(vorticity_tag0, &vert_ents[0], 2, val);
    CHECK_REAL_EQUAL(1.1, val[0], eps);
    CHECK_REAL_EQUAL(1.2, val[1], eps);
    rval = mb.tag_get_data(vorticity_tag1, &vert_ents[0], 2, val);
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
    assert(rval == MB_SUCCESS);
    CHECK_EQUAL((size_t)1920, edges.size());
    CHECK_EQUAL((size_t)1, edges.psize());

    // Check u tag values on two specified edges
    EntityHandle edge_ents[] = {edges[5], edges[6]};
    rval = mb.tag_get_data(u_tag0, &edge_ents[0], 2, val);
    CHECK_REAL_EQUAL(1.113138721544778, val[0], eps);
    CHECK_REAL_EQUAL(-1.113138721930009, val[1], eps);
    rval = mb.tag_get_data(u_tag1, &edge_ents[0], 2, val);
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
    assert(rval == MB_SUCCESS);
    CHECK_EQUAL((size_t)642, cells.size());
#ifdef USE_MPI
    // If MOAB is compiled parallel, sequence size requested are increased
    // by a factor of 1.5, to allow for ghosts. This will introduce a gap
    // between the two face sequences.
    CHECK_EQUAL((size_t)2, cells.psize());
#else
    CHECK_EQUAL((size_t)1, cells.psize());
#endif

    // Check ke tag values on first pentagon and first hexagon
    EntityHandle cell_ents[] = {cells[0], cells[12]};
    rval = mb.tag_get_data(ke_tag0, &cell_ents[0], 2, val);
    CHECK_REAL_EQUAL(15.001, val[0], eps);
    CHECK_REAL_EQUAL(16.013, val[1], eps);
    rval = mb.tag_get_data(ke_tag1, &cell_ents[0], 2, val);
    CHECK_REAL_EQUAL(25.001, val[0], eps);
    CHECK_REAL_EQUAL(26.013, val[1], eps);
  }
}

void test_read_onevar() 
{
  Core moab;
  Interface& mb = moab;
  std::string opts;
  ErrorCode rval = get_options(opts);
  CHECK_ERR(rval);

  opts += std::string(";VARIABLE=ke");
  // Create gather set
  opts += std::string(";GATHER_SET=");
  rval = mb.load_file(example, NULL, opts.c_str());
  CHECK_ERR(rval);

  int procs = 1;
#ifdef USE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  procs = pcomm->proc_config().proc_size();
#endif

  // Make check runs this test in one processor
  if (1 == procs) {
    // Check for proper tags
    Tag ke_tag0, ke_tag1;
    rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("ke1", 1, MB_TYPE_DOUBLE, ke_tag1);
    CHECK_ERR(rval);

    // Get cells (12 pentagons and 630 hexagons)
    Range cells;
    rval = mb.get_entities_by_type(0, MBPOLYGON, cells);
    assert(rval == MB_SUCCESS);
    CHECK_EQUAL((size_t)1284, cells.size()); // Gather set cells included
#ifdef USE_MPI
    // If MOAB is compiled parallel, sequence size requested are increased
    // by a factor of 1.5, to allow for ghosts. This will introduce a gap
    // between the two face sequences.
    CHECK_EQUAL((size_t)4, cells.psize()); // Gather set cells included
#else
    CHECK_EQUAL((size_t)1, cells.psize()); // Gather set cells included
#endif

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

    // Remove gather set cells
    cells = subtract(cells, gather_ents);
    CHECK_EQUAL((size_t)642, cells.size()); // Gather set cells excluded
#ifdef USE_MPI
    // If MOAB is compiled parallel, sequence size requested are increased
    // by a factor of 1.5, to allow for ghosts. This will introduce a gap
    // between the two face sequences.
    CHECK_EQUAL((size_t)2, cells.psize()); // Gather set cells excluded
#else
    CHECK_EQUAL((size_t)1, cells.psize()); // Gather set cells excluded
#endif

    // Check ke tag values on first pentagon and first hexagon
    const double eps = 1e-20;
    double val[2];
    EntityHandle cell_ents[] = {cells[0], cells[12]};
    rval = mb.tag_get_data(ke_tag0, &cell_ents[0], 2, val);
    CHECK_REAL_EQUAL(15.001, val[0], eps);
    CHECK_REAL_EQUAL(16.013, val[1], eps);
    rval = mb.tag_get_data(ke_tag1, &cell_ents[0], 2, val);
    CHECK_REAL_EQUAL(25.001, val[0], eps);
    CHECK_REAL_EQUAL(26.013, val[1], eps);
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
  Tag ke_tag0, ke_tag1;
  rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);
  rval = mb.tag_get_handle("ke1", 1, MB_TYPE_DOUBLE, ke_tag1);
  CHECK_ERR(rval);
}

void test_read_nomesh() 
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

  opts = orig + std::string(";TIMESTEP=0");
  rval = mb.load_file(example, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tags
  Tag ke_tag0, ke_tag1;
  rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
  CHECK_ERR(rval);
  rval = mb.tag_get_handle("ke1", 1, MB_TYPE_DOUBLE, ke_tag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Now read 2nd timestep with nomesh option
  opts = orig + std::string(";TIMESTEP=1;NOMESH");
  rval = mb.load_file(example, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  rval = mb.tag_get_handle("ke1", 1, MB_TYPE_DOUBLE, ke_tag1);
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
  Tag ke_tag0, ke_tag1;
  rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  opts = orig + std::string(";VARIABLE=ke;TIMESTEP=0;NOMESH");
  rval = mb.load_file(example, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tags
  rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
  CHECK_ERR(rval);
  rval = mb.tag_get_handle("ke1", 1, MB_TYPE_DOUBLE, ke_tag1);
  CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);

  // Now read 2nd timestep with nomesh option
  opts = orig + std::string(";VARIABLE=ke;TIMESTEP=1;NOMESH");
  rval = mb.load_file(example, &set, opts.c_str());
  CHECK_ERR(rval);

  // Check for proper tag
  rval = mb.tag_get_handle("ke1", 1, MB_TYPE_DOUBLE, ke_tag1);
  CHECK_ERR(rval);
}

void test_read_no_mixed_elements()
{
  Core moab;
  Interface& mb = moab;
  std::string opts;
  ErrorCode rval = get_options(opts);
  CHECK_ERR(rval);

  opts += std::string(";NO_MIXED_ELEMENTS;VARIABLE=ke");
  rval = mb.load_file(example, NULL, opts.c_str());
  CHECK_ERR(rval);

  int procs = 1;
#ifdef USE_MPI
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  procs = pcomm->proc_config().proc_size();
#endif

  // Make check runs this test in one processor
  if (1 == procs) {
    // Check for proper tags
    Tag ke_tag0, ke_tag1;
    rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle("ke1", 1, MB_TYPE_DOUBLE, ke_tag1);
    CHECK_ERR(rval);

    // Get cells (12 pentagons and 630 hexagons)
    Range cells;
    rval = mb.get_entities_by_type(0, MBPOLYGON, cells);
    assert(rval == MB_SUCCESS);
    CHECK_EQUAL((size_t)642, cells.size());
    // Only one group of cells (each cell is actually represented by a 10-vertex polygon)
    CHECK_EQUAL((size_t)1, cells.psize());

    const double eps = 1e-20;
    double val[2];

    // Check ke tag values on first pentagon and first hexagon
    EntityHandle cell_ents[] = {cells[0], cells[12]};
    rval = mb.tag_get_data(ke_tag0, &cell_ents[0], 2, val);
    CHECK_REAL_EQUAL(15.001, val[0], eps);
    CHECK_REAL_EQUAL(16.013, val[1], eps);
    rval = mb.tag_get_data(ke_tag1, &cell_ents[0], 2, val);
    CHECK_REAL_EQUAL(25.001, val[0], eps);
    CHECK_REAL_EQUAL(26.013, val[1], eps);
  }
}

ErrorCode get_options(std::string &opts) 
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
