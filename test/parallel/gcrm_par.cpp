#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/ProgOptions.hpp"
#include "MBParallelConventions.h"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"

#include <sstream>

using namespace moab;

#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/gcrm_r3.nc";
#endif

void test_read_onevar_trivial();
#if defined(MOAB_HAVE_MPI) && defined(MOAB_HAVE_ZOLTAN)
void test_read_onevar_rcbzoltan();
#endif

void test_read_mesh_parallel_trivial();
#if defined(MOAB_HAVE_MPI) && defined(MOAB_HAVE_ZOLTAN)
void test_read_mesh_parallel_rcbzoltan();
#endif

void test_gather_onevar_on_rank0();
void test_gather_onevar_on_rank1();

void test_multiple_loads_of_same_file();

// Helper functions
void read_one_cell_var(bool rcbzoltan);
void read_mesh_parallel(bool rcbzoltan);
void gather_one_cell_var(int gather_set_rank);
void multiple_loads_of_same_file();

std::string read_options;
const double eps = 1e-6;
const int layers = 3;

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int result = 0;

  result += RUN_TEST(test_read_onevar_trivial);
#if defined(MOAB_HAVE_MPI) && defined(MOAB_HAVE_ZOLTAN)
  result += RUN_TEST(test_read_onevar_rcbzoltan);
#endif

  result += RUN_TEST(test_read_mesh_parallel_trivial);
#if defined(MOAB_HAVE_MPI) && defined(MOAB_HAVE_ZOLTAN)
  result += RUN_TEST(test_read_mesh_parallel_rcbzoltan);
#endif

  result += RUN_TEST(test_gather_onevar_on_rank0);
  result += RUN_TEST(test_gather_onevar_on_rank1);

  result += RUN_TEST(test_multiple_loads_of_same_file);

  MPI_Finalize();
  return result;
}

void test_read_onevar_trivial()
{
  read_one_cell_var(false);
}

void test_read_onevar_rcbzoltan()
{
  read_one_cell_var(true);
}

void test_read_mesh_parallel_trivial()
{
  read_mesh_parallel(false);
}

void test_read_mesh_parallel_rcbzoltan()
{
  read_mesh_parallel(true);
}

void test_gather_onevar_on_rank0()
{
  gather_one_cell_var(0);
}

void test_gather_onevar_on_rank1()
{
  gather_one_cell_var(1);
}

void test_multiple_loads_of_same_file()
{
  multiple_loads_of_same_file();
}

// Helper functions
void read_one_cell_var(bool rcbzoltan)
{
  Core moab;
  Interface& mb = moab;

  read_options = "PARALLEL=READ_PART;PARTITION_METHOD=TRIVIAL;NO_EDGES;VARIABLE=vorticity";
  if (rcbzoltan)
    read_options = "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN;NO_EDGES;VARIABLE=vorticity;DEBUG_IO=1";

  ErrorCode rval = mb.load_file(example, NULL, read_options.c_str());
  CHECK_ERR(rval);

  // Get local edges
  Range local_edges;
  rval = mb.get_entities_by_type(0, MBEDGE, local_edges);
  CHECK_ERR(rval);
  CHECK_EQUAL((size_t)0, local_edges.size());

  // Get local cells
  Range local_cells;
  rval = mb.get_entities_by_type(0, MBPOLYGON, local_cells);
  CHECK_ERR(rval);
  // No mixed elements
  CHECK_EQUAL((size_t)1, local_cells.psize());

  Tag gid_tag;
  rval = mb.tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid_tag, MB_TAG_DENSE);
  CHECK_ERR(rval);

  std::vector<int> gids(local_cells.size());
  rval = mb.tag_get_data(gid_tag, local_cells, &gids[0]);
  Range local_cell_gids;
  std::copy(gids.rbegin(), gids.rend(), range_inserter(local_cell_gids));

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  int procs = pcomm->proc_config().proc_size();
  int rank = pcomm->proc_config().proc_rank();

  // Make check runs this test on two processors
  if (2 == procs) {
    // Check tag for cell variable vorticity at timestep 0
    Tag vorticity_tag0;
    rval = mb.tag_get_handle("vorticity0", layers, MB_TYPE_DOUBLE, vorticity_tag0);
    CHECK_ERR(rval);

    // Check tag for cell variable vorticity at timestep 1
    Tag vorticity_tag1;
    rval = mb.tag_get_handle("vorticity1", layers, MB_TYPE_DOUBLE, vorticity_tag1);
    CHECK_ERR(rval);

    // Get vorticity0 and vorticity1 tag values on 3 local cells
    double vorticity0_val[3 * layers];
    double vorticity1_val[3 * layers];

    if (rcbzoltan) {
      CHECK_EQUAL((size_t)14, local_cell_gids.psize());

      if (0 == rank) {
        CHECK_EQUAL((size_t)319, local_cells.size());
        CHECK_EQUAL((size_t)319, local_cell_gids.size());

        CHECK_EQUAL(3, (int)local_cell_gids[0]);
        CHECK_EQUAL(162, (int)local_cell_gids[159]);
        CHECK_EQUAL(642, (int)local_cell_gids[318]);

        EntityHandle cell_ents[] = {local_cells[0], local_cells[159], local_cells[318]};
        rval = mb.tag_get_data(vorticity_tag0, cell_ents, 3, vorticity0_val);
        CHECK_ERR(rval);

        // Timestep 0
        // Layer 0
        CHECK_REAL_EQUAL(-0.725999, vorticity0_val[0 * layers], eps);
        CHECK_REAL_EQUAL(-1.814997, vorticity0_val[1 * layers], eps);
        CHECK_REAL_EQUAL(-0.554888, vorticity0_val[2 * layers], eps);
        // Layer 1
        CHECK_REAL_EQUAL(-0.725989, vorticity0_val[0 * layers + 1], eps);
        CHECK_REAL_EQUAL(-1.814972, vorticity0_val[1 * layers + 1], eps);
        CHECK_REAL_EQUAL(-0.554881, vorticity0_val[2 * layers + 1], eps);

        rval = mb.tag_get_data(vorticity_tag1, cell_ents, 3, vorticity1_val);
        CHECK_ERR(rval);

        // Timestep 1
        // Layer 0
        CHECK_REAL_EQUAL(-0.706871, vorticity1_val[0 * layers], eps);
        CHECK_REAL_EQUAL(-1.767178, vorticity1_val[1 * layers], eps);
        CHECK_REAL_EQUAL(-0.540269, vorticity1_val[2 * layers], eps);
        // Layer 1
        CHECK_REAL_EQUAL(-0.706861, vorticity1_val[0 * layers + 1], eps);
        CHECK_REAL_EQUAL(-1.767153, vorticity1_val[1 * layers + 1], eps);
        CHECK_REAL_EQUAL(-0.540262, vorticity1_val[2 * layers + 1], eps);
      }
      else if (1 == rank) {
        CHECK_EQUAL((size_t)323, local_cells.size());
        CHECK_EQUAL((size_t)323, local_cell_gids.size());

        CHECK_EQUAL(1, (int)local_cell_gids[0]);
        CHECK_EQUAL(365, (int)local_cell_gids[161]);
        CHECK_EQUAL(557, (int)local_cell_gids[322]);

        EntityHandle cell_ents[] = {local_cells[0], local_cells[161], local_cells[322]};
        rval = mb.tag_get_data(vorticity_tag0, cell_ents, 3, vorticity0_val);
        CHECK_ERR(rval);

        // Timestep 0
        // Layer 0
        CHECK_REAL_EQUAL(3.629994, vorticity0_val[0 * layers], eps);
        CHECK_REAL_EQUAL(-1.173971, vorticity0_val[1 * layers], eps);
        CHECK_REAL_EQUAL(3.526371, vorticity0_val[2 * layers], eps);
        // Layer 1
        CHECK_REAL_EQUAL(3.629944, vorticity0_val[0 * layers + 1], eps);
        CHECK_REAL_EQUAL(-1.173955, vorticity0_val[1 * layers + 1], eps);
        CHECK_REAL_EQUAL(3.526322, vorticity0_val[2 * layers + 1], eps);

        rval = mb.tag_get_data(vorticity_tag1, cell_ents, 3, vorticity1_val);
        CHECK_ERR(rval);

        // Timestep 1
        // Layer 0
        CHECK_REAL_EQUAL(3.534355, vorticity1_val[0 * layers], eps);
        CHECK_REAL_EQUAL(-1.143041, vorticity1_val[1 * layers], eps);
        CHECK_REAL_EQUAL(3.433463, vorticity1_val[2 * layers], eps);
        // Layer 1
        CHECK_REAL_EQUAL(3.534306, vorticity1_val[0 * layers + 1], eps);
        CHECK_REAL_EQUAL(-1.143025, vorticity1_val[1 * layers + 1], eps);
        CHECK_REAL_EQUAL(3.433415, vorticity1_val[2 * layers + 1], eps);
      }
    }
    else {
      CHECK_EQUAL((size_t)321, local_cells.size());
      CHECK_EQUAL((size_t)321, local_cell_gids.size());
      CHECK_EQUAL((size_t)1, local_cell_gids.psize());

      EntityHandle cell_ents[] = {local_cells[0], local_cells[160], local_cells[320]};
      rval = mb.tag_get_data(vorticity_tag0, cell_ents, 3, vorticity0_val);
      CHECK_ERR(rval);

      rval = mb.tag_get_data(vorticity_tag1, cell_ents, 3, vorticity1_val);
      CHECK_ERR(rval);

      if (0 == rank) {
        CHECK_EQUAL(1, (int)local_cell_gids[0]);
        CHECK_EQUAL(161, (int)local_cell_gids[160]);
        CHECK_EQUAL(321, (int)local_cell_gids[320]);

        // Timestep 0
        // Layer 0
        CHECK_REAL_EQUAL(3.629994, vorticity0_val[0 * layers], eps);
        CHECK_REAL_EQUAL(-1.708188, vorticity0_val[1 * layers], eps);
        CHECK_REAL_EQUAL(0.131688, vorticity0_val[2 * layers], eps);
        // Layer 1
        CHECK_REAL_EQUAL(3.629944, vorticity0_val[0 * layers + 1], eps);
        CHECK_REAL_EQUAL(-1.708164, vorticity0_val[1 * layers + 1], eps);
        CHECK_REAL_EQUAL(0.131686, vorticity0_val[2 * layers + 1], eps);

        // Timestep 1
        // Layer 0
        CHECK_REAL_EQUAL(3.534355, vorticity1_val[0 * layers], eps);
        CHECK_REAL_EQUAL(-1.663182, vorticity1_val[1 * layers], eps);
        CHECK_REAL_EQUAL(0.128218, vorticity1_val[2 * layers], eps);
        // Layer 1
        CHECK_REAL_EQUAL(3.534306, vorticity1_val[0 * layers + 1], eps);
        CHECK_REAL_EQUAL(-1.663160, vorticity1_val[1 * layers + 1], eps);
        CHECK_REAL_EQUAL(0.128216, vorticity1_val[2 * layers + 1], eps);
      }
      else if (1 == rank) {
        CHECK_EQUAL(322, (int)local_cell_gids[0]);
        CHECK_EQUAL(482, (int)local_cell_gids[160]);
        CHECK_EQUAL(642, (int)local_cell_gids[320]);

        // Timestep 0
        // Layer 0
        CHECK_REAL_EQUAL(-0.554888, vorticity0_val[0 * layers], eps);
        CHECK_REAL_EQUAL(2.434397, vorticity0_val[1 * layers], eps);
        CHECK_REAL_EQUAL(-0.554888, vorticity0_val[2 * layers], eps);
        // Layer 1
        CHECK_REAL_EQUAL(-0.554881, vorticity0_val[0 * layers + 1], eps);
        CHECK_REAL_EQUAL(2.434363, vorticity0_val[1 * layers + 1], eps);
        CHECK_REAL_EQUAL(-0.554881, vorticity0_val[2 * layers + 1], eps);

        // Timestep 1
        // Layer 0
        CHECK_REAL_EQUAL(-0.540269, vorticity1_val[0 * layers], eps);
        CHECK_REAL_EQUAL(2.370258, vorticity1_val[1 * layers], eps);
        CHECK_REAL_EQUAL(-0.540269, vorticity1_val[2 * layers], eps);
        // Layer 1
        CHECK_REAL_EQUAL(-0.540262, vorticity1_val[0 * layers + 1], eps);
        CHECK_REAL_EQUAL(2.370226, vorticity1_val[1 * layers + 1], eps);
        CHECK_REAL_EQUAL(-0.540262, vorticity1_val[2 * layers + 1], eps);
      }
    }
  }
}

void read_mesh_parallel(bool rcbzoltan)
{
  Core moab;
  Interface& mb = moab;

  read_options = "PARALLEL=READ_PART;PARTITION_METHOD=TRIVIAL;PARALLEL_RESOLVE_SHARED_ENTS;VARIABLE=";
  if (rcbzoltan)
    read_options = "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN;PARALLEL_RESOLVE_SHARED_ENTS;VARIABLE=";

  ErrorCode rval = mb.load_file(example, NULL, read_options.c_str());
  CHECK_ERR(rval);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  int procs = pcomm->proc_config().proc_size();
  int rank = pcomm->proc_config().proc_rank();

  rval = pcomm->check_all_shared_handles();
  CHECK_ERR(rval);

  // Get local vertices
  Range local_verts;
  rval = mb.get_entities_by_type(0, MBVERTEX, local_verts);
  CHECK_ERR(rval);

  int verts_num = local_verts.size();
  if (2 == procs) {
    if (rcbzoltan) {
      if (0 == rank)
        CHECK_EQUAL(684, verts_num);
      else if (1 == rank)
        CHECK_EQUAL(691, verts_num); // Not owned vertices included
    }
    else {
      if (0 == rank)
        CHECK_EQUAL(687, verts_num);
      else if (1 == rank)
        CHECK_EQUAL(688, verts_num); // Not owned vertices included
    }
  }

  rval = pcomm->filter_pstatus(local_verts, PSTATUS_NOT_OWNED, PSTATUS_NOT);
  CHECK_ERR(rval);

  verts_num = local_verts.size();
  if (2 == procs) {
    if (rcbzoltan) {
      if (0 == rank)
        CHECK_EQUAL(684, verts_num);
      else if (1 == rank)
        CHECK_EQUAL(596, verts_num); // Not owned vertices excluded
    }
    else {
      if (0 == rank)
        CHECK_EQUAL(687, verts_num);
      else if (1 == rank)
        CHECK_EQUAL(593, verts_num); // Not owned vertices excluded
    }
  }

  // Get local edges
  Range local_edges;
  rval = mb.get_entities_by_type(0, MBEDGE, local_edges);
  CHECK_ERR(rval);

  int edges_num = local_edges.size();
  if (2 == procs) {
    if (rcbzoltan) {
      if (0 == rank)
        CHECK_EQUAL(1002, edges_num);
      else if (1 == rank)
        CHECK_EQUAL(1013, edges_num); // Not owned edges included
    }
    else {
      if (0 == rank)
        CHECK_EQUAL(1007, edges_num);
      else if (1 == rank)
        CHECK_EQUAL(1008, edges_num); // Not owned edges included
    }
  }

  rval = pcomm->filter_pstatus(local_edges, PSTATUS_NOT_OWNED, PSTATUS_NOT);
  CHECK_ERR(rval);

  edges_num = local_edges.size();
  if (2 == procs) {
    if (rcbzoltan) {
      if (0 == rank)
        CHECK_EQUAL(1002, edges_num);
      else if (1 == rank)
        CHECK_EQUAL(918, edges_num); // Not owned edges excluded
    }
    else {
      if (0 == rank)
        CHECK_EQUAL(1007, edges_num);
      else if (1 == rank)
        CHECK_EQUAL(913, edges_num); // Not owned edges excluded
    }
  }

  // Get local cells
  Range local_cells;
  rval = mb.get_entities_by_type(0, MBPOLYGON, local_cells);
  CHECK_ERR(rval);
  // No mixed elements
  CHECK_EQUAL((size_t)1, local_cells.psize());

  int cells_num = local_cells.size();
  if (2 == procs) {
    if (rcbzoltan) {
      if (0 == rank)
        CHECK_EQUAL(319, cells_num);
      else
        CHECK_EQUAL(323, cells_num);
    }
    else
      CHECK_EQUAL(321, cells_num);
  }

  rval = pcomm->filter_pstatus(local_cells, PSTATUS_NOT_OWNED, PSTATUS_NOT);
  CHECK_ERR(rval);

  cells_num = local_cells.size();
  if (2 == procs) {
    if (rcbzoltan) {
      if (0 == rank)
        CHECK_EQUAL(319, cells_num);
      else
        CHECK_EQUAL(323, cells_num);
    }
    else
      CHECK_EQUAL(321, cells_num);
  }

  std::cout << "proc: " << rank << " verts:" << verts_num << "\n";

  int total_verts_num;
  MPI_Reduce(&verts_num, &total_verts_num, 1, MPI_INT, MPI_SUM, 0, pcomm->proc_config().proc_comm());
  if (0 == rank) {
    std::cout << "total vertices: " << total_verts_num << "\n";
    CHECK_EQUAL(1280, total_verts_num);
  }

  std::cout << "proc: " << rank << " edges:" << edges_num << "\n";

  int total_edges_num;
  MPI_Reduce(&edges_num, &total_edges_num, 1, MPI_INT, MPI_SUM, 0, pcomm->proc_config().proc_comm());
  if (0 == rank) {
    std::cout << "total edges: " << total_edges_num << "\n";
    CHECK_EQUAL(1920, total_edges_num);
  }

  std::cout << "proc: " << rank << " cells:" << cells_num << "\n";

  int total_cells_num;
  MPI_Reduce(&cells_num, &total_cells_num, 1, MPI_INT, MPI_SUM, 0, pcomm->proc_config().proc_comm());
  if (0 == rank) {
    std::cout << "total cells: " << total_cells_num << "\n";
    CHECK_EQUAL(642, total_cells_num);
  }

#ifdef MOAB_HAVE_HDF5_PARALLEL
  std::string write_options("PARALLEL=WRITE_PART;");

  std::string output_file = "test_gcrm";
  if (rcbzoltan)
    output_file += "_rcbzoltan";
  output_file += ".h5m";

  mb.write_file(output_file.c_str(), NULL, write_options.c_str());
#endif
}

void gather_one_cell_var(int gather_set_rank)
{
  Core moab;
  Interface& mb = moab;

  EntityHandle file_set;
  ErrorCode rval = mb.create_meshset(MESHSET_SET, file_set);
  CHECK_ERR(rval);

  read_options = "PARALLEL=READ_PART;PARTITION_METHOD=TRIVIAL;PARALLEL_RESOLVE_SHARED_ENTS";
  std::ostringstream gather_set_option;
  gather_set_option << ";GATHER_SET=" << gather_set_rank;
  read_options += gather_set_option.str();

  rval = mb.load_file(example, &file_set, read_options.c_str());
  CHECK_ERR(rval);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  int procs = pcomm->proc_config().proc_size();
  int rank = pcomm->proc_config().proc_rank();

  // Make sure gather_set_rank is valid
  if (gather_set_rank < 0 || gather_set_rank >= procs)
    return;

  Range cells, cells_owned;
  rval = mb.get_entities_by_type(file_set, MBPOLYGON, cells);
  CHECK_ERR(rval);

  // Get local owned cells
  rval = pcomm->filter_pstatus(cells, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &cells_owned);
  CHECK_ERR(rval);

  EntityHandle gather_set = 0;
  if (gather_set_rank == rank) {
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

  pcomm->gather_data(cells_owned, vorticity_tag0, gid_tag, gather_set, gather_set_rank);

  if (gather_set_rank == rank) {
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
    rval = mb.tag_get_data(vorticity_tag0, &cell_ents[0], 4, vorticity0_val);
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
}

void multiple_loads_of_same_file()
{
  Core moab;
  Interface& mb = moab;

  // Need a file set for nomesh to work right
  EntityHandle file_set;
  ErrorCode rval;
  rval = mb.create_meshset(MESHSET_SET, file_set);
  CHECK_ERR(rval);

  // Read first only header information, no mesh, no variable
  read_options = "PARALLEL=READ_PART;PARTITION;NOMESH;VARIABLE=;PARTITION_METHOD=TRIVIAL";

  rval = mb.load_file(example, &file_set, read_options.c_str());
  CHECK_ERR(rval);

  // Create mesh, no variable
  read_options = "PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION_METHOD=TRIVIAL;VARIABLE=";

  rval = mb.load_file(example, &file_set, read_options.c_str());
  CHECK_ERR(rval);

  // Read variable vorticity at timestep 0, no mesh
  read_options = "PARALLEL=READ_PART;PARTITION;PARTITION_METHOD=TRIVIAL;NOMESH;VARIABLE=vorticity;TIMESTEP=0";

  rval = mb.load_file(example, &file_set, read_options.c_str());
  CHECK_ERR(rval);

  Range local_verts;
  rval = mb.get_entities_by_type(file_set, MBVERTEX, local_verts);
  CHECK_ERR(rval);

  Range local_edges;
  rval = mb.get_entities_by_type(file_set, MBEDGE, local_edges);
  CHECK_ERR(rval);

  Range local_cells;
  rval = mb.get_entities_by_type(file_set, MBPOLYGON, local_cells);
  CHECK_ERR(rval);
  // No mixed elements
  CHECK_EQUAL((size_t)1, local_cells.psize());

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  int procs = pcomm->proc_config().proc_size();
  int rank = pcomm->proc_config().proc_rank();

  // Make check runs this test on two processors
  if (2 == procs) {
    CHECK_EQUAL((size_t)321, local_cells.size());

    // Check tag for cell variable vorticity at timestep 0
    Tag vorticity_tag0;
    rval = mb.tag_get_handle("vorticity0", layers, MB_TYPE_DOUBLE, vorticity_tag0);
    CHECK_ERR(rval);

    // Get vorticity0 tag values on 3 local cells
    double vorticity0_val[3 * layers];
    EntityHandle cell_ents[] = {local_cells[0], local_cells[160], local_cells[320]};
    rval = mb.tag_get_data(vorticity_tag0, cell_ents, 3, vorticity0_val);
    CHECK_ERR(rval);

    if (0 == rank) {
      CHECK_EQUAL((size_t)687, local_verts.size());
      CHECK_EQUAL((size_t)1007, local_edges.size());

      // Layer 0
      CHECK_REAL_EQUAL(3.629994, vorticity0_val[0 * layers], eps);
      CHECK_REAL_EQUAL(-1.708188, vorticity0_val[1 * layers], eps);
      CHECK_REAL_EQUAL(0.131688, vorticity0_val[2 * layers], eps);
      // Layer 1
      CHECK_REAL_EQUAL(3.629944, vorticity0_val[0 * layers + 1], eps);
      CHECK_REAL_EQUAL(-1.708164, vorticity0_val[1 * layers + 1], eps);
      CHECK_REAL_EQUAL(0.131686, vorticity0_val[2 * layers + 1], eps);
    }
    else if (1 == rank) {
      CHECK_EQUAL((size_t)688, local_verts.size());
      CHECK_EQUAL((size_t)1008, local_edges.size());

      // Layer 0
      CHECK_REAL_EQUAL(-0.554888, vorticity0_val[0 * layers], eps);
      CHECK_REAL_EQUAL(2.434397, vorticity0_val[1 * layers], eps);
      CHECK_REAL_EQUAL(-0.554888, vorticity0_val[2 * layers], eps);
      // Layer 1
      CHECK_REAL_EQUAL(-0.554881, vorticity0_val[0 * layers + 1], eps);
      CHECK_REAL_EQUAL(2.434363, vorticity0_val[1 * layers + 1], eps);
      CHECK_REAL_EQUAL(-0.554881, vorticity0_val[2 * layers + 1], eps);
    }
  }
}
