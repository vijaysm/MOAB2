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
static const char example[] = STRINGIFY(MESHDIR) "/io/mpasx1.642.t.2.nc";
#endif

void test_read_onevar_trivial();
void test_read_onevar_trivial_no_mixed_elements();
#if defined(MOAB_HAVE_MPI) && defined(MOAB_HAVE_ZOLTAN)
void test_read_onevar_rcbzoltan();
void test_read_onevar_rcbzoltan_no_mixed_elements();
#endif

void test_read_mesh_parallel_trivial();
void test_read_mesh_parallel_trivial_no_mixed_elements();
#if defined(MOAB_HAVE_MPI) && defined(MOAB_HAVE_ZOLTAN)
void test_read_mesh_parallel_rcbzoltan();
void test_read_mesh_parallel_rcbzoltan_no_mixed_elements();
#endif

void test_gather_onevar_on_rank0();
void test_gather_onevar_on_rank1();

void test_multiple_loads_of_same_file();
void test_multiple_loads_of_same_file_no_mixed_elements();

// Helper functions
void read_one_cell_var(bool rcbzoltan, bool no_mixed_elements);
void read_mesh_parallel(bool rcbzoltan, bool no_mixed_elements);
void gather_one_cell_var(int gather_set_rank);
void multiple_loads_of_same_file(bool no_mixed_elements);

std::string read_options;
const double eps = 1e-20;

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int result = 0;

  result += RUN_TEST(test_read_onevar_trivial);
  result += RUN_TEST(test_read_onevar_trivial_no_mixed_elements);
#if defined(MOAB_HAVE_MPI) && defined(MOAB_HAVE_ZOLTAN)
  result += RUN_TEST(test_read_onevar_rcbzoltan);
  result += RUN_TEST(test_read_onevar_rcbzoltan_no_mixed_elements);
#endif

  result += RUN_TEST(test_read_mesh_parallel_trivial);
  result += RUN_TEST(test_read_mesh_parallel_trivial_no_mixed_elements);
#if defined(MOAB_HAVE_MPI) && defined(MOAB_HAVE_ZOLTAN)
  result += RUN_TEST(test_read_mesh_parallel_rcbzoltan);
  result += RUN_TEST(test_read_mesh_parallel_rcbzoltan_no_mixed_elements);
#endif

  result += RUN_TEST(test_gather_onevar_on_rank0);
  result += RUN_TEST(test_gather_onevar_on_rank1);

  result += RUN_TEST(test_multiple_loads_of_same_file);
  result += RUN_TEST(test_multiple_loads_of_same_file_no_mixed_elements);

  MPI_Finalize();
  return result;
}

void test_read_onevar_trivial()
{
  read_one_cell_var(false, false);
}

void test_read_onevar_trivial_no_mixed_elements()
{
  read_one_cell_var(false, true);
}

void test_read_onevar_rcbzoltan()
{
  read_one_cell_var(true, false);
}

void test_read_onevar_rcbzoltan_no_mixed_elements()
{
  read_one_cell_var(true, true);
}

void test_read_mesh_parallel_trivial()
{
  read_mesh_parallel(false, false);
}

void test_read_mesh_parallel_trivial_no_mixed_elements()
{
  read_mesh_parallel(false, true);
}

void test_read_mesh_parallel_rcbzoltan()
{
  read_mesh_parallel(true, false);
}

void test_read_mesh_parallel_rcbzoltan_no_mixed_elements()
{
  read_mesh_parallel(true, true);
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
  multiple_loads_of_same_file(false);
}

void test_multiple_loads_of_same_file_no_mixed_elements()
{
  multiple_loads_of_same_file(true);
}

// Helper functions
void read_one_cell_var(bool rcbzoltan, bool no_mixed_elements)
{
  Core moab;
  Interface& mb = moab;

  read_options = "PARALLEL=READ_PART;PARTITION_METHOD=TRIVIAL;NO_EDGES;VARIABLE=ke";
  if (rcbzoltan)
    read_options = "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN;NO_EDGES;VARIABLE=ke";

  if (no_mixed_elements)
    read_options += ";NO_MIXED_ELEMENTS";

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
    CHECK_EQUAL((size_t)321, local_cells.size());
    CHECK_EQUAL((size_t)321, local_cell_gids.size());

    // Check tag for cell variable ke at timestep 0
    Tag ke_tag0;
    rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
    CHECK_ERR(rval);

    // Check tag for cell variable ke at timestep 1
    Tag ke_tag1;
    rval = mb.tag_get_handle("ke1", 1, MB_TYPE_DOUBLE, ke_tag1);
    CHECK_ERR(rval);

    // Get ke0 and ke1 tag values on 3 local cells
    EntityHandle cell_ents[] = {local_cells[0], local_cells[160], local_cells[320]};
    double ke0_val[3];
    rval = mb.tag_get_data(ke_tag0, cell_ents, 3, ke0_val);
    CHECK_ERR(rval);
    double ke1_val[3];
    rval = mb.tag_get_data(ke_tag1, cell_ents, 3, ke1_val);
    CHECK_ERR(rval);

    if (rcbzoltan) {
      if (no_mixed_elements)
        CHECK_EQUAL((size_t)1, local_cells.psize());
      else
        CHECK_EQUAL((size_t)2, local_cells.psize());

      CHECK_EQUAL((size_t)37, local_cell_gids.psize());

      if (0 == rank) {
        CHECK_EQUAL(1, (int)local_cell_gids[0]);
        CHECK_EQUAL(284, (int)local_cell_gids[160]);
        CHECK_EQUAL(630, (int)local_cell_gids[320]);

        CHECK_REAL_EQUAL(15.001, ke0_val[0], eps);
        CHECK_REAL_EQUAL(16.284, ke0_val[1], eps);
        CHECK_REAL_EQUAL(16.630, ke0_val[2], eps);
        CHECK_REAL_EQUAL(25.001, ke1_val[0], eps);
        CHECK_REAL_EQUAL(26.284, ke1_val[1], eps);
        CHECK_REAL_EQUAL(26.630, ke1_val[2], eps);
      }
      else if (1 == rank) {
        CHECK_EQUAL(4, (int)local_cell_gids[0]);
        CHECK_EQUAL(341, (int)local_cell_gids[160]);
        CHECK_EQUAL(642, (int)local_cell_gids[320]);

        CHECK_REAL_EQUAL(15.004, ke0_val[0], eps);
        CHECK_REAL_EQUAL(16.341, ke0_val[1], eps);
        CHECK_REAL_EQUAL(16.642, ke0_val[2], eps);
        CHECK_REAL_EQUAL(25.004, ke1_val[0], eps);
        CHECK_REAL_EQUAL(26.341, ke1_val[1], eps);
        CHECK_REAL_EQUAL(26.642, ke1_val[2], eps);
      }
    }
    else {
      CHECK_EQUAL((size_t)1, local_cell_gids.psize());

      if (0 == rank) {
        if (no_mixed_elements)
          CHECK_EQUAL((size_t)1, local_cells.psize());
        else
          CHECK_EQUAL((size_t)2, local_cells.psize());
        CHECK_EQUAL(1, (int)local_cell_gids[0]);
        CHECK_EQUAL(161, (int)local_cell_gids[160]);
        CHECK_EQUAL(321, (int)local_cell_gids[320]);

        CHECK_REAL_EQUAL(15.001, ke0_val[0], eps);
        CHECK_REAL_EQUAL(16.161, ke0_val[1], eps);
        CHECK_REAL_EQUAL(16.321, ke0_val[2], eps);
        CHECK_REAL_EQUAL(25.001, ke1_val[0], eps);
        CHECK_REAL_EQUAL(26.161, ke1_val[1], eps);
        CHECK_REAL_EQUAL(26.321, ke1_val[2], eps);
      }
      else if (1 == rank) {
        CHECK_EQUAL((size_t)1, local_cells.psize());
        CHECK_EQUAL(322, (int)local_cell_gids[0]);
        CHECK_EQUAL(482, (int)local_cell_gids[160]);
        CHECK_EQUAL(642, (int)local_cell_gids[320]);

        CHECK_REAL_EQUAL(16.322, ke0_val[0], eps);
        CHECK_REAL_EQUAL(16.482, ke0_val[1], eps);
        CHECK_REAL_EQUAL(16.642, ke0_val[2], eps);
        CHECK_REAL_EQUAL(26.322, ke1_val[0], eps);
        CHECK_REAL_EQUAL(26.482, ke1_val[1], eps);
        CHECK_REAL_EQUAL(26.642, ke1_val[2], eps);
      }
    }
  }
}

void read_mesh_parallel(bool rcbzoltan, bool no_mixed_elements)
{
  Core moab;
  Interface& mb = moab;

  read_options = "PARALLEL=READ_PART;PARTITION_METHOD=TRIVIAL;PARALLEL_RESOLVE_SHARED_ENTS;VARIABLE=";
  if (rcbzoltan)
    read_options = "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN;PARALLEL_RESOLVE_SHARED_ENTS;VARIABLE=";

  if (no_mixed_elements)
    read_options += ";NO_MIXED_ELEMENTS";

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
        CHECK_EQUAL(685, verts_num);
      else if (1 == rank)
        CHECK_EQUAL(685, verts_num); // Not owned vertices included
    }
    else {
      if (0 == rank)
        CHECK_EQUAL(1120, verts_num);
      else if (1 == rank)
        CHECK_EQUAL(1122, verts_num); // Not owned vertices included
    }
  }

  rval = pcomm->filter_pstatus(local_verts, PSTATUS_NOT_OWNED, PSTATUS_NOT);
  CHECK_ERR(rval);

  verts_num = local_verts.size();
  if (2 == procs) {
    if (rcbzoltan) {
      if (0 == rank)
        CHECK_EQUAL(685, verts_num);
      else if (1 == rank)
        CHECK_EQUAL(595, verts_num); // Not owned vertices excluded
    }
    else {
      if (0 == rank)
        CHECK_EQUAL(1120, verts_num);
      else if (1 == rank)
        CHECK_EQUAL(160, verts_num); // Not owned vertices excluded
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
        CHECK_EQUAL(1005, edges_num);
      else if (1 == rank)
        CHECK_EQUAL(1005, edges_num); // Not owned edges included
    }
    else {
      if (0 == rank)
        CHECK_EQUAL(1438, edges_num);
      else if (1 == rank)
        CHECK_EQUAL(1444, edges_num); // Not owned edges included
    }
  }

  rval = pcomm->filter_pstatus(local_edges, PSTATUS_NOT_OWNED, PSTATUS_NOT);
  CHECK_ERR(rval);

  edges_num = local_edges.size();
  if (2 == procs) {
    if (rcbzoltan) {
      if (0 == rank)
        CHECK_EQUAL(1005, edges_num);
      else if (1 == rank)
        CHECK_EQUAL(915, edges_num); // Not owned edges excluded
    }
    else {
      if (0 == rank)
        CHECK_EQUAL(1438, edges_num);
      else if (1 == rank)
        CHECK_EQUAL(482, edges_num); // Not owned edges excluded
    }
  }

  // Get local cells
  Range local_cells;
  rval = mb.get_entities_by_type(0, MBPOLYGON, local_cells);
  CHECK_ERR(rval);

  int cells_num = local_cells.size();
  if (2 == procs) {
    CHECK_EQUAL(321, cells_num);

    if (rcbzoltan) {
      if (no_mixed_elements)
        CHECK_EQUAL((size_t)1, local_cells.psize());
      else
        CHECK_EQUAL((size_t)2, local_cells.psize());
     }
    else {
      if (0 == rank) {
        if (no_mixed_elements)
          CHECK_EQUAL((size_t)1, local_cells.psize());
        else
          CHECK_EQUAL((size_t)2, local_cells.psize());
      }
      else if (1 == rank)
        CHECK_EQUAL((size_t)1, local_cells.psize());
    }
  }

  rval = pcomm->filter_pstatus(local_cells, PSTATUS_NOT_OWNED, PSTATUS_NOT);
  CHECK_ERR(rval);

  cells_num = local_cells.size();
  if (2 == procs) {
    CHECK_EQUAL(321, cells_num);

    if (rcbzoltan) {
      if (no_mixed_elements)
        CHECK_EQUAL((size_t)1, local_cells.psize());
      else
        CHECK_EQUAL((size_t)2, local_cells.psize());
    }
    else {
      if (0 == rank) {
        if (no_mixed_elements)
          CHECK_EQUAL((size_t)1, local_cells.psize());
        else
          CHECK_EQUAL((size_t)2, local_cells.psize());
      }
      else if (1 == rank)
        CHECK_EQUAL((size_t)1, local_cells.psize());
    }
  }

  std::cout << "proc: " << rank << " verts:" << verts_num << "\n";

  int total_verts_num;
  MPI_Reduce(&verts_num, &total_verts_num, 1, MPI_INTEGER, MPI_SUM, 0, pcomm->proc_config().proc_comm());
  if (0 == rank) {
    std::cout << "total vertices: " << total_verts_num << "\n";
    CHECK_EQUAL(1280, total_verts_num);
  }

  std::cout << "proc: " << rank << " edges:" << edges_num << "\n";

  int total_edges_num;
  MPI_Reduce(&edges_num, &total_edges_num, 1, MPI_INTEGER, MPI_SUM, 0, pcomm->proc_config().proc_comm());
  if (0 == rank) {
    std::cout << "total edges: " << total_edges_num << "\n";
    CHECK_EQUAL(1920, total_edges_num);
  }

  std::cout << "proc: " << rank << " cells:" << cells_num << "\n";

  int total_cells_num;
  MPI_Reduce(&cells_num, &total_cells_num, 1, MPI_INTEGER, MPI_SUM, 0, pcomm->proc_config().proc_comm());
  if (0 == rank) {
    std::cout << "total cells: " << total_cells_num << "\n";
    CHECK_EQUAL(642, total_cells_num);
  }

#ifdef MOAB_HAVE_HDF5_PARALLEL
  std::string write_options("PARALLEL=WRITE_PART;");

  std::string output_file = "test_mpas";
  if (rcbzoltan)
    output_file += "_rcbzoltan";
  if (no_mixed_elements)
    output_file += "_no_mixed_elements";
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

  Tag ke_tag0, gid_tag;
  rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0, MB_TAG_DENSE);
  CHECK_ERR(rval);

  rval = mb.tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid_tag, MB_TAG_DENSE);
  CHECK_ERR(rval);

  pcomm->gather_data(cells_owned, ke_tag0, gid_tag, gather_set, gather_set_rank);

  if (gather_set_rank == rank) {
    // Get gather set cells
    Range gather_set_cells;
    rval = mb.get_entities_by_type(gather_set, MBPOLYGON, gather_set_cells);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)642, gather_set_cells.size());
    CHECK_EQUAL((size_t)2, gather_set_cells.psize());

    // Check ke0 tag values on 4 gather set cells: first pentagon, last pentagon,
    // first hexagon and last hexagon
    EntityHandle cell_ents[] = {gather_set_cells[0], gather_set_cells[11],
                                gather_set_cells[12], gather_set_cells[641]};
    double ke0_val[4];
    rval = mb.tag_get_data(ke_tag0, &cell_ents[0], 4, ke0_val);
    CHECK_ERR(rval);

    CHECK_REAL_EQUAL(15.001, ke0_val[0], eps);
    CHECK_REAL_EQUAL(15.012, ke0_val[1], eps);
    CHECK_REAL_EQUAL(16.013, ke0_val[2], eps);
    CHECK_REAL_EQUAL(16.642, ke0_val[3], eps);
  }
}

void multiple_loads_of_same_file(bool no_mixed_elements)
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
  if (no_mixed_elements)
    read_options += ";NO_MIXED_ELEMENTS";

  rval = mb.load_file(example, &file_set, read_options.c_str());
  CHECK_ERR(rval);

  // Create mesh, no variable
  read_options = "PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION_METHOD=TRIVIAL;VARIABLE=";
  if (no_mixed_elements)
    read_options += ";NO_MIXED_ELEMENTS";

  rval = mb.load_file(example, &file_set, read_options.c_str());
  CHECK_ERR(rval);

  // Read variable ke at timestep 0, no mesh
  read_options = "PARALLEL=READ_PART;PARTITION;PARTITION_METHOD=TRIVIAL;NOMESH;VARIABLE=ke;TIMESTEP=0";
  if (no_mixed_elements)
    read_options += ";NO_MIXED_ELEMENTS";

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

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  int procs = pcomm->proc_config().proc_size();
  int rank = pcomm->proc_config().proc_rank();

  // Make check runs this test on two processors
  if (2 == procs) {
    CHECK_EQUAL((size_t)321, local_cells.size());

    // Check tag for cell variable ke at timestep 0
    Tag ke_tag0;
    rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
    CHECK_ERR(rval);

    // Get ke0 tag values on 3 local cells
    EntityHandle cell_ents[] = {local_cells[0], local_cells[160], local_cells[320]};
    double ke0_val[3];
    rval = mb.tag_get_data(ke_tag0, cell_ents, 3, ke0_val);
    CHECK_ERR(rval);

    if (0 == rank) {
      CHECK_EQUAL((size_t)1120, local_verts.size());
      CHECK_EQUAL((size_t)1438, local_edges.size());
      if (no_mixed_elements)
        CHECK_EQUAL((size_t)1, local_cells.psize());
      else
        CHECK_EQUAL((size_t)2, local_cells.psize());

      CHECK_REAL_EQUAL(15.001, ke0_val[0], eps);
      CHECK_REAL_EQUAL(16.161, ke0_val[1], eps);
      CHECK_REAL_EQUAL(16.321, ke0_val[2], eps);
    }
    else if (1 == rank) {
      CHECK_EQUAL((size_t)1122, local_verts.size());
      CHECK_EQUAL((size_t)1444, local_edges.size());
      CHECK_EQUAL((size_t)1, local_cells.psize());

      CHECK_REAL_EQUAL(16.322, ke0_val[0], eps);
      CHECK_REAL_EQUAL(16.482, ke0_val[1], eps);
      CHECK_REAL_EQUAL(16.642, ke0_val[2], eps);
    }
  }
}
