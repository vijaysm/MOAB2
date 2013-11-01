#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/ProgOptions.hpp"
#include "MBParallelConventions.h"
#include "moab/ReadUtilIface.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/mpasx1.642.t.2.nc";
#endif

void test_read_parallel_mpas_trivial();
void test_read_parallel_mpas_trivial_no_mixed_elements();
void test_read_parallel_mpas_rcbzoltan();
void test_read_parallel_mpas_rcbzoltan_no_mixed_elements();
void test_read_parallel(int num_verts, bool test_nb_nodes, int num_edges, bool test_nb_edges,
                        int num_cells, bool test_nb_cells, bool mixed_elements);
void test_multiple_loads_of_same_file();
void test_multiple_loads_of_same_file_no_mixed_elements();

std::string partition_method;

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int result = 0;

  result += RUN_TEST(test_read_parallel_mpas_trivial);
  result += RUN_TEST(test_read_parallel_mpas_trivial_no_mixed_elements);
  result += RUN_TEST(test_read_parallel_mpas_rcbzoltan);
  result += RUN_TEST(test_read_parallel_mpas_rcbzoltan_no_mixed_elements);
  result += RUN_TEST(test_multiple_loads_of_same_file);
  result += RUN_TEST(test_multiple_loads_of_same_file_no_mixed_elements);

  MPI_Finalize();
  return result;
}

void test_read_parallel_mpas_trivial()
{
  partition_method = std::string(";PARTITION_METHOD=TRIVIAL;PARALLEL_RESOLVE_SHARED_ENTS");
  test_read_parallel(1280, true, 1920, true, 642, true, true);
}

void test_read_parallel_mpas_trivial_no_mixed_elements()
{
  partition_method = std::string(";PARTITION_METHOD=TRIVIAL;PARALLEL_RESOLVE_SHARED_ENTS;NO_MIXED_ELEMENTS");
  test_read_parallel(1280, true, 1920, true, 642, true, false);
}

void test_read_parallel_mpas_rcbzoltan()
{
  partition_method = std::string(";PARTITION_METHOD=RCBZOLTAN;PARALLEL_RESOLVE_SHARED_ENTS");
  test_read_parallel(1280, false, 1920, false, 642, false, true);
}

void test_read_parallel_mpas_rcbzoltan_no_mixed_elements()
{
  partition_method = std::string(";PARTITION_METHOD=RCBZOLTAN;PARALLEL_RESOLVE_SHARED_ENTS;NO_MIXED_ELEMENTS");
  test_read_parallel(1280, false, 1920, false, 642, false, true);
}

void test_read_parallel(int num_verts, bool test_nb_nodes, int num_edges, bool test_nb_edges,
                        int num_cells, bool test_nb_cells, bool mixed_elements)
{
  Core moab;
  Interface& mb = moab;
  EntityHandle file_set;
  ErrorCode rval;
  rval = mb.create_meshset(MESHSET_SET, file_set);
  CHECK_ERR(rval);

  std::string opt = std::string("PARALLEL=READ_PART") + partition_method;
  // Create gather set in processor 0
  opt += std::string(";GATHER_SET=0");
  rval = mb.load_file(example, &file_set, opt.c_str());
  CHECK_ERR(rval);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  int procs = pcomm->proc_config().proc_size();
  int rank = pcomm->proc_config().proc_rank();

  rval = pcomm->check_all_shared_handles();
  CHECK_ERR(rval);

  // Get the total # owned verts
  Range verts;
  rval = mb.get_entities_by_type(0, MBVERTEX, verts);
  CHECK_ERR(rval);

  int my_verts_num = verts.size();
  if (test_nb_nodes && 2 == procs) {
    if (0 == rank)
      CHECK_EQUAL(2400, my_verts_num); // Gather set vertices included
    else if (1 == rank)
      CHECK_EQUAL(1122, my_verts_num); // Not owned vertices included
  }

  rval = pcomm->filter_pstatus(verts, PSTATUS_NOT_OWNED, PSTATUS_NOT);
  CHECK_ERR(rval);

  my_verts_num = verts.size();
  if (test_nb_nodes && 2 == procs) {
    if (0 == rank)
      CHECK_EQUAL(2400, my_verts_num); // Gather set vertices included
    else if (1 == rank)
      CHECK_EQUAL(160, my_verts_num); // Not owned vertices excluded
  }

  // Get the total # owned edges
  Range edges;
  rval = mb.get_entities_by_type(0, MBEDGE, edges);
  CHECK_ERR(rval);

  int my_edges_num = edges.size();
  if (test_nb_edges && 2 == procs) {
    if (0 == rank)
      CHECK_EQUAL(3358, my_edges_num); // Gather set edges included
    else if (1 == rank)
      CHECK_EQUAL(1444, my_edges_num); // Not owned edges included
  }

  rval = pcomm->filter_pstatus(edges, PSTATUS_NOT_OWNED, PSTATUS_NOT);
  CHECK_ERR(rval);

  my_edges_num = edges.size();
  if (test_nb_edges && 2 == procs) {
    if (0 == rank)
      CHECK_EQUAL(3358, my_edges_num); // Gather set edges included
    else if (1 == rank)
      CHECK_EQUAL(482, my_edges_num); // Not owned edges excluded
  }

  // Get the total # owned cells
  Range cells;
  rval = mb.get_entities_by_type(0, MBPOLYGON, cells);
  CHECK_ERR(rval);

  int my_cells_num = cells.size();
  if (test_nb_cells && 2 == procs) {
    if (0 == rank) {
      CHECK_EQUAL(963, my_cells_num); // Gather set cells included
      if (mixed_elements)
        CHECK_EQUAL((size_t)4, cells.psize()); // Gather set cells included
      else
        CHECK_EQUAL((size_t)1, cells.psize()); // Gather set cells included
    }
    else if (1 == rank) {
      CHECK_EQUAL(321, my_cells_num); // Not owned cells included
      CHECK_EQUAL((size_t)1, cells.psize()); // Not owned cells included
    }
  }

   rval = pcomm->filter_pstatus(cells, PSTATUS_NOT_OWNED, PSTATUS_NOT);
   CHECK_ERR(rval);

   my_cells_num = cells.size();
   if (test_nb_cells && 2 == procs) {
     if (0 == rank) {
       CHECK_EQUAL(963, my_cells_num); // Gather set cells included
       if (mixed_elements)
         CHECK_EQUAL((size_t)4, cells.psize()); // Gather set cells included
       else
         CHECK_EQUAL((size_t)1, cells.psize()); // Gather set cells included
     }
     else if (1 == rank) {
       CHECK_EQUAL(321, my_cells_num); // Not owned cells excluded
       CHECK_EQUAL((size_t)1, cells.psize()); // Not owned cells excluded
     }
   }

  if (0 == rank) {
    // Get gather set
    EntityHandle gather_set;
    ReadUtilIface* readUtilIface;
    rval = mb.query_interface(readUtilIface);
    CHECK_ERR(rval);
    rval = readUtilIface->get_gather_set(gather_set);
    CHECK_ERR(rval);

    // Get gather set entities
    Range gather_ents;
    rval = mb.get_entities_by_handle(gather_set, gather_ents);
    CHECK_ERR(rval);

    // Remove gather set vertices in processor 0
    verts = subtract(verts, gather_ents);

    // Remove gather set edges in processor 0
    edges = subtract(edges, gather_ents);

    // Remove gather set cells in processor 0
    cells = subtract(cells, gather_ents);
  }

  my_verts_num = verts.size();
  if (test_nb_nodes && 2 == procs) {
    if (0 == rank)
      CHECK_EQUAL(1120, my_verts_num); // Gather set vertices excluded
    else if (1 == rank)
      CHECK_EQUAL(160, my_verts_num); // Not owned vertices excluded
  }

  std::cout << "proc: " << rank << " verts:" << my_verts_num << "\n";

  int total_verts;
  MPI_Reduce(&my_verts_num, &total_verts, 1, MPI_INTEGER, MPI_SUM, 0, pcomm->proc_config().proc_comm());
  if (0 == rank) {
    std::cout << "total vertices: " << total_verts << "\n";
    if (test_nb_nodes)
      CHECK_EQUAL(total_verts, num_verts);
  }

  my_edges_num = edges.size();
  if (test_nb_edges && 2 == procs) {
    if (0 == rank)
      CHECK_EQUAL(1438, my_edges_num); // Gather set edges excluded
    else if (1 == rank)
      CHECK_EQUAL(482, my_edges_num); // Not owned edges excluded
  }

  std::cout << "proc: " << rank << " edges:" << my_edges_num << "\n";

  int total_edges;
  MPI_Reduce(&my_edges_num, &total_edges, 1, MPI_INTEGER, MPI_SUM, 0, pcomm->proc_config().proc_comm());
  if (0 == rank) {
    std::cout << "total edges: " << total_edges << "\n";
    if (test_nb_edges)
      CHECK_EQUAL(total_edges, num_edges);
  }

  my_cells_num = cells.size();
  if (test_nb_cells && 2 == procs) {
    if (0 == rank) {
      CHECK_EQUAL(321, my_cells_num); // Gather set cells excluded
      if (mixed_elements)
        CHECK_EQUAL((size_t)2, cells.psize()); // Gather set cells excluded
      else
        CHECK_EQUAL((size_t)1, cells.psize()); // Gather set cells excluded
    }
    else if (1 == rank) {
      CHECK_EQUAL(321, my_cells_num); // Not owned cells excluded
      CHECK_EQUAL((size_t)1, cells.psize()); // Not owned cells excluded
    }
  }

   std::cout << "proc: " << rank << " cells:" << my_cells_num << "\n";

   int total_cells;
   MPI_Reduce(&my_cells_num, &total_cells, 1, MPI_INTEGER, MPI_SUM, 0, pcomm->proc_config().proc_comm());
   if (0 == rank) {
     std::cout << "total cells: " << total_cells << "\n";
     if (test_nb_cells)
       CHECK_EQUAL(total_cells, num_cells);
   }

  std::string write_options("PARALLEL=WRITE_PART;");
  mb.write_file("test_mpas.h5m", NULL, write_options.c_str());
}

void test_multiple_loads_of_same_file()
{
  Core moab;
  Interface& mb = moab;
  EntityHandle file_set;
  ErrorCode rval;
  rval = mb.create_meshset(MESHSET_SET, file_set);
  CHECK_ERR(rval);

  // Read first only header information, no mesh, no variable
  std::string opts("PARALLEL=READ_PART;PARTITION;NOMESH;VARIABLE=;PARTITION_METHOD=TRIVIAL");
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Create mesh, no variable
  opts="PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION_METHOD=TRIVIAL;VARIABLE=";
  // Create gather set in processor 1
  opts += std::string(";GATHER_SET=1");
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Read variable ke at timestep 0, no mesh
  opts = "PARALLEL=READ_PART;PARTITION;PARTITION_METHOD=TRIVIAL;NOMESH;VARIABLE=ke;TIMESTEP=0";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  int procs = pcomm->proc_config().proc_size();
  int rank = pcomm->proc_config().proc_rank();

  // Make check runs this test in two processors
  if (2 == procs) {
    Range verts;
    rval = mb.get_entities_by_type(0, MBVERTEX, verts);
    CHECK_ERR(rval);

    Range edges;
    rval = mb.get_entities_by_type(0, MBEDGE, edges);
    CHECK_ERR(rval);

    Range cells;
    rval = mb.get_entities_by_type(0, MBPOLYGON, cells);
    CHECK_ERR(rval);

    int my_verts_num = verts.size();
    int my_edges_num = edges.size();
    int my_cells_num = cells.size();

    if (0 == rank) {
      CHECK_EQUAL(1120, my_verts_num);
      CHECK_EQUAL(1438, my_edges_num);
      CHECK_EQUAL(321, my_cells_num);
      CHECK_EQUAL((size_t)2, cells.psize());

      const double eps = 1e-20;
      double val[2];

      // Check tag for cell variable ke at timestep 0
      Tag ke_tag0;
      rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
      CHECK_ERR(rval);

      // Check ke0 tag values on first pentagon and first hexagon
      EntityHandle cell_ents[] = {cells[0], cells[12]};
      rval = mb.tag_get_data(ke_tag0, &cell_ents[0], 2, val);
      CHECK_REAL_EQUAL(15.001, val[0], eps);
      CHECK_REAL_EQUAL(16.013, val[1], eps);
    }
    else if (1 == rank) {
      CHECK_EQUAL(2402, my_verts_num); // Gather set vertices included; Not owned vertices included
      CHECK_EQUAL(3364, my_edges_num); // Gather set edges included; Not owned edges included
      CHECK_EQUAL(963, my_cells_num); // Gather set cells included; Not owned cells included
      CHECK_EQUAL((size_t)3, cells.psize()); // Gather set cells included; Not owned cells included

      // Get gather set
      EntityHandle gather_set;
      ReadUtilIface* readUtilIface;
      rval = mb.query_interface(readUtilIface);
      CHECK_ERR(rval);
      rval = readUtilIface->get_gather_set(gather_set);
      CHECK_ERR(rval);

      // Get gather set entities
      Range gather_ents;
      rval = mb.get_entities_by_handle(gather_set, gather_ents);
      CHECK_ERR(rval);

      // Remove gather set vertices in processor 1
      verts = subtract(verts, gather_ents);
      my_verts_num = verts.size();
      CHECK_EQUAL(1122, my_verts_num); // Gather set vertices excluded; Not owned vertices included

      // Remove gather set edges in processor 1
      edges = subtract(edges, gather_ents);
      my_edges_num = edges.size();
      CHECK_EQUAL(1444, my_edges_num); // Gather set edges excluded; Not owned edges included

      // Remove gather set cells in processor 1
      cells = subtract(cells, gather_ents);
      my_cells_num = cells.size();
      CHECK_EQUAL(321, my_cells_num); // Gather set cells excluded; Not owned cells included
      CHECK_EQUAL((size_t)1, cells.psize()); // Gather set cells excluded; Not owned cells included
    }
  }
}

void test_multiple_loads_of_same_file_no_mixed_elements()
{
  Core moab;
  Interface& mb = moab;
  EntityHandle file_set;
  ErrorCode rval;
  rval = mb.create_meshset(MESHSET_SET, file_set);
  CHECK_ERR(rval);

  // Read first only header information, no mesh, no variable
  std::string opts("PARALLEL=READ_PART;PARTITION;NOMESH;VARIABLE=;PARTITION_METHOD=TRIVIAL;NO_MIXED_ELEMENTS");
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Create mesh, no variable
  opts="PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION_METHOD=TRIVIAL;NO_MIXED_ELEMENTS;VARIABLE=";
  // Create gather set in processor 1
  opts += std::string(";GATHER_SET=1");
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Read variable ke at timestep 0, no mesh
  opts = "PARALLEL=READ_PART;PARTITION;PARTITION_METHOD=TRIVIAL;NO_MIXED_ELEMENTS;NOMESH;VARIABLE=ke;TIMESTEP=0";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  int procs = pcomm->proc_config().proc_size();
  int rank = pcomm->proc_config().proc_rank();

  // Make check runs this test in two processors
  if (2 == procs) {
    Range verts;
    rval = mb.get_entities_by_type(0, MBVERTEX, verts);
    CHECK_ERR(rval);

    Range edges;
    rval = mb.get_entities_by_type(0, MBEDGE, edges);
    CHECK_ERR(rval);

    Range cells;
    rval = mb.get_entities_by_type(0, MBPOLYGON, cells);
    CHECK_ERR(rval);

    int my_verts_num = verts.size();
    int my_edges_num = edges.size();
    int my_cells_num = cells.size();

    if (0 == rank) {
      CHECK_EQUAL(1120, my_verts_num);
      CHECK_EQUAL(1438, my_edges_num);
      CHECK_EQUAL(321, my_cells_num);
      CHECK_EQUAL((size_t)1, cells.psize());

      const double eps = 1e-20;
      double val[2];

      // Check tag for cell variable ke at timestep 0
      Tag ke_tag0;
      rval = mb.tag_get_handle("ke0", 1, MB_TYPE_DOUBLE, ke_tag0);
      CHECK_ERR(rval);

      // Check ke0 tag values on first pentagon and first hexagon
      EntityHandle cell_ents[] = {cells[0], cells[12]};
      rval = mb.tag_get_data(ke_tag0, &cell_ents[0], 2, val);
      CHECK_REAL_EQUAL(15.001, val[0], eps);
      CHECK_REAL_EQUAL(16.013, val[1], eps);
    }
    else if (1 == rank) {
      CHECK_EQUAL(2402, my_verts_num); // Gather set vertices included; Not owned vertices included
      CHECK_EQUAL(3364, my_edges_num); // Gather set edges included; Not owned edges included
      CHECK_EQUAL(963, my_cells_num); // Gather set cells included; Not owned cells included
      CHECK_EQUAL((size_t)1, cells.psize()); // Gather set cells included; Not owned cells included

      // Get gather set
      EntityHandle gather_set;
      ReadUtilIface* readUtilIface;
      rval = mb.query_interface(readUtilIface);
      CHECK_ERR(rval);
      rval = readUtilIface->get_gather_set(gather_set);
      CHECK_ERR(rval);

      // Get gather set entities
      Range gather_ents;
      rval = mb.get_entities_by_handle(gather_set, gather_ents);
      CHECK_ERR(rval);

      // Remove gather set vertices in processor 1
      verts = subtract(verts, gather_ents);
      my_verts_num = verts.size();
      CHECK_EQUAL(1122, my_verts_num); // Gather set vertices excluded; Not owned vertices included

      // Remove gather set edges in processor 1
      edges = subtract(edges, gather_ents);
      my_edges_num = edges.size();
      CHECK_EQUAL(1444, my_edges_num); // Gather set edges excluded; Not owned edges included

      // Remove gather set cells in processor 1
      cells = subtract(cells, gather_ents);
      my_cells_num = cells.size();
      CHECK_EQUAL(321, my_cells_num); // Gather set cells excluded; Not owned cells included
      CHECK_EQUAL((size_t)1, cells.psize()); // Gather set cells excluded; Not owned cells included
    }
  }
}
