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
void test_read_parallel(int num_verts, bool test_nb_nodes, int num_edges, bool test_nb_edges);

void test_multiple_loads_of_same_file();

std::string partition_method;

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int result = 0;

  result += RUN_TEST(test_read_parallel_mpas_trivial);
  result += RUN_TEST(test_multiple_loads_of_same_file);

  MPI_Finalize();
  return result;
}

void test_read_parallel_mpas_trivial()
{
  partition_method = std::string(";PARTITION_METHOD=TRIVIAL;PARALLEL_RESOLVE_SHARED_ENTS");
  test_read_parallel(1280, true, 1920, true);
}
  
void test_read_parallel(int num_verts, bool test_nb_nodes, int num_edges, bool test_nb_edges)
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
  if (0 == rank)
  {
    std::cout << "total vertices: " << total_verts << "\n";
    if (test_nb_nodes)
      CHECK_EQUAL(total_verts, num_verts);
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

    // Remove gather set edges in processor 0
    edges = subtract(edges, gather_ents);
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
  if (0 == rank)
  {
    std::cout << "total edges: " << total_edges << "\n";
    if (test_nb_edges)
      CHECK_EQUAL(total_edges, num_edges);
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

    int my_verts_num = verts.size();
    if (0 == rank)
      CHECK_EQUAL(1120, my_verts_num);
    else if (1 == rank)
      CHECK_EQUAL(2402, my_verts_num); // Gather set vertices included; Not owned vertices included

    if (1 == rank) {
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
    }

    my_verts_num = verts.size();
    if (0 == rank)
      CHECK_EQUAL(1120, my_verts_num);
    else if (1 == rank)
      CHECK_EQUAL(1122, my_verts_num); // Gather set vertices excluded; Not owned vertices included
  }
}
