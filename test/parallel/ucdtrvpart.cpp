#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/ProgOptions.hpp"
#include "MBParallelConventions.h"
#include "moab/Util.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/homme26x3458.t.3.nc";
#endif

void test_read_parallel_ucd_trivial();
void test_read_parallel_ucd_trivial_spectral();
void test_read_parallel(int num_verts, bool test_nb_nodes);

void test_multiple_loads_of_same_file();

std::string partition_method;

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int result = 0;
  
  result += RUN_TEST(test_read_parallel_ucd_trivial);
  result += RUN_TEST(test_read_parallel_ucd_trivial_spectral);
  result += RUN_TEST(test_multiple_loads_of_same_file);
  
  MPI_Finalize();
  return result;
}

void test_read_parallel_ucd_trivial()
{
  // Disable spectral mesh for the time being, it is not ready yet
  //partition_method = std::string(";PARTITION_METHOD=TRIVIAL_PARTITION;SPECTRAL_MESH;PARALLEL_RESOLVE_SHARED_ENTS");
  partition_method = std::string(";PARTITION_METHOD=TRIVIAL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS");
  test_read_parallel(3458, true);
}
  
void test_read_parallel_ucd_trivial_spectral()
{
  partition_method = std::string(";PARTITION_METHOD=TRIVIAL_PARTITION;SPECTRAL_MESH;PARALLEL_RESOLVE_SHARED_ENTS");
  test_read_parallel(3458, false);
}

void test_read_parallel(int num_verts, bool test_nb_nodes)
{
  Core moab;
  Interface& mb = moab;
  EntityHandle file_set;
  ErrorCode rval;
  rval = mb.create_meshset(MESHSET_SET, file_set);
  CHECK_ERR(rval);

  std::string opt = std::string("PARALLEL=READ_PART") + partition_method;
  rval = mb.load_file(example, &file_set, opt.c_str());
  CHECK_ERR(rval);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);

  rval = pcomm->check_all_shared_handles();
  CHECK_ERR(rval);

  // Get the total # owned verts
  Range verts;
  rval = mb.get_entities_by_type(0, MBVERTEX, verts);
  CHECK_ERR(rval);
  rval = pcomm->filter_pstatus(verts, PSTATUS_NOT_OWNED, PSTATUS_NOT);
  CHECK_ERR(rval);

  int rank = pcomm->proc_config().proc_rank();
  if (0 == rank) {
    // Remove from verts the gather set ents
    EntityHandle gather_set;
    Range gth_ents;
    rval = Util::gather_set_entities(&mb, gather_set, gth_ents);
    CHECK_ERR(rval);
    verts = subtract(verts, gth_ents);
  }

  int my_num = verts.size(), total_verts;
  std::cout << "proc: " << rank << " verts:" << my_num << "\n";
  MPI_Reduce(&my_num, &total_verts, 1, MPI_INTEGER, MPI_SUM, 0, pcomm->proc_config().proc_comm());
  
  if (0 == rank) {
    std::cout << "total vertices: " << total_verts << "\n";
    if (test_nb_nodes)
      CHECK_EQUAL(total_verts, num_verts);
  }

  std::string write_options("PARALLEL=WRITE_PART;");
  mb.write_file("test.h5m", NULL, write_options.c_str());
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
  std::string opts("PARALLEL=READ_PART;PARTITION;NOMESH;VARIABLE=;PARTITION_METHOD=TRIVIAL_PARTITION");
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  opts="PARALLEL=READ_PART;PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION_METHOD=TRIVIAL_PARTITION;VARIABLE=";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  opts = "PARALLEL=READ_PART;PARTITION;PARTITION_METHOD=TRIVIAL_PARTITION;NOMESH;VARIABLE=T;TIMESTEP=0";
  rval = mb.load_file(example, &file_set, opts.c_str());
  CHECK_ERR(rval);

  // Check values of tag T0 (first level) at some strategically chosen places below
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  int procs = pcomm->proc_config().proc_size();

  // Make check runs this test in two processors
  if (2 == procs) {
    Range verts;
    rval = mb.get_entities_by_type(0, MBVERTEX, verts);
    CHECK_ERR(rval);

    int rank = pcomm->proc_config().proc_rank();
    if (0 == rank) {
      // Remove from verts the gather set ents
      EntityHandle gather_set;
      Range gth_ents;
      rval = Util::gather_set_entities(&mb, gather_set, gth_ents);
      CHECK_ERR(rval);
      verts = subtract(verts, gth_ents);
    }
    CHECK_EQUAL((size_t)1825, verts.size());

    Tag Ttag0;
    rval = mb.tag_get_handle("T0", 26, MB_TYPE_DOUBLE, Ttag0, MB_TAG_DENSE);
    CHECK_ERR(rval);

    int count;
    void* Tbuf;
    rval = mb.tag_iterate(Ttag0, verts.begin(), verts.end(), count, Tbuf);
    CHECK_ERR(rval);
    CHECK_EQUAL((size_t)count, verts.size());

    const double eps = 0.0001;
    double* data = (double*) Tbuf;

    if (0 == rank) {
      CHECK_REAL_EQUAL(233.1136, data[0 * 26], eps); // First vert
      CHECK_REAL_EQUAL(237.1977, data[912 * 26], eps); // Median vert
      CHECK_REAL_EQUAL(234.9711, data[1824 * 26], eps); // Last vert
    }
    else if (1 == rank) {
      CHECK_REAL_EQUAL(233.1136, data[0 * 26], eps); // First vert
      CHECK_REAL_EQUAL(231.0446, data[912 * 26], eps); // Median vert
      CHECK_REAL_EQUAL(234.0416, data[1824 * 26], eps); // Last vert
    }
  }
}
