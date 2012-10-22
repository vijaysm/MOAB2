#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/ProgOptions.hpp"
#include "MBParallelConventions.h"

using namespace moab;

#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/homme26x3458.t.3.nc";
#else
static const char example[] = "/io/camEul26x48x96.t3.nc";
#endif

void test_read_parallel_ucd_nodal();
void test_read_parallel(int num_verts);

std::string partition_method;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int result = 0;
  
  result += RUN_TEST(test_read_parallel_ucd_nodal);
  
  MPI_Finalize();
  return result;
}


void test_read_parallel_ucd_nodal()
{
  partition_method = std::string(";PARTITION_METHOD=NODAL_PARTITION");
  test_read_parallel(3458);
}
  
void test_read_parallel(int num_verts)
{
  Core moab;
  Interface& mb = moab;
  EntityHandle file_set;
  ErrorCode rval;
  rval = mb.create_meshset(MESHSET_SET, file_set);
  CHECK_ERR(rval);

  std::string opt = std::string("PARALLEL=READ_PART;") +
      partition_method;
  rval = mb.load_file(example, &file_set, opt.c_str());
  CHECK_ERR(rval);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);

  rval = pcomm->check_all_shared_handles();
  CHECK_ERR(rval);

    // get the total # owned verts
  Range verts;
  rval = mb.get_entities_by_type(0, MBVERTEX, verts);
  CHECK_ERR(rval);
  rval = pcomm->filter_pstatus(verts, PSTATUS_NOT_OWNED, PSTATUS_NOT);
  CHECK_ERR(rval);
  int my_num = verts.size(), total_verts;
  std::cout<<"proc: " << pcomm->proc_config().proc_rank() << " verts:" << my_num << "\n";
  MPI_Reduce(&my_num, &total_verts, 1, MPI_INTEGER, MPI_SUM, 0, pcomm->proc_config().proc_comm());
  
  if (0 == pcomm->proc_config().proc_rank())
  {
    std::cout<<"total vertices: " << total_verts << "\n";
    CHECK_EQUAL(total_verts, num_verts);
  }
  /*std::string write_options("PARALLEL=WRITE_PART;");
  mb.write_file( "test.h5m", NULL, write_options.c_str() );*/

}

