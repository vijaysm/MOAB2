#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/ScdInterface.hpp"
#include "moab/ProgOptions.hpp"
#include "MBParallelConventions.h"

using namespace moab;

#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/cam18x40x48.t2.nc";
#else
static const char example[] = "/io/cam18x40x48.nc";
#endif

void test_read_parallel(int nverts);
void test_read_parallel_alljorkori();
void test_read_parallel_alljkbal();
void test_read_parallel_sqij();
void test_read_parallel_sqjk();
std::string partition_method;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int result = 0;
  
  result += RUN_TEST(test_read_parallel_alljorkori);
  result += RUN_TEST(test_read_parallel_alljkbal);
  result += RUN_TEST(test_read_parallel_sqij);
  result += RUN_TEST(test_read_parallel_sqjk);
  
  MPI_Finalize();
  return result;
}


void test_read_parallel_alljorkori() 
{
  partition_method = std::string(";PARTITION_METHOD=alljorkori");
  test_read_parallel(34560);
}
  
void test_read_parallel_alljkbal() 
{
  partition_method = std::string(";PARTITION_METHOD=alljkbal");
  test_read_parallel(34560);
}
  
void test_read_parallel_sqij() 
{
  partition_method = std::string(";PARTITION_METHOD=sqij");
  test_read_parallel(34560);
}
  
void test_read_parallel_sqjk() 
{
  partition_method = std::string(";PARTITION_METHOD=sqjk");
  test_read_parallel(34560);
}
  
void test_read_parallel(int num_verts)
{
  Core moab;
  Interface& mb = moab;
  EntityHandle file_set;
  ErrorCode rval;
  rval = mb.create_meshset(MESHSET_SET, file_set);
  CHECK_ERR(rval);

  std::string opt = std::string("PARALLEL=READ_PART;PARTITION=;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS") +
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
  MPI_Reduce(&my_num, &total_verts, 1, MPI_INTEGER, MPI_SUM, 0, pcomm->proc_config().proc_comm());
  
  if (0 == pcomm->proc_config().proc_rank()) CHECK_EQUAL(total_verts, num_verts);
}

