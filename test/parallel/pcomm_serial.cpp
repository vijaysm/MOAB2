#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "MBTagConventions.hpp"
#include "moab/Core.hpp"
#include "moab/FileOptions.hpp"
#include "ReadParallel.hpp"
#include "TestUtil.hpp"
#include <vector>

using namespace moab;

void print_usage(char *argv) 
{
  std::cout << "Usage: " << argv << " nprocs filename" << std::endl;
}
  
int main( int argc, char* argv[] )
{
#ifdef MOAB_HAVE_MPI
  MPI_Init( &argc, &argv );
#else
# define MPI_COMM_WORLD 0
#endif

  if (1 < argc && !strcmp(argv[1], "-h")) {
    print_usage(argv[0]);
    return 0;
  }
  
  int nprocs = 2;
  std::string ptag_name("GEOM_DIMENSION");
  std::vector<int> partition_tag_vals;
  std::string filename = TestDir + "/ptest.cub";
  if (argc > 1)
    nprocs = atoi(argv[1]);
  if (argc > 2)
    filename = std::string(argv[2]);
  if (argc > 3) {
    ptag_name = argv[3];
    if (argc > 4) partition_tag_vals.push_back(atoi(argv[4]));
  }
  else partition_tag_vals.push_back(3);

  if (0 == nprocs) {
    print_usage(argv[0]); 
    return 1;
  }
  
  ErrorCode rval;
  Core *moab = new Core[nprocs]();
  std::vector<ParallelComm *> pc(nprocs);
  for (int i = 0; i < nprocs; i++) {
    pc[i] = new ParallelComm(&moab[i], MPI_COMM_WORLD);
    pc[i]->set_rank(i);
    pc[i]->set_size(nprocs);
  }

  std::vector<int> pa_vec;
  pa_vec.push_back(ReadParallel::PA_READ);
  pa_vec.push_back(ReadParallel::PA_GET_FILESET_ENTS);
  pa_vec.push_back(ReadParallel::PA_DELETE_NONLOCAL);
  bool partition_distrib = true;

  FileOptions fopts(NULL);

  const char* fnames = filename.c_str();
  for (int i = 0; i < nprocs; i++) {
    ReadParallel rp(moab+i, pc[i]);
    rval = rp.load_file(&fnames, 1, 0, ReadParallel::POPT_READ_DELETE,
                        ptag_name, 
                        partition_tag_vals, partition_distrib, false, pa_vec, 
                        fopts, NULL, NULL, i, false, -1, -1, -1, -1, 0, 0);
    CHECK_ERR(rval);
  }
  
  rval = ParallelComm::resolve_shared_ents(&pc[0], nprocs, 0, 3);
  CHECK_ERR(rval);

    // exchange interface cells
  rval = ParallelComm::exchange_ghost_cells(&pc[0], nprocs, -1, -1, 0, 0, true);
  CHECK_ERR(rval);
  
    // now 1 layer of hex ghosts
  rval = ParallelComm::exchange_ghost_cells(&pc[0], nprocs, 3, 2, 1, 0, true);
  CHECK_ERR(rval);

    // now 1 layer of hex ghosts with face/edges
  rval = ParallelComm::exchange_ghost_cells(&pc[0], nprocs, 3, 2, 1, 3, true);
  CHECK_ERR(rval);

  for (int i = 0; i < nprocs; i++)
    delete pc[i];

  delete [] moab;
  
#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
