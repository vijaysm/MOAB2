#include "MBParallelComm.hpp"
#include "MBParallelConventions.h"
#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "FileOptions.hpp"
#include "ReadParallel.hpp"
#include "TestUtil.hpp"
#include <vector>

int main( int argc, char* argv[] )
{
#ifdef USE_MPI
  MPI_Init( &argc, &argv );
#endif

  if (1 < argc && !strcmp(argv[1], "-h")) {
    std::cout << "Usage: " << argv[0] << " nprocs filename" << std::endl;
    return 0;
  }
  int nprocs = 2;
  std::string ptag_name("GEOM_DIMENSION");
  std::vector<int> partition_tag_vals;
#ifdef SRCDIR
  const char *fnames[] = {STRINGIFY(SRCDIR) "/ptest.cub"};
#else
  const char *fnames[] = {"./ptest.cub"};
#endif
  if (argc > 1)
    nprocs = atoi(argv[1]);
  if (argc > 2)
    fnames[0] = argv[2];
  if (argc > 3) {
    ptag_name = argv[3];
    if (argc > 4) partition_tag_vals.push_back(atoi(argv[4]));
  }
  else partition_tag_vals.push_back(3);
  
  MBErrorCode rval;
  MBCore *moab = new MBCore[nprocs]();
  std::vector<MBParallelComm *> pc(nprocs);
  for (int i = 0; i < nprocs; i++) {
    pc[i] = new MBParallelComm(&moab[i]);
    pc[i]->set_rank(i);
    pc[i]->set_size(nprocs);
  }

  std::vector<int> pa_vec;
  pa_vec.push_back(ReadParallel::PA_READ);
  pa_vec.push_back(ReadParallel::PA_GET_FILESET_ENTS);
  pa_vec.push_back(ReadParallel::PA_DELETE_NONLOCAL);
  bool partition_distrib = false;

  partition_distrib = true;
  
    //std::string ptag_name("MATERIAL_SET");
    //partition_distrib = true;
  
  FileOptions fopts(NULL);
  
  for (int i = 0; i < nprocs; i++) {
    ReadParallel rp(moab+i, pc[i]);
    MBEntityHandle tmp_set = 0;
    rval = rp.load_file(fnames, 1, tmp_set, ReadParallel::POPT_READ_DELETE,
                        ptag_name, 
                        partition_tag_vals, partition_distrib, false, pa_vec, 
                        fopts, NULL, 0, NULL, i, false, -1, -1, -1, -1, 0);
    CHECK_ERR(rval);
  }
  
  rval = MBParallelComm::resolve_shared_ents(&pc[0], nprocs, 3);
  CHECK_ERR(rval);

    // exchange interface cells
  rval = MBParallelComm::exchange_ghost_cells(&pc[0], nprocs, -1, -1, 0, true);
  CHECK_ERR(rval);
  
    // now 1 layer of hex ghosts
  rval = MBParallelComm::exchange_ghost_cells(&pc[0], nprocs, 3, 0, 1, true);
  CHECK_ERR(rval);

  for (int i = 0; i < nprocs; i++)
    delete pc[i];

  delete [] moab;
  
#ifdef USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
