#include "moab/ParallelComm.hpp"
#include "moab/Core.hpp"
#include "moab_mpi.h"
#include "moab/ParallelMergeMesh.hpp"
#include <iostream>


#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

using namespace moab;

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int nproc, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (nproc!=2)
  {
    std::cout <<" run this test on 2 tasks\n";
    MPI_Finalize();
    return 1;
  }
  const char* filename0 = STRINGIFY(MESHDIR) "/brick1.vtk";
  const char* filename1 = STRINGIFY(MESHDIR) "/brick2.vtk";
  moab::Core *mb = new moab::Core();
  moab::ParallelComm *pc = new moab::ParallelComm(mb, MPI_COMM_WORLD);
  ErrorCode rval = MB_SUCCESS;
  if (0==rank)
    rval = mb->load_file(filename0);
  else
    rval = mb->load_file(filename1);

  if (rval!=MB_SUCCESS)
  {
    std::cout << "fail to load file\n";
    MPI_Finalize();
    return 1;
  }
  ParallelMergeMesh pm(pc, 0.001);
  rval = pm.merge();
  if (rval!=MB_SUCCESS)
  {
    std::cout << "fail to merge in parallel \n";
    MPI_Finalize();
    return 1;
  }
  rval = mb->write_file("testpm.h5m", 0, "PARALLEL=WRITE_PART");
  if (rval!=MB_SUCCESS)
  {
    std::cout << "fail to write output file \n";
    MPI_Finalize();
    return 1;
  }
  MPI_Finalize();
  return 0;
}
