/** test of MBParallelComm functionality
 *
 * To run:
 *
 * mpirun -np <#procs> mbparallelcomm_test
 *
 */

#include "MBParallelComm.hpp"
#include "MBParallelConventions.h"
#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "mpi.h"
#include <iostream>

#define ERROR(a, b) {std::cerr << a << std::endl; return b;}

int main(int argc, char **argv) 
{
    // need to init MPI first, to tell how many procs and rank
  int err = MPI_Init(&argc, &argv);

  int nprocs, rank;
  err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // create MOAB instance based on that
  MBInterface *mbImpl = new MBCore(rank, nprocs);
  if (NULL == mbImpl) return 1;
  
  MBErrorCode result;

    // each interior proc has a vector of N+M vertices, sharing
    // M vertices each with lower- and upper-rank processors, except
    // procs on the end

    // get N, M from command line
  int N, M;
  if (argc < 3) {
    std::cerr << "No arguments passed; assuming N=10, M=2." << std::endl;
    N = 10;
    M = 2;
  }
  else {
    N = atoi(argv[1]);
    M = atoi(argv[2]);
  }
  
  int nverts = N + M;
  if (0 == rank) nverts = N;

    // create my vertices and give them the right global ids
  MBRange my_verts;
  std::vector<double> coords(3*(nverts));
  std::fill(coords.begin(), coords.end(), 0.0);
  result = mbImpl->create_vertices(&coords[0], nverts,
                                   my_verts);
  if (MB_SUCCESS != 0)
    ERROR("Failed to create vertices.", 1);
  
  std::vector<int> global_ids(N+M);
  for (int i = 0; i < nverts; i++)
    global_ids[i] = rank*N - (nverts-N) + i;
  
  int def_val = -1;
  MBTag gid_tag;
  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, 1, MB_TAG_DENSE,
                              MB_TYPE_INTEGER, gid_tag,
                              &def_val, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
    ERROR("Failed to create tag.", 1);
  
  result = mbImpl->tag_set_data(gid_tag, my_verts, &global_ids[0]);
  if (MB_SUCCESS != result) ERROR("Failed to set global_id tag.", 1);
  
    // now figure out what's shared
  MBParallelComm pcomm(mbImpl);
  result = pcomm.resolve_shared_ents(my_verts, 0);
  if (MB_SUCCESS != result) ERROR("Couldn't resolve shared entities.", 1);
  
    // check shared entities
  MBTag sharedproc_tag, sharedprocs_tag;
  result = mbImpl->tag_get_handle(PARALLEL_SHARED_PROC_TAG_NAME, 
                                  sharedproc_tag);
  if (MB_SUCCESS != result) ERROR("Shared processor tag not found.", 1);

  result = mbImpl->tag_get_handle(PARALLEL_SHARED_PROCS_TAG_NAME, 
                                  sharedprocs_tag);
  if (MB_SUCCESS != result) 
    ERROR("Shared processor*s* tag not found.", 1);
  
    // get the tag values
#define MAX_SHARING_PROCS 10
  std::vector<int> shared_proc_tags(MAX_SHARING_PROCS*my_verts.size());
  result = mbImpl->tag_get_data(sharedproc_tag, my_verts, 
                                &shared_proc_tags[0]);
  if (MB_SUCCESS != result) ERROR("Problem getting shared proc tag.", 1);

    // interior procs should have 2*M shared, bdy procs should have M shared
  int nshared = 0;
  for (unsigned int nv = 0; nv < my_verts.size(); nv++)
    if (shared_proc_tags[2*nv] > -1) nshared++;
  
  if ((rank == 0 || rank == nprocs-1) && nshared != (unsigned int) M) {
    std::cerr << "Didn't get correct number of shared vertices on "
              << "processor " << rank << std::endl;
    result = MB_FAILURE;
  }
  
  else if ((rank != 0 && rank != nprocs-1) && nshared != (unsigned int) 2*M) 
  {
    std::cerr << "Didn't get correct number of shared vertices on "
              << "processor " << rank << std::endl;
    result = MB_FAILURE;
  }

    // now check sharedprocs; shouldn't be any 
  MBErrorCode result2 = mbImpl->tag_get_data(sharedprocs_tag, my_verts, 
                                             &shared_proc_tags[0]);
  if (MB_SUCCESS == result2) {
    std::cerr << "Shoudn't get shared proc*s* tag, but did on proc "
              << rank << std::endl;
    result = MB_FAILURE;
  }

  err = MPI_Finalize();

  if (MB_SUCCESS == result)
    std::cerr << "Proc " << rank << ": Success." << std::endl;
    
  return (MB_SUCCESS == result ? 0 : 1);
}
