/* \example HelloMoabPar HelloMoabPar.cpp
 * \brief Read mesh into MOAB in parallel
 * This example shows the simplest way of telling MOAB to read in parallel.
 *
 * This example:
 * 0. Initialize MPI and get the rank and number of processors
 * 1. Process arguments (file name and options for parallel read)
 * 2. Initialize MOAB
 * 3. Load a partitioned file in parallel;
 * 4. Report shared entities and their status
 * 5. Filter owned entities among shared ones on each processor
 * 6. Exchange ghost layers, and repeat the reports
 *
 * To run: mpiexec -np 4 HelloMoabPar
 *
 */

#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab/Core.hpp"
#include <iostream>

using namespace moab;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int nprocs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string filename;
  std::string options;
  if (3 != argc)
  {
    if (rank == 0)
    {
      std::cout << "Usage: " << argv[0] << " <filename>  <options (separated by;)>\n ";
    }
    /* this file has a partition with 4 parts */
    filename = "../MeshFiles/unittest/disk.h5m";
    /*  Options for reading
     *  - read in parallel
     *  - use PARALLEL_PARTITION tag
     *  - resolve shared entities after reading
    */
    options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";
  }
  else
  {
    filename = argv[1];
    options = argv[2];
  }
  if (rank == 0)
    std::cout << "reading file " << filename << "\n  with options:" << options <<
      "\n   on " << nprocs << " processors\n";

  // get MOAB instance and read the file with the specified options
  Interface *mbImpl = new Core;
  if (NULL == mbImpl) return 1;
  ErrorCode rval = mbImpl->load_file(filename.c_str(), 0, options.c_str());
  if (rval != MB_SUCCESS) return 1;

  // get the ParallelComm instance
  ParallelComm* pcomm = ParallelComm::get_pcomm(mbImpl, 0);
  MPI_Comm comm = pcomm->comm();
  if (0 == pcomm) return 1;

  // get all shared entities with other processors
  Range shared_ents;
  rval = pcomm->get_shared_entities(-1, // -1 means all other processors                                    Range &shared_ents,
      shared_ents);
  if (rval != MB_SUCCESS) return 1;
  /* Among shared entities, get those owned by the current processor
   * For this, use a filter operation;
   * Each shared entity is owned by exactly one processor;
   * An entity could be simply-shared (with exactly one other processor) or
   *  multi-shared.
   */
  Range owned_entities;
  rval = pcomm->filter_pstatus(shared_ents, // pass entities that we want to filter
      PSTATUS_NOT_OWNED, // status we are looking for
      PSTATUS_NOT, // operation applied ; so it will return owned entities (!not_owned = owned)
      -1, // this means all processors
      &owned_entities);
  if (rval != MB_SUCCESS) return 1;
  unsigned int nums[4]={0}; // to store the owned entities per dimension
  for (int i=0; i<3; i++)
  {
    nums[i]=(int)owned_entities.num_of_dimension(i);
  }
  int * rbuf;
  if (rank==0)
    rbuf = (int *)malloc(nprocs*4*sizeof(int));
  MPI_Gather( nums, 4, MPI_INT, rbuf, 4, MPI_INT, 0, comm);
  // print the stats gathered:
  if (rank == 0)
  {
    for (int i=0; i<nprocs; i++)
    {
      std::cout << " shared, owned entities on proc " << i << " :" << rbuf[4*i] << " verts, " <<
          rbuf[4*i+1] << " edges, " << rbuf[4*i+2] << " faces\n";
    }

  }

  /*
   * Now exchange 1 layer of ghost elements, using vertices as bridge
   *  we could have done this as part of reading process, by passing an extra option
   *    ";GHOSTS=2.0.1"
   */
  rval = pcomm->exchange_ghost_cells(2, // int ghost_dim,
                                     0, // int bridge_dim,
                                     1, //int num_layers,
                                     0, //int addl_ents,
                                     true); // bool store_remote_handles);
  if (rval != MB_SUCCESS) return 1;

  // repeat the reports, after ghost exchange
  shared_ents.clear();
  owned_entities.clear();
  rval = pcomm->get_shared_entities(-1, // -1 means all other processors                                    Range &shared_ents,
        shared_ents);
  if (rval != MB_SUCCESS) return 1;
  rval = pcomm->filter_pstatus(shared_ents,
        PSTATUS_NOT_OWNED,
        PSTATUS_NOT,
        -1,
        &owned_entities);
  if (rval != MB_SUCCESS)  return 1;

  // find out how many shared entities of each dimension are owned on this processor
  for (int i=0; i<3; i++)
    nums[i]=(int)owned_entities.num_of_dimension(i);

  if (rank==0)
    rbuf = (int *)malloc(nprocs*4*sizeof(int));

  // gather the statistics on processor 0
  MPI_Gather( nums, 4, MPI_INT, rbuf, 4, MPI_INT, 0, comm);
  if (rank == 0)
  {
    std::cout << " \n\n After ghosting: \n";
    for (int i=0; i<nprocs; i++)
    {
      std::cout << " shared, owned entities on proc " << i << " :" << rbuf[4*i] << " verts, " <<
          rbuf[4*i+1] << " edges, " << rbuf[4*i+2] << " faces\n";
    }
    free(rbuf);
  }
  MPI_Finalize();

  return 0;
}
