/** @example HelloParMOAB.cpp \n
 * \brief Read mesh into MOAB and resolve/exchange/report shared and ghosted entities \n
 * <b>To run</b>: mpiexec -np 4 HelloMoabPar [filename]\n
 *
 */

#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab/Core.hpp"
#include <iostream>

using namespace moab;
using namespace std;

string test_file_name = string(MESH_DIR) + string("/64bricks_512hex_256part.h5m");

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  string options;

    // need option handling here for input filename
  if (argc > 1){
    //user has input a mesh file
    test_file_name = argv[1];
  }  

  options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";

  // get MOAB instance and read the file with the specified options
  Interface *mb = new Core;
  if (NULL == mb) return 1;
  // get the ParallelComm instance
  ParallelComm* pcomm = new ParallelComm(mb, MPI_COMM_WORLD);
  int nprocs = pcomm->proc_config().proc_size(), rank = pcomm->proc_config().proc_rank();
  MPI_Comm comm = pcomm->proc_config().proc_comm();

  if (rank == 0)
    cout << "Reading file " << test_file_name << "\n  with options: " << options << endl
         << " on " << nprocs << " processors\n";

  ErrorCode rval = mb->load_file(test_file_name.c_str(), 0, options.c_str());
  if (rval != MB_SUCCESS) return 1;

  Range shared_ents;
    // get entities shared with all other processors
  rval = pcomm->get_shared_entities(-1, shared_ents);
  if (rval != MB_SUCCESS) return 1;

    // filter shared entities with not not_owned, which means owned
  Range owned_entities;
  rval = pcomm->filter_pstatus(shared_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &owned_entities);
  if (rval != MB_SUCCESS) return 1;

  unsigned int nums[4]={0}; // to store the owned entities per dimension
  for (int i=0; i<4; i++) nums[i]=(int)owned_entities.num_of_dimension(i);
  vector<int> rbuf(nprocs*4, 0);
  MPI_Gather( nums, 4, MPI_INT, &rbuf[0], 4, MPI_INT, 0, MPI_COMM_WORLD);
  // print the stats gathered:
  if (rank == 0) {
    for (int i=0; i<nprocs; i++)
      cout << " Shared, owned entities on proc " << i << ": " << rbuf[4*i] << " verts, " <<
          rbuf[4*i+1] << " edges, " << rbuf[4*i+2] << " faces, " << rbuf[4*i+3] << " elements" << endl;
  }

    // Now exchange 1 layer of ghost elements, using vertices as bridge
    // (we could have done this as part of reading process, using the PARALLEL_GHOSTS read option)
  rval = pcomm->exchange_ghost_cells(3, // int ghost_dim,
                                     0, // int bridge_dim,
                                     1, //int num_layers,
                                     0, //int addl_ents,
                                     true); // bool store_remote_handles);
  if (rval != MB_SUCCESS) return 1;

  // repeat the reports, after ghost exchange
  shared_ents.clear();
  owned_entities.clear();
  rval = pcomm->get_shared_entities(-1, shared_ents);
  if (rval != MB_SUCCESS) return 1;
  rval = pcomm->filter_pstatus(shared_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &owned_entities);
  if (rval != MB_SUCCESS)  return 1;

  // find out how many shared entities of each dimension are owned on this processor
  for (int i=0; i<4; i++)
    nums[i]=(int)owned_entities.num_of_dimension(i);

  // gather the statistics on processor 0
  MPI_Gather( nums, 4, MPI_INT, &rbuf[0], 4, MPI_INT, 0, comm);
  if (rank == 0)
  {
    cout << " \n\n After exchanging one ghost layer: \n";
    for (int i=0; i<nprocs; i++)
    {
      cout << " Shared, owned entities on proc " << i << ": " << rbuf[4*i] << " verts, " <<
          rbuf[4*i+1] << " edges, " << rbuf[4*i+2] << " faces, " << rbuf[4*i+3] << " elements" << endl;
    }
  }

  MPI_Finalize();

  return 0;
}
