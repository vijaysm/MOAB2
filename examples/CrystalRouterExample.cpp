/*
 * This example will show one of the building blocks of parallel infrastructure in MOAB
 * More exactly, if we have some homogeneous data to communicate from each processor to a list of other
 * processors, how do we do it?
 *
 * introduce the TupleList and crystal router to MOAB users.
 *
 * This technology is used in resolving shared vertices / sets between partitions
 * It is used in the mbcoupler for sending data (target points) to the proper processor, and communicate
 *   back the results.
 * Also, it is used to communicate departure mesh for intersection in parallel
 *
 *  It is a way of doing  MPI_gatheralltoallv(), when the communication matrix is sparse
 *
 *  It is assumed that every proc needs to communicate only with a few of the other processors.
 *  If every processor needs to communicate with all other, then we will have to use paired isend and irecv, the
 *  communication matrix is full
 *
 *  the example needs to be launched in parallel.
 *  Every proc will build a list of tuples, that will be send to a few procs;
 *
 *  every proc will send 1 tuple, to proc rank + 1 and rank + rank*(size-1)+2 , with value
 *    10000 * send + 100* rank
 *
 *  at the receive, we verify we received
 *    10000 * rank + 100 * from
 *
 *    For some reportrank we also print the tuples.
 *
 *  after routing, we will see if we received, as expected. Should run on at least 2 processors.
 *
 * Note: We do not need a moab instance for this example
 *
 */

/** @example CrystalRouterExample.cpp \n
 * \brief generalized gather scatter using tuples \n
 * <b>To run</b>: mpiexec -np <n> CrystalRouterExample [reportrank] \n
 *
 */
//
#include "moab/ProcConfig.hpp"
#include "moab/TupleList.hpp"
#include <iostream>

using namespace moab;
using namespace std;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int reportrank = 1;
  if (argc>1)
    reportrank = atoi(argv[1]);
  ProcConfig pc(MPI_COMM_WORLD);
  int size = pc.proc_size();
  int rank = pc.proc_rank();

  if (reportrank==rank)
  {
    std::cout << " there are " << size << " procs in example\n";
  }
  // send some data from proc i to i+n/2, also to i +n/2+1 modulo n, wher en is num procs

  gs_data::crystal_data *cd = pc.crystal_router();

  TupleList tl;

  // at most 100 to send
  // we do a preallocate with this; some tuples on some processors might need more memory, to be able
  // to grow locally; 100 is a very large number for this example, considering that each task sends only
  // 2 tuples. Some tasks might receive more tuples though, and in the process, some might grow more than
  // others. By doing these logP sends/receives, we do not grow local memory too much.
  tl.initialize(1, 1, 0, 1, 100);
  tl.enableWriteAccess();
  // form 2 tuples, send to rank+1 and rank+2 (mod size)
  unsigned int n = tl.get_n();
  int sendTo = rank+1;
  sendTo = sendTo%size;
  long intToSend = 100*rank + 10000*sendTo;
  tl.vi_wr[n]= sendTo;
  tl.vl_wr[n]= intToSend;
  tl.vr_wr[n]= 100.*rank;
  tl.inc_n();

  n = tl.get_n();
  sendTo = rank+(rank+1)*rank+2;// just some number relatively different from rank
  sendTo = sendTo%size;
  intToSend = 100*rank + 10000*sendTo;
  tl.vi_wr[n]= sendTo;
  tl.vl_wr[n]= intToSend;
  tl.vr_wr[n]= 1000.*rank;
  tl.inc_n();

  if (reportrank==rank)
  {
    std::cout << "rank " << rank << "\n";
    tl.print(" before sending");
  }

  // all communication happens here:
  ErrorCode rval = cd->gs_transfer(1,tl,0);

  if (MB_SUCCESS!= rval)
  {
    std::cout << "error in tuple transfer\n";
  }

  if (reportrank==rank)
  {
    std::cout << "rank " << rank << "\n";
    tl.print(" after transfer");
  }
  // check that all tuples received have the form 10000* rank + 100*from
  unsigned int received = tl.get_n();
  for (int i=0; i<received; i++)
  {
    int from = tl.vi_rd[i];
    long valrec = tl.vl_rd[i];
    int remainder = valrec -10000*rank -100*from;
    if (remainder != 0 )
      std::cout << " error: tuple " << i << " received at proc rank " << rank << " from proc " << from << " has value " <<
         valrec << " remainder " <<  remainder << "\n";
  }

  MPI_Finalize();

  return 0;
}
