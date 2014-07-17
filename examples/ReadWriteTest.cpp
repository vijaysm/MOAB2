/** @example ReadWriteTest.cpp \n
 * \brief Read mesh into MOAB and write some back \n
 *
 * <b>To run</b>: mpiexec -np 4 ReadWriteTest [input] [output] -O <read_opts> -o <write_opts>\n
 *
 * used for stress test of reader/writer
 *  report times to read and write
 *
 *  example ReadWriteTest ../MeshFiles/io/fv3x46x72.t.3.nc out.nc  \
 *  -O PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;PARALLEL_RESOLVE_SHARED_ENTS;VARIABLE=T,U;  \
 *  -o PARALLEL=WRITE_PART;VARIABLE=T,U
 */

#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab/Core.hpp"
#include <iostream>
#include <time.h>

using namespace moab;
using namespace std;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  MBErrorHandler_Init();

  string options;

  // Need option handling here for input filename
  if (argc < 3)
   return 1;

  char* input_file = argv[1];
  char* output_file = argv[2];
  char* read_opts = NULL;
  char* write_opts = NULL; // Tags to write, separated by commas; it is the name of the tag

  if (argc > 3) {
    int index = 3;
    while (index < argc) {
      if (!strcmp(argv[index], "-O")) // This is for reading options, optional
        read_opts = argv[++index];
      if (!strcmp(argv[index], "-o"))
        write_opts = argv[++index];
      index++;
    }
  }

  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  // Get the ParallelComm instance
  ParallelComm* pcomm = new ParallelComm(mb, MPI_COMM_WORLD);
  int nprocs = pcomm->proc_config().proc_size();
  int rank = pcomm->proc_config().proc_rank();

  EntityHandle set;
  ErrorCode rval = mb->create_meshset(MESHSET_SET, set);CHK_ERR(rval);

  clock_t tt = clock();

  if (0 == rank)
    cout << "Reading file " << input_file << "\n with options: " << read_opts
         << "\n on " << nprocs << " processors\n";

  // Read the file with the specified options
  rval = mb->load_file(input_file, &set, read_opts);CHK_ERR(rval);

  if (0 == rank) {
    cout << "Time: " << (clock() - tt) / (double) CLOCKS_PER_SEC << " seconds" << endl;
    tt = clock();
  }

  // Write back the file with the specified options
  rval = mb->write_file(output_file, 0, write_opts, &set, 1);CHK_ERR(rval);

  if (0 == rank) {
    cout << "Writing file " << output_file << "\n with options: " << write_opts << endl;
    cout << "Time:  " << (clock() - tt) / (double) CLOCKS_PER_SEC << " seconds" << endl;
    tt = clock();
  }

  delete mb;

  MBErrorHandler_Finalize();

  MPI_Finalize();

  return 0;
}
