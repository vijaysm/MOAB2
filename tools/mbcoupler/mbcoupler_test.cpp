#include "MBParallelComm.hpp"
#include "MBCore.hpp"
#include "FileOptions.hpp"
#include "ReadParallel.hpp"
#include "mpi.h"
#include <iostream>
#include <sstream>
#include <assert.h>

#define RRA(a) if (MB_SUCCESS != result) {\
      std::string tmp_str; mbImpl->get_last_error(tmp_str);\
      tmp_str.append("\n"); tmp_str.append(a);\
      dynamic_cast<MBCore*>(mbImpl)->get_error_handler()->set_last_error(tmp_str.c_str()); \
      return result;}

MBErrorCode get_file_options(int argc, char **argv, 
                             std::vector<const char *> &filenames,
                             std::string &opts);

int main(int argc, char **argv) 
{
    // need to init MPI first, to tell how many procs and rank
  int err = MPI_Init(&argc, &argv);

  if (argc < 3) {
    std::cerr << "Usage: ";
    std::cerr << argv[0] << " <nfiles> <fname1> ... <fnamen> [tag_name] [tag_val] [distrib] [with_ghosts]" << std::endl;
    std::cerr << "nfiles        : number of mesh files" << std::endl;
    std::cerr << "fname1..fnamen: mesh files" << std::endl;
    std::cerr << "tag_name      : name of tag used to define partitions [MATERIAL_SET]" << std::endl;
    std::cerr << "tag_val       : tag values denoting partition sets [--]" << std::endl;
    std::cerr << "distrib       : if non-zero, distribute the partition sets with tag_val round-robin" << std::endl;
    std::cerr << "with_ghosts   : if non-zero, after initializing in parallel, also exchange one layer of ghost elements" << std::endl;

    err = MPI_Finalize();
    
    return 1;
  }
  
  int nprocs, rank;
  err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // start time
  double stime, rtime, setime, dtime, ltime;
  if (0 == rank) stime = MPI_Wtime();

    // create MOAB instance based on that
  MBInterface *mbImpl = new MBCore(rank, nprocs);
  if (NULL == mbImpl) return 1;
  
  MBErrorCode result = MB_SUCCESS;

  std::vector<const char *> filenames;
  std::string opts;
  result = get_file_options(argc, argv, filenames, opts);
  

    // read in mesh(es)
  ReadParallel rp(mbImpl);
  MBEntityHandle file_set = 0;
  result = rp.load_file(&filenames[0], filenames.size(), file_set, 
                        FileOptions(opts.c_str()), NULL, 0);
  if (MB_SUCCESS == result) {
    std::cout << "Success." << std::endl;
    err = MPI_Finalize();
    return 0;
  }
  else {
    std::string tmp_str;
    std::cout << "Failure; message:" << std::endl;
    std::cout << mbImpl->get_last_error(tmp_str) << std::endl;
    return 1;
  }
}

MBErrorCode get_file_options(int argc, char **argv, 
                             std::vector<const char *> &filenames,
                             std::string &opts) 
{
  int npos = 1;
  int nfiles = atoi(argv[npos++]);
  
    // get mesh filenames
  filenames.resize(nfiles);
  for (int i = 0; i < nfiles; i++) filenames[i] = argv[npos++];
  
    // get partition information
  const char *tag_name = "MATERIAL_SET";
  int tag_val = -1;
  int distrib = 1;
  int with_ghosts = 0;
  if (npos < argc) tag_name = argv[npos++];
  if (npos < argc) tag_val = strtol(argv[npos++], NULL, 0);
  if (npos < argc) distrib = strtol(argv[npos++], NULL, 0);
  if (npos < argc) with_ghosts = strtol(argv[npos++], NULL, 0);

  std::ostringstream options;
  options << "PARALLEL=BCAST_DELETE;PARTITION=" << tag_name;
  
  if (-1 != tag_val)
    options << ";PARTITION_VAL=" << tag_val;

  if (1 == distrib)
    options << ";PARTITION_DISTRIBUTE";

  options << ";PARALLEL_RESOLVE_SHARED_ENTS";

  if (1 == with_ghosts)
    options << ";PARALLEL_GHOSTS=3.0.1";

  options << ";CPUTIME";
    
  opts = options.str();

  return MB_SUCCESS;
}
