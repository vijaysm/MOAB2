// $Id$
//
//  Simple C++ example of Zoltan library
//  copied and modified for use with MOAB, from Lee Ann Fisk's
//  Zoltan/examples/CPP/zCPPExample1.cpp file by Vitus Leung
//
//  MPICPP - Define this if your C++ bindings for MPI work.
//  NAMESPACES_OK - Define this if your system uses namespaces.
//

#include "MBZoltan.hpp"
#include "moab/Core.hpp"

#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <string>

#define RR                                                            \
    {if (MB_SUCCESS != result) {                                      \
        std::string last_error;                                       \
        mbImpl->get_last_error(last_error);                           \
        std::cerr << "MOAB returned with error: " << std::endl;       \
        std::cerr << last_error << std::endl;                         \
        return 1;                                                     \
      }                                                               \
    }

int main(int argc, char *argv[])
{
  if (argc < 4) {
    std::cout << "Usage: mpirun -np <nprocs> " << argv[0]
              << " <# partitions> "
              << " <mesh_file> <write_out(y/n)> <output_mesh_file> "
              << "[write_as_sets(1=default/0)]"
              << "[write_as_tags(1/0=default)]" 
              << "[partition_dim(default=3)]"
              << "[<method(RCB/RIB/HSFC/Hypergraph(PHG)/PARMETIS/OCTPART)>] "
              << "[<parmetis_method>/<oct_method>]"
              << "[imbalance_tolerance(default=1.10, use 1.03 to get non-empty parts]"
              << std::endl;
    return 1;
  }

  bool write_output = false;
  if (argv[3][0] == 'y') write_output = true;
  else if (argv[3][0] == 'n') write_output = false;
  else {
    std::cout << "Argument 2 must be 'y' or 'n', not '" << argv[3][0] << "'" << std::endl;
    return 1;
  }

    // need to init MPI first, to tell how many procs and rank
  int err = MPI_Init(&argc, &argv);

  int nprocs, rank;
  err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // create MOAB instance based on that
  moab::Interface *mbImpl = new moab::Core();
  if (NULL == mbImpl) return 1;
  
  moab::ErrorCode result = mbImpl->load_mesh(argv[2]); RR;
    
  
  MBZoltan *mbz = new MBZoltan(mbImpl, false, argc, argv);

  int as_sets = 1, as_tags = 0, part_dim = 3;
  if (argc > 5) as_sets = atoi(argv[5]);
  if (argc > 6) as_tags = atoi(argv[6]);

  if (argc > 7) part_dim = atoi(argv[7]);
  
  const char *other_method = NULL, *method = NULL;
  if (argc > 8) method = argv[8];
  if (argc > 9) other_method = argv[9];

  double imbal_tol = 0.0;
  if (argc > 10) {
    imbal_tol = atof(argv[10]);
    if (strcmp(method, "Hypergraph") && strcmp(method, "PHG"))
      std::cerr << "Warning: imbalance tolerance can only be used with Hypergraph or PHG methods; ignoring."
                << std::endl;
  }
  
  int nparts = atoi(argv[1]);
  
  result = mbz->partition_mesh(nparts, method, other_method,
                               as_sets, as_tags, part_dim, imbal_tol); RR;
  
  if (write_output) {
    result = mbImpl->write_mesh(argv[4]); RR;
    std::cout << "Wrote partitioned mesh to " << argv[4] << std::endl;
  }

  delete mbz;
  delete mbImpl;

  err = MPI_Finalize();

  return 0;
}
