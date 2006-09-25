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
#include "MBCore.hpp"

#include <stdlib.h>
#include <unistd.h>
#include <iostream>

#define RR if (MB_SUCCESS != result) return result

int main(int argc, char *argv[])
{
  if (argc < 4) {
    std::cout << "Usage: mpirun -np <nprocs> " << argv[0]
              << " <# partitions> "
              << " <mesh_file> <write_out(y/n)> <output_mesh_file> "
              << "[<method(RCB/RIB/HSFC/Hypergraph(PHG)/PARMETIS/OCTPART)>] "
              << "[<parmetis_method>/<oct_method>]" << std::endl;
    return 1;
  }

  bool write_output = false;
  if (argv[3][0] == 'y') write_output = true;
  else if (argv[3][0] == 'n') write_output = false;
  else {
    std::cout << "Argument 2 must be 'y' or 'n', not '" << argv[3][0] << "'" << std::endl;
    return 1;
  }

  MBInterface *mbImpl = new MBCore();
  
  MBErrorCode result = mbImpl->load_mesh(argv[2]); RR;
  
  MBZoltan *mbz = new MBZoltan(mbImpl);

  const char *other_method = NULL, *method = NULL;
  if (argc > 5) method = argv[5];
  if (argc > 6) other_method = argv[6];
  
  int nparts = atoi(argv[1]);
  
  result = mbz->partition_mesh(nparts, method, other_method); RR;
  
  if (write_output) {
    result = mbImpl->write_mesh(argv[4]); RR;
  }

  delete mbz;
  delete mbImpl;

  return 0;
}
