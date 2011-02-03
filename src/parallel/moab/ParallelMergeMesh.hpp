#ifndef PARALLELMERGEMESH_HPP
#define PARALLELMERGEMESH_HPP

//Almost compiling

//Need to figure out how to not need this.
//Most likely a problem with my installation.
//My installation is installed with --with-mpi and is passing all the make check

#include "moab/Types.hpp"

namespace moab {

class ParallelComm;

class ParallelMergeMesh {
public:
  ParallelMergeMesh(ParallelComm *pc, 
		    const double epsilon);

  ErrorCode merge();

private:
  ParallelComm *myPcomm;
  double myEps;
};

} // namespace moab

#endif
