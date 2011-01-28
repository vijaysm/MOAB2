/* Resolve_Shared_Entities Driver
   Nathan Bertram
   1/17/11
*/

#include "moab/ParallelMergeMesh.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab_mpi.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include "moab/ParallelMergeMesh.hpp"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include <fstream>

int main(int argc, char * argv[])
{
  if(argc != 3){
    std::cerr<<"Usage: ./driver <inputfile> <tolerance>"<<std::endl;
    return 0;
  }

  //Initialize MPI
  int numprocs, id;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  
  //Read in tolerance
  char *end;
  double epsilon = strtod(argv[2], &end);
  if(end != NULL){
    std::cerr<<"Unable to parse epsilon.  Exiting"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    return 1;
  }
  
  //Read in files from input files
  moab::Core *mb = new moab::Core();
  moab::ErrorCode rval;
  std::ifstream file(argv[1]);
  if(file.is_open()){
    std::string line;
    int count = 0;
    while(file.good()){
      getline(file,line);
      if(id==count){
	rval = mb->load_mesh(line.c_str());
	if(rval != moab::MB_SUCCESS){
	  std::cerr<<"Error Opening File "<<file<<std::endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  file.close();
	  return 1;
	}
      }
      count = (count+1)%numprocs;
    }
    file.close();
  }
  else{
    std::cerr<<"Error Opening Input File "<< argv[1]<<std::endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    return 1;
  }

  //Get a pcomm object
  moab::ParallelComm pc(mb); 
 
  //Call the resolve parallel function
  moab::ParallelMergeMesh pm(&pc,epsilon);
  rval = pm.merge();
  if(rval != moab::MB_SUCCESS){
    MPI_Abort(MPI_COMM_WORLD,1);
    return 1;
  }

  //Write out the file
  std::string outfile = "shared_tagged.h5m";
  rval = mb->write_file(outfile.c_str(),0,"PARALLEL=WRITE_PART");
  if(rval != moab::MB_SUCCESS){
    MPI_Abort(MPI_COMM_WORLD,1);
    return 1;
  }

  delete mb;
  MPI_Finalize();

  return 0;
}
