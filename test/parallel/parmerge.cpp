/* Parallel Merge Mesh Test
   Nathan Bertram
   1/31/11
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
#include <sstream>

int main(int argc, char * argv[])
{
  if(argc != 3){
    std::cerr<<"Usage: "<<argv[0]<<" <inputfile> <tolerance>"<<std::endl;
    return 1;
  }

  //Initialize MPI
  int numprocs, myID;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myID);
  
  //Read in tolerance
  double epsilon;
  if(!(std::istringstream(argv[2])>>epsilon)){
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
      if(myID == count && line != ""){
	rval = mb->load_mesh(line.c_str());
	if(rval != moab::MB_SUCCESS){
	  std::cerr<<"Error Opening Mesh File "<< line <<std::endl;
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
  moab::ParallelComm *pc = new moab::ParallelComm(mb); 

  //Call the resolve parallel function
  moab::ParallelMergeMesh pm(pc,epsilon);
  rval = pm.merge();
  if(rval != moab::MB_SUCCESS){
    std::cerr<<"Merge Failed"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    return 1;
  }

  //Write out the file
  std::stringstream np;
  np << numprocs;
  std::string outfile = "/home/nbertram/Desktop/meshes/shared_tagged";
  outfile = outfile + np.str() + ".h5m";
  rval = mb->write_file(outfile.c_str() , 0,"PARALLEL=WRITE_PART");
  if(rval != moab::MB_SUCCESS){
    std::cerr<<"Writing output file failed Code:";
    //Temporary File error info.
    std::cerr<<mb->get_error_string(rval)<<std::endl;
    std::string foo = ""; mb->get_last_error(foo);
    std::cerr<<"File Error: "<<foo<<std::endl;
    return 1;
  }

  //The barrier may be necessary to stop items from being deleted when needed
  //But probably not necessary
  MPI_Barrier(MPI_COMM_WORLD);
  
  delete pc;
  delete mb;
  MPI_Finalize();

  return 0;
}
