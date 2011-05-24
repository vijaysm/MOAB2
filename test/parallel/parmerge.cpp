/* Parallel Merge Mesh Test
   Nathan Bertram
   5/31/11
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

/*  
    Parmerge
    Takes multiple mesh files and merges them into a single output file.
    This is a driver for ParallelMergeMesh
    Does not currently work if #procs > #meshfiles

    <inputfile> is text file containing each mesh file on a line
    i.e.

    /my/path/file1
    /my/path/file2
    ...
    /my/path/fileN

    <outputfile> file is a single file where the entire mesh is written to
    It must be of type ".h5m"

    <tolerance> is the merging tolerance
    
    Typical usage of:
    mpd &
    mpirun -n <#procs> parmerge <inputfile> <outputfile> <tolerance>
*/
int main(int argc, char * argv[])
{
  //Check argument count
  if(argc != 4){
    std::cerr<<"Usage: "<<argv[0]<<" <inputfile> <outputfile> <tolerance>"<<std::endl;
    return 1;
  }
  //Check the output file extension
  std::string outfile(argv[2]);
  if(outfile.compare(outfile.size()-4,4,".h5m")!=0){
    std::cerr<<"Invalid Parallel Output File"<<std::endl;
    return 1;
  }

  //Read in tolerance
  double epsilon;
  if(!(std::istringstream(argv[3])>>epsilon)){
    std::cerr<<"Unable to parse tolerance"<<std::endl;
    return 1;
  }

  //Initialize MPI
  int numprocs, myID;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myID);
  
  //Read in files from input files
  //Round robin distribution of reading meshes
  moab::Core *mb = new moab::Core();
  moab::ErrorCode rval;
  std::ifstream file(argv[1]);
  if(file.is_open()){
    std::string line;
    int count = 0;
    //Read each line
    while(file.good()){
      getline(file,line);
      if(myID == count && line != ""){
	//Read in the file
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
  rval = mb->write_file(outfile.c_str() , 0,"PARALLEL=WRITE_PART");
  if(rval != moab::MB_SUCCESS){
    std::cerr<<"Writing output file failed. Code:";
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
