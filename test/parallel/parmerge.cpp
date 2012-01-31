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
#include "moab/Skinner.hpp"

//Should we print out info on the merged mesh?
#define PrintInfo false
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
    mpirun -n <#procs> parmerge <inputfile> <outputfile> <tolerance> 
*/

//Function to print out info for testing purposes
void print_output(moab::ParallelComm *pc, moab::Core *mb, 
                  int numprocs, int myID, bool perform);

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

  print_output(pc, mb, myID, numprocs, PrintInfo);

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


//This function doesn't normally get called, but is here for debugging
//and verifying that merge is working.  
void print_output(moab::ParallelComm *pc, moab::Core *mb,
                  int myID, int numprocs, bool perform){
  moab::Range ents, skin;
  int o_ct=0, no_ct=0, tmp=0, o_tot=0, no_tot=0;
  if(perform){
    if(myID==0)std::cout<<"------------------------------------------"<<std::endl;
    //Check the count of total vertices
    mb->get_entities_by_dimension(0,0,ents);
    for(moab::Range::iterator rit = ents.begin(); rit != ents.end(); rit++){
      pc->get_owner(*rit, tmp);
      if(tmp==myID){
	o_ct++;
      }
    }
    MPI_Reduce(&o_ct, &o_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(myID==0){
      std::cout<<"There are " << o_tot << " vertices."<<std::endl;
      std::cout<<"------------------------------------------"<<std::endl;
    }
    //Check the count of owned and not owned skin faces.
    //owned-not owned == total skin faces
    moab::Skinner skinner(mb);
    o_ct=0; no_ct=0; o_tot=0; no_tot=0;
    skin.clear();
    ents.clear();
    mb->get_entities_by_dimension(0,3,ents);
    skinner.find_skin(ents, 2, skin);
    for(moab::Range::iterator s_rit = skin.begin(); 
        s_rit != skin.end(); s_rit++){
      pc->get_owner(*s_rit, tmp);
      if(tmp==myID){
	o_ct++;
      }
      else{
	no_ct++;
      }
    }
    MPI_Reduce(&o_ct, &o_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&no_ct, &no_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(myID == 0){
      std::cout<<"There are " << o_tot << " owned skin faces."<<std::endl;
      std::cout<<"There are " << no_tot << " not owned skin faces."<<std::endl;
      std::cout<<"The difference (Global Skin Faces) is " << (o_tot-no_tot) << "." << std::endl;
      std::cout<<"------------------------------------------"<<std::endl;
    }
  }
}
