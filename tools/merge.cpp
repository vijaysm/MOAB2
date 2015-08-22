#include <iostream>
#include <moab/Core.hpp>
#include "moab/MergeMesh.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "moab/ParallelMergeMesh.hpp"
#endif

using namespace moab;

const char usage[] = "[-a|-t|-f|-h] <input filename> <output filename> <tolerance>";

void error( const char* argv0, int rank )
{
  if (rank == 0){
      std::cerr << "Invalid arguments" << std::endl;
      std::cerr << "Usage: " << argv0 << " " << usage << std::endl;
    }
#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif
  exit(0);
}

void help( const char* argv0 ) {
  std::cout << argv0 << " " << usage << std::endl
            << "-a         Default option, will try to merge internal and skin vertices" << std::endl
            << "-f         Use true parallel read (default)" << std::endl
            << "-g         Set debug output level" << std::endl;
  exit(0);
}
int main( int argc, char* argv[] )
{
  int myID = 0, numprocs = 1;
#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &myID );
  MPI_Comm_size( MPI_COMM_WORLD, &numprocs );
#endif
  if (argc<2)
    {
      error(argv[0], myID);
    }

  bool fdo = true, ftag = false, ffile, fall, fhelp = false, first_pass = true;
  double epsilon;
  std::string outfile, inputfile;
  for (int i = 1; i < argc; i++)
    {
      if (!argv[i][0])
        error(argv[0], myID);

      if (fdo && argv[i][0] == '-')
        {
          if (!argv[i][1] || (argv[i][1] != 'H' && argv[i][2]))
            error(argv[0], myID);

          switch ( argv[i][1] )
            {
            case '-': fdo = false;       break;
            case 't': ftag = true;       break;
            case 'f': ffile = true;      break;
            case 'a': fall = true;       break;
            case 'h': fhelp = true;      break;
            default: std::cerr << "[" << myID << "]" << "Invalid option: " << argv[i] << std::endl;
            }
        }
      else if(argv[i]){
          //check if first pass - this is input file name
          if (first_pass){
              inputfile = (std::string)argv[i];
              first_pass = false;
            }
          else if (i < (argc-1)){ // any other string is output file
              outfile = (std::string)argv[i];
            }
          if(i==(argc-1)){ // last variable is always epsilon
              if(!(std::istringstream(argv[i])>>epsilon)){
                  std::cerr<<"Unable to parse tolerance"<<std::endl;
                  error(argv[0], myID);
                  return 1;
                }
            }
        }
    }

  if(fdo == true && ftag == true){

    }
  else if(fdo == true && ffile == true){
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
      //Read in files from input files
      //Round robin distribution of reading meshes
      moab::Core *mb = new moab::Core();
      moab::ErrorCode rval;
      std::ifstream file(inputfile.c_str());
      if(file.is_open()){
          std::string line;
          int count = 0;
          //Read each line
          while(file.good()){
              getline(file,line);
              if(myID == count && line != ""){
                  //Read in the file
                  std::cout << line << std::endl;
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
      moab::ParallelComm *pc = new moab::ParallelComm(mb, MPI_COMM_WORLD);

      //Call the resolve parallel function
      moab::ParallelMergeMesh pm(pc,epsilon);
      rval = pm.merge();
      std::cout <<"Writing ouewtput file failed. Code:" << std::endl;

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
    }
  else if(fdo == true && fall == true){

    }
  else if(fdo == true && fhelp == true ){
      if(myID == 0)
        help(argv[0]);
    }
  else{

    }

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif
  return 1;
}
