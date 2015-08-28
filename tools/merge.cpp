#include <iostream>
#include <moab/Core.hpp>
#include "moab/MergeMesh.hpp"
#include "moab/ProgOptions.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "moab/ParallelMergeMesh.hpp"
#endif

using namespace moab;

const char BRIEF_DESC[] = "Merges mesh files or entities in a mesh file. Use available options as desired.";
std::ostringstream LONG_DESC;

int main( int argc, char* argv[] )
{
#ifdef MOAB_HAVE_MPI
  int myID = 0, numprocs = 1;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &myID );
  MPI_Comm_size( MPI_COMM_WORLD, &numprocs );
#endif
  bool fsimple = true; //default
  bool ffile = false; // parmerge
  bool fall = false; // merge all
  std::string mtag = ""; // tag based merge
  std::string input_file, output_file;
  double merge_tol = 1.0e-4;

  LONG_DESC << "mbmerge tool has the ability to merge nodes in a mesh. For skin-based merge with multiple"
               "files parallel options is also supported." << std::endl
            << "If no method is specified, the default is simple merge. Simple merge case gets all the skins available"
            << " in the mesh file and merges the nodes to obtain a conformal mesh. Options to merge all duplicate nodes"
            << " and merging based on a specific tag on the nodes are also supported."
            << std::endl;

  ProgOptions opts(LONG_DESC.str(), BRIEF_DESC);

  opts.addOpt<void>("all,a", "merge all including interior.", &fall);
  opts.addOpt<std::string>( "mergetag name,t", "merge based on nodes that have a specific tag name assigned", &mtag);
  opts.addOpt<double>("mergetolerance,e", "merge tolerance, default is 1e-4", &merge_tol);
  opts.addOpt<void>("simple,s", "simple merge, merge based on skins provided as in the input mesh (Default)", &fsimple);
  opts.addRequiredArg<std::string>("input_file", "Input file to be merged", &input_file);
  opts.addRequiredArg<std::string>("output_file", "Output mesh file name with extension", &output_file);
#ifdef MOAB_HAVE_MPI
  opts.addOpt<void>("file,f", "files based merge using skin for individual meshes. Input file is a file containing names"
                    "of mesh files to be merged (each on a different line). Only works with parallel build", &ffile);
#endif
  opts.parseCommandLine(argc, argv);

  moab::Core *mb = new moab::Core();
  moab::ErrorCode rval;

  if(mtag != ""){
      rval = mb->load_mesh(input_file.c_str());
      if(rval != moab::MB_SUCCESS){
          std::cerr<<"Error Opening Mesh File "<< input_file << std::endl;
          return 1;
        }
      else{
          std::cout << "Read input mesh file: " << input_file << std::endl;
        }
      int dim = 0;
      moab::Range verts;
      mb->get_entities_by_dimension(0, dim, verts);
      if(rval != moab::MB_SUCCESS){
          std::cerr<< "failed to get entities by dimension" << std::endl;
          return 1;
        }
      Tag  tag_for_merge;
      rval = mb->tag_get_handle(mtag.c_str(), tag_for_merge);
      if(rval != moab::MB_SUCCESS){
          std::cerr<< "unable to get tag: " << mtag << " specified"<< std::endl;
          return 1;
        }
      MergeMesh mm(mb);
      rval = mm.merge_using_integer_tag(verts, tag_for_merge);
      if(rval != moab::MB_SUCCESS){
          std::cerr<< "error in routine merge using integer tag" << std::endl;
          return 1;
        }
      rval = mb->write_file( output_file.c_str());
      if(rval != moab::MB_SUCCESS){
          std::cerr<<"Error Writing Mesh File "<< output_file << std::endl;
          return 1;
        }
      else{
          std::cout << "Wrote output mesh file: " << output_file << std::endl;
        }
    }

  else if(fall == true){

      rval = mb->load_mesh(input_file.c_str());
      if(rval != moab::MB_SUCCESS){
          std::cerr<<"Error Opening Mesh File "<< input_file << std::endl;
          return 1;
        }
      else{
          std::cout << "Read input mesh file: " << input_file << std::endl;
        }
      MergeMesh mm(mb);
      rval = mm.merge_all(0, merge_tol); // root set
      if(rval != moab::MB_SUCCESS){
          std::cerr<< "error in merge_all routine" << std::endl;
          return 1;
        }
      rval = mb->write_file( output_file.c_str());
      if(rval != moab::MB_SUCCESS){
          std::cerr<<"Error Writing Mesh File "<< output_file << std::endl;
          return 1;
        }
      else{
          std::cout << "Wrote output mesh file: " << output_file << std::endl;
        }
    }
  else if(ffile == true){
      /*
          Parmerge
          Takes multiple mesh files and merges them into a single output file.
          This is a driver for ParallelMergeMesh
          Does not currently work if #procs > #meshfiles

          <input_file> is text file containing each mesh file on a line
          i.e.

          /my/path/file1
          /my/path/file2
          ...
          /my/path/fileN

          <output_file> file is a single file where the entire mesh is written to
          It must be of type ".h5m"

          <tolerance> is the merging tolerance

          Typical usage of:
          mpirun -n <#procs> parmerge <input_file> <output_file> <tolerance>
      */
      //Read in files from input files
      //Round robin distribution of reading meshes
#ifdef MOAB_HAVE_MPI
      std::ifstream file(input_file.c_str());
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
                      std::cerr<<"Error Opening Mesh File "<< line << std::endl;
                      MPI_Abort(MPI_COMM_WORLD,1);
                      file.close();
                      return 1;
                    }
                  else{
                      std::cout << "Read input mesh file: " << line << std::endl;
                    }
                }
              count = (count+1)%numprocs;
            }
          file.close();
        }
      else{
          std::cerr<<"Error Opening Input File "<< input_file  << std::endl;
          MPI_Abort(MPI_COMM_WORLD,1);
          return 1;
        }

      //Get a pcomm object
      moab::ParallelComm *pc = new moab::ParallelComm(mb, MPI_COMM_WORLD);

      //Call the resolve parallel function
      moab::ParallelMergeMesh pm(pc,merge_tol);
      rval = pm.merge();
      if(rval != moab::MB_SUCCESS){
          std::cerr<<"Merge Failed"<< std::endl;
          MPI_Abort(MPI_COMM_WORLD,1);
          return 1;
        }

      //Write out the file
      rval = mb->write_file(output_file.c_str() , 0,"PARALLEL=WRITE_PART");
      if(rval != moab::MB_SUCCESS){
          std::cerr<<"Writing output file failed. Code:";
          //Temporary File error info.
          std::cerr<<mb->get_error_string(rval)<< std::endl;
          std::string foo = ""; mb->get_last_error(foo);
          std::cerr<<"File Error: "<<foo<< std::endl;
          return 1;
        }
      else if(myID == 0){
          std::cout << "Wrote output mesh file: " << output_file << std::endl;
        }

      //The barrier may be necessary to stop items from being deleted when needed
      //But probably not necessary
      MPI_Barrier(MPI_COMM_WORLD);

      delete pc;
#endif
    }
  else if(fsimple ==true){
      rval = mb->load_mesh(input_file.c_str());
      if(rval != moab::MB_SUCCESS){
          std::cerr<<"Error Opening Mesh File "<< input_file << std::endl;
          return 1;
        }
      else{
          std::cout << "Read input mesh file: " << input_file << std::endl;
        }
      int dim = 3;
      moab::Range ents;
      mb->get_entities_by_dimension(0, dim, ents);
      if(rval != moab::MB_SUCCESS){
          std::cerr<< "error getting entities by dimension" << std::endl;
          return 1;
        }
      MergeMesh mm(mb);
      rval = mm.merge_entities(ents, merge_tol);
      if(rval != moab::MB_SUCCESS){
          std::cerr<< "error in merge entities routine" << std::endl;
          return 1;
        }
      rval = mb->write_file( output_file.c_str());
      if(rval != moab::MB_SUCCESS){
          std::cerr<<"   Writing Mesh File "<< output_file << std::endl;
          return 1;
        }
      else{
          std::cout << "Wrote output mesh file: " << output_file << std::endl;
        }
    }
  else{
      std::cerr<<" Unhandled option "<< std::endl;
      return 1;
    }

  delete mb;
#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
