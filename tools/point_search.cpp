/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

//Project
#include "moab/Core.hpp"
//Point Locater
#include "moab/point_locater/point_locater.hpp"
#include "moab/point_locater/io.hpp"

//iMesh
#include "imesh/iMesh_extensions.h"
#include "imesh/MBiMesh.hpp"



//STL
#include <string>
#include <vector>
#include <sstream>

namespace io = moab::point_locator::io;

#define PRINT_LAST_ERROR \
    if (MB_SUCCESS != result) {\
      std::string tmp_str;\
      std::cout << "Failure; message:" << std::endl;\
      mbImpl->get_last_error(tmp_str);\
      std::cout << tmp_str << std::endl;\
      MPI_Abort(MPI_COMM_WORLD, result);        \
      return result;\
    }

// Print usage
void print_usage() {
  std::cerr << "Usage: point_search" ;
  std::cerr << "-meshes <source_mesh> <target_mesh> ";
  std::cerr << " -itag <interp_tag> [-gnorm <gnorm_tag>] " ;
  std::cerr << " [-ssnorm <ssnorm_tag> <ssnorm_selection>] [-ropts <roptions>]";
  std::cerr << " [-outfile <out_file> [-wopts <woptions>]]"; 
  std::cerr << " [-dbgout [<dbg_file>]]" << std::endl;
  std::cerr << "    -meshes" << std::endl;
  std::cerr << "        Read in mesh files <source_mesh> and <target_mesh>." 
	    << std::endl;
  std::cerr << "    -itag" << std::endl;
  std::cerr << "        Interpolate tag <interp_tag> from source mesh to target mesh." << std::endl;
  std::cerr << "    -gnorm" << std::endl;
  std::cerr << "        Normalize the value of tag <gnorm_tag> over then entire mesh and save to" << std::endl;
  std::cerr << "        tag \"<gnorm_tag>_normf\" on the mesh set.  Do this for all meshes." << std::endl;
  std::cerr << "    -ssnorm" << std::endl;
  std::cerr << "        Normalize the value of tag <ssnorm_tag> over subsets of a mesh and save to" << std::endl;
  std::cerr << "        tag \"<ssnorm_tag>_normf\" on the Entity Set for each subset.  Subsets are selected" << std::endl;
  std::cerr << "        using criteria in <ssnorm_selection>.  Do this for all meshes." << std::endl;
  std::cerr << "    -ropts" << std::endl;
  std::cerr << "        Read in the mesh files using options in <roptions>." << std::endl;
  std::cerr << "    -outfile" << std::endl;
  std::cerr << "        Write out target mesh to <out_file>." << std::endl;
  std::cerr << "    -wopts" << std::endl;
  std::cerr << "        Write out mesh files using options in <woptions>." << std::endl;
  std::cerr << "    -dbgout" << std::endl;
  std::cerr << "        Write stdout and stderr streams to the file \'<dbg_file>.txt\'." << std::endl;
}
int main(int argc, char* argv[]){
   // need to init MPI first, to tell how many procs and rank
 int err = MPI_Init(&argc, &argv);


 io::File_options<> options;
                                                                                
 moab::ErrorCode result = moab::MB_SUCCESS;
 bool help = false;
 result = io::get_file_options(argc, argv, options);
                                                                                
 if (result != moab::MB_SUCCESS || help) {
   print_usage();
   err = MPI_Finalize();
   return 1;
 }
 
 int nprocs, rank;
 err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
 err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                                                                                
 // redirect stdout and stderr if dbgFile is not null
 if (!options.dbgFile.empty()) {
   std::stringstream dfname;
   dfname << options.dbgFile << rank << ".txt";
   std::freopen(dfname.str().c_str(), "a", stdout);
   std::freopen(dfname.str().c_str(), "a", stderr);
 }
                                                                                
   // start time
 double stime; //, rtime, setime, dtime, ltime;
 if (0 == rank) stime = MPI_Wtime();
                                                                                
 // create MOAB instance based on that
 moab::Interface * mbImpl = new moab::Core();
 if (NULL == mbImpl) return 1;
                                                                                
 // read in mesh(es)
 
 std::vector<moab::ParallelComm*> pcs( options.meshFiles.size());

 iBase_EntitySetHandle *roots = (iBase_EntitySetHandle *) malloc(sizeof(iBase_EntitySetHandle) * options.meshFiles.size());
 iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>( new MBiMesh(mbImpl) );


 for (unsigned int i = 0; i < options.meshFiles.size(); i++) {
   pcs[i] = new moab::ParallelComm( mbImpl, MPI_COMM_WORLD); 
   int index = pcs[i]->get_id();
   std::string newReadopts;
   std::ostringstream extraOpt;
   extraOpt  << ";PARALLEL_COMM=" << index;
   newReadopts = options.readOpts+extraOpt.str();
   iMesh_createEntSet(iMeshInst, 0, &(roots[i]), &err);
   result = mbImpl->load_file( options.meshFiles[i].c_str(), 
			       (moab::EntityHandle *)&roots[i], 
			       newReadopts.c_str() );
   PRINT_LAST_ERROR;
 }
 moab::Range src_elems;
 moab::Point_search locator(*mbImpl, pcs[0], src_elems, 0);
}
