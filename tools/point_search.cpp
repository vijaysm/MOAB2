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
#include "moab/point_locater.hpp"

#include "moab/get_file_options.hpp"

//STL
#include <string>
#include <vector>
#include <sstream>


int main(int argc, char* argv[]){
   // need to init MPI first, to tell how many procs and rank
 int err = MPI_Init(&argc, &argv);
 
 std::vector<const char *> ssTagNames, ssTagValues;
 std::vector<std::string> meshFiles;
 std::string interpTag, gNormTag, ssNormTag, 
	     readOpts, outFile, writeOpts, dbgFile;
                                                                                
 moab::ErrorCode result = moab::MB_SUCCESS;
 bool help = false;
 result = get_file_options(argc, argv, meshFiles, interpTag,
                           gNormTag, ssNormTag, ssTagNames, ssTagValues,
                           readOpts, outFile, writeOpts, dbgFile, help);
                                                                                
 if (result != moab::MB_SUCCESS || help) {
   print_usage();
   err = MPI_Finalize();
   return 1;
 }
 
 int nprocs, rank;
 err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
 err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                                                                                
 // redirect stdout and stderr if dbgFile is not null
 if (!dbgFile.empty()) {
   std::stringstream dfname;
   dfname << dbgFile << rank << ".txt";
   std::freopen(dfname.str().c_str(), "a", stdout);
   std::freopen(dfname.str().c_str(), "a", stderr);
 }
                                                                                
   // start time
 double stime; //, rtime, setime, dtime, ltime;
 if (0 == rank) stime = MPI_Wtime();
                                                                                
 // create MOAB instance based on that
 Interface *mbImpl = new Core();
 if (NULL == mbImpl) return 1;
                                                                                
 // read in mesh(es)
 std::vector<ParallelComm *> pcs(meshFiles.size());                             
 for (unsigned int i = 0; i < meshFiles.size(); i++) {
   pcs[i] = new ParallelComm(mbImpl);
 }
 Range src_elems
 moab::Point_search locator(mbImpl, pcs[0], src_elems, 0);
}
