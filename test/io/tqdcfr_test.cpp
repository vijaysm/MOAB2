/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#include "Tqdcfr.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/FileOptions.hpp"
#include <iostream>
#include <string>

#ifdef USE_MPI
#include "moab_mpi.h"
#endif

#define STRINGIFY_(A) #A
#define STRINGIFY(A) STRINGIFY_(A)

using namespace moab;

int main(int argc, char* argv[])
{
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
    // Check command line arg
  const char* file = STRINGIFY(MESHDIR) "/io/brick_cubit10.2.cub";
  if (argc < 2)
  {
    std::cout << "Usage: tqdcfr <cub_file_name>" << std::endl;
      //exit(1);
  }
  else
    file = argv[1];

  Core *my_impl = new Core();
  Tqdcfr *my_tqd = new Tqdcfr(my_impl);
  FileOptions opts(NULL);
  
  ErrorCode result = my_tqd->load_file(file, 0, opts, 0, 0);

  if (MB_SUCCESS == result)
    std::cout << "Success." << std::endl;
  else {
    std::cout << "load_file returned error:" << std::endl;
    std::string errs;
    result = my_impl->get_last_error(errs);
    if (MB_SUCCESS == result) std::cout << errs << std::endl;
    else std::cout << "(no message)" << std::endl;
  }

  delete my_tqd;
  delete my_impl;
  
    // now check for multiple procs
  my_impl = new Core;
  my_tqd = new Tqdcfr(my_impl);
  
  result = my_tqd->load_file(file, 0, opts, 0, 0);

  if (MB_SUCCESS == result)
    std::cout << "Success." << std::endl;
  else {
    std::cout << "load_file returned error:" << std::endl;
    std::string errstr;
    result = my_impl->get_last_error(errstr);
    if (MB_SUCCESS == result) std::cout << errstr << std::endl;
    else std::cout << "(no message)" << std::endl;
  }

  delete my_tqd;
  delete my_impl;

#ifdef USE_MPI
  int nprocs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // create MOAB instance based on that
  my_impl = new Core ;//(rank, nprocs);
  if (NULL == my_impl) return 1;
  
  std::string options = "PARALLEL=READ_DELETE;PARTITION=MATERIAL_SET;PARTITION_DISTRIBUTE";
  std::cout << "Testing parallel..." << std::endl;
  
  result = my_impl->load_file(file, 0, 
                              options.c_str());

  if (MB_SUCCESS == result)
    std::cout << "Success." << std::endl;
  else {
    std::cout << "load_file returned error:" << std::endl;
    std::string errstr;
    result = my_impl->get_last_error(errstr);
    if (MB_SUCCESS == result) std::cout << errstr << std::endl;
    else std::cout << "(no message)" << std::endl;
  }

  delete my_impl; // done with the read
  MPI_Finalize();
#endif

  return result;
}

