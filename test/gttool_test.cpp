/**
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */
/**
 * \file gttool_test.cpp
 *
 * \brief test geometrize method from GeomTopoTool
 *
 */
#include "moab/Core.hpp"
#include <iostream>

#include <assert.h>
#include <string.h>

#include "moab/GeomTopoTool.hpp"

#define STRINGIFY_(A) #A
#define STRINGIFY(A) STRINGIFY_(A)
#ifdef MESHDIR
std::string TestDir( STRINGIFY(MESHDIR) );
#else
std::string TestDir(".");
#endif

#define PROCESS_ERROR(A, B)  {if (A!=MB_SUCCESS) {  std::cout << B << std::endl; return 1; } }

#define CHECK( STR ) if (rval != MB_SUCCESS) return print_error( STR, rval, __FILE__, __LINE__ )

using namespace moab;

ErrorCode print_error(const char* desc, ErrorCode rval, const char* file,
    int line)
{
  std::cerr << "ERROR: " << desc << std::endl << "  Error code: " << rval
      << std::endl << "  At        : " << file << ':' << line << std::endl;
  return MB_FAILURE; // must always return false or CHECK macro will break
}

ErrorCode geometrize_test(Interface * mb, EntityHandle inputSet);

void handle_error_code(ErrorCode rv, int &number_failed, int &number_successful)
{
  if (rv == MB_SUCCESS) {
    std::cout << "Success";
    number_successful++;
  } else {
    std::cout << "Failure";
    number_failed++;
  }
}

int main(int argc, char *argv[])
{
  std::string filename = TestDir + "/partBed.smf";
  std::string engine_opt;

  if (argc == 1) {
    std::cout << "Using default input file: " << filename << std::endl;
  } else if (argc == 2) {
    filename = argv[1];
  } else {
    std::cerr << "Usage: " << argv[0] << " [geom_filename]" << std::endl;
    return 1;
  }

  int number_tests_successful = 0;
  int number_tests_failed = 0;

  Core mbcore;
  Interface * mb = &mbcore;

  mb->load_file(filename.c_str());

  //   FBEngine * pFacet = new FBEngine(mb, NULL, true);// smooth facetting, no OBB tree passed


  std::cout << "geometrize test: ";
  ErrorCode rval = geometrize_test( mb, 0); // just pass the root set
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";
}
ErrorCode geometrize_test(Interface * mb, EntityHandle inputSet)
{
  GeomTopoTool gtt(mb);
  EntityHandle outSet;
  ErrorCode rval=gtt.geometrize_surface_set(inputSet, outSet);
  if (rval !=MB_SUCCESS)
  {
    std::cout<<"Can't geometrize the set\n";
    return rval;
  }
  std::string ofile("output.h5m");
  std::cout<<"writing output file: " << ofile.c_str() << " ";
  rval=mb->write_file(ofile.c_str(), 0, 0, &outSet, 1);
  if (rval !=MB_SUCCESS)
  {
    std::cout<<"Can't write output file\n";
    return rval;
  }
  return MB_SUCCESS;
}



