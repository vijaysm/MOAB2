/*
 * PolyExoWriter_test.cpp
 *
 *  Created on: Aug 15, 2014
 */


#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "../test/TestUtil.hpp"
#include "PolyExoWriter.hpp"

using namespace moab;

int main(int argc, char* argv[])
{
  // check command line arg
  const char *filename_mesh = STRINGIFY(MESHDIR) "/mbcslam/intx.vtk";
  const char *newFile = "intx.exo";

  if (argc == 3)
  {
    filename_mesh = argv[1];
    newFile = argv[2];
  }
  else
  {
    printf("Usage: %s <mesh_filename>   <newFile>\n", argv[0]);
    if (argc != 1)
      return 1;
    printf("No files specified.  Defaulting to:  %s  %s\n",
        filename_mesh, newFile);
  }

  // read meshes in 2 file sets
  ErrorCode rval = MB_SUCCESS;
  Core moab;
  Interface* mb = &moab;// global
  EntityHandle sf1;
  rval = mb->create_meshset(MESHSET_SET, sf1);
  if (MB_SUCCESS != rval)
    return 1;

  rval=mb->load_file(filename_mesh, &sf1);
  if (MB_SUCCESS != rval)
    return 1;

  PolyExoWriter polyw(mb);

  rval = polyw.write_poly_set(sf1, newFile);

  if (MB_SUCCESS != rval)
   return 1;

  return 0;


}

