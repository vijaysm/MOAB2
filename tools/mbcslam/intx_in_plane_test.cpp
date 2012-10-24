/*
 * intx_in_plane_test.cpp
 *
 *  Created on: Oct 4, 2012
 */
#include <iostream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "Intx2MeshInPlane.hpp"
#include <math.h>

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

using namespace moab;

int main(int argc, char* argv[])
{
  // check command line arg
  const char *filename_mesh1 = STRINGIFY(SRCDIR) "/m1.vtk";
  const char *filename_mesh2 = STRINGIFY(SRCDIR) "/m2.vtk";
  const char *newFile = "intx1.vtk";
  if (argc == 4)
  {
    filename_mesh1 = argv[1];
    filename_mesh2 = argv[2];
    newFile = argv[3];
  }
  else
  {
    printf("Usage: %s <mesh_filename1> <mesh_filename2>  <newFile>\n", argv[0]);
    if (argc != 1)
      return 1;
    printf("No files specified.  Defaulting to: %s  %s  %s\n",
        filename_mesh1, filename_mesh2, newFile);
  }

  // read meshes in 2 file sets
  ErrorCode rval = MB_SUCCESS;
  Core moab;
  Interface* mb = &moab;// global
  EntityHandle sf1, sf2;
  rval = mb->create_meshset(MESHSET_SET, sf1);
  if (MB_SUCCESS != rval)
    return 1;
  rval = mb->create_meshset(MESHSET_SET, sf2);
  if (MB_SUCCESS != rval)
    return 1;
  rval=mb->load_file(filename_mesh1, &sf1);
  if (MB_SUCCESS != rval)
    return 1;
  rval=mb->load_file(filename_mesh2, &sf2);
  if (MB_SUCCESS != rval)
   return 1;


  EntityHandle outputSet;
  rval = mb->create_meshset(MESHSET_SET, outputSet);
  if (MB_SUCCESS != rval)
    return 1;
  Intx2MeshInPlane worker(mb);

  worker.SetErrorTolerance( 1.e-8);

  rval = worker.intersect_meshes(sf1, sf2, outputSet);
  if (MB_SUCCESS != rval)
    return 1;
  rval = mb->write_mesh(newFile, &outputSet, 1);
  if (MB_SUCCESS != rval)
    return 1;
  return 0;


}
