/*
 * intx_on_sphere_test.cpp
 *
 *  Created on: Oct 3, 2012
 *      Author: iulian
 */
#include <iostream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "Intx2MeshOnSphere.hpp"
#include <math.h>

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

using namespace moab;

int main(int argc, char* argv[])
{
  // check command line arg// Euler grid is red, arrival, Lagrangian is blue, departure
  // will will keep the
  const char *filename_mesh1 = STRINGIFY(SRCDIR) "/lagrangeHomme.vtk";
  const char *filename_mesh2 = STRINGIFY(SRCDIR) "/eulerHomme.vtk";
  double R = 6. * sqrt(3.) / 2; // input
  const char *newFile = "intx.vtk";
  if (argc == 5)
  {
    filename_mesh1 = argv[1];
    filename_mesh2 = argv[2];
    R = atof(argv[3]);
    newFile = argv[4];
  }
  else
  {
    printf("Usage: %s <mesh_filename1> <mesh_filename2> <radius>  <newFile>\n",
        argv[0]);
    if (argc != 1)
      return 1;
    printf("No files specified.  Defaulting to: %s  %s  %f %s\n",
        filename_mesh1, filename_mesh2, R, newFile);
  }

  // read meshes in 2 file sets
  ErrorCode rval = MB_SUCCESS;
  Core moab;
  Interface * mb = &moab; // global
  EntityHandle sf1, sf2;
  rval = mb->create_meshset(MESHSET_SET, sf1);
  if (MB_SUCCESS != rval)
    return 1;
  rval = mb->create_meshset(MESHSET_SET, sf2);
  if (MB_SUCCESS != rval)
    return 1;
  rval = mb->load_file(filename_mesh1, &sf1);
  if (MB_SUCCESS != rval)
    return 1;
  rval = mb->load_file(filename_mesh2, &sf2);
  if (MB_SUCCESS != rval)
    return 1;

  EntityHandle outputSet;
  rval = mb->create_meshset(MESHSET_SET, outputSet);
  if (MB_SUCCESS != rval)
    return 1;

  Intx2MeshOnSphere  worker(mb);
  double radius= 6. * sqrt(3.) / 2; // input
  worker.SetErrorTolerance(radius*1.e-8);
  worker.SetRadius(radius);
  rval = worker.intersect_meshes(sf1, sf2, outputSet);
  if (MB_SUCCESS != rval)
    std::cout << " failed to intersect meshes\n";
  rval = mb->write_mesh(newFile, &outputSet, 1);
  if (MB_SUCCESS != rval)
    return 1;
  return 0;

}


