/*
 * spec_visu_test.cpp
 * Will refine a quad mesh to a spectral visu mesh
 * It may even project on the sphere :)
 *
 *  Created on: Oct 21, 2012
 */

#include <stdio.h>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "CslamUtils.hpp"

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

using namespace moab;

int main(int argc, char* argv[])
{
  // check command line arg
  const char *filename_mesh = STRINGIFY(SRCDIR) "/eulerHomme.vtk";
  const char *newFile = "spectral.vtk";
  int NP = 4; // number of nodes in each direction
  double tolerance = 0.0001;
  double R = 6.; // should be input
  if (argc == 6)
  {
    filename_mesh = argv[1];
    NP = atoi(argv[2]);
    newFile = argv[3];
    R = atof(argv[4]);
    tolerance = atof(argv[5]);
  }
  else
  {
    printf("Usage: %s <mesh_filename> <NP>  <newFile> <tolerance>\n", argv[0]);
    if (argc != 1)
      return 1;
  }

  printf("run:  %s  %d  %s %f %f\n",
          filename_mesh, NP, newFile, R, tolerance);
  // read input mesh in a set
  ErrorCode rval = MB_SUCCESS;
  Core moab;
  Interface* mb = &moab;// global
  EntityHandle sf;
  rval = mb->create_meshset(MESHSET_SET, sf);
  if (MB_SUCCESS != rval)
    return 1;

  rval=mb->load_file(filename_mesh, &sf);
  if (MB_SUCCESS != rval)
    return 1;

  Range inputRange;
  // get quads of interest
  rval = mb->get_entities_by_type(sf, MBQUAD, inputRange);
  if (MB_SUCCESS != rval)
    return 1;

  EntityHandle outputSet;
  rval = mb->create_meshset(MESHSET_SET, outputSet);
  if (MB_SUCCESS != rval)
    return 1;

  rval = SpectralVisuMesh( mb,  inputRange, NP, outputSet, tolerance);

  if (MB_SUCCESS != rval)
    return 1;

  rval = ProjectOnSphere(mb, outputSet, R);
  if (MB_SUCCESS != rval)
    return 1;

  rval = mb->write_mesh(newFile, &outputSet, 1);
  if (MB_SUCCESS != rval)
    return 1;
  return 0;
}

