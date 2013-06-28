/*
 * mpas file test
 *
 *  Created on: Feb 12, 2013
 */

// copy from case1 test

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
#include "moab/ProgOptions.hpp"
#include "MBTagConventions.hpp"

#include "CslamUtils.hpp"

#ifdef MESHDIR
std::string TestDir( STRINGIFY(MESHDIR) );
#else
std::string TestDir(".");
#endif

// for M_PI
#include <math.h>

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

using namespace moab;
// some input data
double gtol = 1.e-9; // this is for geometry tolerance

// radius is always 1?
//double CubeSide = 6.; // the above file starts with cube side 6; radius depends on cube side
double t = 0.1, delta_t = 0.1; // check the script

ErrorCode manufacture_lagrange_mesh_on_sphere(Interface * mb,
    EntityHandle euler_set, EntityHandle & lagr_set)
{
  ErrorCode rval = MB_SUCCESS;

  /*
   * get all plys first, then vertices, then move them on the surface of the sphere
   *  radius is 1., always
   *
   */
  double radius = 1.;
  Range polygons;
  rval = mb->get_entities_by_dimension(euler_set, 2, polygons);
  if (MB_SUCCESS != rval)
    return rval;

  Range connecVerts;
  rval = mb->get_connectivity(polygons, connecVerts);
  if (MB_SUCCESS != rval)
    return rval;

  // create new set
  rval = mb->create_meshset(MESHSET_SET, lagr_set);
  if (MB_SUCCESS != rval)
    return rval;

  // get the coordinates of the old mesh, and move it around the sphere according to case 1
  // now put the vertices in the right place....
  //int vix=0; // vertex index in new array

  // first create departure points (vertices in the lagrange mesh)
  // then connect them in quads
  std::map<EntityHandle, EntityHandle> newNodes;
  for (Range::iterator vit = connecVerts.begin(); vit != connecVerts.end();
      vit++)
  {
    EntityHandle oldV = *vit;
    CartVect posi;
    rval = mb->get_coords(&oldV, 1, &(posi[0]));
    if (MB_SUCCESS != rval)
      return rval;
    // cslam utils, case 1
    CartVect newPos;
    departure_point_case1(posi, t, delta_t, newPos);
    newPos = radius * newPos;
    EntityHandle new_vert;
    rval = mb->create_vertex(&(newPos[0]), new_vert);
    if (MB_SUCCESS != rval)
      return rval;
    newNodes[oldV] = new_vert;
  }
  EntityHandle new_conn[MAXEDGES]; // up to 10, for the time being
  for (Range::iterator it = polygons.begin(); it != polygons.end(); it++)
  {
    EntityHandle q = *it;
    int nnodes;
    const EntityHandle * conn4;
    rval = mb->get_connectivity(q, conn4, nnodes);
    if (MB_SUCCESS != rval)
      return rval;

    for (int i = 0; i < nnodes; i++)
    {
      EntityHandle v1 = conn4[i];
      new_conn[i] = newNodes[v1];
    }
    EntityHandle new_poly;
    rval = mb->create_element(MBQUAD, new_conn, nnodes, new_poly);
    if (MB_SUCCESS != rval)
      return rval;
    rval = mb->add_entities(lagr_set, &new_poly, 1);
    if (MB_SUCCESS != rval)
      return rval;
  }

  return rval;
}
int main(int argc, char **argv)
{

  MPI_Init(&argc, &argv);

  std::string fileN= TestDir + "/mpas_p8.h5m";
  const char *filename_mesh1 = fileN.c_str();
  if (argc > 1)
  {
    int index = 1;
    while (index < argc)
    {
      if (!strcmp(argv[index], "-gtol")) // this is for geometry tolerance
      {
        gtol = atof(argv[++index]);
      }
      if (!strcmp(argv[index], "-dt"))
      {
        delta_t = atof(argv[++index]);
      }
      if (!strcmp(argv[index], "-input"))
      {
        filename_mesh1 = argv[++index];
      }
      index++;
    }
  }
  // start copy
  std::string opts = std::string("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION")+
            std::string(";PARALLEL_RESOLVE_SHARED_ENTS");
  Core moab;
  Interface & mb = moab;
  EntityHandle euler_set;
  ErrorCode rval;
  rval = mb.create_meshset(MESHSET_SET, euler_set);
  CHECK_ERR(rval);


  rval = mb.load_file(filename_mesh1, &euler_set, opts.c_str());

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  CHECK_ERR(rval);

  rval = pcomm->check_all_shared_handles();
  CHECK_ERR(rval);
  // end copy
  int rank = pcomm->proc_config().proc_rank();

  if (0==rank)
    std::cout << " case 1: use -gtol " << gtol << " -dt " << delta_t <<
        " -input " << filename_mesh1 << "\n";

  // everybody will get a DP tag, including the non owned entities; so exchange tags is not required for LOC (here)
  EntityHandle lagrange_set;
  rval = manufacture_lagrange_mesh_on_sphere(&mb, euler_set, lagrange_set);
  if (MB_SUCCESS != rval)
    return 1;
  rval = mb.write_file("lagrIni.h5m", 0, 0, &lagrange_set, 1);
 if (MB_SUCCESS != rval)
   std::cout << "can't write lagr set\n";

  rval = enforce_convexity(&mb, lagrange_set);
  if (MB_SUCCESS != rval)
    return 1;

  std::stringstream ste;
  ste<<"lagr0" << rank<<".h5m";
  rval = mb.write_file(ste.str().c_str(), 0, 0, &euler_set, 1);

  if (MB_SUCCESS != rval)
    std::cout << "can't write lagr set\n";

  Intx2MeshOnSphere worker(&mb);

  double radius = 1.; // input

  worker.SetRadius(radius);

  worker.SetErrorTolerance(gtol);

  EntityHandle outputSet;
  rval = mb.create_meshset(MESHSET_SET, outputSet);
  if (MB_SUCCESS != rval)
    return 1;
  rval = worker.intersect_meshes(lagrange_set, euler_set, outputSet);
  if (MB_SUCCESS != rval)
    return 1;

  std::string opts_write("");
  std::stringstream outf;
  outf << "intersect1" << ".h5m";
  rval = mb.write_file(outf.str().c_str(), 0, 0, &outputSet, 1);
  if (MB_SUCCESS != rval)
    std::cout << "can't write output\n";
  double intx_area = area_on_sphere_lHuiller(&mb, outputSet, radius);
  double arrival_area = area_on_sphere_lHuiller(&mb, euler_set, radius);
  std::cout << " Arrival area: " << arrival_area
      << "  intersection area:" << intx_area << " rel error: "
      << fabs((intx_area - arrival_area) / arrival_area) << "\n";

  MPI_Finalize();
  if (MB_SUCCESS != rval)
    return 1;

  return 0;
}

