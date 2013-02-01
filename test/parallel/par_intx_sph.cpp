/*
 * par_intx_sph.cpp
 *  test to trigger intersection on a sphere in parallel
 *  it will start from an eulerian mesh
 *  the mesh is read in parallel; lagrangian mesh is manufactured on the fly (part of the test), in
 *    a different set.
 *
 *  lagrangian mesh will be located on the euler mesh; intersections will be performed on the
 *  euler mesh
 *  Created on: Nov 14, 2012
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
#include "TestUtil.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/ProgOptions.hpp"
#include "MBParallelConventions.h"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"

#include "CslamUtils.hpp"

// for M_PI
#include <math.h>

#ifdef MESHDIR
std::string TestDir( STRINGIFY(MESHDIR) );
#else
std::string TestDir(".");
#endif

using namespace moab;
double EPS1=0.2;
void test_intx_in_parallel();
void test_intx_in_parallel_elem_based();

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  EPS1 = 0.2;
  int result = 0;

  if (argc>1)
  {
    int index=1;
    while (index<argc)
    {
      if (!strcmp( argv[index], "-eps"))
      {
        EPS1=atof(argv[++index]);
      }
      index++;
    }
  }
  //result += RUN_TEST(test_intx_in_parallel);
  result += RUN_TEST(test_intx_in_parallel_elem_based);

  MPI_Finalize();
  return result;
}
// will save the LOC tag on the euler nodes
ErrorCode  manufacture_lagrange_mesh_on_sphere(Interface * mb, EntityHandle euler_set)
{
  ErrorCode rval = MB_SUCCESS;

  /*
   * get all quads first, then vertices, then move them on the surface of the sphere
   *  radius is in, it comes from MeshKit/python/examples/manufHomme.py :
   *  length = 6.
   *  each edge of the cube will be divided using this meshcount
   *  meshcount = 11
   *   circumscribed sphere radius
   *   radius = length * math.sqrt(3) /2
   */
  double radius = 3*sqrt(3.);// our value ....
  Range quads;
  rval = mb->get_entities_by_type(euler_set, MBQUAD, quads);
  CHECK_ERR(rval);

  Range connecVerts;
  rval = mb->get_connectivity(quads, connecVerts);

  // the LOC tag, should be provided by the user?
  Tag tagh = 0;
  std::string tag_name("DP");
  rval = mb->tag_get_handle(tag_name.c_str(), 3, MB_TYPE_DOUBLE, tagh, MB_TAG_DENSE | MB_TAG_CREAT);
  CHECK_ERR(rval);
  void *data; // pointer to the LOC in memory, for each vertex
  int count;

  rval = mb->tag_iterate(tagh, connecVerts.begin(), connecVerts.end(), count, data);
  CHECK_ERR(rval);
  // here we are checking contiguity
  assert(count == (int) connecVerts.size());
  double * ptr_DP=(double*)data;
  // get the coordinates of the old mesh, and move it around the sphere in the same way as in the
  // python script

  // now put the vertices in the right place....
  //int vix=0; // vertex index in new array
  double t=0.1, T=5;// check the script
  double time =0.05;
  double rot= M_PI/10;
  for (Range::iterator vit=connecVerts.begin();vit!=connecVerts.end(); vit++ )
  {
    EntityHandle oldV=*vit;
    CartVect posi;
    rval = mb->get_coords(&oldV, 1, &(posi[0]) );
    CHECK_ERR(rval);
    // do some mumbo jumbo, as in python script
    SphereCoords sphCoord = cart_to_spherical(posi);
    double lat1 = sphCoord.lat-2*M_PI*t/T; // 0.1/5
    double uu = 3*radius/ T * pow(sin(lat1), 2)*sin(2*sphCoord.lon)*cos(M_PI*t/T);
    uu+=2*M_PI*cos(sphCoord.lon)/T;
    double vv = 3*radius/T*(sin(2*lat1))*cos(sphCoord.lon)*cos(M_PI*t/T);
    double vx = -uu*sin(sphCoord.lon)-vv*sin(sphCoord.lat)*cos(sphCoord.lon);
    double vy = -uu*cos(sphCoord.lon)-vv*sin(sphCoord.lat)*sin(sphCoord.lon);
    double vz = vv*cos(sphCoord.lat);
    posi = posi + time * CartVect(vx, vy, vz);
    double x2= posi[0]*cos(rot)-posi[1]*sin(rot);
    double y2= posi[0]*sin(rot) + posi[1]*cos(rot);
    CartVect newPos(x2, y2, posi[2]);
    double len1= newPos.length();
    newPos = radius*newPos/len1;

    ptr_DP[0]=newPos[0];
    ptr_DP[1]=newPos[1];
    ptr_DP[2]=newPos[2];
    ptr_DP+=3; // increment to the next node
  }

  return rval;
}
void test_intx_in_parallel()
{
  std::string opts = std::string("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION")+
      std::string(";PARALLEL_RESOLVE_SHARED_ENTS");
  Core moab;
  Interface & mb = moab;
  EntityHandle euler_set;
  ErrorCode rval;
  rval = mb.create_meshset(MESHSET_SET, euler_set);
  CHECK_ERR(rval);
  std::string example(TestDir + "/Homme_2pt.h5m");

  rval = mb.load_file(example.c_str(), &euler_set, opts.c_str());

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  CHECK_ERR(rval);

  rval = pcomm->check_all_shared_handles();
  CHECK_ERR(rval);

  // everybody will get a LOC tag, including the non owned entities; so exchange tags is not required for LOC (here)
  rval = manufacture_lagrange_mesh_on_sphere(&mb, euler_set);
  CHECK_ERR(rval);
 
  Intx2MeshOnSphere worker(&mb);

  double radius= 3. * sqrt(3.) ; // input
  worker.SetRadius(radius);
  worker.set_box_error(EPS1);//
  worker.SetEntityType(MBQUAD);

  worker.SetErrorTolerance(radius*1.e-8);
  worker.locate_departure_points(euler_set);

  // we need to make sure the covering set is bigger than the euler mesh
  EntityHandle covering_lagr_set;
  rval = mb.create_meshset(MESHSET_SET, covering_lagr_set);
  CHECK_ERR(rval);

  rval = worker.create_departure_mesh(covering_lagr_set);
  CHECK_ERR(rval);
  int rank = pcomm->proc_config().proc_rank();
  std::stringstream ss;
  ss<<"partial" << rank<<".vtk";
  mb.write_file(ss.str().c_str(), 0, 0, &covering_lagr_set, 1);
  EntityHandle outputSet;
  rval = mb.create_meshset(MESHSET_SET, outputSet);
  CHECK_ERR(rval);
  rval = worker.intersect_meshes(covering_lagr_set, euler_set, outputSet);
  CHECK_ERR(rval);

  //std::string opts_write("PARALLEL=WRITE_PART");
  //rval = mb.write_file("manuf.h5m", 0, opts_write.c_str(), &outputSet, 1);
  std::string opts_write("");
  std::stringstream outf;
  outf<<"intersect" << rank<<".h5m";
  rval = mb.write_file(outf.str().c_str(), 0, 0, &outputSet, 1);
  CHECK_ERR(rval);


}
void test_intx_in_parallel_elem_based()
{
  std::string opts = std::string("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION")+
        std::string(";PARALLEL_RESOLVE_SHARED_ENTS");
  Core moab;
  Interface & mb = moab;
  EntityHandle euler_set;
  ErrorCode rval;
  rval = mb.create_meshset(MESHSET_SET, euler_set);
  CHECK_ERR(rval);
  std::string example(TestDir + "/Homme_2pt.h5m");

  rval = mb.load_file(example.c_str(), &euler_set, opts.c_str());

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
  CHECK_ERR(rval);

  rval = pcomm->check_all_shared_handles();
  CHECK_ERR(rval);

  // everybody will get a DP tag, including the non owned entities; so exchange tags is not required for LOC (here)
  rval = manufacture_lagrange_mesh_on_sphere(&mb, euler_set);
  CHECK_ERR(rval);

  int rank = pcomm->proc_config().proc_rank();

  std::stringstream ste;
  ste<<"initial" << rank<<".vtk";
  mb.write_file(ste.str().c_str(), 0, 0, &euler_set, 1);

  Intx2MeshOnSphere worker(&mb);

  double radius= 3. * sqrt(3.) ; // input
  worker.SetRadius(radius);
  worker.set_box_error(EPS1);//
  worker.SetEntityType(MBQUAD);

  worker.SetErrorTolerance(radius*1.e-8);
  //  worker.locate_departure_points(euler_set);

  // we need to make sure the covering set is bigger than the euler mesh
  EntityHandle covering_lagr_set;
  rval = mb.create_meshset(MESHSET_SET, covering_lagr_set);
  CHECK_ERR(rval);

  rval = worker.create_departure_mesh_2nd_alg(euler_set, covering_lagr_set);
  CHECK_ERR(rval);

  std::stringstream ss;
  ss<<"partial" << rank<<".vtk";
  mb.write_file(ss.str().c_str(), 0, 0, &covering_lagr_set, 1);
  EntityHandle outputSet;
  rval = mb.create_meshset(MESHSET_SET, outputSet);
  CHECK_ERR(rval);
  rval = worker.intersect_meshes(covering_lagr_set, euler_set, outputSet);
  CHECK_ERR(rval);

  //std::string opts_write("PARALLEL=WRITE_PART");
  //rval = mb.write_file("manuf.h5m", 0, opts_write.c_str(), &outputSet, 1);
  std::string opts_write("");
  std::stringstream outf;
  outf<<"intersect" << rank<<".h5m";
  rval = mb.write_file(outf.str().c_str(), 0, 0, &outputSet, 1);
  CHECK_ERR(rval);
}
