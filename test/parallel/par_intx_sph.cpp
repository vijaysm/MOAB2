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

void test_intx_in_parallel();

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int result = 0;

  result += RUN_TEST(test_intx_in_parallel);

  MPI_Finalize();
  return result;
}

ErrorCode  manufacture_lagrange_mesh_on_sphere(Interface * mb, EntityHandle euler_set, EntityHandle & out_lagrange_set)
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
  ReadUtilIface* readMeshIface;
  mb->query_interface(readMeshIface);
  double radius = 3*sqrt(3.);// our value ....
  Range quads;
  rval = mb->get_entities_by_type(euler_set, MBQUAD, quads);
  CHECK_ERR(rval);

  Range connecVerts;
  rval = mb->get_connectivity(quads, connecVerts);
  CHECK_ERR(rval);
  assert(connecVerts.psize()==1); // check if indeed we have only one subrange; then it is easy to duplicate mesh
  rval = mb->create_meshset(MESHSET_SET, out_lagrange_set);
  CHECK_ERR(rval);

  EntityHandle *conn_arr;
  EntityHandle start_quad;
  Range new_meshRange;
  // create new quads, and vertices at a manufactured positions
  rval = readMeshIface->get_element_connect(quads.size(), 4,
                                                MBQUAD, 0, start_quad, conn_arr);
  CHECK_ERR(rval);
  new_meshRange.insert(start_quad, start_quad + quads.size() - 1);

  EntityHandle start_vertex;
  std::vector<double*> arrays;
  rval = readMeshIface->get_node_coords(3, connecVerts.size(), 0, start_vertex, arrays);
  CHECK_ERR(rval);
  Range new_vert_range(start_vertex, start_vertex + connecVerts.size() - 1);

  double *xptr = arrays[0], *yptr = arrays[1], *zptr = arrays[2];

  // get the coordinates of the old mesh, and move it around the sphere in the same way as in the
  // python script

  // fill the conn arr with the old one, then modify it to kep the new vertices
  unsigned int shift = start_vertex-connecVerts[0];
  for (Range::iterator qit=quads.begin(); qit!=quads.end(); qit++)
  {
    EntityHandle eh=*qit;
    const EntityHandle * conn4;
    int num_nodes;
    rval = mb->get_connectivity(eh, conn4, num_nodes);
    CHECK_ERR(rval);
    for (int i=0; i<4; i++)
      conn_arr[i]=conn4[i]+shift;
    conn_arr+=4; // advance pointer in array
  }

  // now put the vertices in the right place....
  int vix=0; // vertex index in new array
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
    xptr[vix]=newPos[0];
    yptr[vix]=newPos[1];
    zptr[vix]=newPos[2];
    vix++;
  }
  rval = mb->add_entities(out_lagrange_set, new_meshRange);
  CHECK_ERR(rval);

  // duplicate vertices, quads!
/*
 * x, y, z = mesh.getVtxCoords(node)
  dist1 = math.sqrt( x*x + y*y + z*z )
  ratio = radius/dist1
  x1 = x*ratio
  y1 = y*ratio
  z1 = z*ratio
  mesh.setVtxCoords(node, [x1, y1, z1])
  [r, elev, az] = cart2sph(x,y,z)
  elev1=elev-2*math.pi*t/T
  uu=3*radius/T*math.pow(Sin(elev1),2)*Sin(2*az)*Cos(math.pi*t/T)+2*math.pi*radius*Cos(az)/T
  vv=3*radius/T*(Sin(2*elev1))*Cos(az)*Cos(math.pi*t/T)

  vx=-uu*Sin(az)-vv*Sin(elev)*Cos(az)
  vy= uu*Cos(az)-vv*Sin(elev)*Sin(az)
  vz= vv*Cos(elev)
  vectag[node] = [vx, vy, vz]

  # rotate an extra 360/(4*meshcount) around z, so the deformation will be bigger
rA = math.pi/(4*meshcount)*9.

   x = x + vx*time
  y = y + vy*time
  z = z + vz*time
  x2 = x * Cos(rA) - y*Sin(rA)
  y2 = x * Sin(rA) + y*Cos(rA)
  x = x2
  y = y2


 */
  // give the same gobal ids to the vertices
 
  int dum_val = 0;
  Tag gidTag;
  rval = mb->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gidTag, MB_TAG_DENSE, &dum_val);
  CHECK_ERR(rval);
  void *data;
  int count;
  rval = mb->tag_iterate(gidTag, connecVerts.begin(), connecVerts.end(), count, data);
  CHECK_ERR(rval);
  void * data2;
  int count2;
  rval = mb->tag_iterate(gidTag, new_vert_range.begin(), new_vert_range.end(), count2, data2);
  CHECK_ERR(rval);
  for (int i=0; i<count; i++)
     ((int*)data2)[i]= ((int*)data)[i];
  
  ParallelComm* pcomm = ParallelComm::get_pcomm(mb, 0);
  rval = pcomm->resolve_shared_ents(0, new_meshRange, 2, 0);
  CHECK_ERR(rval);

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

  EntityHandle lagrange_set;
  rval = manufacture_lagrange_mesh_on_sphere(&mb, euler_set, lagrange_set);
  

  CHECK_ERR(rval);
 
  std::string opts_write("PARALLEL=WRITE_PART");
  rval = mb.write_file("manuf.h5m", 0, opts_write.c_str(), &lagrange_set, 1);
  CHECK_ERR(rval);

}
