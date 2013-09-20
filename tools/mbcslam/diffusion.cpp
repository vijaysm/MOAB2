/*
 * diffusion.cpp
 *
 *  Created on: Aug 12, 2013
 *
 */

/* trigger a diffusion computation in serial, first;
  we will start from an unstructured mesh on a sphere;
  do case 1: give some tracers a "slotted cylinder" ; compute the initial concentrations
  along a slotted cylinder, and see how it varies in time;
  use formula from Nair & Lauritzen paper:
  A class of deformational flow test cases for linear transport problems
on the sphere; see CSLAM Utils case1
  scalar fields are defined at page 4-5 in the paper

  steps:
   first create an initial tracer field, and save it as field on a sphere
   (initial time step 0)
   then use velocities (non-divergent, first), to transport the tracer on
   the sphere, according to some flow fields
*/
// copy from intx_mpas test

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
#include "TestUtil.hpp"
#include "moab/ParallelComm.hpp"

#include "CslamUtils.hpp"

// non smooth scalar field
// some input data
double gtol = 1.e-9; // this is for geometry tolerance

double radius = 1.;// in m:  6371220.

int numSteps = 50; // number of times with velocity displayed at points
double T = 5;

int case_number = 1; // 1, 2 (non-divergent) 3 divergent

moab::Tag corrTag;

int field_type = 1 ; // 1 quasi smooth, 2 - smooth, 3 non-smooth,
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
ErrorCode add_field_value(Interface * mb, EntityHandle euler_set, int rank, Tag & tagTracer, Tag & tagElem, Tag & tagArea)
{
  ErrorCode rval = MB_SUCCESS;

  /*
   * get all plys first, then vertices, then move them on the surface of the sphere
   *  radius is 1., most of the time
   *
   */
  Range polygons;
  rval = mb->get_entities_by_dimension(euler_set, 2, polygons);
  if (MB_SUCCESS != rval)
    return rval;

  Range connecVerts;
  rval = mb->get_connectivity(polygons, connecVerts);
  if (MB_SUCCESS != rval)
    return rval;



  void *data; // pointer to the LOC in memory, for each vertex
  int count;

  rval = mb->tag_iterate(tagTracer, connecVerts.begin(), connecVerts.end(), count, data);
  CHECK_ERR(rval);
  // here we are checking contiguity
  assert(count == (int) connecVerts.size());
  double * ptr_DP=(double*)data;
  // lambda is for longitude, theta for latitude
   // param will be: (la1, te1), (la2, te2), b, c; hmax=1, r=1/2
  // nondivergent flow, page 5, case 1, (la1, te1) = (M_PI, M_PI/3)
  //                                    (la2, te2) = (M_PI, -M_PI/3)
  //                 la1,    te1    la2    te2     b     c  hmax  r
  if (field_type==1) // quasi smooth
  {
    double params[] = { M_PI, M_PI/3, M_PI, -M_PI/3, 0.1, 0.9, 1., 0.5};
    for (Range::iterator vit=connecVerts.begin();vit!=connecVerts.end(); vit++ )
    {
      EntityHandle oldV=*vit;
      CartVect posi;
      rval = mb->get_coords(&oldV, 1, &(posi[0]) );
      CHECK_ERR(rval);

      SphereCoords sphCoord = cart_to_spherical(posi);

      ptr_DP[0]=quasi_smooth_field(sphCoord.lon, sphCoord.lat, params);;

      ptr_DP++; // increment to the next node
    }
  }
  else if (2 == field_type) // smooth
  {
    CartVect p1, p2;
    SphereCoords spr;
    spr.R = 1;
    spr.lat = M_PI/3;
    spr.lon= M_PI;
    p1 = spherical_to_cart(spr);
    spr.lat = -M_PI/3;
    p2 = spherical_to_cart(spr);
    //                  x1,    y1,     z1,    x2,   y2,    z2,   h_max, b0
    double params[] = { p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], 1,    5.};
    for (Range::iterator vit=connecVerts.begin();vit!=connecVerts.end(); vit++ )
    {
      EntityHandle oldV=*vit;
      CartVect posi;
      rval = mb->get_coords(&oldV, 1, &(posi[0]) );
      CHECK_ERR(rval);

      SphereCoords sphCoord = cart_to_spherical(posi);

      ptr_DP[0]=smooth_field(sphCoord.lon, sphCoord.lat, params);;

      ptr_DP++; // increment to the next node
    }
  }
  else if (3 == field_type) // slotted
  {
    //                   la1, te1,   la2, te2,       b,   c,   r
    double params[] = { M_PI, M_PI/3, M_PI, -M_PI/3, 0.1, 0.9, 0.5};// no h_max
    for (Range::iterator vit=connecVerts.begin();vit!=connecVerts.end(); vit++ )
    {
      EntityHandle oldV=*vit;
      CartVect posi;
      rval = mb->get_coords(&oldV, 1, &(posi[0]) );
      CHECK_ERR(rval);

      SphereCoords sphCoord = cart_to_spherical(posi);

      ptr_DP[0]=slotted_cylinder_field(sphCoord.lon, sphCoord.lat, params);;

      ptr_DP++; // increment to the next node
    }
  }

  // add average value for quad/polygon (average corners)
  // do some averages


  Range::iterator iter = polygons.begin();
  double total_mass = 0.;
  while (iter != polygons.end())
  {
    rval = mb->tag_iterate(tagElem, iter, polygons.end(), count, data);
    CHECK_ERR(rval);
    double * ptr=(double*)data;

    rval = mb->tag_iterate(tagArea, iter, polygons.end(), count, data);
    CHECK_ERR(rval);
    double * ptrArea=(double*)data;
    for (int i=0; i<count; i++, iter++, ptr++, ptrArea++)
    {
      const moab::EntityHandle * conn = NULL;
      int num_nodes = 0;
      rval = mb->get_connectivity(*iter, conn, num_nodes);
      CHECK_ERR(rval);
      if (num_nodes==0)
        return MB_FAILURE;
      std::vector<double> nodeVals(num_nodes);
      double average=0.;
      rval = mb->tag_get_data(tagTracer, conn, num_nodes, &nodeVals[0] );
      CHECK_ERR(rval);
      for (int j=0; j<num_nodes; j++)
        average+=nodeVals[j];
      average/=num_nodes;
      *ptr = average;

      // now get area
      std::vector<double> coords;
      coords.resize(3*num_nodes);
      rval = mb->get_coords(conn, num_nodes, &coords[0]);
      CHECK_ERR(rval);
      *ptrArea =  area_spherical_polygon_lHuiller (&coords[0], num_nodes, radius);

      // we should have used some
      // total mass:
      total_mass += *ptrArea * average;
    }

  }

  std::stringstream iniPos;
  iniPos<< "Tracer" << rank<<"_"<<0<<  ".vtk";// first time step

  rval = mb->write_file(iniPos.str().c_str(), 0, 0, &euler_set, 1);
  CHECK_ERR(rval);

  std::cout << "initial total mass:" << total_mass << "\n";

  // now we can delete the tags? not yet
  return MB_SUCCESS;
}

ErrorCode compute_velocity_case1(Interface * mb, EntityHandle euler_set, Tag & tagh, int rank, int tStep)
{
  ErrorCode rval = MB_SUCCESS;

  Range polygons;
  rval = mb->get_entities_by_dimension(euler_set, 2, polygons);
  if (MB_SUCCESS != rval)
    return rval;

  Range connecVerts;
  rval = mb->get_connectivity(polygons, connecVerts);
  if (MB_SUCCESS != rval)
    return rval;

  void *data; // pointer to the velo in memory, for each vertex
  int count;

  rval = mb->tag_iterate(tagh, connecVerts.begin(), connecVerts.end(), count, data);
  CHECK_ERR(rval);
  // here we are checking contiguity
  assert(count == (int) connecVerts.size());
  double * ptr_velo=(double*)data;
  // lambda is for longitude, theta for latitude

  for (Range::iterator vit=connecVerts.begin();vit!=connecVerts.end(); vit++ )
  {
    EntityHandle oldV=*vit;
    CartVect posi;
    rval = mb->get_coords(&oldV, 1, &(posi[0]) );
    CHECK_ERR(rval);
    CartVect velo ;
    double t = T * tStep/numSteps; //
    velocity_case1(posi, t, velo);

    ptr_velo[0]= velo[0];
    ptr_velo[1]= velo[1];
    ptr_velo[2]= velo[2];

    // increment to the next node
    ptr_velo+=3;// to next velocity
  }
  std::stringstream velos;
  velos<<"Tracer" << rank<<"_"<<tStep<<  ".vtk";
  rval = mb->write_file(velos.str().c_str(), 0, 0, &euler_set, 1);
  CHECK_ERR(rval);

  return MB_SUCCESS;
}
ErrorCode  create_lagr_mesh(Interface * mb, EntityHandle euler_set, EntityHandle lagr_set)
{
  // create the handle tag for the corresponding element / vertex

  EntityHandle dum = 0;

  ErrorCode rval = mb->tag_get_handle(CORRTAGNAME,
                                           1, MB_TYPE_HANDLE, corrTag,
                                           MB_TAG_DENSE|MB_TAG_CREAT, &dum);
  CHECK_ERR(rval);
  Range polys;
  rval = mb->get_entities_by_dimension(euler_set, 2, polys);
  CHECK_ERR(rval);

  Range connecVerts;
  rval = mb->get_connectivity(polys, connecVerts);
  CHECK_ERR(rval);

  std::map<EntityHandle, EntityHandle> newNodes;
  for (Range::iterator vit = connecVerts.begin(); vit != connecVerts.end(); vit++)
  {
    EntityHandle oldV = *vit;
    CartVect posi;
    rval = mb->get_coords(&oldV, 1, &(posi[0]));
    CHECK_ERR(rval);
    EntityHandle new_vert;
    rval = mb->create_vertex(&(posi[0]), new_vert); // duplicate the position
    CHECK_ERR(rval);
    newNodes[oldV] = new_vert;
    // set also the correspondent tag :)
    rval = mb->tag_set_data(corrTag, &oldV, 1, &new_vert);
    CHECK_ERR(rval);
    // also the other side
    rval = mb->tag_set_data(corrTag, &new_vert, 1, &oldV);
    CHECK_ERR(rval);
  }
  for (Range::iterator it = polys.begin(); it != polys.end(); it++)
  {
    EntityHandle q = *it;
    int nnodes;
    const EntityHandle * conn;
    rval = mb->get_connectivity(q, conn, nnodes);
    CHECK_ERR(rval);
    EntityType typeElem = mb->type_from_handle(q);
    std::vector<EntityHandle> new_conn(nnodes);
    for (int i = 0; i < nnodes; i++)
    {
      EntityHandle v1 = conn[i];
      new_conn[i] = newNodes[v1];
    }
    EntityHandle newElement;
    rval = mb->create_element(typeElem, &new_conn[0], nnodes, newElement);
    CHECK_ERR(rval);
    //set the corresponding tag
    rval = mb->tag_set_data(corrTag, &q, 1, &newElement);
    CHECK_ERR(rval);
    rval = mb->tag_set_data(corrTag, &newElement, 1, &q);
    CHECK_ERR(rval);

    rval = mb->add_entities(lagr_set, &newElement, 1);
    CHECK_ERR(rval);
  }

  return MB_SUCCESS;
}
ErrorCode compute_tracer_case1(Interface * mb, EntityHandle euler_set,
    EntityHandle lagr_set, EntityHandle out_set, Tag & tagElem, int rank,
    int tStep)
{
  ErrorCode rval = MB_SUCCESS;

  if (!corrTag)
    return MB_FAILURE;
  double t = tStep * T / numSteps; // numSteps is global; so is T
  double delta_t = T / numSteps; // this is global too, actually
  Range polys;
  rval = mb->get_entities_by_dimension(euler_set, 2, polys);
  CHECK_ERR(rval);

  Range connecVerts;
  rval = mb->get_connectivity(polys, connecVerts);
  CHECK_ERR(rval);


  // change coordinates of lagr mesh vertices
  for (Range::iterator vit = connecVerts.begin(); vit != connecVerts.end();
      vit++)
  {
    EntityHandle oldV = *vit;
    CartVect posi;
    rval = mb->get_coords(&oldV, 1, &(posi[0]));
    CHECK_ERR(rval);
    // cslam utils, case 1
    CartVect newPos;
    departure_point_case1(posi, t, delta_t, newPos);
    newPos = radius * newPos; // do we need this? the radius should be 1
    EntityHandle new_vert;
    rval = mb->tag_get_data(corrTag, &oldV, 1, &new_vert);
    CHECK_ERR(rval);
    // set the new position for the new vertex
    rval = mb->set_coords(&new_vert, 1, &(newPos[0]));
    CHECK_ERR(rval);
  }

  // so we have now the departure at the previous time
  // intersect the 2 meshes (what about some checking of convexity?) for sufficient
  // small dt, it is not an issue;
  Intx2MeshOnSphere worker(mb);
  worker.SetRadius(radius);

  worker.SetErrorTolerance(gtol);
  // std::cout << "error tolerance epsilon_1=" << gtol << "\n";

  rval = worker.intersect_meshes(lagr_set, euler_set, out_set);
  CHECK_ERR(rval);
  // serially: lagr is the same order as euler;
  // we need to update now the tracer information on each element, based on
  // initial value and areas of each resulting polygons
  rval = worker.update_tracer_data(out_set, tagElem);
  CHECK_ERR(rval);

  std::stringstream newTracer;
  newTracer << "Tracer" << rank << "_" << tStep << ".vtk";
  rval = mb->write_file(newTracer.str().c_str(), 0, 0, &euler_set, 1);
  CHECK_ERR(rval);

  std::stringstream newIntx;
  newIntx << "newIntx" << rank << "_" << tStep << ".vtk";
  rval = mb->write_file(newIntx.str().c_str(), 0, 0, &out_set, 1);
  CHECK_ERR(rval);
  // delete now the polygons and the elements of out_set
  // also, all verts that are not in euler set or lagr_set
  Range allVerts;
  rval = mb->get_entities_by_dimension(0, 0, allVerts);
  CHECK_ERR(rval);

  Range allElems;
  rval = mb->get_entities_by_dimension(0, 2, allElems);
  CHECK_ERR(rval);
  // add to polys range the lagr polys
  rval = mb->get_entities_by_dimension(lagr_set, 2, polys); // do not delete lagr set either, with its vertices
  CHECK_ERR(rval);
 // add to the connecVerts range all verts, from all initial polys
  rval = mb->get_connectivity(polys, connecVerts);
  CHECK_ERR(rval);
  Range todeleteVerts = subtract(allVerts, connecVerts);

  Range todeleteElem = subtract(allElems, polys);

  // empty the out mesh set
  rval = mb->clear_meshset(&out_set, 1);
  CHECK_ERR(rval);

  rval = mb->delete_entities(todeleteElem);
  CHECK_ERR(rval);
  rval = mb->delete_entities(todeleteVerts);
  CHECK_ERR(rval);
  return rval;
}
int main(int argc, char **argv)
{

  MPI_Init(&argc, &argv);

  std::string extra_read_opts;
  // read a homme file, partitioned in 16 so far
  std::string fileN= TestDir + "/HN16.h5m";
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

      if (!strcmp(argv[index], "-input"))
      {
        filename_mesh1 = argv[++index];
      }

      if (!strcmp(argv[index], "-O"))
      {
        extra_read_opts = std::string(argv[++index]);
      }

      if (!strcmp(argv[index], "-f"))
      {
        field_type = atoi(argv[++index]);
      }
      if (!strcmp(argv[index], "-ns"))
      {
        numSteps = atoi(argv[++index]);
      }

      if (!strcmp(argv[index], "-h"))
      {
        std::cout << "usage: -gtol <tol> -input <file> -O <extra_read_opts> \n   "
        <<    "-f <field_type> -h (this help) -ns <numSteps> \n";
        std::cout << " filed type: 1: quasi-smooth; 2: smooth; 3: slotted cylinders (non-smooth)\n";
        return 0;
      }
      index++;
    }
  }
  // start copy
  std::string opts = std::string("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION")+
            std::string(";PARALLEL_RESOLVE_SHARED_ENTS")+extra_read_opts;
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
    std::cout << " case 1: use -gtol " << gtol <<
        " -R " << radius << " -input " << filename_mesh1 <<  " -f " << field_type <<
        " numSteps: " << numSteps << "\n";

  Tag tagTracer = 0;
  std::string tag_name("Tracer");
  rval = mb.tag_get_handle(tag_name.c_str(), 1, MB_TYPE_DOUBLE, tagTracer, MB_TAG_DENSE | MB_TAG_CREAT);
  CHECK_ERR(rval);

  Tag tagElem = 0;
  std::string tag_name2("TracerAverage");
  rval = mb.tag_get_handle(tag_name2.c_str(), 1, MB_TYPE_DOUBLE, tagElem, MB_TAG_DENSE | MB_TAG_CREAT);
  CHECK_ERR(rval);

  Tag tagArea = 0;
  std::string tag_name4("Area");
  rval = mb.tag_get_handle(tag_name4.c_str(), 1, MB_TYPE_DOUBLE, tagArea, MB_TAG_DENSE | MB_TAG_CREAT);
  CHECK_ERR(rval);

  // add a field value, quasi smooth first
  rval = add_field_value(&mb, euler_set, rank, tagTracer, tagElem, tagArea);
  CHECK_ERR(rval);


  // do some velocity fields at some time steps; with animations
  // first delete the



  Tag tagh = 0;
  std::string tag_name3("Case1");
  rval = mb.tag_get_handle(tag_name3.c_str(), 3, MB_TYPE_DOUBLE, tagh, MB_TAG_DENSE | MB_TAG_CREAT);
  CHECK_ERR(rval);
  EntityHandle out_set, lagr_set;
  rval = mb.create_meshset(MESHSET_SET, out_set);
  CHECK_ERR(rval);
  rval = mb.create_meshset(MESHSET_SET, lagr_set);
  CHECK_ERR(rval);
  // copy the initial mesh in the lagrangian set
  // initial vertices will be at the same position as euler;

  rval = create_lagr_mesh(&mb, euler_set, lagr_set);
  CHECK_ERR(rval);
  for (int i=1; i<numSteps+1; i++)
  {
    // time depends on i; t = i*T/numSteps: ( 0, T/numSteps, 2*T/numSteps, ..., T )
    rval = compute_velocity_case1(&mb, euler_set, tagh, rank, i);
    CHECK_ERR(rval);

    // this is to actually compute concentrations, using the current concentrations
    //
    rval = compute_tracer_case1(&mb, euler_set, lagr_set, out_set,
        tagElem, rank, i);

  }
  return 0;
}
