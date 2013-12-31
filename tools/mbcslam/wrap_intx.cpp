/*
 * wrap_intx.cpp
 *  Will implement the intersection method that will be callable from fortran too
 *  will be added to the library mbcslam.a
 *
 *
 *  Created on: Dec 14, 2013
 *      Author: iulian
 */
#include "iMesh.h"
#include "iMeshP.h"
#include "MBiMesh.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "Intx2MeshOnSphere.hpp"

using namespace moab;
double radius = 1.;
double gtol = 1.e-9;
bool debug = true;

#ifdef __cplusplus
extern "C" {
#endif

void update_tracer( iMesh_Instance instance, iBase_EntitySetHandle * opEulerSet, int * ierr)
{
  Range ents;
  moab::Interface * mb =MOABI;
  *ierr =1;
  ErrorCode rval = mb->get_entities_by_dimension(0, 2, ents);// root 0

  ERRORV(rval,  "can't get all 2d entities from root");

  EntityHandle euler_set;
  //EntityHandle lagr_set;

  rval = mb->create_meshset(MESHSET_SET, euler_set);
  ERRORV(rval , "can't create arrival mesh set");

  *opEulerSet = (iBase_EntitySetHandle)euler_set;

  rval = mb->add_entities(euler_set, ents);
  ERRORV(rval , "can't add ents to arrival set");

  Intx2MeshOnSphere worker(mb);
  worker.SetRadius(radius);

  worker.SetErrorTolerance(gtol);

  EntityHandle covering_lagr_set;
  rval = mb->create_meshset(MESHSET_SET, covering_lagr_set);
  ERRORV(rval , "can't create covering set ");

  // we need to update the correlation tag and remote tuples
  rval = worker.create_departure_mesh_2nd_alg(euler_set, covering_lagr_set);
  ERRORV(rval , "can't populate covering set ");

  if (debug)
  {
    rval = mb->write_file("lagr.h5m", 0, 0, &covering_lagr_set, 1  );
    ERRORV(rval , "can't write covering set ");
  }

  EntityHandle outputSet;
  rval = mb->create_meshset(MESHSET_SET, outputSet);
  ERRORV(rval , "can't create output set ");

  rval = worker.intersect_meshes(covering_lagr_set, euler_set, outputSet);
  ERRORV(rval , "can't intersect ");

  if (debug)
  {
    rval = mb->write_file("output.vtk", 0, 0, &outputSet, 1  );
    ERRORV(rval , "can't write covering set ");
  }

  // tagElem is the average computed at each element, from nodal values
  Tag tagElem = 0;
  std::string tag_name2("TracerAverage");
  rval = mb->tag_get_handle(tag_name2.c_str(), 1, MB_TYPE_DOUBLE, tagElem, MB_TAG_DENSE | MB_TAG_CREAT);
  ERRORV(rval , "can't get tracer tag ");

  // area of the euler element is fixed, store it; it is used to recompute the averages at each
  // time step
  Tag tagArea = 0;
  std::string tag_name4("Area");
  rval = mb->tag_get_handle(tag_name4.c_str(), 1, MB_TYPE_DOUBLE, tagArea, MB_TAG_DENSE | MB_TAG_CREAT);
  ERRORV(rval , "can't get area tag");

  rval = worker.update_tracer_data(outputSet, tagElem, tagArea);
  ERRORV(rval , "can't update tracer ");

  // everything can be deleted now from intx data; polygons, etc.

  *ierr = 0;
  return;
}

#ifdef __cplusplus
} // extern "C"
#endif


