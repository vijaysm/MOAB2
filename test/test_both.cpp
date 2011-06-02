#include "iGeom.h"
#include "iMesh.h"
#include "iRel.h"

#include "TestUtil.hpp"

#include <stdlib.h>

void test_both()
{
  int err;
  iGeom_Instance geom;
  iMesh_Instance mesh;
  iRel_Instance rel;
  iRel_PairHandle pair;

  iGeom_newGeom(0, &geom, &err, 0);
  CHECK_ERR(err);
  iMesh_newMesh(0, &mesh, &err, 0);
  CHECK_ERR(err);
  iRel_create(0, &rel, &err, 0);
  CHECK_ERR(err);

  iRel_createPair(rel, geom, iRel_ENTITY, iRel_IGEOM_IFACE, iRel_ACTIVE,
                       mesh, iRel_BOTH,   iRel_IMESH_IFACE, iRel_ACTIVE,
                  &pair, &err);
  CHECK_ERR(err);

  iBase_EntityHandle geom_ent;
  iGeom_createBrick(geom, 2, 2, 2, &geom_ent, &err);
  CHECK_ERR(err);

  double coords[] = {
    0, 0, 0,
    0, 1, 0,
    1, 1, 0,
    1, 0, 0,
  };
  iBase_EntityHandle mesh_ents[4];
  iBase_EntityHandle *mesh_ents_ptr = mesh_ents;
  int mesh_ents_alloc = 4, mesh_ents_size;
  iMesh_createVtxArr(mesh, 4, iBase_INTERLEAVED, coords, 12,
                     &mesh_ents_ptr, &mesh_ents_alloc, &mesh_ents_size, &err);
  CHECK_ERR(err);

  iBase_EntitySetHandle mesh_set;
  iMesh_createEntSet(mesh, false, &mesh_set, &err);
  CHECK_ERR(err);
  iMesh_addEntArrToSet(mesh, mesh_ents, mesh_ents_size, mesh_set, &err);
  CHECK_ERR(err);

  iRel_setEntSetRelation(rel, pair, geom_ent, mesh_set, &err);
  CHECK_ERR(err);

  iBase_EntityHandle *related_ents = NULL;
  int related_ents_alloc = 0, related_ents_size;
  iRel_getEntArrEntArrRelation(rel, pair, mesh_ents, mesh_ents_size, 1,
                               &related_ents, &related_ents_alloc,
                               &related_ents_size, &err);
  CHECK_ERR(err);

  for(int i = 0; i < related_ents_size; i++)
    CHECK_EQUAL(related_ents[i], geom_ent);

  free(related_ents);

  iRel_destroy(rel, &err);
  CHECK_ERR(err);
  iMesh_dtor(mesh, &err);
  CHECK_ERR(err);
  iGeom_dtor(geom, &err);
  CHECK_ERR(err);
}

int main()
{
  int num_fail = 0;

  num_fail += RUN_TEST(test_both);

  return num_fail;
}
