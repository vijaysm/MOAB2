#include "AssocPairC.hpp"

#include "iGeom.h"
#include "iMesh.h"
#include "iRel.h"
#include "Lasso.hpp"

#include <string.h>
#include <stdio.h>

#define PROCESS_GERROR do {                                             \
    if (iBase_SUCCESS != result) {                                      \
      char this_descr[120];                                             \
      iGeom_getDescription((iGeom_Instance)ifaceInstances[iface_no],    \
                           this_descr, 120);                            \
      ERRORR(result, this_descr);                                       \
    }                                                                   \
  } while(false)

#define PROCESS_MERROR do {                                             \
    if (iBase_SUCCESS != result) {                                      \
      char this_descr[120];                                             \
      iMesh_getDescription((iMesh_Instance)ifaceInstances[iface_no],    \
                           this_descr, 120);                            \
      ERRORR(result, this_descr);                                       \
    }                                                                   \
  } while(false)

AssocPairC::AssocPairC(iRel_Instance instance,
                       iBase_Instance iface0,
                       iRel_RelationType ent_or_set0,
                       iRel_IfaceType type0,
                       iBase_Instance iface1,
                       iRel_RelationType ent_or_set1,
                       iRel_IfaceType type1)
  : AssocPair(instance, ent_or_set0, type0, ent_or_set1, type1)
{
  ifaceInstances[0] = iface0;
  ifaceInstances[1] = iface1;

  // finally, create the tags we'll need
  create_tags();
}

AssocPairC::~AssocPairC()
{
  destroy_tags();
}

iBase_Instance AssocPairC::iface_instance(const int iface_no)
{
  return ifaceInstances[iface_no];
}

bool AssocPairC::equivalent(iBase_Instance iface0,
                            iBase_Instance iface1,
                            bool *order_switched)
{
  bool equiv = false;
  if (iface0 == ifaceInstances[0] &&
      iface1 == ifaceInstances[1]) {
    if (NULL != order_switched) *order_switched = false;
    equiv = true;
  }
  else if (iface0 == ifaceInstances[1] &&
           iface1 == ifaceInstances[0]) {
    equiv = true;
    if (NULL != order_switched) *order_switched = true;
  }

  return equiv;
}

bool AssocPairC::contains(iBase_Instance iface)
{
  return (iface == ifaceInstances[0] ||
          iface == ifaceInstances[1]);
}

bool AssocPairC::same_interface(iBase_Instance iface0,
                                iBase_Instance iface1)
{
  return iface0 == iface1;
}

int AssocPairC::get_tags(const int iface_no,
                         iBase_EntityHandle *entities,
                         const int num_entities,
                         iBase_TagHandle tag_handle,
                         void *tag_values,
                         const int tag_size)
{
  int tag_values_alloc = num_entities * tag_size, tag_values_size;
  int result;

  if (iRel_IGEOM_IFACE == iface_type(iface_no)) {
    iGeom_getArrData((iGeom_Instance)ifaceInstances[iface_no],
                     entities, num_entities, tag_handle,
                     &tag_values, &tag_values_alloc,
                     &tag_values_size, &result);
    PROCESS_GERROR;
  }
  else if (iRel_IMESH_IFACE == iface_type(iface_no)) {
    iMesh_getArrData((iMesh_Instance)ifaceInstances[iface_no],
                     entities, num_entities, tag_handle,
                     &tag_values, &tag_values_alloc,
                     &tag_values_size, &result);
    PROCESS_MERROR;
  }
  else
    ERRORR(iBase_NOT_SUPPORTED, "Interface should be geometry or mesh.");

  RETURNR(iBase_SUCCESS);
}

int AssocPairC::get_tags(const int iface_no,
                         iBase_EntitySetHandle *sets,
                         const int num_sets,
                         iBase_TagHandle tag_handle,
                         void *tag_values,
                         const int tag_size)
{
  char *tag_data = static_cast<char*>(tag_values);
  int tag_values_alloc = tag_size, tag_values_size;
  int result;

  if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
    for (int i = 0; i < num_sets; i++) {
      iGeom_getEntSetData((iGeom_Instance)ifaceInstances[iface_no],
                          sets[i], tag_handle,
                          reinterpret_cast<void**>(&tag_data),
                          &tag_values_alloc, &tag_values_size, &result);
      tag_data += tag_size;
      PROCESS_GERROR;
    }
  }
  else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
    for (int i = 0; i < num_sets; i++) {
      iMesh_getEntSetData((iMesh_Instance)ifaceInstances[iface_no],
                          sets[i], tag_handle,
                          reinterpret_cast<void**>(&tag_data),
                          &tag_values_alloc, &tag_values_size, &result);
      tag_data += tag_size;
      PROCESS_MERROR;
    }
  }
  else
    ERRORR(iBase_NOT_SUPPORTED, "Interface should be geometry or mesh.");

  RETURNR(iBase_SUCCESS);
}

int AssocPairC::set_tags(const int iface_no,
                         iBase_EntityHandle *entities,
                         const int num_entities,
                         iBase_TagHandle tag_handle,
                         const void *tag_values,
                         const int tag_size)
{
  const char *tag_data = static_cast<const char*>(tag_values);
  int result;

  if (iRel_IGEOM_IFACE == iface_type(iface_no)) {
    iGeom_setArrData((iGeom_Instance)ifaceInstances[iface_no],
                     entities, num_entities, tag_handle,
                     tag_data, num_entities * tag_size, &result);
    PROCESS_GERROR;
  }
  else if (iRel_IMESH_IFACE == iface_type(iface_no)) {
    iMesh_setArrData((iMesh_Instance)ifaceInstances[iface_no],
                     entities, num_entities, tag_handle,
                     tag_data, num_entities * tag_size, &result);
    PROCESS_MERROR;
  }
  else
    ERRORR(iBase_NOT_SUPPORTED, "Interface should be geometry or mesh.");

  RETURNR(iBase_SUCCESS);
}

int AssocPairC::set_tags(const int iface_no,
                         iBase_EntitySetHandle *sets,
                         const int num_sets,
                         iBase_TagHandle tag_handle,
                         const void *tag_values,
                         const int tag_size)
{
  const char *tag_data = static_cast<const char*>(tag_values);
  int result;

  if (iRel_IGEOM_IFACE == iface_type(iface_no)) {
    for (int i = 0; i < num_sets; i++) {
      iGeom_setEntSetData((iGeom_Instance)ifaceInstances[iface_no],
                          sets[i], tag_handle, tag_data, tag_size, &result);
      tag_data += tag_size;
      PROCESS_GERROR;
    }
  }
  else if (iRel_IMESH_IFACE == iface_type(iface_no)) {
    for (int i = 0; i < num_sets; i++) {
      iMesh_setEntSetData((iMesh_Instance)ifaceInstances[iface_no],
                          sets[i], tag_handle, tag_data, tag_size, &result);
      tag_data += tag_size;
      PROCESS_MERROR;
    }
  }
  else
    ERRORR(iBase_NOT_SUPPORTED, "Interface should be geometry or mesh.");

  RETURNR(iBase_SUCCESS);
}

int AssocPairC::rmv_tags(const int iface_no,
                         iBase_EntityHandle *entities,
                         const int num_entities,
                         iBase_TagHandle tag_handle)
{
  int result;

  if (iRel_IGEOM_IFACE == iface_type(iface_no)) {
    iGeom_rmvArrTag((iGeom_Instance)ifaceInstances[iface_no],
                    entities, num_entities, tag_handle, &result);
    PROCESS_GERROR;
  }
  else if (iRel_IMESH_IFACE == iface_type(iface_no)) {
    iMesh_rmvArrTag((iMesh_Instance)ifaceInstances[iface_no],
                     entities, num_entities, tag_handle, &result);
    PROCESS_MERROR;
  }
  else
    ERRORR(iBase_NOT_SUPPORTED, "Interface should be geometry or mesh.");

  RETURNR(iBase_SUCCESS);
}

int AssocPairC::rmv_tags(const int iface_no,
                         iBase_EntitySetHandle *sets,
                         const int num_sets,
                         iBase_TagHandle tag_handle)
{
  int result;

  if (iRel_IGEOM_IFACE == iface_type(iface_no)) {
    for (int i = 0; i < num_sets; i++) {
      iGeom_rmvEntSetTag((iGeom_Instance)ifaceInstances[iface_no],
                          sets[i], tag_handle, &result);
      PROCESS_GERROR;
    }
  }
  else if (iRel_IMESH_IFACE == iface_type(iface_no)) {
    for (int i = 0; i < num_sets; i++) {
      iMesh_rmvEntSetTag((iMesh_Instance)ifaceInstances[iface_no],
                         sets[i], tag_handle, &result);
      PROCESS_MERROR;
    }
  }
  else
    ERRORR(iBase_NOT_SUPPORTED, "Interface should be geometry or mesh.");

  RETURNR(iBase_SUCCESS);
}

int AssocPairC::get_iterator(const int iface_no,
                             iBase_EntitySetHandle set,
                             iBase_EntityIterator *iter)
{
  int result;

  if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
    iGeom_initEntIter((iGeom_Instance)ifaceInstances[iface_no], set,
                      iBase_ALL_TYPES,
                      iter, &result);
    PROCESS_GERROR;
  }
  else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
    iMesh_initEntIter((iMesh_Instance)ifaceInstances[iface_no], set,
                      iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES,
                      iter, &result);
    PROCESS_MERROR;
  }
  else
    ERRORR(iBase_NOT_SUPPORTED, "Interface should be geometry or mesh.");

  RETURNR(iBase_SUCCESS);
}

int AssocPairC::get_all_entities(const int iface_no,
                                 const int dimension,
                                 iBase_EntityHandle **entities,
                                 int *entities_alloc,
                                 int *entities_size)
{
  int result;

  if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
    iGeom_getEntities((iGeom_Instance)ifaceInstances[iface_no],
                      0, iBase_ALL_TYPES,
                      entities, entities_alloc, entities_size, &result);
    PROCESS_GERROR;
  }
  else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
    iMesh_getEntities((iMesh_Instance)ifaceInstances[iface_no],
                      0, iBase_ALL_TYPES,
                      iMesh_ALL_TOPOLOGIES,
                      entities, entities_alloc, entities_size, &result);
    PROCESS_MERROR;
  }
  else
    ERRORR(iBase_NOT_SUPPORTED, "Interface should be geometry or mesh.");

  RETURNR(iBase_SUCCESS);
}

int AssocPairC::get_all_sets(const int iface_no,
                             iBase_EntitySetHandle **entities,
                             int *entities_alloc,
                             int *entities_size)
{
  int result;

  if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
    iGeom_getEntSets((iGeom_Instance)ifaceInstances[iface_no], 0, 0,
                     reinterpret_cast<iBase_EntitySetHandle**>(entities),
                     entities_alloc, entities_size, &result);
    PROCESS_GERROR;
  }
  else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
    iMesh_getEntSets((iMesh_Instance)ifaceInstances[iface_no], 0, 0,
                     reinterpret_cast<iBase_EntitySetHandle**>(entities),
                     entities_alloc, entities_size, &result);
    PROCESS_MERROR;
  }
  else
    ERRORR(iBase_NOT_SUPPORTED, "Interface should be geometry or mesh.");

  RETURNR(iBase_SUCCESS);
}

int AssocPairC::get_entities(const int iface_no,
                             const int dimension,
                             iBase_EntitySetHandle set_handle,
                             iBase_EntityHandle **entities,
                             int *entities_alloc,
                             int *entities_size)
{
  int result;
  int this_type = (dimension == -1 ? iBase_ALL_TYPES : dimension);

  if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
    iGeom_getEntities((iGeom_Instance)ifaceInstances[iface_no],
                      set_handle, this_type,
                      entities, entities_alloc, entities_size,
                      &result);
    PROCESS_GERROR;
  }
  else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
    iMesh_getEntities((iMesh_Instance)ifaceInstances[iface_no],
                      set_handle, this_type, iMesh_ALL_TOPOLOGIES,
                      entities, entities_alloc, entities_size, &result);
    PROCESS_MERROR;
  }
  else
    ERRORR(iBase_NOT_SUPPORTED, "Interface should be geometry or mesh.");

  RETURNR(iBase_SUCCESS);
}

int AssocPairC::get_ents_dims(const int iface_no,
                              iBase_EntityHandle *entities,
                              int entities_size,
                              int **ent_types,
                              int *ent_types_alloc,
                              int *ent_types_size)
{
  int result;

  if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
    iGeom_getArrType((iGeom_Instance)ifaceInstances[iface_no],
                     entities, entities_size,
                     ent_types, ent_types_alloc, ent_types_size, &result);
    PROCESS_GERROR;
  }
  else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
    iMesh_getEntArrType((iMesh_Instance)ifaceInstances[iface_no],
                        entities, entities_size,
                        ent_types, ent_types_alloc, ent_types_size, &result);
    PROCESS_MERROR;
  }
  else
    ERRORR(iBase_NOT_SUPPORTED, "Interface should be geometry or mesh.");

  RETURNR(iBase_SUCCESS);
}

int AssocPairC::tag_get_handle(const int iface_no,
                               const char *tag_name,
                               const int tag_size_values,
                               const int tag_data_type,
                               const bool create_if_missing,
                               void *default_val,
                               iBase_TagHandle *tag_handle)
{
  int result;

  if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
    iGeom_getTagHandle((iGeom_Instance)ifaceInstances[iface_no],
                       tag_name, tag_handle, &result, strlen(tag_name));
    if (result == iBase_TAG_NOT_FOUND && create_if_missing) {
      iGeom_createTag((iGeom_Instance)ifaceInstances[iface_no],
                      tag_name, tag_size_values, tag_data_type, tag_handle,
                      &result, strlen(tag_name));
    }

    PROCESS_GERROR;
  }
  else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
    iMesh_getTagHandle((iMesh_Instance)ifaceInstances[iface_no], tag_name,
                       tag_handle, &result, strlen(tag_name));
    if (result == iBase_TAG_NOT_FOUND && create_if_missing) {
      iMesh_createTag((iMesh_Instance)ifaceInstances[iface_no],
                      tag_name, tag_size_values, tag_data_type, tag_handle,
                      &result, strlen(tag_name));
    }

    PROCESS_MERROR;
  }
  else
    ERRORR(iBase_NOT_SUPPORTED, "Interface should be geometry or mesh.");

  RETURNR(iBase_SUCCESS);
}

int AssocPairC::tag_destroy(const int iface_no,
                            iBase_TagHandle tag_handle)
{
  int result;

  if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
    iGeom_destroyTag((iGeom_Instance)ifaceInstances[iface_no],
                     tag_handle, true, &result);
    PROCESS_GERROR;
  }
  else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
    iMesh_destroyTag((iMesh_Instance)ifaceInstances[iface_no],
                     tag_handle, true, &result);
    PROCESS_MERROR;
  }

  RETURNR(iBase_SUCCESS);
}
