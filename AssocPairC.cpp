#include "AssocPairC.hpp"

#include "iGeom.h"
#include "iMesh.h"
#include "iRel.h"

#include <string.h>
#include <stdio.h>

#define RETURN(a) {iRel_LAST_ERROR.error_type = a; return a;}
#define iRel_processError(a, b) {sprintf(iRel_LAST_ERROR.description, "%s", b); iRel_LAST_ERROR.error_type = a;}

#define PROCESS_GERROR  if (iBase_SUCCESS != result) {\
          char this_descr[120];\
          int err;\
          iGeom_getDescription((iGeom_Instance)ifaceInstances[iface_no], this_descr, \
                               &err, 120);\
          iRel_processError(result, this_descr);}
#define PROCESS_MERROR  if (iBase_SUCCESS != result) {\
          char this_descr[120];\
          int err;\
          iMesh_getDescription((iMesh_Instance)ifaceInstances[iface_no], this_descr, \
                               &err, 120);\
          iRel_processError(result, this_descr);}


AssocPairC::AssocPairC(iBase_Instance iface0,
                       RelationType ent_or_set0,
                       IfaceType type0,
                       iBase_Instance iface1,
                       RelationType ent_or_set1,
                       IfaceType type1)
  : AssocPair(ent_or_set0, type0, ent_or_set1, type1)
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
  
  if (iRel_IGEOM_IFACE != iface_type(iface_no) &&
      iRel_IMESH_IFACE != iface_type(iface_no)) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  if (iRel_IGEOM_IFACE == iface_type(iface_no)) {
    iGeom_getArrData((iGeom_Instance)ifaceInstances[iface_no],
                     entities, num_entities, tag_handle,
                     reinterpret_cast<char**>(&tag_values), &tag_values_alloc, 
                     &tag_values_size, &result);
    PROCESS_GERROR;
  }
  else if (iRel_IMESH_IFACE == iface_type(iface_no)) {
    iMesh_getArrData((iMesh_Instance)ifaceInstances[iface_no],
                     entities, num_entities, tag_handle,
                     reinterpret_cast<char**>(&tag_values), &tag_values_alloc,
                     &tag_values_size, &result);
    PROCESS_MERROR;
  }

  return iRel_LAST_ERROR.error_type;
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
  
  if (iRel_IGEOM_IFACE != iface_type(iface_no) &&
      iRel_IMESH_IFACE != iface_type(iface_no)) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
    for (int i = 0; i < num_sets; i++) {
      iGeom_getEntSetData((iGeom_Instance)ifaceInstances[iface_no],
                          sets[i], tag_handle, &tag_data, &tag_values_alloc,
                          &tag_values_size, &result);
      tag_data += tag_size;
      PROCESS_GERROR;
    }
  }
  else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
    for (int i = 0; i < num_sets; i++) {
      iMesh_getEntSetData((iMesh_Instance)ifaceInstances[iface_no],
                          sets[i], tag_handle, &tag_data, &tag_values_alloc,
                          &tag_values_size, &result);
      tag_data += tag_size;
      PROCESS_MERROR;
    }
  }

  return iRel_LAST_ERROR.error_type;
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

  if (iRel_IGEOM_IFACE != iface_type(iface_no) &&
      iRel_IMESH_IFACE != iface_type(iface_no)) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
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

  return iRel_LAST_ERROR.error_type;
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

  if (iRel_IGEOM_IFACE != iface_type(iface_no) &&
      iRel_IMESH_IFACE != iface_type(iface_no)) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
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

  return iRel_LAST_ERROR.error_type;
}

int AssocPairC::get_all_entities(const int iface_no,
                                 const int dimension,
                                 iBase_EntityHandle **entities,
                                 int *entities_alloc,
                                 int *entities_size) 
{
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  
  int result;

  if (iRel_IGEOM_IFACE != iface_type(iface_no) &&
      iRel_IMESH_IFACE != iface_type(iface_no)) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }

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
  
  RETURN(iRel_LAST_ERROR.error_type);
}

int AssocPairC::get_all_sets(const int iface_no,
                                 iBase_EntitySetHandle **entities,
                                 int *entities_alloc,
                                 int *entities_size) 
{
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  
  int result;

  if (iRel_IGEOM_IFACE != iface_type(iface_no) &&
      iRel_IMESH_IFACE != iface_type(iface_no)) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }

  if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
    iGeom_getEntSets((iGeom_Instance)ifaceInstances[iface_no], 0, 
                       0, reinterpret_cast<iBase_EntitySetHandle**>(entities), 
                       entities_alloc, entities_size, &result);
      
    PROCESS_GERROR;
  }
  else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
    iMesh_getEntSets((iMesh_Instance)ifaceInstances[iface_no], 
                       0, 0, reinterpret_cast<iBase_EntitySetHandle**>(entities), 
                       entities_alloc, entities_size, &result);

    PROCESS_MERROR;
  }
  
  RETURN(iRel_LAST_ERROR.error_type);
}

int AssocPairC::get_entities(const int iface_no,
                             const int dimension,
                             iBase_EntitySetHandle set_handle,
                             iBase_EntityHandle **entities,
                             int *entities_alloc,
                             int *entities_size) 
{
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  
  int result;

  if (iRel_IGEOM_IFACE != iface_type(iface_no) &&
      iRel_IMESH_IFACE != iface_type(iface_no)) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }

  if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
    int this_type =
      (-1 == dimension ? iBase_ALL_TYPES : dimension);
    
    iGeom_getEntities((iGeom_Instance)ifaceInstances[iface_no], 
                      set_handle, this_type,
                      entities, entities_alloc, entities_size,
                      &result);
    PROCESS_GERROR;
  }
  else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
    int this_type =
      (-1 == dimension ? iBase_ALL_TYPES : dimension);

    iMesh_getEntities((iMesh_Instance)ifaceInstances[iface_no], 
                      set_handle, this_type,
                      iMesh_ALL_TOPOLOGIES, 
                      entities, entities_alloc, entities_size, &result);
    PROCESS_MERROR;
  }
  
  RETURN(iRel_LAST_ERROR.error_type);
}

int AssocPairC::get_ents_dims(const int iface_no,
                              iBase_EntityHandle *entities,
                              int entities_size,
                              int **ent_types,
                              int *ent_types_alloc,
                              int *ent_types_size) 
{
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  
  int result;

  if (iRel_IGEOM_IFACE != iface_type(iface_no) &&
      iRel_IMESH_IFACE != iface_type(iface_no)) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }

  if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
    iGeom_getArrType((iGeom_Instance)ifaceInstances[iface_no], 
                     entities, entities_size, 
                     ent_types, ent_types_alloc, ent_types_size, 
                     &result);
    PROCESS_GERROR;
  }
  else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
    iMesh_getEntArrType((iMesh_Instance)ifaceInstances[iface_no], 
                        entities, entities_size, 
                        ent_types, 
                        ent_types_alloc, ent_types_size, &result);
    PROCESS_MERROR;
  }
  
  RETURN(iRel_LAST_ERROR.error_type);
}

iBase_TagHandle AssocPairC::tag_get_handle(const int iface_no,
                                           const char *tag_name,
                                           const int tag_size_values,
                                           const int tag_data_type,
                                           const bool create_if_missing,
                                           void *default_val) 
{
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;

  int result;
  
  if (iRel_IGEOM_IFACE != iface_type(iface_no) &&
      iRel_IMESH_IFACE != iface_type(iface_no)) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    return 0;
  }

  iBase_TagHandle this_tag = 0;
  if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
    iGeom_getTagHandle((iGeom_Instance)ifaceInstances[iface_no],
                       tag_name, &this_tag, &result, strlen(tag_name));
    PROCESS_GERROR;
  }
  else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
    iMesh_getTagHandle((iMesh_Instance)ifaceInstances[iface_no],
                       tag_name, &this_tag, &result, strlen(tag_name));
    PROCESS_MERROR;
  }
  
  if (0 != this_tag || !create_if_missing) return this_tag;
  
  if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
    iGeom_createTag((iGeom_Instance)ifaceInstances[iface_no],
                    tag_name, tag_size_values,
                    tag_data_type, &this_tag,
                    &result, strlen(tag_name));
    PROCESS_GERROR;
  }
  else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
    iMesh_createTag((iMesh_Instance)ifaceInstances[iface_no],
                    tag_name, tag_size_values,
                    tag_data_type, &this_tag, &result, strlen(tag_name));
    PROCESS_MERROR;
  }
  
  return this_tag;
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

  return iBase_SUCCESS;
}
