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
          iGeom_getDescription((iGeom_Instance)ifaceInstances[iface_no], this_descr, \
                               &result, 120);\
          iRel_processError(result, this_descr);}
#define PROCESS_MERROR  if (iBase_SUCCESS != result) {\
          char this_descr[120];\
          iMesh_getDescription((iMesh_Instance)ifaceInstances[iface_no], this_descr, \
                               &result, 120);\
          iRel_processError(result, this_descr);}


AssocPairC::AssocPairC(iBase_Instance iface0,
                       const int ent_or_set0,
                       const IfaceType type0,
                       iBase_Instance iface1,
                       const int ent_or_set1,
                       const IfaceType type1,
                       Lasso *lasso) 
    : AssocPair(ent_or_set0, ent_or_set1, lasso)
{
  if (type0 < type1) {
    ifaceInstances[0] = iface0;
    ifaceInstances[1] = iface1;
    ifaceTypes[0] = type0;
    ifaceTypes[1] = type1;    
  }
  else {
    ifaceInstances[0] = iface1;
    ifaceInstances[1] = iface0;
    ifaceTypes[0] = type1;
    ifaceTypes[1] = type0;
  }

    // finally, create the tags we'll need
  create_tags();
}

AssocPairC::~AssocPairC() 
{
    // need to destroy tags for this assoc pair
  for (int iface_no = 0; iface_no < 2; iface_no++) {
    int result;
    if (iface_type(iface_no) == iRel_IGEOM_IFACE) {
      iGeom_destroyTag((iGeom_Instance)ifaceInstances[iface_no],
                       assocTags[iface_no], true, &result);
      PROCESS_GERROR;
    }
    else if (iface_type(iface_no) == iRel_IMESH_IFACE) {
      iMesh_destroyTag((iMesh_Instance)ifaceInstances[iface_no],
                       assocTags[iface_no], true, &result);
      PROCESS_MERROR;
    }
    assocTags[iface_no] = 0;
  }
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

int AssocPairC::get_int_tags(const int iface_no,
                             iBase_EntityHandle *entities,
                             const int num_entities,
                             iBase_TagHandle tag_handle,
                             int *tag_values) 
{
  int tag_values_alloc = num_entities, tag_values_size;
  int result;
  
  if (iRel_IGEOM_IFACE != ifaceTypes[iface_no] &&
      iRel_IMESH_IFACE != ifaceTypes[iface_no]) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  if (iRel_IGEOM_IFACE == ifaceTypes[iface_no]) {
    iGeom_getIntArrData((iGeom_Instance)ifaceInstances[iface_no],
                        entities,
                        num_entities, tag_handle,
                        &tag_values, &tag_values_alloc, 
                        &tag_values_size, &result);
    PROCESS_GERROR;
  }
  else if (iRel_IMESH_IFACE == ifaceTypes[iface_no]) {
    iMesh_getIntArrData((iMesh_Instance)ifaceInstances[iface_no],
                        entities,
                        num_entities, tag_handle,
                        &tag_values, &tag_values_alloc,
                        &tag_values_size, &result);
    PROCESS_MERROR;
  }

  return iRel_LAST_ERROR.error_type;
}

int AssocPairC::get_int_tags(const int iface_no,
                             iBase_EntitySetHandle *entities,
                             const int num_entities,
                             iBase_TagHandle tag_handle,
                             int *tag_values) 
{
  int tag_values_alloc = num_entities, tag_values_size;
  int result;
  
  if (iRel_IGEOM_IFACE != ifaceTypes[iface_no] &&
      iRel_IMESH_IFACE != ifaceTypes[iface_no]) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  if (iRel_IGEOM_IFACE == ifaceTypes[iface_no]) {
    for (int i = 0; i < num_entities; i++) {
      iGeom_getEntSetIntData((iGeom_Instance)ifaceInstances[iface_no],
                             entities[i],
                             tag_handle, tag_values+i, &result);
      PROCESS_GERROR;
    }
  }
  else if (iRel_IMESH_IFACE == ifaceTypes[iface_no]) {
    for (int i = 0; i < num_entities; i++) {
      iMesh_getEntSetIntData((iMesh_Instance)ifaceInstances[iface_no],
                             entities[i],
                             tag_handle, tag_values+i, &result);
      PROCESS_MERROR;
    }
  }

  return iRel_LAST_ERROR.error_type;
}
 
int AssocPairC::get_eh_tags(const int iface_no,
                            iBase_EntityHandle *entities,
                            const int num_entities,
                            iBase_TagHandle tag_handle,
                            iBase_EntityHandle *tag_values) 
{
  int tag_values_alloc = num_entities, tag_values_size;
  int result;
  
  if (iRel_IGEOM_IFACE != ifaceTypes[iface_no] &&
      iRel_IMESH_IFACE != ifaceTypes[iface_no]) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  if (iRel_IGEOM_IFACE == ifaceTypes[iface_no]) {
    iGeom_getEHArrData((iGeom_Instance)ifaceInstances[iface_no],
                       entities,
                       num_entities, tag_handle,
                       &tag_values, &tag_values_alloc,
                       &tag_values_size, &result);
    PROCESS_GERROR;
  }
  else if (iRel_IMESH_IFACE == ifaceTypes[iface_no]) {
    iMesh_getEHArrData((iMesh_Instance)ifaceInstances[iface_no],
                       entities,
                       num_entities, tag_handle,
                       &tag_values, &tag_values_alloc, 
                       &tag_values_size, &result);
    PROCESS_MERROR;
  }

  return iRel_LAST_ERROR.error_type;
}
 
int AssocPairC::get_eh_tags(const int iface_no,
                            iBase_EntitySetHandle *entities,
                            const int num_entities,
                            iBase_TagHandle tag_handle,
                            iBase_EntityHandle *tag_values) 
{
  int tag_values_alloc = num_entities, tag_values_size;
  int result;
  
  if (iRel_IGEOM_IFACE != ifaceTypes[iface_no] &&
      iRel_IMESH_IFACE != ifaceTypes[iface_no]) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  if (iRel_IGEOM_IFACE == ifaceTypes[iface_no]) {
    for (int i = 0; i < num_entities; i++) {
      iGeom_getEntSetEHData((iGeom_Instance)ifaceInstances[iface_no],
                            entities[i],
                            tag_handle, tag_values+i, &result);
      PROCESS_GERROR;
    }
  }
  else if (iRel_IMESH_IFACE == ifaceTypes[iface_no]) {
    for (int i = 0; i < num_entities; i++) {
      iMesh_getEntSetEHData((iMesh_Instance)ifaceInstances[iface_no],
                            entities[i],
                            tag_handle, tag_values+i, &result);
      PROCESS_MERROR;
    }
  }

  return iRel_LAST_ERROR.error_type;
}
  
int AssocPairC::set_int_tags(const int iface_no,
                             iBase_EntityHandle *entities,
                             const int num_entities,
                             iBase_TagHandle tag_handle,
                             int *tag_values)
{
  int result;

  if (iRel_IGEOM_IFACE != ifaceTypes[iface_no] &&
      iRel_IMESH_IFACE != ifaceTypes[iface_no]) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  if (iRel_IGEOM_IFACE == ifaceTypes[iface_no]) {
    iGeom_setIntArrData((iGeom_Instance)ifaceInstances[iface_no],
                        entities,
                        num_entities, tag_handle,
                        tag_values, num_entities, &result);
    PROCESS_GERROR;
  }
  else if (iRel_IMESH_IFACE == ifaceTypes[iface_no]) {
    iMesh_setIntArrData((iMesh_Instance)ifaceInstances[iface_no],
                        entities,
                        num_entities, tag_handle,
                        tag_values, num_entities, &result);
    PROCESS_MERROR;
  }

  return iRel_LAST_ERROR.error_type;
}
   
int AssocPairC::set_int_tags(const int iface_no,
                             iBase_EntitySetHandle *entities,
                             const int num_entities,
                             iBase_TagHandle tag_handle,
                             int *tag_values)
{
  int result;

  if (iRel_IGEOM_IFACE != ifaceTypes[iface_no] &&
      iRel_IMESH_IFACE != ifaceTypes[iface_no]) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  if (iRel_IGEOM_IFACE == ifaceTypes[iface_no]) {
    for (int i = 0; i < num_entities; i++) {
      iGeom_setEntSetIntData((iGeom_Instance)ifaceInstances[iface_no],
                             entities[i],
                             tag_handle, tag_values[i], &result);
      PROCESS_GERROR;
    }
  }
  else if (iRel_IMESH_IFACE == ifaceTypes[iface_no]) {
    for (int i = 0; i < num_entities; i++) {
      iMesh_setEntSetIntData((iMesh_Instance)ifaceInstances[iface_no],
                             entities[i],
                             tag_handle, tag_values[i], &result);
      PROCESS_MERROR;
    }
  }

  return iRel_LAST_ERROR.error_type;
}
 
int AssocPairC::set_eh_tags(const int iface_no,
                            iBase_EntityHandle *entities,
                            const int num_entities,
                            iBase_TagHandle tag_handle,
                            iBase_EntityHandle *tag_values)
{
  int result;

  if (iRel_IGEOM_IFACE != ifaceTypes[iface_no] &&
      iRel_IMESH_IFACE != ifaceTypes[iface_no]) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  if (iRel_IGEOM_IFACE == ifaceTypes[iface_no]) {
    iGeom_setEHArrData((iGeom_Instance)ifaceInstances[iface_no],
                       entities,
                       num_entities, tag_handle,
                       tag_values, num_entities, &result);
    PROCESS_GERROR;
  }
  else if (iRel_IMESH_IFACE == ifaceTypes[iface_no]) {
    iMesh_setEHArrData((iMesh_Instance)ifaceInstances[iface_no],
                       entities,
                       num_entities, tag_handle,
                       tag_values, num_entities, &result);
    PROCESS_MERROR;
  }

  return iRel_LAST_ERROR.error_type;
}
 
int AssocPairC::set_eh_tags(const int iface_no,
                            iBase_EntitySetHandle *entities,
                            const int num_entities,
                            iBase_TagHandle tag_handle,
                            iBase_EntityHandle *tag_values)
{
  int result;

  if (iRel_IGEOM_IFACE != ifaceTypes[iface_no] &&
      iRel_IMESH_IFACE != ifaceTypes[iface_no]) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                      "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }
  
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  if (iRel_IGEOM_IFACE == ifaceTypes[iface_no]) {
    for (int i = 0; i < num_entities; i++) {
      iGeom_setEntSetEHData((iGeom_Instance)ifaceInstances[iface_no],
                            entities[i],
                            tag_handle, tag_values[i], &result);
      PROCESS_GERROR;
    }
  }
  else if (iRel_IMESH_IFACE == ifaceTypes[iface_no]) {
    for (int i = 0; i < num_entities; i++) {
      iMesh_setEntSetEHData((iMesh_Instance)ifaceInstances[iface_no],
                            entities[i],
                            tag_handle, tag_values[i], &result);
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

  if (iRel_IGEOM_IFACE != ifaceTypes[iface_no] &&
      iRel_IMESH_IFACE != ifaceTypes[iface_no]) {
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

  if (iRel_IGEOM_IFACE != ifaceTypes[iface_no] &&
      iRel_IMESH_IFACE != ifaceTypes[iface_no]) {
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

  if (iRel_IGEOM_IFACE != ifaceTypes[iface_no] &&
      iRel_IMESH_IFACE != ifaceTypes[iface_no]) {
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

  if (iRel_IGEOM_IFACE != ifaceTypes[iface_no] &&
      iRel_IMESH_IFACE != ifaceTypes[iface_no]) {
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
  
  if (iRel_IGEOM_IFACE != ifaceTypes[iface_no] &&
      iRel_IMESH_IFACE != ifaceTypes[iface_no]) {
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

int AssocPairC::create_mesh_vertex(const double x,
                                   const double y,
                                   const double z,
                                   iBase_EntityHandle &vertex) 
{
  int iface_no;
  if (ifaceTypes[0] == iRel_IMESH_IFACE) iface_no = 0;
  else if (ifaceTypes[1] == iRel_IMESH_IFACE) iface_no = 1;
  else {
    iRel_processError(iBase_INVALID_ARGUMENT, 
                      "One of the interfaces must be mesh.");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  
  int result;
  iMesh_createVtx((iMesh_Instance)ifaceInstances[iface_no],
                  x, y, z, &vertex, &result);
  RETURN(result);
}

int AssocPairC::create_mesh_entity(iMesh_EntityTopology ent_topology,
                                   iBase_EntityHandle *lower_handles,
                                   int lower_handles_size,
                                   iBase_EntityHandle &new_ent,
                                   int &status)
{
  int iface_no;
  if (ifaceTypes[0] == iRel_IMESH_IFACE) iface_no = 0;
  else if (ifaceTypes[1] == iRel_IMESH_IFACE) iface_no = 1;
  else {
    iRel_processError(iBase_INVALID_ARGUMENT, 
                      "One of the interfaces must be mesh.");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  
  int result;
  iMesh_createEnt((iMesh_Instance)ifaceInstances[iface_no],
                  ent_topology, lower_handles, 
                  lower_handles_size, 
                  &new_ent, &status, &result);
  RETURN(result);
}

int AssocPairC::set_mesh_coords(iBase_EntityHandle *verts,
                                int verts_size,
                                double *coords,
                                iBase_StorageOrder order) 
{
  int iface_no;
  if (iRel_IMESH_IFACE == ifaceTypes[0]) iface_no = 0;
  else if (iRel_IMESH_IFACE == ifaceTypes[1]) iface_no = 1;
  else {
    iRel_processError(iBase_INVALID_ARGUMENT, 
                      "One of the interfaces must be mesh.");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  
  int result;
  iMesh_setVtxArrCoords((iMesh_Instance)ifaceInstances[iface_no],
                        verts, verts_size,
                        order, coords, 3*verts_size, &result);
  RETURN(result);
}

int AssocPairC::get_mesh_coords(iBase_EntityHandle *verts,
                                int verts_size,
                                double **coords,
                                int *coords_alloc,
                                int *coords_size,
                                iBase_StorageOrder order) 
{
  int iface_no;
  if (iRel_IMESH_IFACE == ifaceTypes[0]) iface_no = 0;
  else if (iRel_IMESH_IFACE == ifaceTypes[1]) iface_no = 1;
  else {
    iRel_processError(iBase_INVALID_ARGUMENT, 
                      "One of the interfaces must be mesh.");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  
  int result;
  iMesh_getVtxArrCoords((iMesh_Instance)ifaceInstances[iface_no],
                        verts, verts_size, order,
                        coords, coords_alloc,
                        coords_size, &result);
  RETURN(result);
}

int AssocPairC::get_closest_pt(iBase_EntityHandle *gents, 
                               int gents_size,
                               iBase_StorageOrder storage_order,
                               double *near_coords, 
                               int near_coords_size,
                               double **on_coords,
                               int *on_coords_alloc, 
                               int *on_coords_size) 
{
  int iface_no;
  if (iRel_IGEOM_IFACE == ifaceTypes[0]) iface_no = 0;
  else if (iRel_IGEOM_IFACE == ifaceTypes[1]) iface_no = 1;
  else {
    iRel_processError(iBase_INVALID_ARGUMENT, 
                      "One of the interfaces must be geom.");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  
  int result;
  iGeom_getArrClosestPt((iGeom_Instance)ifaceInstances[iface_no],
                        gents, gents_size,
                        storage_order,
                        near_coords, near_coords_size,
                        on_coords, on_coords_alloc,
                        on_coords_size, &result);
  RETURN(result);
}

