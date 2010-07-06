#include "sidlArray.h"
#include "iBase.hh"
#include "iGeom_CGM.h"
#include "iGeom_EntityType_IOR.h"

#include "iBase.hh"
#include "iGeom.hh"
#include "iMesh.hh"
#include "iRel.hh"

#include "AssocPairSIDL.hpp"

#include <iostream>

#define RETURN(a) {iRel_LAST_ERROR.error_type = a; return a;}
#define iRel_processError(a, b) {sprintf(iRel_LAST_ERROR.description, b); iRel_LAST_ERROR.error_type = (int)a;}

AssocPairSIDL::AssocPairSIDL(::sidl::BaseInterface iface0,
                             const int ent_or_set0,
                             ::sidl::BaseInterface iface1,
                             const int ent_or_set1,
                             Lasso *lasso) 
  : AssocPair(ent_or_set0, get_type(iface0), ent_or_set1, get_type(iface1), 
              lasso)
{
  // create the tags we'll need
  create_tags();
}

AssocPairSIDL::~AssocPairSIDL() 
{}

IfaceType AssocPairSIDL::get_type(::sidl::BaseInterface iface)
{
  iGeom::Geometry geom = iface;
  if (!geom._is_nil()) return iGeom_IFACE;

  iMesh::Mesh mesh = iface;
  if (!mesh._is_nil()) return iMesh_IFACE;
  
  iRel::Associate assoc = iface;
  if (!assoc._is_nil()) return iRel_IFACE;

  iBase::Error err;
  err.set(iBase::ErrorType_INVALID_ARGUMENT,
          "Couldn't find proper type of interface.");
  throw err;
}

bool AssocPairSIDL::same_interface(iBase_Instance iface0, 
                                   iBase_Instance iface1) 
{
  ::sidl::BaseInterface *iface0_tmp = (::sidl::BaseInterface*) iface0, 
      *iface1_tmp = (::sidl::BaseInterface*) iface1;
  return (iface0_tmp == iface1_tmp);
}

bool AssocPairSIDL::equivalent(iBase_Instance iface0, iBase_Instance iface1,
                               bool *order_switched) 
{
  bool equiv = false;
  ::sidl::BaseInterface *iface0_tmp = (::sidl::BaseInterface*) iface0, 
      *iface1_tmp = (::sidl::BaseInterface*) iface1;
  if (iface0_tmp->isSame(ifaceInstances[0]) && 
      iface1_tmp->isSame(ifaceInstances[1])) {
    if (NULL != order_switched) *order_switched = false;
    equiv = true;
  }
  else if (iface0_tmp->isSame(ifaceInstances[1]) && 
           iface1_tmp->isSame(ifaceInstances[0])) {
    equiv = true;
    if (NULL != order_switched) *order_switched = true;
  }

  return equiv;
}
  
bool AssocPairSIDL::contains(iBase_Instance iface) 
{
  ::sidl::BaseInterface *iface0_tmp = (::sidl::BaseInterface*) iface;
  return (iface0_tmp->isSame(ifaceInstances[0]) || 
          iface0_tmp->isSame(ifaceInstances[1]));
}

int AssocPairSIDL::get_tags(const int iface_no,
                            iBase_EntityHandle *entities,
                            const int num_entities,
                            iBase_TagHandle tag_handle,
                            void *tag_values,
                            const int tag_size) 
{
  int result = iBase_SUCCESS;
  
  CAST_iBase_INTERFACE(ifaceInstances[iface_no], iface_atag, ArrTag,
                       iBase_FAILURE);
  sidl::array<char> tag_values_tmp = 
    convert_to_sidl_vector(tag_values, num_entities*tag_size);
    
  try {
    iface_atag.getArrData(convert_to_sidl_vector(entities, num_entities),
                          num_entities, tag_handle, tag_values_tmp,
                          num_entities);
  }
  catch (iBase::Error err) {
    if (iBase::ErrorType_SUCCESS != err.getErrorType() &&
        iBase::ErrorType_TAG_NOT_FOUND != err.getErrorType()) {
      iRel_LAST_ERROR.error_type = iBase_FAILURE;
      sprintf(iRel_LAST_ERROR.description, "%s", err.getDescription().c_str());
      result = iBase_FAILURE;
    }
  }

  RETURN(iBase_SUCCESS);
}

int AssocPairSIDL::get_tags(const int iface_no,
                            iBase_EntitySetHandle *entities,
                            const int num_entities,
                            iBase_TagHandle tag_handle,
                            void *tag_values,
                            const int tag_size) 
{
  int result = iBase_SUCCESS;
  
  CAST_iBase_INTERFACE(ifaceInstances[iface_no], iface_stag, SetTag,
                         iBase_FAILURE);
  try {
    for (int i = 0; i < num_entities; i++) {
      // TODO: this is certainly wrong
      tag_values[i] = iface_stag.getEntSetData(entities[i], tag_handle);
    }
  }
  catch (iBase::Error err) {
    if (iBase::ErrorType_SUCCESS != err.getErrorType() &&
        iBase::ErrorType_TAG_NOT_FOUND != err.getErrorType()) {
      iRel_LAST_ERROR.error_type = iBase_FAILURE;
      sprintf(iRel_LAST_ERROR.description, "%s", err.getDescription().c_str());
      result = iBase_FAILURE;
    }
  }

  RETURN(iBase_SUCCESS);
}
  
int AssocPairSIDL::set_tags(const int iface_no,
                            iBase_EntityHandle *entities,
                            const int num_entities,
                            iBase_TagHandle tag_handle,
                            const void *tag_values,
                            const int tag_size)
{
  try {
    CAST_iBase_INTERFACE(ifaceInstances[iface_no], iface_atag, ArrTag,
                         iBase_FAILURE);
    sidl::array<char> tag_values_tmp = 
      convert_to_sidl_vector(tag_values, num_entities * tag_size);

    iface_atag.setArrData(convert_to_sidl_vector(entities, num_entities),
                          num_entities, tag_handle,
                          tag_values_tmp, num_entities);
  }
  catch (iBase::Error err) {
    iRel_LAST_ERROR.error_type = iBase_FAILURE;
    sprintf(iRel_LAST_ERROR.description, "%s", err.getDescription().c_str());
    RETURN(iBase_FAILURE);
  }

  RETURN(iBase_SUCCESS);
}
  
int AssocPairSIDL::set_tags(const int iface_no,
                            iBase_EntitySetHandle *sets,
                            const int num_sets,
                            iBase_TagHandle tag_handle,
                            const void *tag_values,
                            const int tag_size)
{
  try {
    CAST_iBase_INTERFACE(ifaceInstances[iface_no], iface_stag, SetTag,
                         iBase_FAILURE);
    for (int i = 0; i < num_sets; i++) {
      // TODO: this is certainly wrong
      iface_stag.setEntSetData(sets[i], tag_handle, tag_values[i]);
    }
  }
  catch (iBase::Error err) {
    iRel_LAST_ERROR.error_type = iBase_FAILURE;
    sprintf(iRel_LAST_ERROR.description, "%s", err.getDescription().c_str());
    RETURN(iBase_FAILURE);
  }

  RETURN(iBase_SUCCESS);
}
  
int AssocPairSIDL::get_all_entities(const int iface_no,
                                    const int dimension,
                                    iBase_EntityHandle **entities,
                                    int *entities_alloc,
                                    int *entities_size) 
{
  sidl::array<iBase_EntityHandle> entities_tmp;
  int entities_tmp_size = 0;

  if (iface_type(iface_no) == iGeom_IFACE) {
    CAST_iGeom_INTERFACE(ifaceInstances[iface_no], geom_topo, Topology,
                         iBase_FAILURE);
    geom_topo.getEntities(0, iGeom::EntityType_ALL_TYPES, 
                          entities_tmp, entities_tmp_size);
  }
  else if (iface_type(iface_no) == iMesh_IFACE) {
    CAST_iMesh_INTERFACE(ifaceInstances[iface_no], mesh, Mesh,
                         iBase_FAILURE);
    mesh.getEntities(0, iMesh::EntityType_ALL_TYPES, 
                     iMesh::EntityTopology_ALL_TOPOLOGIES, 
                     entities_tmp, entities_tmp_size);
  }
  else RETURN(iBase_FAILURE);

  CHECK_SIZE_VOID(iBase_EntityHandle, entities, entities_alloc,
                  entities_tmp_size);
  iBase_EntityHandle *eh = ARRAY_PTR(entities_tmp, iBase_EntityHandle);
  std::copy(eh, eh+entities_tmp_size, *entities);
  *entities_size = entities_tmp_size;

  RETURN(iBase_SUCCESS);
}
int AssocPairSIDL::get_all_entities(const int iface_no,
                                    iBase_EntitySetHandle **entities,
                                    int *entities_alloc,
                                    int *entities_size) 
{
  sidl::array<iBase_EntityHandle> entities_tmp;
  int entities_tmp_size = 0;

    CAST_iBase_INTERFACE(ifaceInstances[iface_no], eset, EntSet,
                         iBase_FAILURE);
    eset.getEntSets(0, 0, entities_tmp, entities_tmp_size);

  CHECK_SIZE_VOID(iBase_EntityHandle, entities, entities_alloc,
                  entities_tmp_size);
  iBase_EntityHandle *eh = ARRAY_PTR(entities_tmp, iBase_EntityHandle);
  std::copy(eh, eh+entities_tmp_size, *entities);
  *entities_size = entities_tmp_size;

  RETURN(iBase_SUCCESS);
}

int AssocPairSIDL::get_entities(const int iface_no,
                                const int dimension,
                                iBase_EntityHandle set_handle,
                                iBase_EntityHandle **entities,
                                int *entities_allocated,
                                int *entities_size) 
{
  sidl::array<iBase_EntityHandle> entities_tmp;
  int entities_tmp_size = 0;
  
  if (iface_type(iface_no) == iGeom_IFACE) {
    iGeom::EntityType this_type =
      (-1 == dimension ? iGeom::EntityType_ALL_TYPES : 
       (iGeom::EntityType)dimension);

    CAST_iGeom_INTERFACE(ifaceInstances[iface_no], geom_topo, Topology,
                         iBase_FAILURE);
    geom_topo.getEntities(set_handle, this_type,
                          entities_tmp, entities_tmp_size);
  }
  else if (iface_type(iface_no) == iMesh_IFACE) {
    iMesh::EntityType this_type =
      (-1 == dimension ? iMesh::EntityType_ALL_TYPES : 
       (iMesh::EntityType)dimension);

    CAST_iMesh_INTERFACE(ifaceInstances[iface_no], mesh, Mesh,
                         iBase_FAILURE);
    mesh.getEntities(set_handle, this_type,
                     iMesh::EntityTopology_ALL_TOPOLOGIES, 
                     entities_tmp, entities_tmp_size);
  }
  else RETURN(iBase_FAILURE);

  CHECK_SIZE_VOID(iBase_EntityHandle, entities, entities_allocated,
                  entities_tmp_size);
  iBase_EntityHandle *eh = ARRAY_PTR(entities_tmp, iBase_EntityHandle);
  std::copy(eh, eh+entities_tmp_size, *entities);
  *entities_size = entities_tmp_size;

  RETURN(iBase_SUCCESS);
}


int AssocPairSIDL::get_ents_dims(const int iface_no,
                                 iBase_EntityHandle *entities,
                                 int entities_size,
                                 int **ent_types,
                                 int *ent_types_alloc,
                                 int *ent_types_size) 
{
  iRel_LAST_ERROR.error_type = iBase_SUCCESS;
  
  if (iGeom_IFACE != iface_type(iface_no) &&
      iMesh_IFACE != iface_type(iface_no)) {
    iRel_processError(iBase_NOT_SUPPORTED, 
                       "Interface should be geometry or mesh.");
    RETURN(iBase_NOT_SUPPORTED);
  }

  if (iface_type(iface_no) == iGeom_IFACE) {
    try {
      CAST_iGeom_INTERFACE(ifaceInstances[iface_no], geom_topo, Topology,
                           iBase_FAILURE);
      sidl::array<iGeom::EntityType> ent_types_tmp;
      int ent_types_tmp_size = 0;
      geom_topo.getArrType(convert_to_sidl_vector(entities, 
                                                  entities_size), 
                           entities_size, 
                           ent_types_tmp, ent_types_tmp_size);
      CHECK_SIZE_VOID(int, ent_types, ent_types_alloc,
                      ent_types_tmp_size);
      for (int i = 0; i < ent_types_tmp_size; i++)
        (*ent_types)[i] = (int) ent_types_tmp.get(i);
      *ent_types_size = ent_types_tmp_size;
    }
    catch (iBase::Error err) {
      iRel_processError((int)err.getErrorType(),
                         err.getDescription().c_str());
      RETURN(iRel_LAST_ERROR.error_type);
    }
  }
  else if (iface_type(iface_no) == iMesh_IFACE) {
    try {
      CAST_iMesh_INTERFACE(ifaceInstances[iface_no], mesh_arr, Arr,
                           iBase_FAILURE);
      sidl::array<iMesh::EntityType> ent_types_tmp;
      int ent_types_tmp_size = 0;
      mesh_arr.getEntArrType(convert_to_sidl_vector(entities, 
                                                    entities_size), 
                             entities_size, 
                             ent_types_tmp, ent_types_tmp_size);
      CHECK_SIZE_VOID(int, ent_types, ent_types_alloc,
                      ent_types_tmp_size);
      for (int i = 0; i < ent_types_tmp_size; i++)
        (*ent_types)[i] = (int) ent_types_tmp.get(i);
      *ent_types_size = ent_types_tmp_size;
    }
    catch (iBase::Error err) {
      iRel_processError((int)err.getErrorType(),
                         err.getDescription().c_str());
      RETURN(iRel_LAST_ERROR.error_type);
    }
  }
  
  RETURN(iBase_SUCCESS);
}

iBase_TagHandle AssocPairSIDL::tag_get_handle(const int iface_no,
                                              const char *tag_name,
                                              const int tag_size_values,
                                              const int tag_data_type,
                                              const bool create_if_missing,
                                              void *default_val) 
{
  CAST_iBase_INTERFACE(ifaceInstances[iface_no], iface_tag, Tag,
                       NULL);
  iBase_TagHandle this_tag = 0;
  try {
    this_tag = iface_tag.getTagHandle(tag_name);
  }
  catch (iBase::Error err) {
    if (iBase::ErrorType_TAG_NOT_FOUND != err.getErrorType()) return 0;
  }
  
  if (0 != this_tag || !create_if_missing) return this_tag;
  
  try {
    iface_tag.createTag(tag_name, tag_size_values, 
                        (iBase::TagValueType)tag_data_type,
                        this_tag);
  }
  catch (iBase::Error err) {
    iRel_processError(err.getErrorType(), err.getDescription().c_str());
    return 0;
  }
  
  return this_tag;
}
