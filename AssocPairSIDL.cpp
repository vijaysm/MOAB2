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
    : AssocPair(ent_or_set0, ent_or_set1, lasso)
{
    // try casting to figure out which type it is
  if (is_type(iface0, iMesh_IFACE)) {
      // if we get here, it's mesh
    ifaceTypes[0] = iMesh_IFACE;
    ifaceInstances[0] = iface0;
  }
  else if (is_type(iface0, iGeom_IFACE)) {
      // if we get here, it's geom
    ifaceTypes[0] = iGeom_IFACE;
    ifaceInstances[0] = iface0;
  }
  else if (is_type(iface0, iRel_IFACE)) {
      // if we get here, it's associate
    ifaceTypes[0] = iRel_IFACE;
    ifaceInstances[0] = iface0;
  }
  else {
    iBase::Error err;
    err.set(iBase::ErrorType_INVALID_ARGUMENT,
            "Couldn't find proper type of interface for createAssociation.");
    throw err;
  }

    // try casting to figure out which type it is
  if (is_type(iface1, iMesh_IFACE)) {
      // if we get here, it's mesh
    ifaceTypes[1] = iMesh_IFACE;
    ifaceInstances[1] = iface1;
  }
  else if (is_type(iface1, iGeom_IFACE)) {
      // if we get here, it's geom
    ifaceTypes[1] = iGeom_IFACE;
    ifaceInstances[1] = iface1;
  }
  else if (is_type(iface1, iRel_IFACE)) {
      // if we get here, it's associate
    ifaceTypes[1] = iRel_IFACE;
    ifaceInstances[1] = iface1;
  }
  else {
    iBase::Error err;
    err.set(iBase::ErrorType_INVALID_ARGUMENT,
            "Couldn't find proper type of interface for createAssociation.");
    throw err;
  }


    // sort by iface type
#define SWITCH(this_type, arg1, arg2) \
  {this_type arg_tmp = arg1; arg1 = arg2; arg2 = arg_tmp;}

  if (ifaceTypes[0] > ifaceTypes[1]) {
    SWITCH(IfaceType, ifaceTypes[0], ifaceTypes[1]);
    SWITCH(::sidl::BaseInterface, ifaceInstances[0], ifaceInstances[1]);
    SWITCH(int, entOrSet[0], entOrSet[1]);
  }

    // finally, create the tags we'll need
  create_tags();
}

AssocPairSIDL::~AssocPairSIDL() 
{}

bool AssocPairSIDL::is_type(::sidl::BaseInterface iface,
                            IfaceType this_type) 
{
  if (iGeom_IFACE == this_type) {
    iGeom::Geometry geom = iface;
    if (geom._is_nil()) return false;
    else return true;
  }
  
  if (iMesh_IFACE == this_type) {
    iMesh::Mesh mesh = iface;
    if (mesh._is_nil()) return false;
    else return true;
  }
  
  if (iRel_IFACE == this_type) {
    iRel::Associate assoc = iface;
    if (assoc._is_nil()) return false;
    else return true;
  }

  return false;
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

int AssocPairSIDL::get_int_tags(const int iface_no,
                                            iBase_EntityHandle *entities,
                                            int num_entities,
                                            iBase_TagHandle tag_handle,
                                            bool are_sets,
                                            int *tag_values) 
{
  int result = iBase_SUCCESS;
  
  if (are_sets) {
    CAST_iBase_INTERFACE(ifaceInstances[iface_no], iface_stag, SetTag,
                         iBase_FAILURE);
    for (int i = 0; i < num_entities; i++) {
      try {
        tag_values[i] = iface_stag.getEntSetIntData(entities[i],
                                                    tag_handle);
      }
      catch (iBase::Error err) {
        if (iBase::ErrorType_SUCCESS != err.getErrorType() &&
            iBase::ErrorType_TAG_NOT_FOUND != err.getErrorType()) {
          iRel_LAST_ERROR.error_type = iBase_FAILURE;
          sprintf(iRel_LAST_ERROR.description, "%s", err.getDescription().c_str());
          result = iBase_FAILURE;
        }
      }
    }
  }
  else {
    CAST_iBase_INTERFACE(ifaceInstances[iface_no], iface_atag, ArrTag,
                         iBase_FAILURE);
    sidl::array<int> tag_values_tmp = 
      convert_to_sidl_vector(tag_values, num_entities);
    
    try {
      iface_atag.getIntArrData(convert_to_sidl_vector(entities, num_entities),
                               num_entities, tag_handle,
                               tag_values_tmp, num_entities);
    }
    catch (iBase::Error err) {
      if (iBase::ErrorType_SUCCESS != err.getErrorType() &&
          iBase::ErrorType_TAG_NOT_FOUND != err.getErrorType()) {
        iRel_LAST_ERROR.error_type = iBase_FAILURE;
        sprintf(iRel_LAST_ERROR.description, "%s", err.getDescription().c_str());
        result = iBase_FAILURE;
      }
    }
    
  }

  RETURN(iBase_SUCCESS);
}
  
int AssocPairSIDL::get_eh_tags(const int iface_no,
                                           iBase_EntityHandle *entities,
                                           int num_entities,
                                           iBase_TagHandle tag_handle,
                                           bool are_sets,
                                           iBase_EntityHandle *tag_values) 
{
  int result = iBase_SUCCESS;
  
  if (are_sets) {
    CAST_iBase_INTERFACE(ifaceInstances[iface_no], iface_stag, SetTag,
                         iBase_FAILURE);
    for (int i = 0; i < num_entities; i++) {
      try {
        tag_values[i] = iface_stag.getEntSetEHData(entities[i],
                                                   tag_handle);
      }
      catch (iBase::Error err) {
        if (iBase::ErrorType_SUCCESS != err.getErrorType() &&
            iBase::ErrorType_TAG_NOT_FOUND != err.getErrorType()) {
          iRel_LAST_ERROR.error_type = iBase_FAILURE;
          sprintf(iRel_LAST_ERROR.description, "%s", err.getDescription().c_str());
          result = iBase_FAILURE;
        }
      }
    }
  }
  else {
    CAST_iBase_INTERFACE(ifaceInstances[iface_no], iface_atag, ArrTag,
                         iBase_FAILURE);
    sidl::array<iBase_EntityHandle> tag_values_tmp = 
      convert_to_sidl_vector(tag_values, num_entities);

    try {
      iface_atag.getEHArrData(convert_to_sidl_vector(entities, num_entities),
                              num_entities, tag_handle,
                              tag_values_tmp, num_entities);
    }
    catch (iBase::Error err) {
      if (iBase::ErrorType_SUCCESS != err.getErrorType() &&
          iBase::ErrorType_TAG_NOT_FOUND != err.getErrorType()) {
        iRel_LAST_ERROR.error_type = iBase_FAILURE;
        sprintf(iRel_LAST_ERROR.description, "%s", err.getDescription().c_str());
        result = iBase_FAILURE;
      }
    }
  }

  RETURN(result);
}

  
int AssocPairSIDL::set_int_tags(const int iface_no,
                                            iBase_EntityHandle *entities,
                                            int num_entities,
                                            iBase_TagHandle tag_handle,
                                            bool are_sets,
                                            int *tag_values)
{
  try {
    if (are_sets) {
      CAST_iBase_INTERFACE(ifaceInstances[iface_no], iface_stag, SetTag,
                           iBase_FAILURE);
      for (int i = 0; i < num_entities; i++) {
        tag_values[i] = iface_stag.getEntSetIntData(entities[i],
                                                    tag_handle);
      }
    }
    else {
      CAST_iBase_INTERFACE(ifaceInstances[iface_no], iface_atag, ArrTag,
                           iBase_FAILURE);
      sidl::array<int> tag_values_tmp = 
        convert_to_sidl_vector(tag_values, num_entities);
    
      iface_atag.setIntArrData(convert_to_sidl_vector(entities, num_entities),
                               num_entities, tag_handle,
                               tag_values_tmp, num_entities);
    }
  }
  catch (iBase::Error err) {
    iRel_LAST_ERROR.error_type = iBase_FAILURE;
    sprintf(iRel_LAST_ERROR.description, "%s", err.getDescription().c_str());
    RETURN(iBase_FAILURE);
  }

  RETURN(iBase_SUCCESS);
}
  
int AssocPairSIDL::set_eh_tags(const int iface_no,
                                           iBase_EntityHandle *entities,
                                           int num_entities,
                                           iBase_TagHandle tag_handle,
                                           bool are_sets,
                                           iBase_EntityHandle *tag_values)
{
  try {
    if (are_sets) {
      CAST_iBase_INTERFACE(ifaceInstances[iface_no], iface_stag, SetTag,
                           iBase_FAILURE);
      for (int i = 0; i < num_entities; i++) {
        iface_stag.setEntSetEHData(entities[i],
                                   tag_handle,
                                   tag_values[i]);
      }
    }
    else {
      CAST_iBase_INTERFACE(ifaceInstances[iface_no], iface_atag, ArrTag,
                           iBase_FAILURE);
      sidl::array<iBase_EntityHandle> tag_values_tmp = 
        convert_to_sidl_vector(tag_values, num_entities);
    
      iface_atag.setEHArrData(convert_to_sidl_vector(entities, num_entities),
                              num_entities, tag_handle,
                              tag_values_tmp, num_entities);
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
                                                const bool are_sets,
                                                iBase_EntityHandle **entities,
                                                int *entities_alloc,
                                                int *entities_size) 
{
  sidl::array<iBase_EntityHandle> entities_tmp;
  int entities_tmp_size = 0;

  if (are_sets) {
    CAST_iBase_INTERFACE(ifaceInstances[iface_no], eset, EntSet,
                         iBase_FAILURE);
    eset.getEntSets(0, 0, entities_tmp, entities_tmp_size);
  }
  else if (iface_type(iface_no) == iGeom_IFACE) {
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
  
  if (iGeom_IFACE != ifaceTypes[iface_no] &&
      iMesh_IFACE != ifaceTypes[iface_no]) {
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

int AssocPairSIDL::create_mesh_vertex(const double x,
                                                  const double y,
                                                  const double z,
                                                  iBase_EntityHandle &vertex) 
{
  int iface_no;
  if (ifaceTypes[0] == iMesh_IFACE) iface_no = 0;
  else if (ifaceTypes[1] == iMesh_IFACE) iface_no = 1;
  else {
    iRel_processError(iBase_INVALID_ARGUMENT, 
                       "One of the interfaces must be mesh.");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  
  try {
    CAST_iMesh_INTERFACE(ifaceInstances[iface_no], mesh_mod, Modify,
                         iBase_FAILURE);
    mesh_mod.createVtx(x, y, z, vertex);
  }
  catch (iBase::Error err) {
    iRel_processError(err.getErrorType(), err.getDescription().c_str());
    RETURN((int)err.getErrorType());
  }
  
  RETURN(iBase_SUCCESS);
}

int AssocPairSIDL::create_mesh_entity(iMesh_EntityTopology ent_topology,
                                                  iBase_EntityHandle *lower_handles,
                                                  int lower_handles_size,
                                                  iBase_EntityHandle &new_ent,
                                      int &status)
{
  int iface_no;
  if (ifaceTypes[0] == iMesh_IFACE) iface_no = 0;
  else if (ifaceTypes[1] == iMesh_IFACE) iface_no = 1;
  else {
    iRel_processError(iBase_INVALID_ARGUMENT, 
                       "One of the interfaces must be mesh.");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  
  try {
    CAST_iMesh_INTERFACE(ifaceInstances[iface_no], mesh_mod, Modify,
                         iBase_FAILURE);
    mesh_mod.createEnt((iMesh::EntityTopology)ent_topology, 
                       convert_to_sidl_vector(lower_handles, 
                                              lower_handles_size), 
                       lower_handles_size,
                       new_ent, (iMesh::CreationStatus&)status);
  }
  catch (iBase::Error err) {
    iRel_processError(err.getErrorType(), err.getDescription().c_str());
    RETURN((int)err.getErrorType());
  }
  
  RETURN(iBase_SUCCESS);
}

int AssocPairSIDL::set_mesh_coords(iBase_EntityHandle *verts,
                                               int verts_size,
                                               double *coords,
                                               iBase_StorageOrder order) 
{
  int iface_no;
  if (iMesh_IFACE == ifaceTypes[0]) iface_no = 0;
  else if (iMesh_IFACE == ifaceTypes[1]) iface_no = 1;
  else {
    iRel_processError(iBase_INVALID_ARGUMENT, 
                       "One of the interfaces must be mesh.");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  
  try {
    CAST_iMesh_INTERFACE(ifaceInstances[iface_no], mesh_amod, ArrMod,
                         iBase_FAILURE);
    mesh_amod.setVtxArrCoords(convert_to_sidl_vector(verts, verts_size), 
                              verts_size, (iMesh::StorageOrder)order, 
                              convert_to_sidl_vector(coords, 3*verts_size), 
                              3*verts_size);
  }
  catch (iBase::Error err) {
    iRel_processError(err.getErrorType(), err.getDescription().c_str());
    RETURN((int)err.getErrorType());
  }
  
  RETURN(iBase_SUCCESS);
}

int AssocPairSIDL::get_mesh_coords(iBase_EntityHandle *verts,
                                               int verts_size,
                                               double **coords,
                                               int *coords_alloc,
                                               int *coords_size,
                                               int *order) 
{
  int iface_no;
  if (iMesh_IFACE == ifaceTypes[0]) iface_no = 0;
  else if (iMesh_IFACE == ifaceTypes[1]) iface_no = 1;
  else {
    iRel_processError(iBase_INVALID_ARGUMENT, 
                       "One of the interfaces must be mesh.");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  CHECK_SIZE_VOID(double, coords, coords_alloc, 3*verts_size);
  *coords_size = 3*verts_size;
  
  try {
    CAST_iMesh_INTERFACE(ifaceInstances[iface_no], mesh_mod, Modify,
                         iBase_FAILURE);
    sidl::array<double> tmp_coords = convert_to_sidl_vector(*coords, *coords_size);
    mesh_mod.getVtxArrCoords(convert_to_sidl_vector(verts, verts_size), 
                             verts_size, (iMesh::StorageOrder&)*order, 
                             tmp_coords, 
                             *coords_size);
  }
  catch (iBase::Error err) {
    iRel_processError(err.getErrorType(), err.getDescription().c_str());
    RETURN((int)err.getErrorType());
  }
  
  RETURN(iBase_SUCCESS);
}

int AssocPairSIDL::get_closest_pt(iBase_EntityHandle *gents, 
                                              int gents_size,
                                              iBase_StorageOrder storage_order,
                                              double *near_coords, 
                                              int near_coords_size,
                                              double **on_coords,
                                              int *on_coords_alloc, 
                                              int *on_coords_size) 
{
  int iface_no;
  if (iGeom_IFACE == ifaceTypes[0]) iface_no = 0;
  else if (iGeom_IFACE == ifaceTypes[1]) iface_no = 1;
  else {
    iRel_processError(iBase_INVALID_ARGUMENT, 
                       "One of the interfaces must be geom.");
    RETURN(iBase_INVALID_ARGUMENT);
  }
  CHECK_SIZE_VOID(double, on_coords, on_coords_alloc, 3*gents_size);
  *on_coords_size = 3*gents_size;

  try {
    CAST_iGeom_INTERFACE(ifaceInstances[iface_no], geom_shape, Shape,
                         iBase_FAILURE);
    sidl::array<double> tmp_on_coords = 
      convert_to_sidl_vector(*on_coords, *on_coords_size);
    
    geom_shape.getArrClosestPt(convert_to_sidl_vector(gents, gents_size),
                               gents_size,
                               (iBase::StorageOrder)storage_order,
                               convert_to_sidl_vector(near_coords, near_coords_size),
                               near_coords_size,
                               tmp_on_coords,
                               *on_coords_size);
  }
  catch (iBase::Error err) {
    iRel_processError(err.getErrorType(), err.getDescription().c_str());
    RETURN((int)err.getErrorType());
  }
  
  RETURN(iBase_SUCCESS);
}
