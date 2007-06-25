#ifndef __iRel_LASSO_HPP__
#define __iRel_LASSO_HPP__

#include "iGeom.h"
#include "iMesh.h"
#include "iRel_protos.h"

#ifdef __cplusplus

extern "C" 
{
#endif

  typedef void* iRel_Instance;
  typedef void* iRel_RelationHandle;

  enum IfaceType 
  {iRel_IBASE_IFACE, 
   iRel_IGEOM_IFACE, 
   iRel_IMESH_IFACE, 
   iRel_IFIELD_IFACE, 
   iRel_IREL_IFACE};

  extern struct iBase_Error iRel_LAST_ERROR;

  void iRel_dtor(iRel_Instance instance);

  void iRel_createAssociation (
    iRel_Instance instance,
    iBase_Instance iface1,
    const int ent_or_set1,
    const int iface_type1,
    iBase_Instance iface2,
    const int ent_or_set2,
    const int iface_type2,
    iRel_RelationHandle *rel,
    int *ierr);
  
  void iRel_destroyAssociation (
    iRel_Instance instance, 
    iRel_RelationHandle rel,
    int ierr);

  void iRel_getAssociatedInterfaces (
    iRel_Instance instance,
    iBase_Instance iface,
    iBase_Instance **interfaces,
    int *interfaces_allocated,
    int *interfaces_size,
    int ierr);

  void iRel_setEntEntAssociation (
    iRel_Instance instance,
    iRel_RelationHandle rel,
    iBase_EntityHandle ent1,
    int is_set1,
    iBase_EntityHandle ent2,
    int is_set2,
    int *ierr);

  void iRel_setEntArrAssociation (
    iRel_Instance instance,
    iRel_RelationHandle rel,    
    iBase_EntityHandle ent1,
    int is_set1,
    int switch_order,
    iBase_EntityHandle *ent_array_2,
    int num_entities,
    int is_set2,
    int *ierr);

  void iRel_setArrAssociation (
    iRel_Instance instance,
    iRel_RelationHandle rel,    
    iBase_EntityHandle *ent_array_1,
    int num_ent1,
    int is_set1,
    iBase_EntityHandle *ent_array_2,
    int num_ent2,
    int is_set2,
    int *ierr);

  void iRel_getEntEntAssociation (
    iRel_Instance instance,
    iRel_RelationHandle rel,    
    iBase_EntityHandle ent1,
    int is_set1,
    int switch_order,
    iBase_EntityHandle *ent2,
    int *is_set2,
    int *ierr);

  void iRel_getEntArrAssociation (
    iRel_Instance instance,
    iRel_RelationHandle rel,    
    iBase_EntityHandle ent1,
    int is_set1,
    int return_sets,
    int switch_order,
    iBase_EntityHandle **ent_array_2,
    int *ent_array_2_allocated,
    int *ent_array_2_size,
    int *ierr);

  void iRel_getArrAssociation (
    iRel_Instance instance,
    iRel_RelationHandle rel,    
    iBase_EntityHandle *ent_array_1,
    int ent_array_1_size,
    int is_set1,
    int return_sets,
    int switch_order,
    iBase_EntityHandle **ent_array_2,
    int *ent_array_2_allocated,
    int *ent_array_2_size,
    int **offset,
    int *offset_allocated,
    int *offset_size,
    int *ierr);

  void iRel_createVtxAndAssociate (
    iRel_Instance instance,
    double x,
    double y,
    double z,
    iBase_EntityHandle associatedGeomEnt,
    iBase_EntityHandle *new_entity_handle,
    int *ierr);

  void iRel_createEntAndAssociate (
    iRel_Instance instance,
    int new_entity_topology,
    iBase_EntityHandle *lower_order_entity_handles,
    int lower_order_entity_handles_size,
    iBase_EntityHandle associatedGeomEnt,
    iBase_EntityHandle *new_entity_handle,
    int *creation_status,
    int *ierr);

  void iRel_createVtxArrAndAssociate (
    iRel_Instance instance,
    int num_verts,
    int storage_order,
    double *new_coords,
    int new_coords_size,
    iBase_EntityHandle *associatedGeomEnts,
    iBase_EntityHandle **new_vertex_handles,
    int *new_vertex_handles_allocated,
    int *new_vertex_handles_size,
    int *ierr);

  void iRel_createEntArrAndAssociate (
    iRel_Instance instance,
    int new_entity_topology,
    iBase_EntityHandle *lower_order_entity_handles,
    int lower_order_entity_handles_size,
    int *offsets,
    int offsets_size,
    iBase_EntityHandle *associatedGeomEnts,
    iBase_EntityHandle **new_entity_handles,
    int *new_entity_handles_allocated,
    int *new_entity_handles_size,
    int **status,
    int *status_allocated,
    int *status_size,
    int *ierr);

  void iRel_inferAllAssociations (
    iRel_Instance instance,
    iRel_RelationHandle rel,
    int *ierr);

  void iRel_inferEntAssociations (
    iRel_Instance instance,
    iRel_RelationHandle rel,    
    iBase_EntityHandle entity,
    int is_set,
    int iface_no,
    int *ierr);

  void iRel_inferArrAssociations (
    iRel_Instance instance,
    iRel_RelationHandle rel,    
    iBase_EntityHandle *entities,
    int entities_size,
    int is_set,
    int iface_no,
    int *ierr);

  void iRel_moveTo(iRel_Instance instance,
                   iGeom_Instance geom, iMesh_Instance mesh,
                   iBase_EntityHandle geom_entity_handle,
                   int *ierr);

  void iRel_newAssoc(const char *options,
                     iRel_Instance *instance,
                     int *ierr,
                     const int options_len);
  
#ifdef __cplusplus
} // extern "C"
#endif

#endif // #ifndef __iRel_LASSO_HPP__

