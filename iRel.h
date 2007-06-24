#ifndef __iRel_LASSO_HPP__
#define __iRel_LASSO_HPP__

#include "iGeom.h"
#include "iMesh.h"

#ifdef __cplusplus

extern "C" 
{
#endif

  typedef void* iRel_Instance;

  enum IfaceType 
  {iRel_IBASE_IFACE, 
   iRel_IGEOM_IFACE, 
   iRel_IMESH_IFACE, 
   iRel_IFIELD_IFACE, 
   iRel_IREL_IFACE};

  extern struct iBase_Error iRel_LAST_ERROR;

  int
  iRel_dtor(iRel_Instance instance);

  int 
  iRel_createAssociation (
    iRel_Instance instance,
    iBase_Instance iface1,
    const int ent_or_set1,
    enum IfaceType type1,
    iBase_Instance iface2,
    const int ent_or_set2,
    enum IfaceType type2);
  
  int iRel_destroyAssociation (
    iRel_Instance instance,
    iBase_Instance iface1,
    iBase_Instance iface2
    );

  int iRel_getAssociatedInterfaces (
    iRel_Instance instance,
    iBase_Instance iface,
    iBase_Instance **interfaces,
    int *interfaces_allocated,
    int *interfaces_size
    );

  int iRel_setEntEntAssociation (
    iRel_Instance instance,
    iBase_Instance iface1,
    iBase_Instance iface2,
    iBase_EntityHandle ent1,
    int is_set1,
    iBase_EntityHandle ent2,
    int is_set2);

  int iRel_setEntArrAssociation (
    iRel_Instance instance,
    iBase_Instance iface1,
    iBase_Instance iface2,
    iBase_EntityHandle ent1,
    int is_set1,
    iBase_EntityHandle *ent_array_2,
    int num_entities,
    int is_set2);

  int iRel_setArrAssociation (
    iRel_Instance instance,
    iBase_Instance iface1,
    iBase_Instance iface2,
    iBase_EntityHandle *ent_array_1,
    int num_ent1,
    int is_set1,
    iBase_EntityHandle *ent_array_2,
    int num_ent2,
    int is_set2);

  int iRel_getEntEntAssociation (
    iRel_Instance instance,
    iBase_Instance iface1,
    iBase_Instance iface2,
    iBase_EntityHandle ent1,
    int is_set1,
    iBase_EntityHandle *ent2,
    int *is_set2);

  int iRel_getEntArrAssociation (
    iRel_Instance instance,
    iBase_Instance iface1,
    iBase_Instance iface2,
    iBase_EntityHandle ent1,
    int is_set1,
    int return_sets,
    iBase_EntityHandle **ent_array_2,
    int *ent_array_2_allocated,
    int *ent_array_2_size);

  int iRel_getArrAssociation (
    iRel_Instance instance,
    iBase_Instance iface1,
    iBase_Instance iface2,
    iBase_EntityHandle *ent_array_1,
    int ent_array_1_size,
    int is_set1,
    int return_sets,
    iBase_EntityHandle **ent_array_2,
    int *ent_array_2_allocated,
    int *ent_array_2_size,
    int **offset,
    int *offset_allocated,
    int *offset_size);

  int iRel_createVtxAndAssociate (
    iRel_Instance instance,
    double x,
    double y,
    double z,
    iBase_EntityHandle associatedGeomEnt,
    iBase_EntityHandle *new_entity_handle
    );

  int iRel_createEntAndAssociate (
    iRel_Instance instance,
    enum iMesh_EntityTopology new_entity_topology,
    iBase_EntityHandle *lower_order_entity_handles,
    int lower_order_entity_handles_size,
    iBase_EntityHandle associatedGeomEnt,
    iBase_EntityHandle *new_entity_handle,
    enum iBase_CreationStatus *status
    );

  int iRel_createVtxArrAndAssociate (
    iRel_Instance instance,
    int num_verts,
    enum iBase_StorageOrder storage_order,
    double *new_coords,
    int new_coords_size,
    iBase_EntityHandle *associatedGeomEnts,
    iBase_EntityHandle **new_vertex_handles,
    int *new_vertex_handles_allocated,
    int *new_vertex_handles_size
    );

  int iRel_createEntArrAndAssociate (
    iRel_Instance instance,
    enum iMesh_EntityTopology new_entity_topology,
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
    int *status_size
    );

  int iRel_inferAllAssociations (
    iRel_Instance instance,
    iBase_Instance iface1,
    iBase_Instance iface2
    );

  int iRel_inferEntAssociations (
    iRel_Instance instance,
    iBase_Instance iface1,
    iBase_Instance iface2,
    iBase_EntityHandle entity,
    int is_set);

  int iRel_inferArrAssociations (
    iRel_Instance instance,
    iBase_Instance iface1,
    iBase_Instance iface2,
    iBase_EntityHandle *entities,
    int entities_size,
    int is_set);

  int 
  iRel_moveTo(iRel_Instance instance,
              iGeom_Instance geom, iMesh_Instance mesh,
              iBase_EntityHandle geom_entity_handle);

  iRel_Instance iRel_newAssoc(const char **options,
                                const int num_options);
  
#ifdef __cplusplus
} // extern "C"
#endif

#endif // #ifndef __iRel_LASSO_HPP__

