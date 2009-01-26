#ifndef IMESH_REC_CBIND_H__
#define IMESH_REC_CBIND_H__

#include "iMesh.h"
#include "iMesh_extensions_protos.h"

#ifdef __cplusplus
extern "C" {
#endif

    /**\brief  Get entities of specific type and/or topology in set or instance, recursive
     *
     * Get entities of specific type and/or topology in set or instance.  If recursive
     * is passed in non-zero, includes entities in owned sets.  All 
     * entities of a given type or topology are requested by specifying
     * iBase_ALL_TOPOLOGIES or iBase_ALL_TYPES, respectively.  Specified type
     * or topology must be a value in the iBase_EntityType or iMesh_EntityTopology
     * enumeration, respectively.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param entity_type Type of entities being requested
     * \param entity_topology Topology of entities being requested
     * \param recursive If non-zero, gets entities in owned sets too
     * \param *entity_handles Pointer to array of entity handles returned 
     *        from function
     * \param *entity_handles_allocated Pointer to allocated size of 
     *        entity_handles array
     * \param *entity_handles_size Pointer to occupied size of entity_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntitiesRec(iMesh_Instance instance,
                            /*in*/ const iBase_EntitySetHandle entity_set_handle,
                            /*in*/ const int entity_type,
                            /*in*/ const int entity_topology,
                            /*in*/ const int recursive,
                            /*out*/ iBase_EntityHandle** entity_handles,
                            /*out*/ int* entity_handles_allocated,
                            /*out*/ int* entity_handles_size,
                            /*out*/ int *err);

    /**\brief  Get the number of entities with the specified type in the instance or set, recursive
     *
     * Get the number of entities with the specified type in the instance 
     * or set.  If recursive is passed in non-zero, includes entities in owned sets.  
     * If entity set handle is zero, return information for instance, 
     * otherwise for set.  Value of entity type must be from the
     * iBase_EntityType enumeration.  If iBase_ALL_TYPES is specified,
     * total number of entities (excluding entity sets) is returned.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param entity_type Type of entity requested
     * \param recursive If non-zero, includes entities in owned sets too
     * \param num_type Pointer to number of entities, returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNumOfTypeRec(iMesh_Instance instance,
                             /*in*/ const iBase_EntitySetHandle entity_set_handle,
                             /*in*/ const int entity_type,
                             /*in*/ const int recursive,
                             /*out*/ int *num_type, 
                             /*out*/ int *err);

    /**\brief  Get the number of entities with the specified topology in the instance or set
     *
     * Get the number of entities with the specified topology in the instance 
     * or set.  If recursive is passed in non-zero, includes entities in owned sets.  
     * If entity set handle is zero, return information for instance,
     * otherwise for set.  Value of entity topology must be from the
     * iMesh_EntityTopology enumeration.  If iMesh_ALL_TOPOLOGIES is specified,
     * total number of entities (excluding entity sets) is returned.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param entity_topology Topology of entity requested
     * \param recursive If non-zero, includes entities in owned sets too
     * \param num_topo Pointer to number of entities, returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNumOfTopoRec(iMesh_Instance instance,
                             /*in*/ const iBase_EntitySetHandle entity_set_handle,
                             /*in*/ const int entity_topology,
                             /*in*/ const int recursive,
                             /*out*/ int *num_topo, 
                             /*out*/ int *err);


    /**\brief  Get entities with specified type, topology, tag(s) and (optionally) tag value(s)
     *
     * Get entities with the specified type, topology, tag(s), and optionally tag value(s).
     * If tag values pointer is input as zero, entities with specified tag(s) are returned,
     * regardless of their value.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param entity_type Type of entities being requested
     * \param entity_topology Topology of entities being requested
     * \param tag_handles Array of tag handles
     * \param tag_vals Array of tag values (zero if values not requested)
     * \param num_tags_vals Number of tags and optionally values
     * \param recursive If non-zero, gets entities in owned sets too
     * \param *entity_handles Pointer to array of entity handles returned 
     *        from function
     * \param *entity_handles_allocated Pointer to allocated size of 
     *        entity_handles array
     * \param *entity_handles_size Pointer to occupied size of entity_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntsByTagsRec(iMesh_Instance instance,
                              /*in*/ const iBase_EntitySetHandle entity_set_handle,
                              /*in*/ const int entity_type,
                              /*in*/ const int entity_topology,
                              /*in*/ const iBase_TagHandle *tag_handles,
                              /*in*/ const char * const *tag_vals,
                              /*in*/ const int num_tags_vals,
                              /*in*/ const int recursive,
                              /*out*/ iBase_EntityHandle** entity_handles,
                              /*out*/ int* entity_handles_allocated,
                              /*out*/ int* entity_handles_size,
                              /*out*/ int *err);

    /**\brief  Get entity sets with specified type, topology, tag(s) and (optionally) tag value(s)
     *
     * Get entity sets with the specified type, topology, tag(s), and optionally tag value(s).
     * If tag values pointer is input as zero, entities with specified tag(s) are returned,
     * regardless of their value.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param tag_handles Array of tag handles
     * \param tag_vals Array of tag values (zero if values not requested)
     * \param num_tags_vals Number of tags and optionally values
     * \param recursive If non-zero, gets entities in owned sets too
     * \param *set_handles Pointer to array of entity handles returned 
     *        from function
     * \param *set_handles_allocated Pointer to allocated size of 
     *        set_handles array
     * \param *set_handles_size Pointer to occupied size of entity_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntSetsByTagsRec(iMesh_Instance instance,
                                 /*in*/ const iBase_EntitySetHandle entity_set_handle,
                                 /*in*/ const iBase_TagHandle *tag_handles,
                                 /*in*/ const char * const *tag_vals,
                                 /*in*/ const int num_tags_vals,
                                 /*in*/ const int recursive,
                                 /*out*/ iBase_EntitySetHandle** set_handles,
                                 /*out*/ int* set_handles_allocated,
                                 /*out*/ int* set_handles_size,
                                 /*out*/ int *err);

    /**\brief Get MBCN type corresponding to iMesh topology value
     *
     * Get MBCN type corresponding to iMesh topology value.  Required for input
     * to MBCN canonical numbering functions, which are written in terms of 
     * MBCN entity types.  Returns -1 for type if entity topology is out of
     * bounds, or MBMAXTYPE if no corresponding MBCN type exists.
     * \param imesh_entity_topology iMesh_EntityTopology value
     * \param mbcn_type MBEntityType corresponding to entity topology
     */
  void iMesh_MBCNType(/*in*/ const int imesh_entity_topology,
                      /*out*/ int *mbcn_type);
    

#ifdef __cplusplus
}
#endif

#endif
