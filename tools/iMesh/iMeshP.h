#ifndef IMESHP_H__
#define IMESHP_H__

  /** \mainpage The ITAPS Parallel Mesh Interface iMesh
   *
   */

#ifndef ITAPS
#define ITAPS
#endif

#include "iMesh_extensions.h"
#include "iMeshP_protos.h"

#ifdef __cplusplus
extern "C" {
#endif

    /**\brief  Type used to store a part handle
     *
     * Type used to store a part handle
     */
  typedef void* iMeshP_PartHandle;

    /**\brief  Type used to store a partition handle
     *
     * Type used to store a partition handle
     */
  typedef void* iMeshP_PartitionHandle;

  const int iMeshP_INTERNAL = 0;
  const int iMeshP_BOUNDARY = 1;
  const int iMeshP_GHOST = 2;
  
// Given a partition handle, return the total global number of parts 
// in the partition.
//  COMMUNICATION:  None++.
  void iMeshP_getNumParts(iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          int *num_global_part, 
                          int *err); 

//  Map from parts to processes:  
//  Given a partition handle and a part handle, return the rank of the 
//  process that owns the part.
//  The part_handle may be local or remote.
//  COMMUNICATION:  None++.
  void iMeshP_getRankOfPart(iMesh_Instance instance,
                            const iMeshP_PartitionHandle partition_handle,
                            /*in*/    const iMeshP_PartHandle part_handle,
                            /*out*/   int *rank,
                            int *err); 
  void iMeshP_getRankOfPartArr(iMesh_Instance instance,
                               const iMeshP_PartitionHandle partition_handle,
                               /*in*/    const iMeshP_PartHandle *part_handle,
                               /*in*/    const int part_handle_size,
                               /*inout*/ int **rank, 
                               /*inout*/ int *rank_allocated, 
                               int *err); 

//  Map from processes to parts:  
//  Given a partition handle and a process rank,
//  return the part handles owned by the process.
//  COMMUNICATION:  None++.
  void iMeshP_getPartsOnRank(iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             /*in*/    const int rank,
                             /*inout*/ iMeshP_PartHandle **part_handles, 
                             /*inout*/ int *part_handles_allocated, 
                             /*out*/   int *part_handles_size, 
                             int *err); 
  void iMeshP_getPartsArrOnRank(iMesh_Instance instance,
                                const iMeshP_PartitionHandle partition_handle,
                                /*in*/    const int *rank,
                                /*in*/    const int rank_size,
                                /*inout*/ iMeshP_PartHandle **part_handles, 
                                /*inout*/ int *part_handles_allocated, 
                                /*out*/   int *part_handles_size, 
                                int *err); 

//  Provide global mesh information about a partition.  
//  Note that these functions may require communication and, thus, 
//  would have to be called by all processes in the partition handle.
//  Given a mesh instance and partition handle, return the
//  total number of entities with given type or topology in the partition.
//  COMMUNICATION:  Collective.
  void iMeshP_getNumOfTypeAll(iMesh_Instance instance,
                              const iMeshP_PartitionHandle partition_handle,
                              const iBase_EntitySetHandle entity_set_handle,
                              const int entity_type, int *num_type, int *err);
  void iMeshP_getNumOfTopoAll(iMesh_Instance instance,
                              const iMeshP_PartitionHandle partition_handle,
                              const iBase_EntitySetHandle entity_set_handle,
                              const int entity_topology, int *num_topo, int *err);


// Given a partition handle and a part handle, test whether the part
// handle refers to a local or remote part, and whether it is valid.
  void iMeshP_testPart(iMesh_Instance instance,
                       /* in */ iMeshP_PartitionHandle partition_handle,
                       /* in */ iMeshP_PartHandle part_handle,
                       /* out */ int *part_status, /* LOCAL, REMOTE, INVALID */
                       int *err);
 
//  Provide part information about an entity: Given an entity and a 
//  partition handle, return the part handle of the part that owns the entity.
//  Return an error code if an entity is not in the partition.
//  COMMUNICATION:  None++.
  void iMeshP_getEntOwnerPart(iMesh_Instance instance,
                              /*in*/  const iMeshP_PartitionHandle partition_handle, 
                              /*in*/  iBase_EntityHandle entity_handle,
                              /*out*/  iMeshP_PartHandle *part_handle,
                              int* err); 
  void iMeshP_getEntOwnerPartArr(iMesh_Instance instance,
                                 /*in*/  const iMeshP_PartitionHandle partition_handle, 
                                 /*in*/  iBase_EntityHandle *entity_handles,
                                 /*in*/  int entity_handles_size,
                                 /*inout*/  iMeshP_PartHandle **part_handles,
                                 /*inout*/  int *part_handles_allocated,
                                 int* err); 
  
//  Provide entity categorization within part.
//  Given a partition handle, a part handle, and an entity handle, return a
//  flag indicating whether the entity is owned by the part.
//  If part_handle is remote, an error is returned.
//  COMMUNICATION:  None.
  void iMeshP_isEntOwner(iMesh_Instance instance,
                         /*in*/ const iMeshP_PartHandle part_handle, 
                         /*in*/ const iBase_EntityHandle entity_handle, 
                         /*out*/ int* is_owner, 
                         int *err); 
  void iMeshP_isEntOwnerArr(iMesh_Instance instance,
                            /*in*/ const iMeshP_PartHandle part_handle, 
                            /*in*/ const iBase_EntityHandle *entity_handles, 
                            /*in*/ const int entity_handles_size, 
                            /*inout*/ int** is_owner, 
                            /*inout*/ int* is_owner_allocated, 
                            /*inout*/ int* is_owner_size, 
                            int *err); 

//  Given a partition handle, a part handle, and an entity handle, return a
//  flag indicating whether the entity is strictly internal, on a 
//  boundary, or a ghost.
//  If part_handle is remote, an error is returned.
//  COMMUNICATION:  None.
  void iMeshP_getEntStatus(iMesh_Instance instance,
                           /*in*/ const iMeshP_PartHandle part_handle, 
                           /*in*/ const iBase_EntityHandle entity_handle, 
                           /*out*/ int* par_status, // Values=INTERNAL,BOUNDARY,GHOST
                           int *err); 
  void iMeshP_getEntStatusArr(iMesh_Instance instance,
                              /*in*/ const iMeshP_PartHandle part_handle, 
                              /*in*/ const iBase_EntityHandle *entity_handles, 
                              /*in*/ const int entity_handles_size, 
                              /*inout*/ int** par_status, // Values=INTERNAL,BOUNDARY,GHOST
                              /*inout*/ int* par_status_allocated, 
                              /*inout*/ int* par_status_size, 
                              int *err); 

//  Provide information about copies of entities.  
//  All these functions should work on the local process as well as 
//  remote processes; entity handles returned are likely but not 
//  necessarily remote. 
//  Given a partition handle and an entity handle, return the number 
//  of copies of the entity in the partition.
//  COMMUNICATION:  None++.
  void iMeshP_getNumCopies(iMesh_Instance instance,
                           /*in*/ const iMeshP_PartitionHandle partition_handle, 
                           /*in*/ const iBase_EntityHandle entity_handle, 
                           /*out*/ int *num_copies_ent,
                           int *err); 

//  Given a partition handle and an entity handle, return (remote) entity
//  handles and part handles of all copies of the entity.
//  COMMUNICATION:  None++.
  void iMeshP_getCopies(iMesh_Instance instance,
                        /*in*/ const iMeshP_PartitionHandle partition_handle, 
                        /*in*/ const iBase_EntityHandle entity_handle, 
                        /*inout*/ iMeshP_PartHandle **part_handles, 
                        /*inout*/ int *part_handles_allocated, 
                        /*out*/ int *part_handles_size, 
                        /*inout*/ iBase_EntityHandle **copies_entity_handles, 
                        /*inout*/ int *copies_entity_handles_allocated, 
                        /*inout*/ int *copies_entity_handles_size, 
                        int *err); 

#ifdef __cplusplus
}
#endif

#endif
