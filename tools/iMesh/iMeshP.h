#ifndef iMeshP_H
#define iMeshP_H

#include "iMesh.h"
#include "iMeshP_protos.h"
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ITAPS_DECLARE_HANDLE( NAME ) typedef struct NAME##DummyStruct *NAME

ITAPS_DECLARE_HANDLE( iMeshP_PartitionHandle );
ITAPS_DECLARE_HANDLE( iMeshP_PartHandle );
ITAPS_DECLARE_HANDLE( iMeshP_RequestHandle );
ITAPS_DECLARE_HANDLE( iMeshP_Status );

typedef unsigned iMeshP_Part;

enum iMeshP_EntStatus 
{
  iMeshP_INTERNAL,
  iMeshP_BOUNDARY,
  iMeshP_GHOST
};


/*
------------------------------------------------
Background and Terminology:

Terminology -- Abstract Data Model:
-  The term "mesh" refers to an abstraction in the data model; 
   it does not imply a serial or parallel distribution.
-  The term "partition" refers to an assignment of a set of entities to 
   subsets; like a "mesh," it does not imply a serial or parallel 
   implementation.
-  An application may use one or more meshes.  
-  Parititions can create subsets of entities from one or more meshes.
-  Meshes can be subdivided by one or more partitions.
-  Partitions contain parts.  Parts contain the subsets of entities in the
   partition.

Terminology -- Parallelism
-  A "process" can be thought of as an MPI process. The
   number of processes can be considered to be the result of MPI_Comm_size.
   The rank of a process can be thought of as the result of MPI_Comm_rank.
   We will think in terms of processes rather than processors.  Initial
   implementations of the parallel interface will likely use MPI terminology
   directly; future implementations may accommodate other communication 
   paradigms and libraries.
-  Partitions have communicators associated with them.  These communicators
   can be thought of as MPI communicators.  
-  "Global" operations are operations performed with respect to a 
   partition's communicator.
-  "Local" operations are operations performed with respect to a part or
   a mesh instance within a process.
-  Part A "neighbors" Part B if Part A has copies of entities owned by Part B
   and/or if Part B has copies of entities owned by Part A.
-  Each function description includes its communication requirements.  The
   options are described here:
   +  COMMUNICATION:  Collective -- the function must be called by all 
      processes in the partition's communicator.
   +  COMMUNICATION:  Point-to-Point -- communication is used, but the 
      communication is from one process to only one other process.  The
      receiving process must issue an appropriate receive call to receive 
      the message.
   +  COMMUNICATION:  None -- the function does not use communication; only
      local operations are performed.
   +  COMMUNICATION:  None++ -- no communication is done; the values
      are precomputed by iMeshP_syncPartitionAll or 
      iMeshP_completeParallelMeshPar.
                  
Terminology -- Interfaces
-  Each process has one or more "mesh instances."  A mesh instance can be
   thought of as a mesh database.  An implementation should support the 
   existence of more than one mesh instance per process (e.g., it should 
   always associate mesh data with a mesh instance).  However, we expect 
   applications would most often use only one mesh instance per process.
-  There is one root set per mesh instance.
-  Each process may have one or more partition handles.
-  A partition partitions entities in one mesh instance.  
-  Entities in a mesh instance can be partitioned by one or more partitions.  
   Mesh instances know which partitions they contain.
-  Parts are uniquely identified globally by part IDs of type iMeshP_Part.
   Local parts can also be accessed by part handles.  
   Functions accepting part handles operate correctly on only local 
   parts (parts on the calling process); they will return an error 
   for remote (off-process) parts.  


Assumptions:
-  Each part is wholly contained within a process.  
-  A process may have zero, one or multiple parts.
-  Many iMesh functions that accept EntitySet handles
   are also useful in the context of part handles.  
   These functions will be reinterpreted so that they can accept either an
   EntitySet handle, or a part handle.  See Carl's Oct 15 email for more detail.
-  Generation and management of global IDs for entities 
   will be provided as a service above the parallel iMesh interface.
   Uniqueness of global IDs is managed at the partition level.
-  For each entity that is copied onto remote parts, the owning part knows 
   both the remote part ID and remote entity handle of all copies.
-  All parts with copies of a boundary entity know the remote part ID 
   and remote entity handle of all copies of the entity.  
-  All parts with copies of any entity know the part ID and
   entity handle corresponding to the owner of the entity.
-  Functions that return entity information for a part, set or mesh 
   instance return the information for all entities (including copies and
   ghosts) in that part, set or mesh instance.  Applications can check 
   whether an entity is owned or a ghost using iMeshP_isEntOwner or
   iMeshP_getEntStatus.
*/

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*                          Partition Functionality                       */
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

/**  Create a partition handle:  Given a mesh instance to contain
 *  the partition and a communicator, return a partition
 *  handle.  May have different creation routines for different 
 *  communication systems; once created, the application shouldn't have to
 *  worry about the communication system again.
 *  For serial use, *communicator may be MPI_COMM_SELF or communicator may
 *  be NULL.
 *  COMMUNICATION:  Collective.*/
void iMeshP_createPartitionAll(iMesh_Instance instance,
                    /*in*/  MPI_Comm communicator,
                    /*out*/ iMeshP_PartitionHandle *partition_handle,
                            int *err);
 
/**  Destroy a partition handle: Given a partition handle, destroy it and 
 *  invalidate it.
 *  COMMUNICATION:  Collective.*/
void iMeshP_destroyPartitionAll(iMesh_Instance instance,
                    /*inout*/ iMeshP_PartitionHandle *partition_handle,
                              int *err);

/**  Given a partition handle, return its communicator.
 *  COMMUNICATION:  None*/
void iMeshP_getPartitionComm(iMesh_Instance instance,
                    /*in*/  iMeshP_PartitionHandle partition_handle,
                    /*out*/ MPI_Comm *communicator,
                            int *err);
    
/**  Update a partition after parts have been added.
 *  Gives the implementation an opportunity to locally store info
 *  about the partition so that queries on the partition can be 
 *  performed without synchronous communication. 
 *  Values that are precomputed by syncPartitionAll:
 *  +  The total number of parts in a partition.
 *  +  Mapping between part IDs and processes.
 *  +  Updated remote entity handle information.
 *  COMMUNICATION:  Collective.*/
void iMeshP_syncPartitionAll(iMesh_Instance instance,
                       iMeshP_PartitionHandle partition_handle,
                       int *err); 

/**  Given a mesh instance, return number of partition handles and 
 *  all partition handles associated with the mesh instance.
 *  COMMUNICATION:  None.*/
void iMeshP_getNumPartitions(iMesh_Instance instance,
                    /*out*/ int *num_partitions,
                            int *err);

void iMeshP_getPartitions(iMesh_Instance instance,
             /*inout*/ iMeshP_PartitionHandle **partition_handle,
             /*inout*/ int *partition_handle_allocated, 
             /*out*/   int *partition_handle_size, 
                       int *err); 

/** Given a partition handle, return the total global number of parts 
 * in the partition.
 *  COMMUNICATION:  None++.*/
void iMeshP_getNumGlobalParts(iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          int *num_global_part, 
                          int *err); 

/** Given a partition handle, return the number of local (on-process) parts 
 * in the partition.
 *  COMMUNICATION:  None++.*/
void iMeshP_getNumLocalParts(iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          int *num_local_part, 
                          int *err); 

/** Given a partition handle, return the number of local parts 
 * in the partition as well as the partition handles for the local parts.
 *  COMMUNICATION:  None.*/
void iMeshP_getLocalParts(iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          int *num_local_part, 
                          iMeshP_PartHandle *part_handles,
                          int *part_handles_allocated,
                          int *part_handles_size,
                          int *err); 

/**  Map from parts to processes:  
 *  Given a partition handle and a part ID, return the rank of the 
 *  process that owns the part.
 *  The part may be local or remote.
 *  COMMUNICATION:  None++.*/
void iMeshP_getRankOfPart(iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
             /*in*/    const iMeshP_Part part_id,
             /*out*/   int *rank,
                       int *err); 
void iMeshP_getRankOfPartArr(iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
             /*in*/    const iMeshP_Part *part_ids,
             /*in*/    const int part_ids_size,
             /*inout*/ int **rank, 
             /*inout*/ int *rank_allocated, 
                       int *rank_size,
                       int *err); 

/**  Provide global mesh information about a partition.  
 *  Note that these functions may require communication and, thus, 
 *  would have to be called by all processes in the partition handle.
 *  Given a mesh instance and partition handle, return the
 *  total number of entities with given type or topology in the partition.
 *  COMMUNICATION:  Collective.*/
void iMeshP_getNumOfTypeAll(iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
                       const iBase_EntitySetHandle entity_set_handle,
                       const int entity_type, int *num_type, int *err);
void iMeshP_getNumOfTopoAll(iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
                       const iBase_EntitySetHandle entity_set_handle,
                       const int entity_topology, int *num_topo, int *err);


/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*                        Part Functionality                              */
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

/**  Given a partition handle, create a part and add it to the
 *  partition on the process invoking the creation.  Return part handle.
 *  COMMUNICATION:  None.*/
void iMeshP_createPart(iMesh_Instance instance,
                   /*in*/  iMeshP_PartitionHandle partition_handle,
                   /*out*/ iMeshP_PartHandle *part_handle,
                           int *err);
 
/** Given a partition handle and a part handle, remove the part
 * from the partition, destroy the part, and invalidate the part handle.
 * If part_handle is remote, an error is returned.
 *  COMMUNICATION:  None.*/
void iMeshP_destroyPart(iMesh_Instance instance,
                    /*in*/ iMeshP_PartitionHandle partition_handle,
                 /*inout*/ iMeshP_PartHandle *part_handle,
                           int *err);

/********************************
 * Map between part handles and part IDs.
 * Given a partition handle and a local part handle, return the part ID.
 * If the part handle is not a valid part handle for a local part,
 * return an error.
 *  COMMUNICATION:  None.*/

void iMeshP_getPartIdFromPartHandle(iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
             /*in*/    const iMeshP_PartHandle part_handle,
             /*out*/   iMeshP_Part *part_id,
                       int *err);

void iMeshP_getPartIdsFromPartHandlesArr(iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
             /*in*/    const iMeshP_PartHandle *part_handles,
             /*in*/    const int part_handles_size,
             /*inout*/ iMeshP_Part **part_ids,
             /*inout*/ int *part_ids_allocated,
             /*out*/   int *part_ids_size,
                       int *err);

/** Given a partition handle and a part ID, return the part handle 
 * if the part is local; otherwise, return error code.
 *  COMMUNICATION:  None.*/

void iMeshP_getPartHandleFromPartId(iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
             /*in*/    iMeshP_Part part_id,
             /*out*/   iMeshP_PartHandle *part_handle,
                       int *err);

void iMeshP_getPartHandlesFromPartsIdsArr(iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
             /*in*/    const iMeshP_Part *part_ids,
             /*in*/    const int part_ids_size,
             /*inout*/ iMeshP_PartHandle **part_handles,
             /*inout*/ int *part_handles_allocated,
             /*out*/   int *part_handles_size,
                       int *err);

/**  Add/remove on-process entities to/from on-process part:  Given
 *  a partition handle, an entity handle, and a part handle, add/remove the
 *  entity to/from the part.  
 *  On-process add/remove can be accomplished by functions 
 *  that add/remove entities to/from EntitySets.
 *  The part handle is passed in the entity set argument of
 *  the appropriate entity set function.
 *  The part handle must be local.
 *  COMMUNICATION:  None.
 *  iMeshP_addEntToPart --> iMesh_addEntToSet
 *  iMeshP_rmvEntFromPart --> iMesh_rmvEntFromSet
 *  iMeshP_addEntArrToPart --> iMesh_addEntArrToSet
 * iMeshP_rmvEntArrFromPart --> iMesh_rmvEntArrFromSet
 */


/**  Identify parts that neighbor a given part.
 *  Given a partition handle, a part handle, and an entity type, 
 *  return the number of parts in the partition that neighbor the given part
 *  (i.e., that (1) have copies of entities of the given entity type owned by 
 *  the given part or (2) own entities of the given entity type that are 
 *  copied on the given part).
 *  If the part_handle is invalid, an error is returned.
 *  COMMUNICATION:  None++.*/
void iMeshP_getNumPartNbors(iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
             /*in*/    const iMeshP_PartHandle part_handle,
             /*in*/    int entity_type,
             /*out*/   int *num_part_nbors,
                       int *err); 
void iMeshP_getNumPartNborsArr(iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
             /*in*/    const iMeshP_PartHandle *part_handles,
             /*in*/    const int part_handles_size,
             /*in*/    int entity_type,
             /*inout*/ int **num_part_nbors,
             /*inout*/ int *num_part_nbors_allocated,
             /*out*/   int *num_part_nbors_size,
                       int *err); 

/**  Given a partition handle, a part handle, and an entity type, 
 *  return the number of neighbor part IDs and the part IDs of 
 *  all parts in the partition that neighbor the given part
 *  (i.e., that (1) have copies of entities of the given entity type owned by 
 *  the given part or (2) own entities of the given entity type that are 
 *  copied on the given part).
 *  If part_handle is invalid, an error is returned.
 *  COMMUNICATION:  None++.*/
void iMeshP_getPartNbors(iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
             /*in*/    const iMeshP_PartHandle part_handle,
             /*in*/    int entity_type,
             /*out*/   int *num_part_nbors,
             /*inout*/ iMeshP_Part **nbor_part_ids,
             /*inout*/ int *nbor_part_ids_allocated,
             /*out*/   int *nbor_part_ids_size,
                       int *err); 
void iMeshP_getPartNborsArr(iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
             /*in*/    const iMeshP_PartHandle *part_handles,
             /*in*/    const int part_handles_size,
             /*in*/    int entity_type,
             /*inout*/ int **num_part_nbors,
             /*inout*/ int *num_part_nbors_allocated,
             /*out*/   int *num_part_nbors_size,
             /*inout*/ iMeshP_Part **nbor_part_ids,
             /*inout*/ int *nbor_part_ids_allocated,
             /*out*/   int *nbor_part_ids_size,
                       int *err); 

/**  Provide part boundary info:
 *  Note: Allow optional specification of desired entity types and 
 *  topologies; allow target part handle = iMeshP_ALL_PARTS to 
 *  count/include all qualifying boundary entities of part.
 *  Given a partition handle, a part handle, entity type and topology, and a
 *  target part ID, return the number of qualifying entities on 
 *  the part boundary shared with the target part.  
 *  If part_handle is invalid, an error is returned.
 *  COMMUNICATION:  None.*/
void iMeshP_getNumPartBdryEnts(iMesh_Instance instance,
                      /*in*/  const iMeshP_PartitionHandle partition_handle,
                      /*in*/  const iMeshP_PartHandle part_handle, 
                      /*in*/  const int entity_type, 
                      /*in*/  const int entity_topology, 
                      /*in*/  const iMeshP_Part target_part_id, 
                      /*out*/ int *num_entities, 
                              int *err); 

/**  Given a partition handle, a part handle, entity type and topology, and a 
 *  target part ID, return an array of entity handles for all 
 *  qualifying entities along the part boundary shared with the target 
 *  part.
 *  If part_handle is invalid, an error is returned.
 *  COMMUNICATION:  None.*/
void iMeshP_getPartBdryEnts(iMesh_Instance instance,
                      /*in*/  const iMeshP_PartitionHandle partition_handle,
                      /*in*/  const iMeshP_PartHandle part_handle, 
                      /*in*/  const int entity_type, 
                      /*in*/  const int entity_topology, 
                      /*in*/  const iMeshP_Part target_part_id, 
                   /*inout*/  iBase_EntityHandle **entity_handles,
                   /*inout*/  int *entity_handles_allocated,
                      /*out*/ int *entity_handles_size, 
                              int *err); 

/**  Boundary iterators:  Given a partition handle, a part handle, and a 
 *  target part ID, return an iterator over all entities along
 *  the part boundary shared with the target part.  
 *  Functionality for getNext, reset, and end is 
 *  provided through the regular iMesh iterator functions.
 *  If part_handle is invalid, an error is returned.
 *  COMMUNICATION:  None.*/
void iMeshP_initPartBdryEntIter(iMesh_Instance instance,
                      /*in*/  const iMeshP_PartitionHandle partition_handle,
                      /*in*/  const iMeshP_PartHandle part_handle, 
                      /*in*/  const int entity_type, 
                      /*in*/  const int entity_topology, 
                      /*in*/  const iMeshP_Part nbor_part_id, 
                   /*inout*/  iMesh_EntityIterator* entity_iterator, 
                              int* err); 
void iMeshP_initPartBdryEntArrIter(iMesh_Instance instance,
                      /*in*/  const iMeshP_PartitionHandle partition_handle,
                      /*in*/  const iMeshP_PartHandle part_handle, 
                      /*in*/  const int entity_type, 
                      /*in*/  const int entity_topology, 
                      /*in*/  const int array_size, 
                      /*in*/  const iMeshP_Part nbor_part_id, 
                   /*inout*/  iMesh_EntityIterator* entity_iterator, 
                              int* err); 


/**  Provide entity information about a part and entity set. 
 *  Given an entity set handle 
 *  and a part handle, return the number of entities of given type/topo
 *  that are in both the part and the entity set.
 *  If part_handle is invalid, an error is returned.
 *  COMMUNICATION:  None.*/
void iMeshP_getNumOfType(iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
                       const iMeshP_PartHandle part_handle,
                       const iBase_EntitySetHandle entity_set_handle,
                       const int entity_type, 
                       int *num_type, 
                       int *err);
void iMeshP_getNumOfTopo(iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
                       const iMeshP_PartHandle part_handle,
                       const iBase_EntitySetHandle entity_set_handle,
                       const int entity_topology, 
                       int *num_topo, 
                       int *err);

/**  Given an entity set handle 
 *  and a part handle, return vertex information for vertices
 *  that are in both the part and the entity set.
 *  If part_handle is invalid, an error is returned.
 *  COMMUNICATION:  None.*/
void iMeshP_getAllVtxCoords(iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             const iMeshP_PartHandle part_handle,
                             const iBase_EntitySetHandle entity_set_handle,
                             double** coordinates,
                             int* coordinates_allocated,
                             int* coordinates_size,
                             int** in_entity_set,
                             int* in_entity_set_allocated,
                             int* in_entity_set_size,
                             int* storage_order, 
                             int *err);
void iMeshP_getVtxCoordIndex(iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             const iMeshP_PartHandle part_handle,
                             const iBase_EntitySetHandle entity_set_handle,
                             const int requested_entity_type,
                             const int requested_entity_topology,
                             const int entity_adjacency_type,
                             int** offset,
                             int* offset_allocated,
                             int* offset_size,
                             int** index,
                             int* index_allocated,
                             int* index_size,
                             int** entity_topologies,
                             int* entity_topologies_allocated,
                             int* entity_topologies_size, 
                             int *err);

/**  Given an entity set handle 
 *  and a part handle, return entity handles for entities
 *  that are in both the part and the entity set.
 *  If part_handle is invalid, an error is returned.
 *  COMMUNICATION:  None.*/
void iMeshP_getEntities(iMesh_Instance instance,
                        const iMeshP_PartitionHandle partition_handle,
                        const iMeshP_PartHandle part_handle,
                        const iBase_EntitySetHandle entity_set_handle,
                        const int entity_type,
                        const int entity_topology,
                        iBase_EntityHandle** entity_handles,
                        int* entity_handles_allocated,
                        int* entity_handles_size,
                        int *err);

/**  Given an entity set handle 
 *  and a part handle, return entity adjacencies for entities
 *  that are in both the part and the entity set.
 *  If part_handle is invalid, an error is returned.
 *  COMMUNICATION:  None.*/
void iMeshP_getAdjEntities(iMesh_Instance instance,
                           const iMeshP_PartitionHandle partition_handle,
                           const iMeshP_PartHandle part_handle,
                           const iBase_EntityHandle entity_set_handle,
                           const int entity_type_requestor,
                           const int entity_topology_requestor,
                           const int entity_type_requested,
                           iBase_EntityHandle** adj_entity_handles,
                           int* adj_entity_handles_allocated,
                           int* adj_entity_handles_size,
                           int** offset,
                           int* offset_allocated,
                           int* offset_size,
                           int** in_entity_set,
                           int* in_entity_set_allocated,
                           int* in_entity_set_size, 
                           int *err);

/**  Provide entity and array iterators for entities within both a given part 
 *  and a given entity set.  
 *  Functionality for getNext, reset, and end is 
 *  provided through the regular iMesh iterator functions.
 *  If part_handle is invalid, an error is returned.
 *  COMMUNICATION:  None.*/
void iMeshP_initEntIter(iMesh_Instance instance,
                        const iMeshP_PartitionHandle partition_handle,
                        const iMeshP_PartHandle part_handle,
                        const iBase_EntitySetHandle entity_set_handle,
                        const int requested_entity_type,
                        const int requested_entity_topology,
                        iMesh_EntityIterator* entity_iterator,
                        int *err);
void iMeshP_initEntArrIter(iMesh_Instance instance,
                           const iMeshP_PartitionHandle partition_handle,
                           const iMeshP_PartHandle part_handle,
                           const iBase_EntitySetHandle entity_set_handle,
                           const int requested_entity_type,
                           const int requested_entity_topology,
                           const int requested_array_size,
                           iMesh_EntityArrIterator* entArr_iterator,
                           int *err);

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*                           Entity Functionality                         */
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

/**  Provide part information about an entity: Given an entity and a 
 *  partition handle, return the part ID of the part that owns the entity.
 *  Return an error code if an entity is not in the partition.
 *  COMMUNICATION:  None++.*/
void iMeshP_getEntOwnerPart(iMesh_Instance instance,
                      /*in*/  const iMeshP_PartitionHandle partition_handle, 
                      /*in*/  const iBase_EntityHandle entity_handle,
                     /*out*/  iMeshP_Part *part_id,
                              int* err); 
void iMeshP_getEntOwnerPartArr(iMesh_Instance instance,
                      /*in*/  const iMeshP_PartitionHandle partition_handle, 
                      /*in*/  const iBase_EntityHandle *entity_handles,
                      /*in*/  const int entity_handles_size,
                   /*inout*/  iMeshP_Part **part_ids,
                   /*inout*/  int *part_ids_allocated,
                     /*out*/  int *part_ids_size,
                              int* err); 
  
/**  Provide entity categorization within part.
 *  Given a partition handle, a part handle, and an entity handle, return a
 *  flag indicating whether the entity is owned by the part.
 *  If part_handle is invalid, an error is returned.
 *  COMMUNICATION:  None.*/
void iMeshP_isEntOwner(iMesh_Instance instance,
                    /*in*/ const iMeshP_PartitionHandle partition_handle, 
                    /*in*/ const iMeshP_PartHandle part_handle, 
                    /*in*/ const iBase_EntityHandle entity_handle, 
                   /*out*/ int* is_owner, 
                           int *err); 
void iMeshP_isEntOwnerArr(iMesh_Instance instance,
                    /*in*/ const iMeshP_PartitionHandle partition_handle, 
                    /*in*/ const iMeshP_PartHandle part_handle, 
                    /*in*/ const iBase_EntityHandle *entity_handles, 
                    /*in*/ const int entity_handles_size, 
		 /*inout*/ int** is_owner, 
                 /*inout*/ int* is_owner_allocated, 
                 /*inout*/ int* is_owner_size, 
                           int *err); 

/**  Given a partition handle, a part handle, and an entity handle, return a
 *  flag indicating whether the entity is strictly internal, on a 
 *  boundary, or a ghost.
 *  If part_handle is invalid, an error is returned.
 *  COMMUNICATION:  None.*/
void iMeshP_getEntStatus(iMesh_Instance instance,
                    /*in*/ const iMeshP_PartitionHandle partition_handle, 
                    /*in*/ const iMeshP_PartHandle part_handle, 
                    /*in*/ const iBase_EntityHandle entity_handle, 
                   /*out*/ int* par_status, /* enum iMeshP_EntStatus */
                           int *err); 
void iMeshP_getEntStatusArr(iMesh_Instance instance,
                    /*in*/ const iMeshP_PartitionHandle partition_handle, 
                    /*in*/ const iMeshP_PartHandle part_handle, 
                    /*in*/ const iBase_EntityHandle *entity_handles, 
                    /*in*/ const int entity_handles_size, 
		 /*inout*/ int** par_status, /* enum iMeshP_EntStatus */
                 /*inout*/ int* par_status_allocated, 
                 /*inout*/ int* par_status_size, 
                           int *err); 

/**  Provide information about copies of entities.  
 *  All these functions should work on the local process as well as 
 *  remote processes; entity handles returned are likely but not 
 *  necessarily remote. 
 *  Given a partition handle and an entity handle, return the number 
 *  of copies of the entity in the partition.
 *  COMMUNICATION:  None++.*/
void iMeshP_getNumCopies(iMesh_Instance instance,
                    /*in*/ const iMeshP_PartitionHandle partition_handle, 
                    /*in*/ const iBase_EntityHandle entity_handle, 
                   /*out*/ int *num_copies_ent,
                           int *err); 

/**  Given a partition handle and an entity handle, return the part IDs
 *  of copies of the entity in the partition. (Requested by Onkar.)
 *  COMMUNICATION:  None++.*/
void iMeshP_getCopyParts(iMesh_Instance instance,
                    /*in*/ const iMeshP_PartitionHandle partition_handle, 
                    /*in*/ const iBase_EntityHandle entity_handle, 
                 /*inout*/ iMeshP_Part **part_ids, 
                 /*inout*/ int *part_ids_allocated, 
                   /*out*/ int *part_ids_size, 
                           int *err); 

/**  Given a partition handle and an entity handle, return (remote) entity
 *  handles and part IDs of all copies of the entity.
 *  COMMUNICATION:  None++.*/
void iMeshP_getCopies(iMesh_Instance instance,
                    /*in*/ const iMeshP_PartitionHandle partition_handle, 
                    /*in*/ const iBase_EntityHandle entity_handle, 
                 /*inout*/ iMeshP_Part **part_ids, 
                 /*inout*/ int *part_ids_allocated, 
                   /*out*/ int *part_ids_size, 
                 /*inout*/ iBase_EntityHandle **copies_entity_handles, 
                 /*inout*/ int *copies_entity_handles_allocated, 
                           int *copies_entity_handles_size,
                           int *err); 

/**  Given a partition handle, an entity handle and a part ID, return the
 *  (possibly remote) entity handle of the entity in the specified part.
 *  Return an error if the entity does not exist in the specified part.
 *  COMMUNICATION:  None++.*/
void iMeshP_getCopyOnPart(iMesh_Instance instance,
                    /*in*/ const iMeshP_PartitionHandle partition_handle, 
                    /*in*/ const iBase_EntityHandle entity_handle, 
                    /*in*/ const iMeshP_Part* part_id, 
                   /*out*/ iBase_EntityHandle* copy_entity_handle, 
                           int *err); 


/**  Given a partition handle and an entity handle, return the entity
 *  handle and part ID from the owner of the entity (e.g., the 
 *  entity handle of the copy that has right-to-modify).
 *  COMMUNICATION:  None++.*/
void iMeshP_getOwnerCopy(iMesh_Instance instance,
                    /*in*/ const iMeshP_PartitionHandle partition_handle, 
                    /*in*/ const iBase_EntityHandle entity_handle, 
                   /*out*/ iMeshP_Part *owner_part_id, 
                   /*out*/ iBase_EntityHandle *owner_entity_handle, 
                           int *err); 

/*************************************************
 *******         COMMUNICATION          **********
 *************************************************/
/***********************************
 * Non-blocking calls for off-processor mesh-modification return a request 
 * that indicates whether or not the operation has completed.
 * The following functions receive and process requests.*/

/**\brief  Wait for specified request to complete
 *
 * Wait for the specified request to complete.
 * \param req Request object to wait on
 * \param stat Status of communication request
 * \param err Error returned from function
 *  COMMUNICATION:  Blocking point-to-point.
 */
void iMeshP_Wait(iMesh_Instance instance,
      const iMeshP_PartitionHandle partition_handle,
             /*in*/ iMeshP_RequestHandle req,
            /*out*/ iMeshP_Status *stat,
            /*out*/ int *err);

/**\brief  Wait for any of the specified requests to complete
 *
 * Wait for any of the specified requests to complete.
 * \param req Request objects to wait on
 * \param req_size Number of request objects
 * \param index Index of request in req which was completed
 * \param stat Status of communication request
 * \param err Error returned from function
 *  COMMUNICATION:  Blocking point-to-point.
 */
void iMeshP_WaitAny(iMesh_Instance instance,
         const iMeshP_PartitionHandle partition_handle,
               /*in*/ iMeshP_RequestHandle *req,
               /*in*/ int req_size,
               /*out*/ int *index,
               /*out*/ iMeshP_Status *stat,
               /*out*/ int *err);

/**\brief  Wait for all of the specified requests to complete
 *
 * Wait for all of the specified requests to complete.
 * \param req Request objects to wait on
 * \param req_size Number of request objects
 * \param stat Status of communication request
 * \param err Error returned from function
 *  COMMUNICATION:  Blocking point-to-point.
 */
void iMeshP_WaitAll(iMesh_Instance instance,
         const iMeshP_PartitionHandle partition_handle,
               /*in*/ iMeshP_RequestHandle *req,
               /*in*/ int req_size,
               /*out*/ iMeshP_Status *stat,
               /*out*/ int *err);


/**\brief  Wait for specified request to complete
 *
 * Wait for the specified request to complete.  Returns entities
 * for which information was \em received
 * \param req Request object to wait on
 * \param entities Entities for which information was received
 * \param entities_alloc Allocated size of entities vector
 * \param entities_size Occupied size of entities vector
 * \param err Error returned from function
 *  COMMUNICATION:  Blocking point-to-point.
 */
void iMeshP_WaitEnt(iMesh_Instance instance,
                 const iMeshP_PartitionHandle partition_handle,
                  /*in*/ iMeshP_RequestHandle req,
               /*inout*/ iBase_EntityHandle **out_entities,
               /*inout*/ int *out_entities_alloc,
               /*inout*/ int *out_entities_size,
               /*out*/ int *err);

/**\brief Test for whether specified request has completed
 *
 * Test for whether specified request has completed
 * \param req Request objects to wait on
 * \param flag Returned true if request completed
 * \param stat Status of communication request
 * \param err Error returned from function
 *  COMMUNICATION:  Point-to-point; non-blocking.
 */
void iMeshP_Test(iMesh_Instance instance,
      const iMeshP_PartitionHandle partition_handle,
             /*in*/ iMeshP_RequestHandle req,
            /*out*/ int *flag,
            /*out*/ iMeshP_Status *stat,
            /*out*/ int *err);

/** Poll for requests.  The internals of this function are going to have
 * to cover a lot of ground.  The array in the return is there as a
 * placeholder to tell the application that something interesting / useful
 * has been done to a handle.  This might indicate successful in-migration,
 * a recent change in vertex location, or successful completion of handle
 * matching. 
 * 
 * Returns an array of requests that have been handled.  If
 * the array has a size allocated already, then the implementation stops
 * working when it has generated that many completed requests, even if there
 * are more messages waiting.  The syntax on this call should perhaps be
 * modified somewhat to make it more compatible with an analogous MPI call.
 * 
 * Should also perhaps have its name changed to handleMessages or some such
 * thing, since the function should not only poll for unexpected
 * messages, but also handle expected ones; these is no reason to
 * separately call two different functions to handle non-blocking receives
 * for different types of messages.
 * COMMUNICATION:  non-blocking; point-to-point.
 */
void iMeshP_pollForRequests(iMesh_Instance instance,
           iMeshP_PartitionHandle partition_handle,
           iMeshP_RequestHandle **requests_completed,
			    int *requests_allocated,
			    int *requests_size,
			    int *err);

/********************************************************************
 *******     Requests for off-processor mesh modification      ******
 ********************************************************************/

/*********************************
 *  Add entities to on-process and/or off-process parts:  
 *  Collective-communication version:  
 *  iMeshP_exchEntArrToPartsPar is a collective, non-blocking operation
 *  to be called by all processes in the partition's communicator.  
 *  An iMeshP_RequestHandle is returned; any of the iMeshP_Wait functions can be
 *  used to block until the request is completed.
 *  COMMUNICATION:  Collective.  Non-blocking.*/
void iMeshP_exchEntArrToPartsAll(iMesh_Instance instance,
                       /*in*/  const iMeshP_PartitionHandle partition_handle,
                       /*in*/  const iBase_EntityHandle *entity_handles,
                       /*in*/  const int entity_handles_size,
                       /*in*/  const iMeshP_Part *target_part_ids,
                       /*in*/  int command_code,  /* e.g., MIGRATE,COPY */
                       /*in*/  int update_ghost,  /* e.g., YES,NO */
                      /*out*/  iMeshP_RequestHandle *request,
                               int *err);

/************************************************************************/

/** Request in-migration to a given part of an entity and its upward adjacencies.
 * This is a pull migration.  The
 * entity must be on the part bdry and is identified by local handle, and
 * the implementation handles the rest.  This operation will require multiple
 * rounds of communication, and at some times certain entities may be
 * locked (unavailable for local modification) while info about their
 * remote copies is still in question.
 * It's worth mentioning that implementations will have to take care to
 * migrate tags and parallel set membership as well as adjacency info.
 * CHANGES: Returns a request handle.  Assumes you wouldn't be asking if
 * you didn't need the upward adjacencies as well.
 * COMMUNICATION:  point-to-point, non-blocking, pull. */
void iMeshP_migrateEntity(iMesh_Instance instance, 
			  const iMeshP_PartitionHandle partition_handle,
                          const iMeshP_PartHandle part_handle,
		          const iBase_EntityHandle local_entity_handle,
			  iMeshP_RequestHandle *request,
		          int *err);

/** Update vertex coordinates.  One could argue that we could overload
 * the setVtxCoords function to do this, and maybe we should.  But that
 * obfuscates when communication could occur.  The communication here is
 * push-and-forget.
 * Because this is push-and-forget, no request handle -should- be generated.
 * COMMUNICATION:  point-to-point, non-blocking, push-and-forget.*/
void iMeshP_updateVtxCoords(iMesh_Instance instance, 
			    const iMeshP_PartitionHandle partition_handle,
		            const iBase_EntityHandle local_vertex_handle,
			    int *err);

/** Replace entities.  This refers to changes on the part bdry where the
 * application/service is responsible for ensuring that things are done
 * identically on both sides and that the args are passed in an order that
 * can be matched.  (Specifically, matching new entities should appear in
 * the same order in the call array.)  Communication here could be a
 * two-way push-and-forget, or some variant on push-and-confirm.
 * CHANGES: At Onkar's suggestion, added an offset array (similar to array
 * adjacency requests) so that a single call can easily handle coordination
 * with multiple entities on part-boundary.
 * COMMUNICATION:  point-to-point, non-blocking, push-and-forget. */
void iMeshP_replaceOnPartBdry(iMesh_Instance instance, 
		       const iMeshP_PartitionHandle partition_handle,
		       const iBase_EntityHandle *old_entities,
		       const int old_entities_size,
		       const iBase_EntityHandle *new_entities,
		       const int new_entities_size,
		       const int *offset,
		       const int offset_size,
		       int *err);


/** The ability to create and delete copies is likely
 * to be useful, even though for common topologically-based cases, copies
 * maintainance can (and IMO should) be done automagically, either at
 * migration time or during iMeshP_syncMeshAll.  
 * Communication here is
 * push-and-confirm for creation (so that the original knows ID's of the
 * ghosts), push-and-forget for deletion (in the latter case, no request
 * handle is needed).  I'm assuming here that the closure of a new ghost
 * will be pushed automatically as part of the underlying communication,
 * and that the remote part will clean up the closure as appropriate during
 * deletion.
 * COMMUNICATION:  point-to-point, non-blocking, push.*/
void iMeshP_addGhostOf(iMesh_Instance instance, 
	/* in */       const iMeshP_PartitionHandle partition_handle,
	/* in */       const iMeshP_Part target_part_id,
	/* in */       const iBase_EntityHandle entity_to_copy,
	/* out */      iMeshP_RequestHandle *request,
		       int *err);

void iMeshP_rmvGhostOf(iMesh_Instance instance, 
	/* in */       const iMeshP_PartitionHandle partition_handle,
	/* in */       const iMeshP_Part target_part_id,
	/* in */       const iBase_EntityHandle copy_to_purge,
		       int *err);

/** Done with mesh modification.  This is a blocking call, to get
 * everything up-to-date and back in synch.  Essentially, waits for all
 * message traffic to clear, as well as (possibly) rebuilding a bunch of
 * ghost info that was allowed to go obsolete.*/
void iMeshP_syncMeshAll(iMesh_Instance instance, 
		       const iMeshP_PartitionHandle partition_handle,
		       int *err);
		            
/****************************************************************************/
/* Functions to send Tag data from owning entities to copies.*/

/**\brief  Synchronously send tag data for entity type/topo
 * and part
 *
 * Send tag information for shared entities of specified type/
 * topology to specified part.  
 * The tag data is "pushed" from the
 * owner entities to all copies.
 * This version operates on all
 * shared entities of specified type/topology (or all
 * types/topologies if iMesh_ALL_TYPE or iMesh_ALL_TOPOLOGIES are
 * given).  This function assumes tag handles given on various
 * calling parts are consistent, i.e. they have the same name,
 * data type, size, etc.  This call blocks until communication is
 * completed
 * \param tag Tag handle to exchange
 * \param entity_type Tag data exchanged only for this entity type
 * \param entity_topo Tag data exchanged only for this entity topology
 * \param err Error returned from function
 */
void iMeshP_pushTags(iMesh_Instance instance,
          const iMeshP_PartitionHandle partition_handle,
                /*in*/ iBase_TagHandle tag, 
                /*in*/ int entity_type, 
                /*in*/ int entity_topo, 
                /*out*/ int *err);

/**\brief  Synchronously send tag data for entities and part
 *
 * Send tag information for specified entities, which must
 * already be shared.  
 * The tag data is "pushed" from the
 * owner entities to all copies.
 * This function assumes tag handles given on various
 * calling parts are consistent, i.e. they have the same name,
 * data type, size, etc.  This call blocks until communication is
 * completed
 * \param tag Tag handle to exchange
 * \param entities Owned entities for which to send data
 * \param entities_size Number of entities
 * \param err Error returned from function
 */
void iMeshP_pushTagsEnt(iMesh_Instance instance,
                 const iMeshP_PartitionHandle partition_handle,
                   /*in*/ iBase_TagHandle tag, 
                   /*in*/ const iBase_EntityHandle *entities,
                   /*in*/ int entities_size,
                   /*out*/ int *err);


/**\brief  Asynchronously send tag data for entity type/topo
 * and part
 *
 * Send tag information for shared entities of specified type/
 * topology to specified part.  
 * The tag data is "pushed" from the
 * owner entities to all copies.
 * This version operates on all
 * shared entities of specified type/topology (or all
 * types/topologies if iMesh_ALL_TYPE or iMesh_ALL_TOPOLOGIES are
 * given).  This function assumes tag handles given on various
 * calling parts are consistent, i.e. they have the same name,
 * data type, size, etc.  Applications can call iMeshP_Wait or
 * iMeshP_WaitEnt to block until associated call completes
 * \param tag Tag handle to exchange
 * \param entity_type Tag data exchanged only for this entity type
 * \param entity_topo Tag data exchanged only for this entity topology
 * \param iMeshP_RequestHandle Request object used in call to
 *    iMeshP_Wait or iMeshP_WaitEnt
 * \param err Error returned from function
 */
void iMeshP_iPushTags(iMesh_Instance instance,
                 const iMeshP_PartitionHandle partition_handle,
                 /*in*/ iBase_TagHandle tag, 
                 /*in*/ int entity_type, 
                 /*in*/ int entity_topo, 
                /*out*/ iMeshP_RequestHandle *req,
                /*out*/ int *err);

/**\brief  Asynchronously send tag data for entities and part
 *
 * Sedn tag information for specified entities, which must
 * already be shared.  
 * The tag data is "pushed" from the
 * owner entities to all copies.
 * This function assumes tag handles given on various
 * calling parts are consistent, i.e. they have the same name,
 * data type, size, etc.  Applications can call iMeshP_Wait or
 * iMeshP_WaitEnt to block until associated call completes
 * \param tag Tag handle to exchange
 * \param entities Owned entities for which to send data
 * \param entities_size Number of entities
 * \param iMeshP_RequestHandle Request object used in call to
 *    iMeshP_Wait or iMeshP_WaitEnt
 * \param err Error returned from function
 */
void iMeshP_iPushTagsEnt(iMesh_Instance instance,
                     const iMeshP_PartitionHandle partition_handle,
                    /*in*/ iBase_TagHandle tag, 
                    /*in*/ const iBase_EntityHandle *entities,
                    /*in*/ int entities_size,
                   /*out*/ iMeshP_RequestHandle *req,
                   /*out*/ int *err);

/*  GHOST ENTITY SUPPORT */
/*
The triplet describing a ghosting "rule" (ghost dim, bridge dim, #
layers) will be stored in the partition and re-established by the end
of iMeshP_syncPartitionAll and iMeshP_syncMeshAll.  Implementations can
choose to keep ghosting consistent throughout mesh modification, but ghosts
are not required to be consistent until the end of these two functions.

The number of layers specified is with respect to the global mesh;
that is, ghosting may extend beyond a single neighboring processor if the
number of layers is high.

iMeshP_createGhostEntities is a preprocessing function.  It is
cumulative; that is, multiple calls may add more ghosts, not eliminate
previous ghosts.  
*/

/* \brief Create ghost entities between parts/processors.
  * Ghost entities are specified similar to 2nd-order adjacencies, i.e.
  * through a "bridge" dimension.  Number of layers is measured from
  * the inter-processor interfaces.  For example, to get 2 layers of hex
  * entities in the ghost layer, measured from faces on the interface,
  * use ghost_dim=3, bridge_dim=2, and num_layers=2.
  *
  * Ghost information is cached on the partition.
  *
  * \param instance iMesh instance
  * \param partition_handle Partition on which to create ghosts
  * \param ghost_dim Dimension of entities ghosted
  * \param bridge_dim Dimension through which bridge adjacencies are found
  * \param num_layers Number of layers of ghost entities
  * \param err Error returned from function
  * \param include_copies  Create ghosts of non-owned part boundary entities?
  *         (1/0)
  * COMMUNICATION:  Collective.  Blocking.
  */
void iMeshP_createGhostEnts(/*in*/ iMesh_Instance instance,
                       /*in*/ iMeshP_PartitionHandle partition_handle,
                       /*in*/ int ghost_dim,
                       /*in*/ int bridge_dim,
                       /*in*/ int num_layers,
                       /*in*/ int include_copies,
                       /*out*/ int *err);

/* \brief Delete all ghost entities between parts/processors.
  * \param instance iMesh instance
  * \param partition_handle Partition on which ghost entities are deleted
  * \param err Error returned from function
  * COMMUNICATION:  None.
  */
void iMeshP_deleteGhostEnts(/*in*/ iMesh_Instance instance,
                       /*in*/ iMeshP_PartitionHandle partition_handle,
                       /*out*/ int *err);


/* \brief Return information about all ghosting on a partition.
  * \param instance iMesh instance
  * \param partition_handle Partition on which ghost entities are deleted
  * \param num_ghost_rules  Number of ghosting rules currently registered 
  * \param ghost_dim Dimension of ghost entities
  * \param bridge_dim Dimension of bridge entities
  * \param num_layers Number of layers of ghost entities
  * \param err Error returned from function
  * COMMUNICATION:  None.
  */
void iMeshP_ghostEntInfo(/*in*/ iMesh_Instance instance,
                    /*in*/ iMeshP_PartitionHandle partition_handle,
                    /*out*/ int *num_ghost_rules,
                    /*inout*/ int *ghost_rules_allocated, 
                    /*inout*/ int *ghost_rules_size, 
                    /*out*/ int **ghost_dim,
                    /*out*/ int **bridge_dim,
                    /*out*/ int **num_layers,
                    /*out*/ int *err);

/********************************************************************
 *******    FILE I/O                                           ****** 
 ********************************************************************/
/* \brief iMeshP file I/O closely aligns with iMesh file I/O.  The major
 * change is the addition of a iMeshP_PartitionHandle argument to both
 * iMeshP_load and iMeshP_save, enabling I/O from parallel processes.
 * For now, individual implementations will support different sets of
 * options; Tim and Ken will work to unify the options by SC08.
 */

/* \brief iMeshP_load: function to populate a mesh instance and 
 * a partition by reading data from files.
 * The application creates both the mesh instance and the partition
 * handle before calling iMeshP_load.  iMeshP_load then reads the
 * specified file, inserts entities into the mesh instance, constructs
 * parts within the partition, and inserts entities into the parts.
 * Optional capabilities of iMeshP_load include computing an initial
 * partition (e.g., if a serial mesh file without part assignments is read)
 * and creating ghost entities as requested by the application.
 * Options allow n>=1 files on p processes.
 */

void iMeshP_load(iMesh_Instance instance,
                 const iMeshP_PartitionHandle partition,
                 const iBase_EntitySetHandle entity_set_handle,
                 const char *name, 
                 const char *options,
                 int *err, 
                 int name_len, 
                 int options_len);

/* \brief iMeshP_save: function to write data from a mesh instance and 
 * a partition to files.
 * iMeshP_save writes mesh and partition data to the
 * specified file.
 * Options allow n>=1 files on p processes.
 */
void iMeshP_save(iMesh_Instance instance,
                 const iMeshP_PartitionHandle partition,
                 const iBase_EntitySetHandle entity_set_handle,
                 const char *name, 
                 const char *options,
                 int *err, 
                 const int name_len, 
                 int options_len);


/*
------------------------------------------------
Major Items left to do:
-  Support for multiple partitions.
   We discussed designating a given partition as 
   the "active" partition; i.e., the partition that is actually used in 
   the distribution of mesh data in distributed memory.  We were concerned 
   that when multiple partitions were used, multiple copies of mesh 
   entities would be needed to fully support multiple partitions at the 
   same time.  Designating one partition as "active" would store data 
   with respect to only one partition.
-  File I/O support.
   Need common set of options to allow interoperability.
   Support single files, N << P files on P processes, and P files.
   Support reading and writing partition information.
   Support initial parallel partitioning of serial file data.
   Support storing mapping of parts to processes in files.

------------------------------------------------
Minor Items left to do:
-  Determine which capabilities need both "getNumX" and "getX" functions.
   That is, when would an application need "getNumX" to allocate memory
   for "getX" or separately from "getX".  When could we use only "getX"
   and return numX as a byproduct.

-  Determine with functions need "Ent" and "EntArr" versions, or whether
   we should adopt only the more general "EntArr" version.

-  Determine whether to revise iMeshP_createPartition to make it less MPI 
   specific.  We don't want to require applications to link with MPI if the 
   implementation doesn't require it.  We may define an ITAPS_Comm name 
   typedef'ed appropriately.

-  iMeshP_getOwnerCopy could be achieved by calling iMeshP_getOwnerPart
   followed by iMeshP_getCopyOnPart.  Do we want to keep iMeshP_getOwnerCopy?

-  Need function to receive tag data from part-boundary entities in owner.
   Possible options:  return the tag data values received directly, or 
   include a mathematical operation (similar to MPI_SUM). 9/15/08

------------------------------------------------
Comments and resolved questions:  

- Applications will build partitions by (1) creating a partition handle
  on each process to be included in the partition; (2) adding parts to 
  the partition handle within the process; (3) populating the parts with 
  entities, and (4) calling iMeshP_syncPartitionAll to allow the 
  implementation to compute global data for the partition.

- For now, we will not include an iterator over local (to a
  process) parts within a partition.  If it is needed, it can be added
  later.

- We will not provide capability to move entire parts to new
  processes;  instead, the new process must create the part in its
  partition handle and then receive (perhaps in bulk) the entities to 
  populate the part.  In other words, parts can be added to only a local 
  partition handle.

- Currently, iMesh doesn't have the functionality to get entities or 
  entity sets by type and tag in serial.  Should it?  
  Many people said it would be useful; others said it could be costly
  (in parallel) or numerically difficult (for floating point values).
  This issue is an iMesh issue, not a parallel interface issue, so
  for this document, the issue is resolved.  The resolution:  If 
  iMesh adopts this capability, we will add it to the
  parallel interface.

- We will not include functions that return all entities with 
  given characteristics within a partition; the memory use of these
  functions can be large.  Instead, we will return entity information
  with respect to parts and/or mesh instances.  If the user wants such
  info, he should go through the mechanics of gathering it himself so
  that he is painfully aware of how much memory he is allocating.
  Removed the following global queries:
  + All tag names over the partition;
  + All entities in this partition having a given type, tag and/or 
    tag name.
  + All entity sets in this partition having a given
    type, tag and/or tag name.

- We will not include functions that return information about each
  part and/or process in a partition.  Such functions limit memory 
  scalability for large numbers of parts.  If the user wants such
  info, he should go through the mechanics of gathering it himself so
  that he is painfully aware of how much memory he is allocating.
  Removed the following global queries:
  + The number of entities in each part of the partition;
  + The number of entity sets in each part of the partition;
  + The number of entities with given type, tag, and/or
    tag name in each part of the partition;
  + The number of entity sets with given type, tag, 
    and/or tag name in each part of the partition;
  + All tag names in each part of the partition;

- For functions that replace a set handle with a part handle, return
  all appropriate entities in a part, whether they are owned or are 
  copies.  The application can test for ownership if needed.

- Part assignments computed with respect to a set of 
  entities induce part assignments to adjacent entities in an
  implementation-dependent fashion.  That is, if a partition is computed
  with respect to regions, queries about ownership of faces and vertices
  are valid.

------------------------------------------------
Discussed but unresolved questions:

- We discussed adding functions that give
  hints to an implementation about which data mappings the application 
  will use, allowing the implementation to pre-compute them if it chooses 
  to.  The example discussed was mapping between entities and parts, but 
  other examples in iMesh may also exist.

- We discussed adding an iterator over entities 
  with given type/topology in a set or part.  We have non-iterator 
  functionality, but not an iterator.  
  KDD:  Is this true?  What is iMesh_initEntIter (and its analogous
  KDD:  iMeshP_initEntIter)?

- We discussed storing in a partition 
  information about which "objects" were used in computing the partition.  
  These objects can be single entities or groups of entities.
  KDD:  Perhaps this capability should be part of the load-balancing service.

- We discussed designating a given partition as 
  the "active" partition; i.e., the partition that is actually used in 
  the distribution of mesh data in distributed memory.  We were concerned 
  that when multiple partitions were used, multiple copies of mesh 
  entities would be needed to fully support multiple partitions at the 
  same time.  Designating one partition as "active" would store data 
  with respect to only one partition.

------------------------------------------------
Not-yet-discussed, unresolved questions

Entity questions:
- From Carl:  "getTag*Operate: Again, we haven't got this in serial.  Does 
  the existence of such operations imply that we expect to implement 
  fields as tags? (Because that wasn't what I was assuming about field 
  implementations at all, personally...)  Note that I'm not opposed to 
  this sort of global reduction operation, I just wonder whether it'll see 
  use outside of field-like situations.  If not, then it should be in 
  parallel fields, not parallel mesh, and usage for 
  fields-implemented-as-tags should be handled there."
*/

/**********************************/
/* NOTES FROM BOOTCAMP MARCH 2008 */
/**********************************/
/*
-  Changed getPartRank to getRankOfPart.  (Carl)
-  Made sure iMeshP_getNumOfTypeAll and iMeshP_getNumOfTopoAll were
documented as collective operations.  (Carl)
-  Changed suffix "Par" to "All".  (Lori)
-  Added iMeshP_testPart() to test status of part handle, returning
LOCAL, REMOTE, or INVALID.  (Mark M, Lori).
6/25/08:  Removed this function since part handles can no longer be remote.
If an application wants to test the validity of a part handle, it can try
to compute its Part ID.
-  Changed iMeshP_addCopyOf and iMeshP_rmvCopyOf back to
iMeshP_addGhostOf and iMeshP_rmvGhostOf.  If we wanted to use these
functions for adding boundary copies, we'd have to include a list of
already existing remote copies in the arguments, as well as
communicate with parts already owning copies to let them know a ghost
copy has been made.  Actually, this raises an interesting question:
does a boundary copy need to know about all ghost copies of it?
-  Change getEntParStatus to getEntStatus.  (Lori)
-  Changed sendEntArrToPartsPar to exchEntArrToPartsAll.  (Lori,Tim)


Parts and Processes:
-  Martin argued for consecutive unique Part IDs in addition to or
instead of Part handles.  He will send use cases.   If we decide to
add them back to the interface, we could compute them in
iMeshP_syncPartitionAll rather than in iMeshP_createPart.  That is, an
application couldn't access them until after iMeshP_syncPartitionAll.
6/25/08:  On follow-up, Martin couldn't recall why having consecutive
PartIDs was necessary.  While we all agree they are conceptually nice,
they are difficult to implement and not really necessary.  Part IDs will
be globally unique but not necessarily consecutive.
-  Are part handles globally unique?  They probably need to be
globally unique in order for them to be useful as remote part
handles.  Also, does the process rank need to be encoded in the part
handle in order to map from parts to processes for communication?
6/25/08:  DECIDED:  We will have globally unique part IDs.  Part handles
will be valid for only local parts.  Accessing remote parts must be done
via Part IDs.
-  If in iMeshP_syncPartitionAll, we computed a mapping from part
handles to integers 0,..., k-1, we could store only ranges of
integers to achieve the part-to-process and process-to-parts mappings;
this would require O(P) storage per process for P processes.
6/5/08:  DECIDED:  Do not need getPartOnRank or getNumPartOnRank.  These
functions were troublesome due to their storage or communication requirements.
We decided to remove them.
-  Alternatively, the mapping of all parts to processes can be stored
in O(k) total memory, distributed across processors (e.g., a
distributed data directory) but interrogating the directory requires
communication.  
6/5/08:  See note above.
-  iMeshP_getPartsOnRank created discussion and needs to be resolved.
IMeshP_getPartsOnRank would likely require either O(k) storage per
process for k parts or communication.  For other points, please see
Mark M's 3/12/08 email.  
6/5/08:  See note above.

CreateEnt:
-  Carl asked if we should have a version of createEnt that accepts a
part handle.  Should this function be used only for creating owned
entities?   How do you envision creating part boundary entities when a
parallel mesh is initially loaded?  

Ghost entities:
-  We currently have a mechanism only for pushing ghosts onto other
parts.  Will we want a mechanism for pulling them, too?  (E.g., a
part says, "I want ghosts for this entity.")

PartNbor functions:
-  Did we agree to remove the entity type from these functions?  That
is, do we want to return the part IDs for all parts that have
any copies?  The entity type was intended to allow us to get the part
IDs for all parts that have copies of a given type (perhaps
ALL_TYPES).  

Functions handling both Parts and Entity Sets:
-  Tim said these function names (e.g., iMeshP_getNumOfType,
iMeshP_getAllVtxCoord) are too close to existing iMesh function
names, even though the argument lists would be different.  He agreed
to email suggestions for better names.

Copies:
-  Functions getNumCopies, getCopies, getCopyParts, and getCopyOnPart
have different behavior for ghost and part-boundary copies.  Ghosts
will return only itself and its owner in getCopies; part-boundary
entities will return copies on other parts as well.
-  Tim envisions applications (e.g., FETI methods) updating tag data
in their copies that they would like to accumulate back to the
owner.  Xiaolin said that he writes in his ghosts, but doesn't send
those values back to the owner.  Currently, we have the ability 
to send tag data only from owners to ghosts.  Tim will look at this issue
and propose a solution.

Communication:
-  Although we should think in terms of parts, communication really
occurs with respect to processes.  We need to make sure all our
communication routines make sense with respect to both processes and
parts, and perhaps, revise their descriptions.  Also, if we send to
parts, the implementation really does need the mapping of parts to
processes. 

Entity Owner/Status Queries:
-  Should we combine functions getEntOwnerPart and getEntStatus into
one function?  Or should we combine functions getOwnerCopy and
getEntOwner into one function?  Or should we remove getOwnerCopy and
make applications call getOwnerPart followed by getCopyOnPart?

Reducing storage:
-  Mark Miller proposed allowing the user to specify the amount of 
copying done by the implementation, depending on applications' needs.
For example, for a static viz application, every boundary entity may not
need to know where all its copies are, so the implementation would not
have to store them.  Can the implementations accept a flag advising them how
much copying is needed?  If so, the implementations could choose to 
optimize storage or ignore the flag.
*/



////////////////////////////////////////////////////////////
//SVN File Information
//
//  $Author:$
//  $Date:$
//  $Revision:$
////////////////////////////////////////////////////////////


#ifdef __cplusplus
} // extern "C" 
#endif

#endif /* defined(iMeshP_h) */

