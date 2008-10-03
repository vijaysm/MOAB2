#include "iMeshP.h"
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBCN.hpp"
#include "MeshTopoUtil.hpp"
#include "FileOptions.hpp"
#include "MBParallelComm.hpp"
#include "MBParallelConventions.h"

#ifdef USE_MPI    
#include "mpi.h"
#endif

#define MBI reinterpret_cast<MBInterface*>(instance)

#define RETURN(a) {iMesh_LAST_ERROR.error_type = a; *err = a;return;}
#define iMesh_processError(a, b) {sprintf(iMesh_LAST_ERROR.description, "%s", b); iMesh_LAST_ERROR.error_type = a; *err = a;}

// TAG_CHECK_SIZE is like CHECK_SIZE except it checks for and makes the allocated memory
// size a multiple of sizeof(void*), and the pointer is assumed to be type char*
#define TAG_CHECK_SIZE(array, allocated, size)  \
  if (0 != allocated && NULL != array && allocated < (size)) {\
    iMesh_processError(iBase_MEMORY_ALLOCATION_FAILED, \
          "Allocated array not large enough to hold returned contents.");\
    RETURN(iBase_MEMORY_ALLOCATION_FAILED);\
  }\
  if (NULL == array || allocated == 0) {\
    allocated=(size); \
    if (allocated%sizeof(void*) != 0) allocated=((size)/sizeof(void*)+1)*sizeof(void*);\
    array = (char*)malloc(allocated); \
    if (NULL == array) {iMesh_processError(iBase_MEMORY_ALLOCATION_FAILED, \
          "Couldn't allocate array.");RETURN(iBase_MEMORY_ALLOCATION_FAILED); }\
  }
#define HANDLE_ARRAY_PTR(array) reinterpret_cast<MBEntityHandle*>(array)
#define CONST_HANDLE_ARRAY_PTR(array) reinterpret_cast<const MBEntityHandle*>(array)
#define TAG_HANDLE(handle) reinterpret_cast<MBTag>(handle)
#define CONST_TAG_HANDLE(handle) static_cast<const MBTag>(handle)
#define ENTITY_HANDLE(handle) reinterpret_cast<MBEntityHandle>(handle)
#define CONST_ENTITY_HANDLE(handle) reinterpret_cast<const MBEntityHandle>(handle)
#define RANGE_ITERATOR(it) reinterpret_cast<RangeIterator*>(it)

const iBase_ErrorType iBase_ERROR_MAP[] = 
{
  iBase_SUCCESS, // MB_SUCCESS = 0,
  iBase_INVALID_ENTITY_HANDLE, // MB_INDEX_OUT_OF_RANGE,
  iBase_INVALID_ENTITY_TYPE, // MB_TYPE_OUT_OF_RANGE,
  iBase_MEMORY_ALLOCATION_FAILED, // MB_MEMORY_ALLOCATION_FAILED,
  iBase_INVALID_ENTITY_HANDLE, // MB_ENTITY_NOT_FOUND,
  iBase_NOT_SUPPORTED, // MB_MULTIPLE_ENTITIES_FOUND,
  iBase_TAG_NOT_FOUND, // MB_TAG_NOT_FOUND,
  iBase_FILE_NOT_FOUND, // MB_FILE_DOES_NOT_EXIST,
  iBase_FILE_WRITE_ERROR, // MB_FILE_WRITE_ERROR,
  iBase_NOT_SUPPORTED, // MB_NOT_IMPLEMENTED,
  iBase_TAG_ALREADY_EXISTS, // MB_ALREADY_ALLOCATED,
  iBase_FAILURE // MB_FAILURE};
};

extern iBase_Error iMesh_LAST_ERROR;

template <typename ArrType> inline bool
allocate_itaps_array( ArrType*& array, int& allocated, int& size, int requested )
{
  size = requested;
  if (allocated) {
    return (allocated >= requested);
  }
  else {
    array = (ArrType*)malloc( requested * sizeof(ArrType) );
    allocated = requested;
    return (array != 0);
  }
}

inline int allocate_itaps_array_failed()
{
  strcpy( iMesh_LAST_ERROR.description, 
          "Insufficient allocated array size or insufficient "
          "memory to allocate array." );
  iMesh_LAST_ERROR.error_type = iBase_MEMORY_ALLOCATION_FAILED;
  return iBase_MEMORY_ALLOCATION_FAILED;
}

#define ALLOCATE_ARRAY( NAME, SIZE ) \
  if (!allocate_itaps_array( *NAME, *NAME##_allocated, *NAME##_size, (SIZE) )) \
    RETURN(allocate_itaps_array_failed())

#ifdef __cplusplus
extern "C" {
#endif

void iMeshP_createPartitionAll( iMesh_Instance /*instance*/,
                        /*in*/  MPI_Comm /*communicator*/,
                        /*out*/ iMeshP_PartitionHandle */*partition_handle*/,
                                int *err )
{
  RETURN(iBase_NOT_SUPPORTED);
}

//  Given a partition handle, a part handle, and an entity handle, return a
//  flag indicating whether the entity is strictly internal, on a 
//  boundary, or a ghost.
//  If part_handle is remote, an error is returned.
//  COMMUNICATION:  None.
  void iMeshP_getEntStatus(iMesh_Instance instance,
                           /*in*/ const iMeshP_PartitionHandle,
                           /*in*/ const iMeshP_PartHandle, 
                           /*in*/ const iBase_EntityHandle entity_handle, 
                           /*out*/ int* par_status, // Values=INTERNAL,BOUNDARY,GHOST
                           int *err) 
  {
    MBParallelComm pc(MBI);
    unsigned char pstatus;
    MBErrorCode result = MBI->tag_get_data(pc.pstatus_tag(), 
                                           CONST_HANDLE_ARRAY_PTR(&entity_handle), 1, 
                                           &pstatus);
    if (MB_SUCCESS != result) RETURN(result);
    *par_status = iMeshP_INTERNAL;
    if (pstatus & PSTATUS_GHOST) *par_status = iMeshP_GHOST;
    else if (pstatus & PSTATUS_INTERFACE) *par_status = iMeshP_BOUNDARY;

    RETURN(iBase_SUCCESS);
  }
  
  void iMeshP_getEntStatusArr(iMesh_Instance instance,
                              /*in*/ const iMeshP_PartitionHandle,
                              /*in*/ const iMeshP_PartHandle part_handle, 
                              /*in*/ const iBase_EntityHandle *entity_handles, 
                              /*in*/ const int entity_handles_size, 
                              /*inout*/ int** par_status, // Values=INTERNAL,BOUNDARY,GHOST
                              /*inout*/ int* par_status_allocated, 
                              /*inout*/ int* par_status_size, 
                              int *err) 
  {
    MBParallelComm pc(MBI);
    std::vector<unsigned char> pstatus(entity_handles_size);
    MBErrorCode result = MBI->tag_get_data(pc.pstatus_tag(), 
                                           CONST_HANDLE_ARRAY_PTR(entity_handles), 
                                           entity_handles_size,
                                           &pstatus[0]);
    if (MB_SUCCESS != result) RETURN(result);

    ALLOCATE_ARRAY( par_status, entity_handles_size );
  
    for (int i = 0; i < entity_handles_size; i++) {
      if (!pstatus[i]) (*par_status)[i] = iMeshP_INTERNAL;
      else if (pstatus[i] & PSTATUS_GHOST) (*par_status)[i] = iMeshP_GHOST;
      else if (pstatus[i] & PSTATUS_INTERFACE) (*par_status)[i] = iMeshP_BOUNDARY;
    }

    RETURN(iBase_SUCCESS);
  }

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
                             int *err) 
  {
    MBParallelComm *pc = MBParallelComm::get_pcomm(MBI, 0);
    if (!pc) RETURN(iBase_ERROR_MAP[MB_FAILURE]);

    MBRange part_sets;
  
    ALLOCATE_ARRAY( part_handles, pc->partition_sets().size() );
    MBRange::iterator rit;
    int i;
    for (i = 0, rit = pc->partition_sets().begin(); 
         rit != pc->partition_sets().end(); rit++, i++)
      (*part_handles)[i] = reinterpret_cast<iMeshP_PartHandle>(*rit);
  
    RETURN(iBase_SUCCESS);
  }
    
  void iMeshP_getPartsArrOnRank(iMesh_Instance instance,
                                const iMeshP_PartitionHandle partition_handle,
                                /*in*/    const int *rank,
                                /*in*/    const int rank_size,
                                /*inout*/ iMeshP_PartHandle **part_handles, 
                                /*inout*/ int *part_handles_allocated, 
                                /*out*/   int *part_handles_size, 
                                int *err) 
  {
    MBParallelComm *pc = MBParallelComm::get_pcomm(MBI, 0);
    if (!pc) RETURN(iBase_ERROR_MAP[MB_FAILURE]);

    if (rank[0] != pc->proc_config().proc_rank() || rank_size > 1) {
      RETURN(iBase_ERROR_MAP[MB_NOT_IMPLEMENTED]);
    }
  
    iMeshP_getPartsOnRank(instance, partition_handle, rank[0],
                          part_handles, part_handles_allocated, part_handles_size,
                          err);
  }

//  Provide entity categorization within part.
//  Given a partition handle, a part handle, and an entity handle, return a
//  flag indicating whether the entity is owned by the part.
//  If part_handle is remote, an error is returned.
//  COMMUNICATION:  None.
  void iMeshP_isEntOwner(iMesh_Instance instance,
                         /*in*/ const iMeshP_PartitionHandle,
                         /*in*/ const iMeshP_PartHandle part_handle, 
                         /*in*/ const iBase_EntityHandle entity_handle, 
                         /*out*/ int* is_owner, 
                         int *err) 
  {
    MBParallelComm pc(MBI);
    unsigned char pstatus;
    MBErrorCode result = MBI->tag_get_data(pc.pstatus_tag(), 
                                           CONST_HANDLE_ARRAY_PTR(&entity_handle), 1, 
                                           &pstatus);
    if (MB_SUCCESS != result) RETURN(result);

    if (pstatus & PSTATUS_NOT_OWNED) *is_owner = 0;
    else *is_owner = 1;

    RETURN(iBase_SUCCESS);
  }

  void iMeshP_isEntOwnerArr(iMesh_Instance instance,
                            /*in*/ const iMeshP_PartitionHandle,
                            /*in*/ const iMeshP_PartHandle part_handle, 
                            /*in*/ const iBase_EntityHandle *entity_handles, 
                            /*in*/ const int entity_handles_size, 
                            /*inout*/ int** is_owner, 
                            /*inout*/ int* is_owner_allocated, 
                            /*inout*/ int* is_owner_size, 
                            int *err)
  {
    MBParallelComm pc(MBI);
    std::vector<unsigned char> pstatus(entity_handles_size);
    MBErrorCode result = MBI->tag_get_data(pc.pstatus_tag(), 
                                           CONST_HANDLE_ARRAY_PTR(entity_handles), 
                                           entity_handles_size,
                                           &pstatus[0]);
    if (MB_SUCCESS != result) RETURN(result);

    ALLOCATE_ARRAY( is_owner, entity_handles_size );
  
    for (int i = 0; i < entity_handles_size; i++) {
      if (pstatus[i] & PSTATUS_NOT_OWNED) (*is_owner)[i] = 0;
      else (*is_owner)[i] = 1;
    }

    RETURN(iBase_SUCCESS);
  }

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
                           int *err) 
  {
    MBParallelComm pc(MBI);
    int shared_proc;
    MBErrorCode result = MBI->tag_get_data(pc.sharedp_tag(), 
                                           CONST_HANDLE_ARRAY_PTR(&entity_handle), 1,
                                           &shared_proc);
    if (MB_SUCCESS != result) {
      RETURN(iBase_ERROR_MAP[result]);
    }
  
    if (-1 != shared_proc) {
      *num_copies_ent = 1;
      RETURN(iBase_SUCCESS);
    }
  
    static int tag_size = 0;
    if (!tag_size) {
      result = MBI->tag_get_size(pc.sharedps_tag(), tag_size);
      if (MB_SUCCESS != result) {
        RETURN(iBase_ERROR_MAP[result]);
      }
    }
    static std::vector<int> shared_procs(tag_size);

    result = MBI->tag_get_data(pc.sharedps_tag(), 
                               CONST_HANDLE_ARRAY_PTR(&entity_handle), 1,
                               &shared_procs[0]);
    if (MB_SUCCESS != result) {
      RETURN(iBase_ERROR_MAP[result]);
    }
  
    *num_copies_ent = std::find(shared_procs.begin(), shared_procs.end(), -1) - shared_procs.begin();

    RETURN(iBase_SUCCESS);
  }

//  Given a partition handle and an entity handle, return (remote) entity
//  handles and part handles of all copies of the entity.
//  COMMUNICATION:  None++.
  void iMeshP_getCopies(iMesh_Instance instance,
                        /*in*/ const iMeshP_PartitionHandle partition_handle, 
                        /*in*/ const iBase_EntityHandle entity_handle, 
                        /*inout*/ iMeshP_Part **part_ids, 
                        /*inout*/ int *part_ids_allocated, 
                        /*out*/ int *part_ids_size, 
                        /*inout*/ iBase_EntityHandle **copies_entity_handles, 
                        /*inout*/ int *copies_entity_handles_allocated, 
                        /*inout*/ int *copies_entity_handles_size, 
                        int *err)
  {
    MBParallelComm pc(MBI);
    MBEntityHandle shared_handle;
    MBErrorCode result = MBI->tag_get_data(pc.sharedh_tag(), 
                                           CONST_HANDLE_ARRAY_PTR(&entity_handle), 1,
                                           &shared_handle);
    if (MB_SUCCESS != result) {
      RETURN(iBase_ERROR_MAP[result]);
    }
  
    std::vector<iMeshP_PartHandle> part_handles;
    if (0 != shared_handle) {
      part_handles.resize( 1 );
      ALLOCATE_ARRAY(copies_entity_handles, 1);
      part_handles[0] = 0;
      (*copies_entity_handles)[0] = reinterpret_cast<iBase_EntityHandle>(shared_handle);
    }
    else {
  
      static int tag_size = 0;
      if (!tag_size) {
        result = MBI->tag_get_size(pc.sharedhs_tag(), tag_size);
        if (MB_SUCCESS != result) {
          RETURN(iBase_ERROR_MAP[result]);
        }
      }
      static std::vector<MBEntityHandle> shared_handles(tag_size);

      result = MBI->tag_get_data(pc.sharedhs_tag(), 
                                 CONST_HANDLE_ARRAY_PTR(&entity_handle), 1,
                                 &shared_handles[0]);
      if (MB_SUCCESS != result) {
        RETURN(iBase_ERROR_MAP[result]);
      }

      int index = std::find(shared_handles.begin(), shared_handles.end(), -1) - shared_handles.begin();

      part_handles.resize(index+1, 0);
      ALLOCATE_ARRAY( copies_entity_handles, index+1 );
      std::copy(&shared_handles[0], &shared_handles[index], 
                HANDLE_ARRAY_PTR(*copies_entity_handles));
    }
    
    ALLOCATE_ARRAY(part_ids, part_handles.size());
    iMeshP_getPartIdsFromPartHandlesArr( instance, partition_handle, 
                                         &part_handles[0], part_handles.size(),
                                         part_ids, part_ids_allocated,
                                         part_ids_size, err );
    *part_ids_size = part_handles.size();
    RETURN(*err);
  }
  

#ifdef __cplusplus
} // extern "C"
#endif
