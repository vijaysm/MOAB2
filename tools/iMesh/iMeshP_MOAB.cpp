#include "iMeshP.h"
#include "iMesh_MOAB.hpp"
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBCN.hpp"
#include "MeshTopoUtil.hpp"
#include "FileOptions.hpp"
#include "MBParallelComm.hpp"
#include "MBParallelConventions.h"

#define IS_BUILDING_MB
#include "MBInternals.hpp"
#undef IS_BUILDING_MB

#include <assert.h>

#ifdef USE_MPI    
#include "mpi.h"
#endif

/********************* Error Handling **************************/

#define FIXME printf("Warning: function has incomplete implementation: %s\n", __func__ )


/********************* Part IDs **************************/

static inline
unsigned PART_ID( unsigned rank, unsigned local_id ) {
  const unsigned rank_bits = sizeof(unsigned)*4 - 1;
  const unsigned id_bits = sizeof(unsigned)*4;
  assert(local_id < (1 << id_bits));
  assert(rank < (1 << rank_bits));
  return (rank << id_bits) | local_id;
}

static inline 
unsigned RANK( unsigned part_id ) {
  const unsigned id_bits = sizeof(unsigned)*4;
  return (part_id >> id_bits);
}

static inline
unsigned LOCAL_ID( unsigned part_id ) {
  const unsigned id_bits = sizeof(unsigned)*4;
  const unsigned mask = ~((1 << id_bits) - 1);
  return (part_id & mask);
}


/******** Type-safe casting between MOAB and ITAPS types *********/

#ifndef TEMPLATE_FUNC_SPECIALIZATION
// if no template specializtion, disable some type checking
template <typename T, typename S> inline
T itaps_cast( S handle )
{ 
  assert(sizeof(S) >= sizeof(T)); 
  return reinterpret_cast<T>(handle);
}
#else

// basic template method : only works to cast to equivalent types (no-op)
template <typename T, typename S> inline
T itaps_cast( S h )
{ return h; }
// verify size and do reinterpret cast
template <typename T> inline T itaps_cast_internal_( MBEntityHandle h )
{
  assert(sizeof(T) >= sizeof(MBEntityHandle));
  return reinterpret_cast<T>(h);
}
// verify size and do reinterpret cast
template <typename T> inline MBEntityHandle* itaps_cast_ptr_( T* h )
{
  assert(sizeof(T) == sizeof(MBEntityHandle));
  return reinterpret_cast<MBEntityHandle*>(h);
}
// verify size and do reinterpret cast
template <typename T> inline const MBEntityHandle* itaps_cast_const_ptr_( const T* h )
{
  assert(sizeof(T) == sizeof(MBEntityHandle));
  return reinterpret_cast<const MBEntityHandle*>(h);
}
// verify set-type handle before cast
template <typename T> inline T itaps_set_cast_( MBEntityHandle h )
{
  assert(TYPE_FROM_HANDLE(h) == MBENTITYSET);
  return itaps_cast_internal_<T>(h);
}

// define conversion routines between itaps handle and MBEntityHandle types
#define DECLARE_ALLOWED_ITAPS_CONVERSION( ITAPS_HANDLE_TYPE ) \
  template <> inline \
  ITAPS_HANDLE_TYPE \
  itaps_cast<ITAPS_HANDLE_TYPE,MBEntityHandle>( MBEntityHandle h ) \
  { return itaps_cast_internal_<ITAPS_HANDLE_TYPE>(h); } \
  \
  template <> inline \
  MBEntityHandle \
  itaps_cast<MBEntityHandle,ITAPS_HANDLE_TYPE>( ITAPS_HANDLE_TYPE handle ) \
  { return reinterpret_cast<MBEntityHandle>(handle); } \
  \
  template <> inline \
  MBEntityHandle* \
  itaps_cast<MBEntityHandle*,ITAPS_HANDLE_TYPE*>( ITAPS_HANDLE_TYPE* ptr ) \
  { return itaps_cast_ptr_(ptr); } \
  \
  template <> inline \
  const MBEntityHandle* \
  itaps_cast<const MBEntityHandle*,const ITAPS_HANDLE_TYPE*>( const ITAPS_HANDLE_TYPE* ptr ) \
  { return itaps_cast_const_ptr_(ptr); }


// define conversion routines between itaps handle and MBEntityHandle types
// but limit to MBEntityHandle for MBENTITYSET type.
#define DECLARE_ALLOWED_ITAPS_SET_CONVERSION( ITAPS_HANDLE_TYPE ) \
  template <> inline \
  ITAPS_HANDLE_TYPE \
  itaps_cast<ITAPS_HANDLE_TYPE,MBEntityHandle>( MBEntityHandle h ) \
  { return itaps_set_cast_<ITAPS_HANDLE_TYPE>(h); } \
  \
  template <> inline \
  MBEntityHandle \
  itaps_cast<MBEntityHandle,ITAPS_HANDLE_TYPE>( ITAPS_HANDLE_TYPE handle ) \
  { return reinterpret_cast<MBEntityHandle>(handle); } \
  \
  template <> inline \
  MBEntityHandle* \
  itaps_cast<MBEntityHandle*,ITAPS_HANDLE_TYPE*>( ITAPS_HANDLE_TYPE* ptr ) \
  { return itaps_cast_ptr_(ptr); } \
  \
  template <> inline \
  const MBEntityHandle* \
  itaps_cast<const MBEntityHandle*,const ITAPS_HANDLE_TYPE*>( const ITAPS_HANDLE_TYPE* ptr ) \
  { return itaps_cast_const_ptr_(ptr); }

DECLARE_ALLOWED_ITAPS_SET_CONVERSION( iMeshP_PartitionHandle )
DECLARE_ALLOWED_ITAPS_SET_CONVERSION( iMeshP_PartHandle )
//DECLARE_ALLOWED_ITAPS_SET_CONVERSION( iBase_EntitySetHandle )
DECLARE_ALLOWED_ITAPS_SET_CONVERSION( iBase_EntityHandle )

#endif

// Need a different function name for MBTag because (currently)
// both MBTag and iBase_EntityHandle are void**.
iBase_TagHandle itaps_tag_cast( MBTag t )
{ 
  assert(sizeof(iBase_TagHandle) >= sizeof(MBTag)); 
  return reinterpret_cast<iBase_TagHandle>(t);
}

/********************* ITAPS arrays **************************/

// Access this method using ALLOCATE_ARRAY macro, rather than callind directly.
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

// For use by ALLOCATE_ARRAY macro
inline int allocate_itaps_array_failed()
{
  strcpy( iMesh_LAST_ERROR.description, 
          "Insufficient allocated array size or insufficient "
          "memory to allocate array." );
  iMesh_LAST_ERROR.error_type = iBase_MEMORY_ALLOCATION_FAILED;
  return iBase_MEMORY_ALLOCATION_FAILED;
}

// If ITAPS array is NULL, allocate it to the requested size, otherwise
// verify that it can hold at least SIZE values.
#define ALLOCATE_ARRAY( NAME, SIZE ) \
  if (!allocate_itaps_array( *NAME, *NAME##_allocated, *NAME##_size, (SIZE) )) \
    RETURN(allocate_itaps_array_failed())

// Handle returning MBRange in ITAPS array (do ALLOCATE_ARRAY and copy).
#define MBRANGE_TO_ITAPS_ARRAY( RANGE, NAME ) do { \
  ALLOCATE_ARRAY( NAME, (RANGE).size() ); \
  std::copy( (RANGE).begin(), (RANGE).end(), itaps_cast<MBEntityHandle*>(*(NAME)) ); \
  } while (false)


/********************* iMeshP API **************************/

#ifdef __cplusplus
extern "C" {
#endif

void iMeshP_createPartitionAll( iMesh_Instance instance,
                        /*in*/  MPI_Comm communicator,
                        /*out*/ iMeshP_PartitionHandle *partition_handle,
                                int *err )
{
  partition_handle = 0;

  MBTag prtn_tag;
  MBErrorCode rval = MBI->tag_create( PARALLEL_PARITIONING_TAG_NAME, 
                                      sizeof(int), 
                                      MB_TAG_SPARSE,
                                      MB_TYPE_INTEGER, 
                                      prtn_tag, 
                                      0, 
                                      true ); CHKERR(rval);
  
  MBEntityHandle handle;
  rval = MBI->create_meshset( MESHSET_SET, handle ); CHKERR(rval);
  MBParallelComm* pcomm = MBParallelComm::get_pcomm( MBI, handle, &communicator );
  if (!pcomm) {
    MBI->delete_entities( &handle, 1 );
    RETURN(MB_FAILURE);
  }
  
  *partition_handle = itaps_cast<iMeshP_PartitionHandle>(handle);
  RETURN (iBase_SUCCESS);
}

void iMeshP_destroyPartitionAll( iMesh_Instance instance,
                                 iMeshP_PartitionHandle partition_handle,
                                 int *err)
{
  MBEntityHandle handle = itaps_cast<MBEntityHandle>(partition_handle);
  MBParallelComm* pcomm = MBParallelComm::get_pcomm( MBI, handle );
  if (pcomm)
    delete pcomm;
  MBErrorCode rval = MBI->delete_entities( &handle, 1 ); CHKERR(rval);
  RETURN (iBase_SUCCESS);
}

void iMeshP_getPartitionComm( iMesh_Instance instance,
                              iMeshP_PartitionHandle partition_handle,
                              MPI_Comm* communicator_out,
                              int* err )
{
  MBEntityHandle handle = itaps_cast<MBEntityHandle>(partition_handle);
  MBParallelComm* pcomm = MBParallelComm::get_pcomm( MBI, handle );
  if (pcomm)
    RETURN (iBase_FAILURE);
  *communicator_out = pcomm->proc_config().proc_rank();
  RETURN (iBase_SUCCESS);
}

void iMeshP_syncPartitionAll( iMesh_Instance /*instance*/,
                              iMeshP_PartitionHandle /*partition_handle*/,
                              int* err )
{
  RETURN (iBase_SUCCESS);
}

void iMeshP_getNumPartitions( iMesh_Instance instance,
                              int* num_partitions_out,
                              int* err )
{
  MBTag prtn_tag;
  MBErrorCode rval = MBI->tag_get_handle( PARALLEL_PARITIONING_TAG_NAME, prtn_tag );
  if (MB_TAG_NOT_FOUND == rval) {
    *num_partitions_out = 0;
    RETURN (iBase_SUCCESS);
  }
  CHKERR(rval);

  rval = MBI->get_number_entities_by_type_and_tag( 0, MBENTITYSET, &prtn_tag, 0, 1,
                                                   *num_partitions_out );
  CHKERR(rval);
  RETURN (iBase_SUCCESS);
}

void iMeshP_getPartitions( iMesh_Instance instance,
                           iMeshP_PartitionHandle **partition_handle,
                           int *partition_handle_allocated, 
                           int *partition_handle_size, 
                           int *err )
{
  MBRange range;
  MBTag prtn_tag;
  MBErrorCode rval = MBI->tag_get_handle( PARALLEL_PARITIONING_TAG_NAME, prtn_tag );
  if (MB_SUCCESS == rval)
    rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, &prtn_tag, 0, 1,
                                              range );
  else if (MB_TAG_NOT_FOUND == rval)
    rval = MB_SUCCESS;
  
  CHKERR(rval);
  MBRANGE_TO_ITAPS_ARRAY( range, partition_handle );
  RETURN (iBase_SUCCESS );
}

void iMeshP_getNumGobalParts( iMesh_Instance instance,
                              const iMeshP_PartitionHandle partition_handle,
                              int *num_global_part, 
                              int *err )
{
  FIXME;
  RETURN (iBase_NOT_SUPPORTED);
}

void iMeshP_getNumLocalParts(iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          int *num_local_part, 
                          int *err)
{
  MBEntityHandle prtn = itaps_cast<MBEntityHandle>(partition_handle);
  MBParallelComm* pcomm = MBParallelComm::get_pcomm( MBI, prtn );
  if (!pcomm)
    RETURN (iBase_FAILURE);
  
  MBRange parts;
  MBErrorCode rval = pcomm->get_partition_sets( prtn, parts );
  CHKERR(rval);
  
  *num_local_part = parts.size();
  RETURN (iBase_SUCCESS);
}

void iMeshP_getLocalParts( iMesh_Instance instance,
                           const iMeshP_PartitionHandle partition_handle,
                           int *num_local_part, 
                           iMeshP_PartHandle **part_handles,
                           int *part_handles_allocated,
                           int *part_handles_size,
                           int *err )
{
  MBEntityHandle prtn = itaps_cast<MBEntityHandle>(partition_handle);
  MBParallelComm* pcomm = MBParallelComm::get_pcomm( MBI, prtn );
  if (!pcomm)
    RETURN (iBase_FAILURE);
  
  MBRange parts;
  MBErrorCode rval = pcomm->get_partition_sets( prtn, parts );
  CHKERR(rval);
  
  MBRANGE_TO_ITAPS_ARRAY( parts, part_handles );
  RETURN (iBase_SUCCESS);
}

void iMeshP_getRankOfPart( iMesh_Instance instance,
                           const iMeshP_PartitionHandle partition_handle,
                           const iMeshP_Part part_id,
                           int *rank,
                           int *err )
{
  int junk1 = 1, junk2 = 1;
  iMeshP_getRankOfPartArr( instance, partition_handle, &part_id,
                           1, &rank, &junk1, &junk2, err );
}

void iMeshP_getRankOfPartArr( iMesh_Instance instance,
                              const iMeshP_PartitionHandle partition_handle,
                              const iMeshP_Part *part_ids,
                              const int part_ids_size,
                              int **rank, 
                              int *rank_allocated, 
                              int *rank_size,
                              int *err )
{
  FIXME; // need to handle handles to "remote parts" ?

  MBEntityHandle prtn = itaps_cast<MBEntityHandle>(partition_handle);
  MBParallelComm* pcomm = MBParallelComm::get_pcomm( MBI, prtn );
  if (!pcomm)
    RETURN (iBase_FAILURE);
  
  ALLOCATE_ARRAY( rank, part_ids_size );
  std::fill( *rank, (*rank) + part_ids_size, pcomm->proc_config().proc_rank() );
  RETURN (iBase_SUCCESS);
}

void iMeshP_getNumOfTypeAll( iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             const iBase_EntitySetHandle entity_set_handle,
                             const int entity_type, 
                             int *num_type, 
                             int *err )
{
  FIXME; // need to prune out entities that are remote

  MBEntityHandle prtn = itaps_cast<MBEntityHandle>(partition_handle);
  MBParallelComm* pcomm = MBParallelComm::get_pcomm( MBI, prtn );
  if (!pcomm)
    RETURN (iBase_FAILURE);
  
  iMesh_getNumOfType( instance, entity_set_handle, entity_type, num_type, err );
  int vals[2] = { *num_type, *err }, sums[2];
  int ierr = MPI_Allreduce( vals, sums, 2, MPI_INT, MPI_SUM, pcomm->proc_config().proc_comm() );
  assert(iBase_SUCCESS == 0);
  if (ierr || sums[1])
    RETURN (iBase_FAILURE);
  
  *num_type = sums[0];
  RETURN (iBase_SUCCESS);
}

void iMeshP_getNumOfTopoAll( iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             const iBase_EntitySetHandle entity_set_handle,
                             const int entity_topology, 
                             int *num_topo, 
                             int *err )
{
  FIXME; // need to prune out entities that are remote

  MBEntityHandle prtn = itaps_cast<MBEntityHandle>(partition_handle);
  MBParallelComm* pcomm = MBParallelComm::get_pcomm( MBI, prtn );
  if (!pcomm)
    RETURN (iBase_FAILURE);
  
  iMesh_getNumOfTopo( instance, entity_set_handle, entity_topology, num_topo, err );
  int vals[2] = { *num_topo, *err }, sums[2];
  int ierr = MPI_Allreduce( vals, sums, 2, MPI_INT, MPI_SUM, pcomm->proc_config().proc_comm() );
  assert(iBase_SUCCESS == 0);
  if (ierr || sums[1])
    RETURN (iBase_FAILURE);
  
  *num_topo = sums[0];
  RETURN (iBase_SUCCESS);
}

void iMeshP_createPart( iMesh_Instance instance,
                        iMeshP_PartitionHandle partition_handle,
                        iMeshP_PartHandle *part_handle,
                        int *err )
{
  MBEntityHandle h, p = itaps_cast<MBEntityHandle>(partition_handle);
  MBErrorCode rval;
  MBParallelComm* pcomm = MBParallelComm::get_pcomm( MBI, p );
  if (!pcomm)
    RETURN (iBase_FAILURE);
  
  const MBTag tag = pcomm->partition_tag();

  MBRange parts;
  rval = MBI->get_entities_by_handle( p, parts ); CHKERR(rval);
  std::vector<int> part_ids( parts.size() );
  rval = MBI->tag_get_data( tag, parts, &part_ids[0] ); CHKERR(rval);
  int max_loc_id = -1, loc_id;
  for (size_t i = 0; i < part_ids.size(); ++i) {
    loc_id = LOCAL_ID( part_ids[i] );
    if (loc_id > max_loc_id)
      max_loc_id = loc_id;
  }
  int new_id = PART_ID( pcomm->proc_config().proc_rank(), max_loc_id + 1 );
  

  rval = MBI->create_meshset( MESHSET_SET, h ); CHKERR(rval);
  if (MB_SUCCESS == rval)
    rval = MBI->add_entities( p, &h, 1 );
  if (MB_SUCCESS == rval)
    rval = MBI->tag_set_data( tag, &h, 1, &new_id );
  if (MB_SUCCESS != rval) {
    MBI->delete_entities( &h, 1 );
    CHKERR(rval);
  }
  
  *part_handle = itaps_cast<iMeshP_PartHandle>(h);
  RETURN (iBase_SUCCESS);
}

void iMeshP_destroyPart( iMesh_Instance instance,
                         iMeshP_PartitionHandle partition_handle,
                         iMeshP_PartHandle part_handle,
                         int *err )
{
  MBErrorCode rval;
  MBEntityHandle h = itaps_cast<MBEntityHandle>(part_handle), 
                 p = itaps_cast<MBEntityHandle>(partition_handle);

  
  rval = MBI->remove_entities( p, &h, 1 ); CHKERR(rval);
  rval = MBI->delete_entities( &h, 1 ); CHKERR(rval);
}

void iMeshP_getPartIdFromPartHandle( iMesh_Instance instance,
                                     const iMeshP_PartitionHandle partition_handle,
                                     const iMeshP_PartHandle part_handle,
                                     iMeshP_Part *part_id,
                                     int *err )
{
  int junk1 = 1, junk2 = 1;
  iMeshP_getPartIdsFromPartHandlesArr( instance, 
                                       partition_handle,
                                       &part_handle,
                                       1,
                                       &part_id,
                                       &junk1,
                                       &junk2,
                                       err );
}

void iMeshP_getPartIdsFromPartHandlesArr( iMesh_Instance instance,
                                          const iMeshP_PartitionHandle partition_handle,
                                          const iMeshP_PartHandle *part_handles,
                                          const int part_handles_size,
                                          iMeshP_Part **part_ids,
                                          int *part_ids_allocated,
                                          int *part_ids_size,
                                          int *err )
{
  MBEntityHandle p = itaps_cast<MBEntityHandle>(partition_handle);
  MBParallelComm* pcomm = MBParallelComm::get_pcomm( MBI, p );
  if (!pcomm)
    RETURN (iBase_FAILURE);
  
  const MBTag tag = pcomm->partition_tag();
  
  ALLOCATE_ARRAY( part_ids, part_handles_size );
  
  int* array;
  std::vector<int> tmp_storage;
  if (sizeof(iMeshP_Part) == sizeof(int)) {
    array = reinterpret_cast<int*>(*part_ids);
  }
  else {
    tmp_storage.resize(part_handles_size);
    array = &tmp_storage[0];
  }
  
  MBErrorCode rval = MBI->tag_get_data( tag, 
                                        itaps_cast<const MBEntityHandle*>(part_handles),
                                        part_handles_size,
                                        array ); CHKERR(rval);
                                        
  if (sizeof(iMeshP_Part) != sizeof(int))
    std::copy( array, array + part_handles_size, *part_ids );
  
  RETURN (iBase_SUCCESS);
}

void iMeshP_getNumPartNbors( iMesh_Instance instance,
                             iMeshP_PartitionHandle partition_handle,
                             iMeshP_PartHandle part_handle,
                             int entity_type,
                             int *num_part_nbors,
                             int *err )
{
  int junk1 = 1, junk2 = 1;
  iMeshP_getNumPartNborsArr( instance, partition_handle,
                             &part_handle, 1, entity_type,
                             &num_part_nbors, &junk1, &junk2,
                             err );
}

void iMeshP_getNumPartNborsArr( iMesh_Instance instance,
                                const iMeshP_PartitionHandle partition_handle,
                                const iMeshP_PartHandle *part_handles,
                                const int part_handles_size,
                                int entity_type,
                                int **num_part_nbors,
                                int *num_part_nbors_allocated,
                                int *num_part_nbors_size,
                                int *err )
{
  FIXME;
  RETURN( iBase_NOT_SUPPORTED );
}


void iMeshP_getPartNbors( iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          const iMeshP_PartHandle part_handle,
                          int entity_type,
                          int *num_part_nbors,
                          iMeshP_Part **nbor_part_ids,
                          int *nbor_part_ids_allocated,
                          int *nbor_part_ids_size,
                          int *err )
{
  int junk1 = 1, junk2 = 1;
  iMeshP_getPartNborsArr( instance, partition_handle, 
                          &part_handle, 1, entity_type, 
                          &num_part_nbors, &junk1, &junk2,
                          nbor_part_ids, nbor_part_ids_allocated, 
                          nbor_part_ids_size, err );
}

void iMeshP_getPartNborsArr( iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             const iMeshP_PartHandle *part_handles,
                             const int part_handles_size,
                             int entity_type,
                             int **num_part_nbors,
                             int *num_part_nbors_allocated,
                             int *num_part_nbors_size,
                             iMeshP_Part **nbor_part_ids,
                             int *nbor_part_ids_allocated,
                             int *nbor_part_ids_size,
                             int *err ) 
{
  FIXME;
  RETURN(iBase_NOT_SUPPORTED);
}

void iMeshP_getNumPartBdryEnts( iMesh_Instance instance,
                                const iMeshP_PartitionHandle partition_handle,
                                const iMeshP_PartHandle part_handle, 
                                const int entity_type, 
                                const int entity_topology, 
                                const iMeshP_Part target_part_id, 
                                int *num_entities, 
                                int *err )
{
  FIXME;
  RETURN(iBase_NOT_SUPPORTED);
}

void iMeshP_getPartBdryEnts( iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             const iMeshP_PartHandle part_handle,
                             const int entity_type,
                             const int entity_topology,
                             const iMeshP_Part target_part_id,
                             iBase_EntityHandle **entity_handles,
                             int *entity_handles_allocated,
                             int *entity_handles_size,
                             int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_initPartBdryEntIter( iMesh_Instance instance,
                                 const iMeshP_PartitionHandle partition_handle,
                                 const iMeshP_PartHandle part_handle,
                                 const int entity_type,
                                 const int entity_topology,
                                 const iMeshP_Part nbor_part_id,
                                 iMesh_EntityIterator* entity_iterator,
                                 int* err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_initPartBdryEntArrIter( iMesh_Instance instance,
                                    const iMeshP_PartitionHandle partition_handle,
                                    const iMeshP_PartHandle part_handle,
                                    const int entity_type,
                                    const int entity_topology,
                                    const int array_size,
                                    const iMeshP_Part nbor_part_id,
                                    iMesh_EntityIterator* entity_iterator,
                                    int* err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }



static void set_intersection_query( iMesh_Instance instance,
                                    iMeshP_PartHandle set1,
                                    iBase_EntitySetHandle set2,
                                    int type,
                                    int topo,
                                    MBRange& result,
                                    int* err )
{
  MBErrorCode rval;
  MBRange r1, r2;
  MBEntityHandle h1 = itaps_cast<MBEntityHandle>(set1);
  MBEntityHandle h2 = itaps_cast<MBEntityHandle>(set2);
  
  if (topo != iMesh_ALL_TOPOLOGIES) {
    if ((unsigned)topo > sizeof(mb_topology_table)/sizeof(mb_topology_table[0]))
      RETURN (iBase_INVALID_ENTITY_TYPE);
    MBEntityType t = mb_topology_table[topo];
    if (t == MBMAXTYPE) {
      result.clear();
      RETURN (iBase_SUCCESS); // if unsupported topo, then there are zero of them
    }
    rval = MBI->get_entities_by_type( h1, t, r1 ); CHKERR(rval);
    rval = MBI->get_entities_by_type( h2, t, r2 ); CHKERR(rval);
  }
  else if (type != iBase_ALL_TYPES) {
    if (type < 0 || type > 3)
      RETURN (iBase_INVALID_ENTITY_TYPE);
    rval = MBI->get_entities_by_dimension( h1, type, r1 ); CHKERR(rval);
    rval = MBI->get_entities_by_dimension( h2, type, r2 ); CHKERR(rval);
  }
  else {
    rval = MBI->get_entities_by_handle( h1, r1 ); CHKERR(rval);
    rval = MBI->get_entities_by_handle( h2, r2 ); CHKERR(rval);
  }
  
  result = r1.intersect( r2 );
  RETURN (iBase_SUCCESS);
}
  

void iMeshP_getNumOfType( iMesh_Instance instance,
                          const iMeshP_PartitionHandle ,
                          const iMeshP_PartHandle part_handle,
                          const iBase_EntitySetHandle entity_set_handle,
                          const int entity_type,
                          int *num_type,
                          int *err )
{
  MBRange r;
  set_intersection_query( instance, part_handle, entity_set_handle, 
                          entity_type, iMesh_ALL_TOPOLOGIES, r, err );
  *num_type = r.size();
}

void iMeshP_getNumOfTopo( iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          const iMeshP_PartHandle part_handle,
                          const iBase_EntitySetHandle entity_set_handle,
                          const int entity_topology,
                          int *num_topo,
                          int *err )
{
  MBRange r;
  set_intersection_query( instance, part_handle, entity_set_handle, 
                          iBase_ALL_TYPES, entity_topology, r, err );
  *num_topo = r.size();
}

void iMeshP_getAllVtxCoords( iMesh_Instance instance,
                             const iMeshP_PartitionHandle ,
                             const iMeshP_PartHandle part_handle,
                             const iBase_EntitySetHandle entity_set_handle,
                             double** coordinates,
                             int* coordinates_allocated,
                             int* coordinates_size,
                             int** in_entity_set,
                             int* in_entity_set_allocated,
                             int* in_entity_set_size,
                             int* storage_order,
                             int *err )
{
  int dim;
  MBErrorCode rval = MBI->get_dimension( dim ); CHKERR(rval);

  MBRange r;
  set_intersection_query( instance, part_handle, entity_set_handle, 
                          iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES, r, err );
  if (*err != iBase_SUCCESS)
    return;
    
  MBRange v;
  rval = MBI->get_adjacencies( r, 0, false, v, MBInterface::UNION );
  CHKERR(rval);
    
  const size_t n = v.size();
  ALLOCATE_ARRAY( coordinates, dim*n );
  ALLOCATE_ARRAY( in_entity_set, n );
  std::fill( *in_entity_set, (*in_entity_set)+n, 1 );

  switch (*storage_order) {
    case iBase_UNDETERMINED:
      *storage_order = iBase_INTERLEAVED;
      // NO BREAK -- fall through
    case iBase_INTERLEAVED:
      rval = MBI->get_coords( v, *coordinates );
      break;
    case iBase_BLOCKED:
      rval = MBI->get_coords( v,
                               *coordinates,
                              (*coordinates)+n,
                              (*coordinates)+2*n );
      break;
    default:
      RETURN (iBase_INVALID_ARGUMENT);
  }
  
  CHKERR(rval);
  RETURN (iBase_SUCCESS);
}

void iMeshP_getVtxCoordIndex( iMesh_Instance instance,
                              const iMeshP_PartitionHandle ,
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
                              int *err )
{
  MBRange::iterator i;
  MBRange r;
  set_intersection_query( instance, part_handle, entity_set_handle, 
                          iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES, r, err );
  if (*err != iBase_SUCCESS)
    return;
    
  MBRange v;
  MBErrorCode rval = MBI->get_adjacencies( r, 0, false, v, MBInterface::UNION );
  CHKERR(rval);

  r.clear();
  set_intersection_query( instance, part_handle, entity_set_handle, 
                          requested_entity_type, requested_entity_topology, 
                          r, err );
  if (*err != iBase_SUCCESS)
    return;
    
  std::vector<MBEntityHandle> sorted(v.size());
  std::copy( v.begin(), v.end(), sorted.begin() );
  
  if (!r.all_of_dimension(entity_adjacency_type)) {
    MBRange r2;
    rval =  MBI->get_adjacencies( r, entity_adjacency_type, false, r2, MBInterface::UNION );
    CHKERR(rval);
    r.swap(r2);
  }
  
    // count size of connectivity data
  size_t cn = 0;
  int num;
  const MBEntityHandle* conn;
  std::vector<MBEntityHandle>::iterator s;
  std::vector<MBEntityHandle> store;
  for (i = r.begin(); i != r.end(); ++i) {
    rval = MBI->get_connectivity( *i, conn, num, true, &store );
    CHKERR(rval);
    cn += num;
  }

  const size_t n = r.size();
  ALLOCATE_ARRAY( offset, n );
  ALLOCATE_ARRAY( index, cn );
  ALLOCATE_ARRAY( entity_topologies, n );
    
    // populate output arrays
  int* offset_iterator = *offset;
  int* index_iterator = *index;
  int* topo_iterator = *entity_topologies;
  for (i = r.begin(); i != r.end(); ++i) {
    *offset_iterator = index_iterator - *index;
    ++offset_iterator;
    
    rval = MBI->get_connectivity( *i, conn, num, true, &store );
    CHKERR(rval);
    for (int j = 0; j < num; ++j) {
      s = std::lower_bound( sorted.begin(), sorted.end(), conn[j] );
      if (s == sorted.end() || *s != conn[j]) 
        RETURN(iBase_FAILURE);
      *index_iterator = s - sorted.begin();
      ++index_iterator;
    }
    
    *topo_iterator = tstt_topology_table[ TYPE_FROM_HANDLE(*i) ];
    ++topo_iterator;
  }
}

void iMeshP_getEntities( iMesh_Instance instance,
                         const iMeshP_PartitionHandle ,
                         const iMeshP_PartHandle part_handle,
                         const iBase_EntitySetHandle entity_set_handle,
                         const int entity_type,
                         const int entity_topology,
                         iBase_EntityHandle** entity_handles,
                         int* entity_handles_allocated,
                         int* entity_handles_size,
                         int *err )
{
  MBRange r;
  set_intersection_query( instance, part_handle, entity_set_handle, 
                          entity_type, entity_topology, r, err );
  if (iBase_SUCCESS != *err)
    return;
  
  MBRANGE_TO_ITAPS_ARRAY( r, entity_handles );
  RETURN(iBase_SUCCESS);
}

void iMeshP_getAdjEntities( iMesh_Instance instance,
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
                            int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_initEntIter( iMesh_Instance instance,
                         const iMeshP_PartitionHandle partition_handle,
                         const iMeshP_PartHandle part_handle,
                         const iBase_EntitySetHandle entity_set_handle,
                         const int requested_entity_type,
                         const int requested_entity_topology,
                         iMesh_EntityIterator* entity_iterator,
                         int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_initEntArrIter( iMesh_Instance instance,
                            const iMeshP_PartitionHandle partition_handle,
                            const iMeshP_PartHandle part_handle,
                            const iBase_EntitySetHandle entity_set_handle,
                            const int requested_entity_type,
                            const int requested_entity_topology,
                            const int requested_array_size,
                            iMesh_EntityArrIterator* entArr_iterator,
                            int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_getEntOwnerPart( iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             const iBase_EntityHandle entity_handle,
                             iMeshP_Part *part_id,
                             int* err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_getEntOwnerPartArr( iMesh_Instance instance,
                                const iMeshP_PartitionHandle partition_handle,
                                const iBase_EntityHandle *entity_handles,
                                const int entity_handles_size,
                                iMeshP_Part **part_ids,
                                int *part_ids_allocated,
                                int *part_ids_size,
                                int* err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }
  
void iMeshP_isEntOwner( iMesh_Instance instance,
                        const iMeshP_PartitionHandle partition_handle,
                        const iMeshP_PartHandle part_handle,
                        const iBase_EntityHandle entity_handle,
                        int* is_owner,
                        int *err )
{
  int junk1 = 1, junk2 = 1;
  iMeshP_isEntOwnerArr( instance, partition_handle,
                        part_handle, &entity_handle, 1,
                        &is_owner, &junk1, &junk2,
                        err );
}

void iMeshP_isEntOwnerArr( iMesh_Instance instance,
                           const iMeshP_PartitionHandle partition_handle,
                           const iMeshP_PartHandle part_handle,
                           const iBase_EntityHandle *entity_handles,
                           const int entity_handles_size,
                           int** is_owner,
                           int* is_owner_allocated,
                           int* is_owner_size,
                           int *err )
{
  MBEntityHandle p = itaps_cast<MBEntityHandle>(partition_handle);
  MBParallelComm* pcomm = MBParallelComm::get_pcomm( MBI, p );
  if (!pcomm)
    RETURN (iBase_FAILURE);
  
  const MBEntityHandle* handles = itaps_cast<const MBEntityHandle*>(entity_handles);
  std::vector<unsigned char> pstatus(entity_handles_size);
  MBErrorCode result = MBI->tag_get_data( pcomm->pstatus_tag(), 
                                          handles, 
                                          entity_handles_size,
                                          &pstatus[0] ); CHKERR(result);
  
  MBEntityHandle part = itaps_cast<MBEntityHandle>(part_handle);
  MBRange part_ents;
  result = MBI->get_entities_by_handle( part, part_ents );
  CHKERR(result);

  ALLOCATE_ARRAY( is_owner, entity_handles_size );
  for (int i = 0; i < entity_handles_size; i++) {
    if (pstatus[i] & PSTATUS_NOT_OWNED) 
      (*is_owner)[i] = 0;
    else if (part_ents.find(handles[i]) == part_ents.end())
      (*is_owner)[i] = 0;
    else 
      (*is_owner)[i] = 1;
  }

  RETURN(iBase_SUCCESS);
}

void iMeshP_getEntStatus(iMesh_Instance instance,
                         /*in*/ const iMeshP_PartitionHandle partition_handle,
                         /*in*/ const iMeshP_PartHandle part_handle, 
                         /*in*/ const iBase_EntityHandle entity_handle, 
                         /*out*/ int* par_status, // Values=INTERNAL,BOUNDARY,GHOST
                         int *err) 
{
  int junk1 = 1, junk2 = 1;
  iMeshP_getEntStatusArr( instance, partition_handle,
                          part_handle, &entity_handle, 1,
                          &par_status, &junk1, &junk2,
                          err );
}

void iMeshP_getEntStatusArr(iMesh_Instance instance,
                            /*in*/ const iMeshP_PartitionHandle partition_handle,
                            /*in*/ const iMeshP_PartHandle part_handle, 
                            /*in*/ const iBase_EntityHandle *entity_handles, 
                            /*in*/ const int entity_handles_size, 
                            /*inout*/ int** par_status, // Values=INTERNAL,BOUNDARY,GHOST
                            /*inout*/ int* par_status_allocated, 
                            /*inout*/ int* par_status_size, 
                            int *err) 
{
  MBEntityHandle p = itaps_cast<MBEntityHandle>(partition_handle);
  MBParallelComm* pcomm = MBParallelComm::get_pcomm( MBI, p );
  if (!pcomm)
    RETURN (iBase_FAILURE);

  std::vector<unsigned char> pstatus(entity_handles_size);
  MBErrorCode result = MBI->tag_get_data(pcomm->pstatus_tag(), 
                                         itaps_cast<const MBEntityHandle*>(entity_handles), 
                                         entity_handles_size,
                                         &pstatus[0]); CHKERR(result);

  ALLOCATE_ARRAY( par_status, entity_handles_size );
  for (int i = 0; i < entity_handles_size; i++) {
    if (!pstatus[i]) 
      (*par_status)[i] = iMeshP_INTERNAL;
    else if (pstatus[i] & PSTATUS_GHOST) 
      (*par_status)[i] = iMeshP_GHOST;
    else if (pstatus[i] & PSTATUS_INTERFACE)
      (*par_status)[i] = iMeshP_BOUNDARY;
  }

  RETURN(iBase_SUCCESS);
}

void iMeshP_getNumCopies( iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          const iBase_EntityHandle entity_handle,
                          int *num_copies_ent,
                          int *err )
{
  MBEntityHandle p = itaps_cast<MBEntityHandle>(partition_handle);
  MBParallelComm* pcomm = MBParallelComm::get_pcomm( MBI, p );
  if (!pcomm)
    RETURN (iBase_FAILURE);
  
  int shared_proc;
  MBErrorCode result = MBI->tag_get_data(pcomm->sharedp_tag(), 
                                         itaps_cast<const MBEntityHandle*>(&entity_handle), 1,
                                         &shared_proc); CHKERR(result);

  if (-1 != shared_proc) {
    *num_copies_ent = 1;
    RETURN(iBase_SUCCESS);
  }

  const void* data_ptr = 0;
  int data_size = 0;
  result = MBI->tag_get_data( pcomm->sharedps_tag(), 
                              itaps_cast<const MBEntityHandle*>(&entity_handle), 
                              1,
                              &data_ptr,
                              &data_size ); CHKERR(result);

  const int* shared_procs = reinterpret_cast<const int*>(data_ptr);
  int tag_size = data_size / sizeof(int);
  *num_copies_ent = std::find(shared_procs, shared_procs + tag_size, -1) - shared_procs;

  RETURN(iBase_SUCCESS);
}

void iMeshP_getCopyParts( iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          const iBase_EntityHandle entity_handle,
                          iMeshP_Part **part_ids,
                          int *part_ids_allocated,
                          int *part_ids_size,
                          int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

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
                                           itaps_cast<const MBEntityHandle*>(&entity_handle), 1,
                                           &shared_handle);
    if (MB_SUCCESS != result) {
      RETURN(iBase_ERROR_MAP[result]);
    }
  
    std::vector<iMeshP_PartHandle> part_handles;
    if (0 != shared_handle) {
      part_handles.resize( 1 );
      ALLOCATE_ARRAY(copies_entity_handles, 1);
      part_handles[0] = 0;
      (*copies_entity_handles)[0] = itaps_cast<iBase_EntityHandle>(shared_handle);
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
                                 itaps_cast<const MBEntityHandle*>(&entity_handle), 1,
                                 &shared_handles[0]);
      if (MB_SUCCESS != result) {
        RETURN(iBase_ERROR_MAP[result]);
      }

      int index = std::find(shared_handles.begin(), shared_handles.end(), -1) - shared_handles.begin();

      part_handles.resize(index+1, 0);
      ALLOCATE_ARRAY( copies_entity_handles, index+1 );
      std::copy(&shared_handles[0], &shared_handles[index], 
                itaps_cast<MBEntityHandle*>(*copies_entity_handles));
    }
    
    ALLOCATE_ARRAY(part_ids, part_handles.size());
    iMeshP_getPartIdsFromPartHandlesArr( instance, partition_handle, 
                                         &part_handles[0], part_handles.size(),
                                         part_ids, part_ids_allocated,
                                         part_ids_size, err );
    *part_ids_size = part_handles.size();
    RETURN(*err);
  }

void iMeshP_getCopyOnPart( iMesh_Instance instance,
                           const iMeshP_PartitionHandle partition_handle,
                           const iBase_EntityHandle entity_handle,
                           const iMeshP_Part* part_id,
                           iBase_EntityHandle* copy_entity_handle,
                           int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }


void iMeshP_getOwnerCopy( iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          const iBase_EntityHandle entity_handle,
                          iMeshP_Part *owner_part_id,
                          iBase_EntityHandle *owner_entity_handle,
                          int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_Wait( iMesh_Instance instance,
                  const iMeshP_PartitionHandle partition_handle,
                  iMeshP_RequestHandle req,
                  iMeshP_Status *stat,
                  int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_WaitAny( iMesh_Instance instance,
                     const iMeshP_PartitionHandle partition_handle,
                     iMeshP_RequestHandle *req,
                     int req_size,
                     int *index,
                     iMeshP_Status *stat,
                     int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_WaitAll( iMesh_Instance instance,
                     const iMeshP_PartitionHandle partition_handle,
                     iMeshP_RequestHandle *req,
                     int req_size,
                     iMeshP_Status *stat,
                     int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_WaitEnt( iMesh_Instance instance,
                     const iMeshP_PartitionHandle partition_handle,
                     iMeshP_RequestHandle req,
                     iBase_EntityHandle **out_entities,
                     int *out_entities_alloc,
                     int *out_entities_size,
                     int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_Test( iMesh_Instance instance,
                  const iMeshP_PartitionHandle partition_handle,
                  iMeshP_RequestHandle req,
                  int *flag,
                  iMeshP_Status *stat,
                  int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_pollForRequests( iMesh_Instance instance,
                             iMeshP_PartitionHandle partition_handle,
                             iMeshP_RequestHandle **requests_completed,
                             int *requests_allocated,
                             int *requests_size,
                             int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_exchEntArrToPartsAll( iMesh_Instance instance,
                                  const iMeshP_PartitionHandle partition_handle,
                                  const iBase_EntityHandle *entity_handles,
                                  const int entity_handles_size,
                                  const iMeshP_Part *target_part_ids,
                                  int command_code,
                                  int update_ghost,
                                  iMeshP_RequestHandle *request,
                                  int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_migrateEntity( iMesh_Instance instance,
                           const iMeshP_PartitionHandle partition_handle,
                           const iMeshP_PartHandle part_handle,
                           const iBase_EntityHandle local_entity_handle,
                           iMeshP_RequestHandle *request,
                           int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_updateVtxCoords( iMesh_Instance instance,
                             const iMeshP_PartitionHandle partition_handle,
                             const iBase_EntityHandle local_vertex_handle,
                             int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_replaceOnPartBdry( iMesh_Instance instance,
                               const iMeshP_PartitionHandle partition_handle,
                               const iBase_EntityHandle *old_entities,
                               const int old_entities_size,
                               const iBase_EntityHandle *new_entities,
                               const int new_entities_size,
                               const int *offset,
                               const int offset_size,
                               int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_addGhostOf( iMesh_Instance instance,
                        const iMeshP_PartitionHandle partition_handle,
                        const iMeshP_Part target_part_id,
                        const iBase_EntityHandle entity_to_copy,
                        iMeshP_RequestHandle *request,
                        int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_rmvGhostOf( iMesh_Instance instance,
                        const iMeshP_PartitionHandle partition_handle,
                        const iMeshP_Part target_part_id,
                        const iBase_EntityHandle copy_to_purge,
                        int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_syncMeshAll( iMesh_Instance instance,
                         const iMeshP_PartitionHandle partition_handle,
                         int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }
		            
void iMeshP_pushTags( iMesh_Instance instance,
                      const iMeshP_PartitionHandle partition_handle,
                      iBase_TagHandle tag,
                      int entity_type,
                      int entity_topo,
                      int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_pushTagsEnt( iMesh_Instance instance,
                         const iMeshP_PartitionHandle partition_handle,
                         iBase_TagHandle tag,
                         const iBase_EntityHandle *entities,
                         int entities_size,
                         int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_iPushTags( iMesh_Instance instance,
                       const iMeshP_PartitionHandle partition_handle,
                       iBase_TagHandle tag,
                       int entity_type,
                       int entity_topo,
                       iMeshP_RequestHandle *req,
                       int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_iPushTagsEnt( iMesh_Instance instance,
                          const iMeshP_PartitionHandle partition_handle,
                          iBase_TagHandle tag,
                          const iBase_EntityHandle *entities,
                          int entities_size,
                          iMeshP_RequestHandle *req,
                          int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_createGhostEnts( iMesh_Instance instance,
                             iMeshP_PartitionHandle partition_handle,
                             int ghost_dim,
                             int bridge_dim,
                             int num_layers,
                             int include_copies,
                             int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_deleteGhostEnts( iMesh_Instance instance,
                             iMeshP_PartitionHandle partition_handle,
                             int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_ghostEntInfo( iMesh_Instance instance,
                          iMeshP_PartitionHandle partition_handle,
                          int *num_ghost_rules,
                          int *ghost_rules_allocated,
                          int *ghost_rules_size,
                          int **ghost_dim,
                          int **bridge_dim,
                          int **num_layers,
                          int *err )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_load( iMesh_Instance instance,
                  const iMeshP_PartitionHandle partition,
                  const iBase_EntitySetHandle entity_set_handle,
                  const char *name,
                  const char *options,
                  int *err,
                  int name_len,
                  int options_len )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }

void iMeshP_save( iMesh_Instance instance,
                  const iMeshP_PartitionHandle partition,
                  const iBase_EntitySetHandle entity_set_handle,
                  const char *name,
                  const char *options,
                  int *err,
                  const int name_len,
                  int options_len )
{ FIXME; RETURN(iBase_NOT_SUPPORTED); }




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
    MBEntityHandle p = itaps_cast<MBEntityHandle>(partition_handle);
    MBParallelComm *pc = MBParallelComm::get_pcomm(MBI, p);
    if (!pc) RETURN(iBase_ERROR_MAP[MB_FAILURE]);

    MBRange part_sets;
  
    ALLOCATE_ARRAY( part_handles, pc->partition_sets().size() );
    MBRange::iterator rit;
    int i;
    for (i = 0, rit = pc->partition_sets().begin(); 
         rit != pc->partition_sets().end(); rit++, i++)
      (*part_handles)[i] = itaps_cast<iMeshP_PartHandle>(*rit);
  
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
    MBEntityHandle p = itaps_cast<MBEntityHandle>(partition_handle);
    MBParallelComm *pc = MBParallelComm::get_pcomm(MBI, p);
    if (!pc) RETURN(iBase_ERROR_MAP[MB_FAILURE]);

    if (rank[0] != (int)pc->proc_config().proc_rank() || rank_size > 1) {
      RETURN(iBase_ERROR_MAP[MB_NOT_IMPLEMENTED]);
    }
  
    iMeshP_getPartsOnRank(instance, partition_handle, rank[0],
                          part_handles, part_handles_allocated, part_handles_size,
                          err);
  }

#ifdef __cplusplus
} // extern "C"
#endif
