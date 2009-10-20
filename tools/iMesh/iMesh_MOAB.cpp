#include "iMesh_extensions.h"
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBCN.hpp"
#include "MeshTopoUtil.hpp"
#include "FileOptions.hpp"
#include "iMesh_MOAB.hpp"
#define IS_BUILDING_MB
#include "MBInternals.hpp"
#undef IS_BUILDING_MB

#ifdef USE_MPI    
#include "mpi.h"
#endif

#include <iostream>
#include <cassert>
#include <cstring>
#include <stdarg.h>
#include <stdio.h>
#define MIN(a,b) (a < b ? a : b)

MBiMesh::MBiMesh()
    : haveDeletedEntities(false)
{
  memset(AdjTable, 0, 16*sizeof(int));
  for (int i = 0; i < 4; i++) AdjTable[4*i] = AdjTable[i] = 1;
  AdjTable[15] = 1;
}

MBiMesh::~MBiMesh() {}

MBErrorCode MBiMesh::delete_mesh() {
  haveDeletedEntities = true;
  return MBCore::delete_mesh();
}

MBErrorCode MBiMesh::delete_entities( const MBEntityHandle* a, const int n )
{
  if (n > 0)
    haveDeletedEntities = true;
  return MBCore::delete_entities( a, n );
}

MBErrorCode MBiMesh::delete_entities( const MBRange& r )
{
  if (!r.empty())
    haveDeletedEntities = true;
  return MBCore::delete_entities( r );
}

MBErrorCode create_int_ents(MBInterface *instance,
                            MBRange &from_ents,
                            MBEntityHandle in_set = 0);

#define CHECK_SIZE(array, allocated, size, type, retval)  \
  if (0 != allocated && NULL != array && allocated < (size)) {\
    iMesh_processError(iBase_MEMORY_ALLOCATION_FAILED, \
          "Allocated array not large enough to hold returned contents.");\
    RETURN(iBase_MEMORY_ALLOCATION_FAILED);\
  }\
  if ((size) && ((allocated) == 0 || NULL == (array))) {\
    array = (type*)malloc((size)*sizeof(type));\
    allocated=(size);\
    if (NULL == array) {iMesh_processError(iBase_MEMORY_ALLOCATION_FAILED, \
          "Couldn't allocate array.");RETURN(iBase_MEMORY_ALLOCATION_FAILED); }\
  }
// TAG_CHECK_SIZE is like CHECK_SIZE except it checks for and makes the allocated memory
// size a multiple of sizeof(void*), and the pointer is assumed to be type char*
#define TAG_CHECK_SIZE(array, allocated, size)  \
  if (0 != allocated && NULL != array && allocated < (size)) {\
    iMesh_processError(iBase_MEMORY_ALLOCATION_FAILED, \
          "Allocated array not large enough to hold returned contents.");\
    RETURN(iBase_MEMORY_ALLOCATION_FAILED);\
  }\
  if ((size) && (NULL == (array) || (allocated) == 0)) {\
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
#define CAST_TO_VOID(ptr) reinterpret_cast<void*>(ptr)

// map from MB's entity type to TSTT's entity topology
const iMesh_EntityTopology tstt_topology_table[] =
{
  iMesh_POINT,          // MBVERTEX
  iMesh_LINE_SEGMENT,   // MBEDGE
  iMesh_TRIANGLE,       // MBTRI
  iMesh_QUADRILATERAL,  // MBQUAD
  iMesh_POLYGON,        // MBPOLYGON
  iMesh_TETRAHEDRON,    // MBTET
  iMesh_PYRAMID,        // MBPYRAMID
  iMesh_PRISM,          // MBPRISM
  iMesh_ALL_TOPOLOGIES, // MBKNIFE
  iMesh_HEXAHEDRON,     // MBHEX
  iMesh_POLYHEDRON,     // MBPOLYHEDRON
  iMesh_ALL_TOPOLOGIES, // MBENTITYSET
  iMesh_ALL_TOPOLOGIES, // MBMAXTYPE
};

// map from MB's entity type to TSTT's entity type
const iBase_EntityType tstt_type_table[] =
{
  iBase_VERTEX,         // MBVERTEX
  iBase_EDGE,           // MBEDGE
  iBase_FACE,           // MBTRI
  iBase_FACE,           // MBQUAD
  iBase_FACE,           // MBPOLYGON
  iBase_REGION,         // MBTET
  iBase_REGION,         // MBPYRAMID
  iBase_REGION,         // MBPRISM
  iBase_REGION,         // MBKNIFE
  iBase_REGION,         // MBHEX
  iBase_REGION,         // MBPOLYHEDRON
  iBase_ALL_TYPES, // MBENTITYSET
  iBase_ALL_TYPES  // MBMAXTYPE
};

// map to MB's entity type from TSTT's entity topology
const MBEntityType mb_topology_table[] =
{
  MBVERTEX,
  MBEDGE,
  MBPOLYGON,
  MBTRI,
  MBQUAD,
  MBPOLYHEDRON,
  MBTET,
  MBHEX,
  MBPRISM,
  MBPYRAMID,
  MBMAXTYPE,
  MBMAXTYPE
};

// map from TSTT's tag types to MOAB's
const MBDataType mb_data_type_table[] = 
{
  MB_TYPE_INTEGER,
  MB_TYPE_DOUBLE ,
  MB_TYPE_HANDLE ,
  MB_TYPE_OPAQUE
};

// map from MOAB's tag types to tstt's
const iBase_TagValueType tstt_data_type_table[] = 
{
  iBase_BYTES,
  iBase_INTEGER,
  iBase_DOUBLE,
  iBase_BYTES,
  iBase_ENTITY_HANDLE
};

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
  iBase_FAILURE, // MB_VARIABLE_DATA_LENGTH,
  iBase_FAILURE, // MB_INVALID_SIZE,
  iBase_NOT_SUPPORTED, // MB_UNSUPPORTED_OPERATION,
  iBase_FAILURE // MB_FAILURE};
};

struct RangeIterator
{
  MBRange iteratorRange;
  MBRange::iterator currentPos;
  int requestedSize;
};

iMesh_EntityIterator create_itaps_iterator( MBRange& range, int array_size )
{
  RangeIterator* iter = new RangeIterator;
  if (iter) {
    iter->iteratorRange.swap(range);
    iter->currentPos = iter->iteratorRange.begin();
    iter->requestedSize = array_size;
  }
  return reinterpret_cast<iMesh_EntityIterator>(iter);
}

static inline void
iMesh_processError( iBase_ErrorType code, const char* desc ) 
{
  std::strncpy( iMesh_LAST_ERROR.description, desc,
                sizeof(iMesh_LAST_ERROR.description) );
  iMesh_LAST_ERROR.error_type = code;
}

#ifdef __cplusplus
extern "C" {
#endif

  void eatwhitespace(std::string &this_string);

  MBErrorCode iMesh_tag_set_vertices(iMesh_Instance instance,
                                     MBEntityHandle in_set, 
                                     const int req_dimension, 
                                     const MBEntityType req_type,
                                     MBTag &tag, MBRange &req_entities, 
                                     int &num_verts, int *err);

  iBase_Error iMesh_LAST_ERROR;

  void iMesh_getErrorType(iMesh_Instance instance, 
                          int *error_type, int *err) 
  {
    *error_type = iMesh_LAST_ERROR.error_type;
    RETURN(iBase_SUCCESS);
  }
  
  void iMesh_getDescription(iMesh_Instance instance, 
                            char *descr, int *err, int descr_len)
  {
    unsigned int len = MIN(strlen(iMesh_LAST_ERROR.description), ((unsigned int) descr_len));
    strncpy(descr, iMesh_LAST_ERROR.description, len);
    descr[len] = '\0';
    RETURN(iBase_SUCCESS);
  }

  void iMesh_setError(iMesh_Instance instance,
                      int err_type, char *descr, int *err, int descr_len) 
  {
    iMesh_LAST_ERROR.error_type = static_cast<iBase_ErrorType>(err_type);
    strncpy(iMesh_LAST_ERROR.description, descr, descr_len);
    RETURN(iBase_SUCCESS);
  }

  void iMesh_getError(iMesh_Instance instance,
                      int *err_type, char *descr, int *err, int descr_len) 
  {
    *err_type = iMesh_LAST_ERROR.error_type;
    strncpy(descr, iMesh_LAST_ERROR.description,
            MIN(strlen(iMesh_LAST_ERROR.description), (unsigned int)descr_len));
    RETURN(iBase_SUCCESS);
  }

  void iMesh_newMesh(const char *options, 
                     iMesh_Instance *instance, int *err, int options_len) 
  {
    std::string tmp_options(options, options_len);
    eatwhitespace(tmp_options);
    FileOptions opts((std::string(";") + tmp_options).c_str());

    MBInterface* core;

    MBErrorCode result = opts.get_null_option("PARALLEL");
    if (MB_SUCCESS == result) {
#ifdef USE_MPI    
      int flag = 1;
      int retval = MPI_Initialized(&flag);
      if (MPI_SUCCESS != retval || !flag) {
        int argc = 0;
        char **argv = NULL;
    
          // mpi not initialized yet - initialize here
        retval = MPI_Init(&argc, &argv);
      }
      core = new MBiMesh;
#else
        //mError->set_last_error( "PARALLEL option not valid, this instance"
        //                        " compiled for serial execution.\n" );
      *err = MB_NOT_IMPLEMENTED;
      return;
#endif
    }
    else core = new MBiMesh;

    *instance = reinterpret_cast<iMesh_Instance>(core);
    if (0 == *instance) {
      iMesh_processError(iBase_FAILURE, "Failed to instantiate mesh instance.");
      RETURN(iBase_FAILURE);
    }
  
    RETURN(iBase_SUCCESS);
  }

  void iMesh_dtor(iMesh_Instance instance, int *err) 
  {
    delete MBI;
    RETURN(iBase_SUCCESS);
  }
   
  void iMesh_load(iMesh_Instance instance,
                  const iBase_EntitySetHandle handle,
                  const char *name, const char *options, 
                  int *err, int name_len, int options_len) 
  {
      // get filename, option & null-terminate
    std::string tmp_filename(name, name_len), 
      tmp_options(options, options_len);

    eatwhitespace(tmp_filename);
    eatwhitespace(tmp_options);
    std::string opts = ";"; opts += tmp_options;
  
    MBEntityHandle file_set;
  
    MBErrorCode result = MBI->load_file(tmp_filename.c_str(), file_set, opts.c_str());

    if (MB_SUCCESS != result) {
      std::string msg("iMesh_load:ERROR loading a mesh, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      RETURN(iBase_ERROR_MAP[result]);
    }

    if (handle) {
      MBRange set_ents;
      result = MBI->get_entities_by_handle(file_set, set_ents);
      if (MB_SUCCESS != result) RETURN(iBase_ERROR_MAP[result]);
      result = MBI->add_entities(ENTITY_HANDLE(handle), set_ents);
      if (MB_SUCCESS != result) RETURN(iBase_ERROR_MAP[result]);
    }

      // create interior edges/faces if requested
    if (MBimesh->AdjTable[5] || MBimesh->AdjTable[10]) {
      MBRange set_ents;
      result = MBI->get_entities_by_handle(file_set, set_ents, true);
      if (MB_SUCCESS != result) RETURN(iBase_ERROR_MAP[result]);
      result = create_int_ents(MBI, set_ents, file_set);
      if (MB_SUCCESS != result) RETURN(iBase_ERROR_MAP[result]);
    }

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_save(iMesh_Instance instance,
                  const iBase_EntitySetHandle handle,
                  const char *name, const char *options, 
                  int *err, const int name_len, int options_len) 
  {
      // get filename & attempt to NULL-terminate
    std::string tmp_filename( name, name_len ), tmp_options(options, options_len);

    eatwhitespace(tmp_filename);
    eatwhitespace(tmp_options);
  
    MBEntityHandle set = ENTITY_HANDLE(handle);
    MBErrorCode result = MBI->write_file(tmp_filename.c_str(), NULL, tmp_options.c_str(),
                                         &set, 1);

    if (MB_SUCCESS != result) {
      std::string msg("iMesh_save:ERROR saving a mesh, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    }
  
    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_getRootSet(iMesh_Instance instance,
                        iBase_EntitySetHandle *root_set, int *err) 
  {
    *root_set = 0;
      //return CAST_TO_VOID(MBI->get_root_set());
    RETURN(iBase_SUCCESS);

  }

  void iMesh_getGeometricDimension(iMesh_Instance instance,
                                   int *geom_dim, int *err)
  {
    MBI->get_dimension(*geom_dim);
    RETURN(iBase_SUCCESS);
  }

  void iMesh_setGeometricDimension(iMesh_Instance instance,
                                   int geom_dim, int *err)
  {
    MBErrorCode rval = MBI->set_dimension(geom_dim);
    RETURN((MB_SUCCESS == rval ? iBase_SUCCESS : iBase_FAILURE));
  }

  void iMesh_getDfltStorage(iMesh_Instance instance,
                            int *order, int *err)
  {
    *order = iBase_BLOCKED;
    RETURN(iBase_SUCCESS);
  }
  
  void iMesh_getAdjTable (iMesh_Instance instance,
                          int** adjacency_table,
                          /*inout*/ int* adjacency_table_allocated, 
                          /*out*/ int* adjacency_table_size, int *err)
  {
    *adjacency_table_size = 16;
    CHECK_SIZE(*adjacency_table, *adjacency_table_allocated,
               *adjacency_table_size, int, iBase_MEMORY_ALLOCATION_FAILED);
    memcpy(*adjacency_table, MBimesh->AdjTable, 16*sizeof(int));
    RETURN(iBase_ERROR_MAP[MB_SUCCESS]);
  }

  void iMesh_setAdjTable (iMesh_Instance instance,
                          int* adj_table,
                          /*inout*/ int adj_table_size, 
                          int *err)
  {
    if (16 != adj_table_size) {
      RETURN(iBase_INVALID_ARGUMENT);
    }

    memcpy(MBimesh->AdjTable, adj_table, 16*sizeof(int));
    RETURN(iBase_SUCCESS);
  }

  void iMesh_getNumOfType(iMesh_Instance instance,
                          /*in*/ const iBase_EntitySetHandle entity_set_handle,
                          /*in*/ const int entity_type,
                          int *num_type, int *err)
  {
    iMesh_getNumOfTypeRec(instance, entity_set_handle, entity_type, false,
                          num_type, err);
  }
  
  void iMesh_getNumOfTopo(iMesh_Instance instance,
                          /*in*/ const iBase_EntitySetHandle entity_set_handle,
                          /*in*/ const int entity_topology,
                          int *num_topo, int *err)
  {
    iMesh_getNumOfTopoRec(instance, entity_set_handle, entity_topology,
                          false, num_topo, err);
  }

  void iMesh_areEHValid( iMesh_Instance instance,
                         int doReset,
                         int* areHandlesInvarient, 
                         int* err )
  {
    MBiMesh* mbi = dynamic_cast<MBiMesh*>(MBI);
    if (!mbi) {
      RETURN(iBase_FAILURE);
    }
  
    *areHandlesInvarient = !mbi->have_deleted_ents( !!doReset );
    RETURN(iBase_SUCCESS);
  }

  void iMesh_getEntities(iMesh_Instance instance,
                         /*in*/ const iBase_EntitySetHandle entity_set_handle,
                         /*in*/ const int entity_type,
                         /*in*/ const int entity_topology,
                         /*inout*/ iBase_EntityHandle** entity_handles,
                         /*inout*/ int* entity_handles_allocated,
                         /*out*/ int* entity_handles_size, 
                         int *err) 
  {
    iMesh_getEntitiesRec(instance, entity_set_handle, entity_type,
                         entity_topology, false,
                         entity_handles, entity_handles_allocated, entity_handles_size,
                         err);
  }

  void iMesh_getVtxArrCoords (iMesh_Instance instance,
                              /*in*/ const iBase_EntityHandle* vertex_handles,
                              /*in*/ const int vertex_handles_size,
                              /*inout*/ int storage_order,
                              /*inout*/ double** coords,
                              /*inout*/ int* coords_allocated,
                              /*out*/ int* coords_size, int *err) 
  {

      // make sure we can hold them all
    CHECK_SIZE(*coords, *coords_allocated, 3*vertex_handles_size, double, iBase_MEMORY_ALLOCATION_FAILED);
  
      // now get all the coordinates
      // coords will come back interleaved by default
    MBErrorCode result;
    if (storage_order == iBase_INTERLEAVED) {
      result = MBI->get_coords(CONST_HANDLE_ARRAY_PTR(vertex_handles), 
                               vertex_handles_size, *coords);
    }
  
    else {
      std::vector<double> dum_coords(3*vertex_handles_size);
      result = MBI->get_coords(CONST_HANDLE_ARRAY_PTR(vertex_handles), 
                               vertex_handles_size,
                               &dum_coords[0]);
      if (MB_SUCCESS == result) {
        int i;
        double *x = *coords;
        double *y = *coords+vertex_handles_size;
        double *z = *coords+2*vertex_handles_size;
        std::vector<double>::const_iterator c_iter = dum_coords.begin();
        for (i = 0; i < vertex_handles_size; i++) {
          *x = *c_iter; ++x; ++c_iter;
          *y = *c_iter; ++y; ++c_iter;
          *z = *c_iter; ++z; ++c_iter;
        }
      }
    }
  
    if (MB_SUCCESS != result) {
      std::string msg("iMesh_getVtxArrCoords: problem getting vertex coords, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      RETURN(iBase_ERROR_MAP[result]);
    }

    *coords_size = 3*vertex_handles_size;
  
    RETURN(iBase_SUCCESS);
  }

/**
 * Method:  initEntArrIter[]
 */
  void iMesh_initEntArrIter (iMesh_Instance instance,
                             /*in*/ const iBase_EntitySetHandle entity_set_handle,
                             /*in*/ const int requested_entity_type,
                             /*in*/ const int requested_entity_topology,
                             /*in*/ const int requested_array_size,
                             /*out*/ iMesh_EntityArrIterator* entArr_iterator,
                             int *err) 
  {
    MBEntityType req_type = mb_topology_table[requested_entity_topology];
    int req_dimension = (req_type == MBMAXTYPE ? (int) requested_entity_type : -1);
    MBRange range;
  
    MBErrorCode result;
    if (requested_entity_type == iBase_ALL_TYPES &&
        requested_entity_topology == iMesh_ALL_TOPOLOGIES)
      result = MBI->get_entities_by_handle( ENTITY_HANDLE(entity_set_handle),
                                            range );
    else if (requested_entity_topology == iMesh_SEPTAHEDRON)
      result = MB_SUCCESS; // never any septahedrons because MOAB doesn't support them
    else if (MBMAXTYPE != req_type)
      result = MBI->get_entities_by_type(ENTITY_HANDLE(entity_set_handle),
                                         req_type,
                                         range);
    else
      result = MBI->get_entities_by_dimension(ENTITY_HANDLE(entity_set_handle),
                                              req_dimension,
                                              range);

    if (result != MB_SUCCESS) {
      std::string msg("iMesh_initEntArrIter: ERROR getting entities of proper type or topology, "
                      "with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      iMesh_LAST_ERROR.error_type = iBase_ERROR_MAP[result];
      RETURN(iMesh_LAST_ERROR.error_type);
    }

    *entArr_iterator = reinterpret_cast<iMesh_EntityArrIterator>
                       (create_itaps_iterator( range, requested_array_size ));
    RETURN(iBase_SUCCESS);
  }

/**
 * Method:  getEntArrNextIter[]
 */
  void iMesh_getNextEntArrIter (iMesh_Instance instance,
                                /*in*/ iMesh_EntityArrIterator entArr_iterator,
                                /*inout*/ iBase_EntityHandle** entity_handles,
                                /*inout*/ int* entity_handles_allocated,
                                /*out*/ int* entity_handles_size,
                                int *has_data, int *err) 
  {
    RangeIterator *this_it = RANGE_ITERATOR(entArr_iterator);

      // check the size of the destination array
    int expected_size = (this_it->requestedSize < (int)this_it->iteratorRange.size() ? 
                         this_it->requestedSize : this_it->iteratorRange.size());
    CHECK_SIZE(*entity_handles, *entity_handles_allocated, expected_size,
               iBase_EntityHandle, false);
  
    int i = 0;
    while (i < this_it->requestedSize && this_it->currentPos != this_it->iteratorRange.end())
    {
      if (MBimesh->is_valid(*this_it->currentPos)) 
        (*entity_handles)[i++] = (iBase_EntityHandle)*(this_it->currentPos);
      ++(this_it->currentPos);
    }
  
    *has_data = (i!=0);
    *entity_handles_size = i;
    RETURN(iBase_SUCCESS);
  }

/**
 * Method:  resetEntArrIter[]
 */
  void iMesh_resetEntArrIter (iMesh_Instance instance,
                              /*in*/ iMesh_EntityArrIterator entArr_iterator, int *err) 
  {
    RangeIterator *this_it = RANGE_ITERATOR(entArr_iterator);

    this_it->currentPos = this_it->iteratorRange.begin();

    RETURN(iBase_SUCCESS);
  }

  void iMesh_endEntArrIter (iMesh_Instance instance,
                            /*in*/ iMesh_EntityArrIterator entArr_iterator, int *err) 
  {
    RangeIterator *this_it = RANGE_ITERATOR(entArr_iterator);

    this_it->currentPos = this_it->iteratorRange.end();

    RETURN(iBase_SUCCESS);
  }

  void iMesh_getEntArrTopo(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*inout*/ int** topology,
                           /*inout*/ int* topology_allocated,
                           /*out*/ int* topology_size, int *err) 
  {
      // go through each entity and look up its type
    CHECK_SIZE(*topology, *topology_allocated, entity_handles_size, 
               int, iBase_MEMORY_ALLOCATION_FAILED);

    for (int i = 0; i < entity_handles_size; i++)
      (*topology)[i] = 
        tstt_topology_table[MBI->type_from_handle(ENTITY_HANDLE(entity_handles[i]))];

    *topology_size = entity_handles_size;

    RETURN(iBase_SUCCESS);
  }
  
  void iMesh_getEntArrType(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*inout*/ int** etype,
                           /*inout*/ int* etype_allocated,
                           /*out*/ int* etype_size, int *err) 
  {
      // go through each entity and look up its type
    CHECK_SIZE(*etype, *etype_allocated, entity_handles_size, 
               int, iBase_MEMORY_ALLOCATION_FAILED);

    for (int i = 0; i < entity_handles_size; i++)
      (*etype)[i] = 
        tstt_type_table[MBI->type_from_handle(ENTITY_HANDLE(entity_handles[i]))];

    *etype_size = entity_handles_size;

    RETURN(iBase_SUCCESS);
  }

  void iMesh_getEntArrAdj(iMesh_Instance instance,
                          /*in*/ const iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size,
                          /*in*/ const int entity_type_requested,
                          /*inout*/ iBase_EntityHandle** adjacentEntityHandles,
                          /*inout*/ int* adjacentEntityHandles_allocated,
                          /*out*/ int* adjacentEntityHandles_size,
                          /*inout*/ int** offset,
                          /*inout*/ int* offset_allocated,
                          /*out*/ int* offset_size, int *err) 
  {
    MBErrorCode result = MB_SUCCESS;

    CHECK_SIZE(*offset, *offset_allocated, entity_handles_size+1, 
               int, iBase_MEMORY_ALLOCATION_FAILED);
  
    const MBEntityHandle* entity_iter = (const MBEntityHandle*)entity_handles;
    const MBEntityHandle* const entity_end = entity_iter + entity_handles_size;
    int* off_iter = *offset;
    int prev_off = 0;
  
    std::vector<MBEntityHandle> conn_storage;
    std::vector<MBEntityHandle> adj_ents;
    const MBEntityHandle *connect;
    int num_connect;
    
    MBEntityHandle* array; // ptr to working array of result handles
    int array_alloc;       // allocated size of 'array'
    const bool allocated_array = !*adjacentEntityHandles_allocated || !*adjacentEntityHandles;
    if (allocated_array) {
      array = 0;
      array_alloc = 0;
    }
    else {
      array = reinterpret_cast<MBEntityHandle*>(*adjacentEntityHandles);
      array_alloc = *adjacentEntityHandles_allocated;
    }
    
    for ( ; entity_iter != entity_end; ++entity_iter)
    {
      *off_iter = prev_off;
      off_iter++;
      
      if (iBase_VERTEX == entity_type_requested &&
          TYPE_FROM_HANDLE(*entity_iter) != MBPOLYHEDRON) {
        result = MBI->get_connectivity(*entity_iter, connect, num_connect, false, &conn_storage);
        if (MB_SUCCESS != result) {
          iMesh_processError(iBase_ERROR_MAP[result], "iMesh_getEntArrAdj: trouble getting adjacency list.");
          RETURN(iBase_ERROR_MAP[result]);
        }
      }
      else if (iBase_ALL_TYPES == entity_type_requested) {
        adj_ents.clear();
        for (int dim = 0; dim < 4; ++dim) {
          if (MBCN::Dimension(TYPE_FROM_HANDLE(*entity_iter)) == dim)
            continue;
          result = MBI->get_adjacencies( entity_iter, 1, dim, false, adj_ents );
          if (MB_SUCCESS != result) {
            iMesh_processError(iBase_ERROR_MAP[result], "iMesh_getEntArrAdj: trouble getting adjacency list.");
            RETURN(iBase_ERROR_MAP[result]);
          }
        }
        connect = &adj_ents[0];
        num_connect = adj_ents.size();
      }
      else {
        adj_ents.clear();
        result = MBI->get_adjacencies( entity_iter, 1, 
                                       entity_type_requested, false, adj_ents );
        if (MB_SUCCESS != result) {
          iMesh_processError(iBase_ERROR_MAP[result], "iMesh_getEntArrAdj: trouble getting adjacency list.");
          RETURN(iBase_ERROR_MAP[result]);
        }
        connect = &adj_ents[0];
        num_connect = adj_ents.size();
      }
      
      if (prev_off + num_connect <= array_alloc) {
        std::copy(connect, connect+num_connect, array+prev_off);
      }
      else if (allocated_array) {
          // if array is not allocated yet, guess at initial size
          // as the number of adjacencies for the first entity times
          // the number of input entities.  This will result in a single
          // exact allocation if one input entity or typical queries 
          // such as connectivity of a non-mixed mesh or regions adjacent
          // to faces.
        if (!array_alloc) 
          array_alloc = entity_handles_size * num_connect;
        else
          array_alloc = std::max(array_alloc*2, prev_off+num_connect);  
        array = (MBEntityHandle*)realloc( array, array_alloc*sizeof(MBEntityHandle) );
        if (!array) {
          RETURN(iBase_MEMORY_ALLOCATION_FAILED);
        }
        std::copy(connect, connect+num_connect, array+prev_off);
      }
      // else do nothing.  Will catch error later when comparing
      //  occupied to allocated sizes.  Continue here because
      //  must pass back required size.
      
      prev_off += num_connect;
    }
    *off_iter = prev_off;
    *offset_size = entity_handles_size+1;
    *adjacentEntityHandles_size = prev_off;

    if (*adjacentEntityHandles_size > array_alloc) {
      RETURN(iBase_BAD_ARRAY_SIZE);
    }
    else if (allocated_array) {
      *adjacentEntityHandles = reinterpret_cast<iBase_EntityHandle*>(array);
      *adjacentEntityHandles_allocated = array_alloc;
    }

    RETURN(iBase_SUCCESS);
  }

  void iMesh_getEntArr2ndAdj( iMesh_Instance instance,
                              iBase_EntityHandle const* entity_handles,
                              int entity_handles_size,
                              int order_adjacent_key,
                              int requested_entity_type,
                              iBase_EntityHandle** adj_entity_handles,
                              int* adj_entity_handles_allocated,
                              int* adj_entity_handles_size,
                              int** offset,
                              int* offset_allocated,
                              int* offset_size,
                              int* err )
  {
    MBErrorCode result = MB_SUCCESS;

    CHECK_SIZE(*offset, *offset_allocated, entity_handles_size+1, 
               int, iBase_MEMORY_ALLOCATION_FAILED);
  
    const MBEntityHandle* entity_iter = (const MBEntityHandle*)entity_handles;
    const MBEntityHandle* const entity_end = entity_iter + entity_handles_size;
    int* off_iter = *offset;
    int prev_off = 0;
  
    std::vector<MBEntityHandle> all_adj_ents;
    MeshTopoUtil mtu(MBI);
    
    for ( ; entity_iter != entity_end; ++entity_iter)
    {
      *off_iter = prev_off;
      off_iter++;
      MBRange adj_ents;

      result = mtu.get_bridge_adjacencies( *entity_iter,
                                           order_adjacent_key,
                                           requested_entity_type, adj_ents );
      if (MB_SUCCESS != result) {
        iMesh_processError(iBase_ERROR_MAP[result], "iMesh_getEntArrAdj: trouble getting adjacency list.");
        RETURN(iBase_ERROR_MAP[result]);
      }

      std::copy(adj_ents.begin(), adj_ents.end(), std::back_inserter(all_adj_ents));
      prev_off += adj_ents.size();
    }
    *off_iter = prev_off;

    CHECK_SIZE(*adj_entity_handles, *adj_entity_handles_allocated, 
               (int)all_adj_ents.size(), 
               iBase_EntityHandle, iBase_MEMORY_ALLOCATION_FAILED);
    memcpy(*adj_entity_handles, &all_adj_ents[0], 
           sizeof(MBEntityHandle)*all_adj_ents.size() );

    *adj_entity_handles_size = all_adj_ents.size();
    *offset_size = entity_handles_size+1;
    RETURN(iBase_SUCCESS);
  }

  void iMesh_getAdjEntIndices(iMesh_Instance instance,
                      /*in*/    iBase_EntitySetHandle entity_set_handle,
                      /*in*/    int entity_type_requestor,
                      /*in*/    int entity_topology_requestor,
                      /*in*/    int entity_type_requested,
                      /*inout*/ iBase_EntityHandle** entity_handles,
                      /*inout*/ int* entity_handles_allocated,
                      /*out*/   int* entity_handles_size,
                      /*inout*/ iBase_EntityHandle** adj_entity_handles,
                      /*inout*/ int* adj_entity_handles_allocated,
                      /*out*/   int* adj_entity_handles_size,
                      /*inout*/ int** adj_entity_indices,
                      /*inout*/ int* adj_entity_indices_allocated,
                      /*out*/   int* adj_entity_indices_size,
                      /*inout*/ int** offset,
                      /*inout*/ int* offset_allocated,
                      /*out*/   int* offset_size,
                      /*out*/   int *err)
  {
    const int allocated_entity_handles = (*entity_handles_allocated == 0);
    const int allocated_indices = (*adj_entity_indices_allocated == 0);
    const int allocated_offset = (*offset_allocated == 0);

    // get source entities
    iMesh_getEntities( instance, 
                       entity_set_handle,
                       entity_type_requestor, 
                       entity_topology_requestor,
                       entity_handles,
                       entity_handles_allocated,
                       entity_handles_size,
                       err );
    if (iBase_SUCCESS != *err)
      return;
    
    // get adjacencies
    iBase_EntityHandle* all_adj_handles = 0;
    int size = 0, alloc = 0;
    iMesh_getEntArrAdj( instance,
                        *entity_handles, *entity_handles_size,
                        entity_type_requested,
                        &all_adj_handles, &alloc, &size,
                        offset, offset_allocated, offset_size,
                        err );
    if (*err != iBase_SUCCESS) {
      if (allocated_entity_handles) {
        free( *entity_handles );
        *entity_handles = 0;
        *entity_handles_allocated = 0;
      }
      return;
    }
    
    // allocate or check size of adj_entity_indices
    *adj_entity_indices_size = size;
    if (allocated_indices) {
      *adj_entity_indices = (int*)malloc(sizeof(iBase_EntityHandle)*size);
      if (!*adj_entity_indices) 
        *err = iBase_MEMORY_ALLOCATION_FAILED;
      else
        *adj_entity_indices_allocated = size;
    }
    else if (*adj_entity_indices_allocated < size) {
      *err = iBase_BAD_ARRAY_DIMENSION;
    }
    if (iBase_SUCCESS != *err) {
      free( all_adj_handles );
      if (allocated_entity_handles) {
        free( *entity_handles );
        *entity_handles = 0;
        *entity_handles_allocated = 0;
      }
      if (allocated_offset) {
        free( *offset );
        *offset = 0;
        *offset_allocated = 0;
      }
      return;
    }
    
    // Now create an array of unique sorted handles from all_adj_handles.
    // We need to create a copy because we still need all_adj_handles.  We
    // will eventually need to copy the resulting unique list into 
    // adj_entity_handles, so if adj_entity_handles is already allocated and
    // of sufficient size, use it rather than allocating another temporary.
    iBase_EntityHandle* unique_adj = 0;
    if (*adj_entity_handles_allocated >= size) {
      unique_adj = *adj_entity_handles;
    }
    else {
      unique_adj = (iBase_EntityHandle*)malloc(sizeof(iBase_EntityHandle) * size);
    }
    std::copy( all_adj_handles, all_adj_handles+size, unique_adj );
    std::sort( unique_adj, unique_adj + size );
    *adj_entity_handles_size = std::unique( unique_adj, unique_adj + size ) - unique_adj;
    
    // If we created a temporary array for unique_adj rather than using
    // already allocated space in adj_entity_handles, allocate adj_entity_handles
    // and copy the unique handle list into it
    if (*adj_entity_handles != unique_adj) {
      if (!*adj_entity_handles_allocated) {
        *adj_entity_handles = (iBase_EntityHandle*)malloc(
                              sizeof(iBase_EntityHandle) * *adj_entity_handles_size);
        if (!*adj_entity_handles)
          *err = iBase_MEMORY_ALLOCATION_FAILED;
        else
          *adj_entity_handles_allocated = *adj_entity_handles_size;
      }
      else if (*adj_entity_handles_allocated < *adj_entity_handles_size) 
        *err = iBase_BAD_ARRAY_DIMENSION;
      if (iBase_SUCCESS != *err) {
        free( unique_adj );
        free( all_adj_handles );
        if (allocated_entity_handles) {
          free( *entity_handles );
          *entity_handles = 0;
          *entity_handles_allocated = 0;
        }
        if (allocated_offset) {
          free( *offset );
          *offset = 0;
          *offset_allocated = 0;
        }
        if (allocated_indices) {
          free( *adj_entity_indices );
          *adj_entity_indices = 0;
          *adj_entity_indices_allocated = 0;
        }
        return;
      }

      std::copy( unique_adj, unique_adj + *adj_entity_handles_size, *adj_entity_handles );
      free( unique_adj );
      unique_adj = *adj_entity_handles;
    }
    
    // convert from adjacency list to indices into unique_adj
    for (int i = 0; i < *adj_entity_indices_size; ++i)
      (*adj_entity_indices)[i] = std::lower_bound( unique_adj, 
        unique_adj + *adj_entity_handles_size, all_adj_handles[i] ) - unique_adj;
    free( all_adj_handles );
  }


  void iMesh_createEntSet(iMesh_Instance instance,
                          /*in*/ const int isList,
                          /*out*/ iBase_EntitySetHandle* entity_set_created, int *err) 
  {
      // create the entity set
    MBEntityHandle meshset;
    MBErrorCode result;

    if (isList)
      result = MBI->create_meshset(MESHSET_ORDERED, meshset);
    else
      result = MBI->create_meshset(MESHSET_SET, meshset);
  
    if (MB_SUCCESS != result) {
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_createEntSet: ERROR creating a entityset instance");
      RETURN(iBase_ERROR_MAP[result]);
    }
  
      // return EntitySet_Handle
    *entity_set_created = (iBase_EntitySetHandle)meshset;
    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_destroyEntSet (iMesh_Instance instance,
                            /*in*/ iBase_EntitySetHandle entity_set, int *err) 
  {
    MBEntityHandle set = ENTITY_HANDLE(entity_set);
    MBErrorCode result = MBI->delete_entities(&set, 1);
    if (MB_SUCCESS != result)
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_destroyEntSet: couldn't delete the set.");

    RETURN(iBase_SUCCESS);
  }

  void iMesh_isList (iMesh_Instance instance,
                     /*in*/ const iBase_EntitySetHandle entity_set,
                     int *is_list, int *err) 
  {
    unsigned int options;
    MBErrorCode result = MBI->get_meshset_options(ENTITY_HANDLE(entity_set), options);
    if (MB_SUCCESS != result)
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_isList: couldn't query set.");
    if (options & MESHSET_ORDERED)
      *is_list = true;
    else *is_list = false;
  
    RETURN(iBase_SUCCESS);
  }

  void iMesh_getNumEntSets(iMesh_Instance instance,
                           /*in*/ const iBase_EntitySetHandle entity_set_handle,
                           /*in*/ const int num_hops,
                           int *num_sets, int *err) 
  {
    MBErrorCode rval = MBI->num_contained_meshsets( ENTITY_HANDLE(entity_set_handle),
                                                    num_sets,
                                                    std::max(0,num_hops) );
    if (rval != MB_SUCCESS) {
      std::string msg("iMesh_entitysetGetNumberEntitySets:ERROR getting number of entitysets "
                      "in EntitySet, with error type: ");
      msg += MBI->get_error_string(rval);
      iMesh_processError(iBase_ERROR_MAP[rval], msg.c_str());
    }

    RETURN(iBase_ERROR_MAP[rval]);
  } 

  void iMesh_getEntSets(iMesh_Instance instance,
                        /*in*/ const iBase_EntitySetHandle entity_set_handle,
                        /*in*/ const int num_hops,
                        /*inout*/ iBase_EntitySetHandle** contained_entset_handles,
                        /*inout*/ int* contained_entset_handles_allocated,
                        /*inout*/ int* contained_entset_handles_size, int *err) 
  {
    std::vector<MBEntityHandle> sets;
    MBErrorCode rval = MBI->get_contained_meshsets( ENTITY_HANDLE(entity_set_handle),
                                                    sets, 
                                                    std::max( num_hops, 0 ) );
    if (MB_SUCCESS != rval) {
      iMesh_processError(iBase_ERROR_MAP[rval], "iMesh_entitysetGetEntitySets: problem getting entities by type.");
      RETURN(iBase_ERROR_MAP[rval]);
    }
    CHECK_SIZE(*contained_entset_handles, *contained_entset_handles_allocated,
               (int)sets.size(), iBase_EntitySetHandle, iBase_MEMORY_ALLOCATION_FAILED);

    std::copy( sets.begin(), sets.end(), (MBEntityHandle*)*contained_entset_handles );
    *contained_entset_handles_size = sets.size();
    RETURN(iBase_SUCCESS);
  }

  void iMesh_addEntArrToSet(iMesh_Instance instance,
                            /*in*/ const iBase_EntityHandle* entity_handles,
                            /*in*/ int entity_handles_size,
                            /*in*/ iBase_EntitySetHandle entity_set, 
                            int *err)
  {
    const MBEntityHandle *ents = CONST_HANDLE_ARRAY_PTR(entity_handles);
    MBErrorCode result = MBI->add_entities(ENTITY_HANDLE(entity_set),
                                           ents, entity_handles_size);

    if (result != MB_SUCCESS) {
      std::string msg("iMesh_addEntArrToSet:ERROR adding entities in EntitySet, "
                      "with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    }
  
    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_addEntToSet(iMesh_Instance instance,
                         /*in*/ iBase_EntityHandle entity_handle,
                         /*in*/ iBase_EntitySetHandle entity_set, int *err)
  {
    iMesh_addEntArrToSet(instance, &entity_handle, 1, entity_set, err);
  }

  void iMesh_rmvEntArrFromSet(iMesh_Instance instance,
                              /*in*/ const iBase_EntityHandle* entity_handles,
                              /*in*/ int entity_handles_size,
                              /*in*/ iBase_EntitySetHandle entity_set, int *err)
  {
    const MBEntityHandle *ents = CONST_HANDLE_ARRAY_PTR(entity_handles);

    MBErrorCode result = MBI->remove_entities
      (ENTITY_HANDLE(entity_set), ents, entity_handles_size);
  
    if (result != MB_SUCCESS) {
      std::string msg("iMesh_rmvEntArrFromSet:ERROR removing entities in EntitySet, "
                      "with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    }
  
    RETURN(iBase_ERROR_MAP[result]);
  }
  
  void iMesh_rmvEntFromSet(iMesh_Instance instance,
                           /*in*/ iBase_EntityHandle entity_handle,
                           /*in*/ iBase_EntitySetHandle entity_set, 
                           int *err)
  {
    iMesh_rmvEntArrFromSet(instance, &entity_handle, 1, entity_set, err);
  }
  
  void iMesh_addEntSet(iMesh_Instance instance,
                       /*in*/ iBase_EntitySetHandle entity_set_to_add,
                       /*in*/ iBase_EntitySetHandle entity_set_handle,
                       int *err)
  {
    MBEntityHandle to_add = ENTITY_HANDLE(entity_set_to_add);
    MBErrorCode result = MBI->add_entities(ENTITY_HANDLE(entity_set_handle), &to_add, 1);

    if (result != MB_SUCCESS) {
      std::string msg("iMesh_addEntSet:ERROR adding entitysets, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    }
  

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_rmvEntSet(iMesh_Instance instance,
                       /*in*/ iBase_EntitySetHandle entity_set_to_remove,
                       /*in*/ iBase_EntitySetHandle entity_set_handle, 
                       int *err)
  {
    MBEntityHandle to_remove = ENTITY_HANDLE(entity_set_to_remove);
    MBErrorCode result = MBI->remove_entities
      (ENTITY_HANDLE(entity_set_handle), &to_remove, 1);
  
    if (result != MB_SUCCESS) {
      std::string msg("iMesh_rmvEntSet:ERROR removing entitysets in EntitySet, "
                      "with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    }
  
    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_isEntContained (iMesh_Instance instance,
                             /*in*/ iBase_EntitySetHandle containing_entity_set,
                             /*in*/ iBase_EntityHandle contained_entity,
                             int *is_contained, int *err) 
  {
    int junk1 = 1, junk2 = 1;
    iMesh_isEntArrContained( instance, containing_entity_set, 
                             &contained_entity, 1, &is_contained, 
                             &junk1, &junk2, err );
  }

  void iMesh_isEntArrContained( iMesh_Instance instance,
                            /*in*/ iBase_EntitySetHandle containing_set,
                            /*in*/ const iBase_EntityHandle* entity_handles,
                            /*in*/ int num_entity_handles,
                         /*inout*/ int** is_contained,
                         /*inout*/ int* is_contained_allocated,
                           /*out*/ int* is_contained_size,
                           /*out*/ int* err )

  {
    MBEntityHandle set = ENTITY_HANDLE(containing_set);
    CHECK_SIZE(*is_contained, *is_contained_allocated,
               (int)num_entity_handles, int, iBase_MEMORY_ALLOCATION_FAILED);
    *is_contained_size = num_entity_handles;
    
    if (containing_set) {
      for (int i = 0; i < num_entity_handles; ++i) {
        MBEntityHandle h = ENTITY_HANDLE(entity_handles[i]);
        (*is_contained)[i] = MBI->contains_entities( set, &h, 1 );
      }
    }
    else {
      std::fill( *is_contained, (*is_contained)+num_entity_handles, 1 );
    }
    RETURN(iBase_SUCCESS);
  }

  void iMesh_isEntSetContained (iMesh_Instance instance,
                                /*in*/ const iBase_EntitySetHandle containing_entity_set,
                                /*in*/ const iBase_EntitySetHandle contained_entity_set,
                                int *is_contained, int *err) 
  {
    iMesh_isEntContained(instance, containing_entity_set, 
                         reinterpret_cast<iBase_EntityHandle>(contained_entity_set),
                         is_contained, err);
  }

  void iMesh_addPrntChld(iMesh_Instance instance,
                         /*inout*/ iBase_EntitySetHandle parent_entity_set,
                         /*inout*/ iBase_EntitySetHandle child_entity_set, 
                         int *err) 
  {
    MBErrorCode result = MBI->add_parent_child
      (ENTITY_HANDLE(parent_entity_set),
       ENTITY_HANDLE(child_entity_set));

    if (result != MB_SUCCESS) {
      std::string msg("MB Mesh::addPrntChld: ERROR addParentChild failed, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    }
  

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_rmvPrntChld(iMesh_Instance instance,
                         /*inout*/ iBase_EntitySetHandle parent_entity_set,
                         /*inout*/ iBase_EntitySetHandle child_entity_set, 
                         int *err)
  {
    MBErrorCode result = MBI->remove_parent_child
      (ENTITY_HANDLE(parent_entity_set),
       ENTITY_HANDLE(child_entity_set));
  
    if (result != MB_SUCCESS) {
      std::string msg("iMesh_rmvPrntChld: ERROR RemoveParentChild failed, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    }

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_isChildOf(iMesh_Instance instance,
                       /*in*/ const iBase_EntitySetHandle parent_entity_set,
                       /*in*/ const iBase_EntitySetHandle child_entity_set,
                       int *is_child, int *err)
  {
    std::vector<MBEntityHandle> children;

    MBErrorCode result = MBI->get_child_meshsets
      (ENTITY_HANDLE(parent_entity_set), children);

    if (result != MB_SUCCESS) {
      std::string msg("iMesh_isChildOf: ERROR IsParentChildRelated failed, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      *is_child = false;
      RETURN(iBase_ERROR_MAP[result]);
    }
  

  
    if (std::find(children.begin(), children.end(), ENTITY_HANDLE(child_entity_set))
        != children.end())
      *is_child = true;

    else
      *is_child = false;

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_getNumChld(iMesh_Instance instance,
                        /*in*/ const iBase_EntitySetHandle entity_set,
                        /*in*/ const int num_hops,
                        int *num_child, int *err)
  {
    *num_child = 0;
    MBErrorCode result = MBI->num_child_meshsets
      (ENTITY_HANDLE(entity_set), num_child, num_hops);

    if (result != MB_SUCCESS) {
      std::string msg("iMesh_getNumChld: ERROR GetNumChildren failed, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      *num_child = 0;
      RETURN(iBase_ERROR_MAP[result]);
    }

    RETURN(iBase_SUCCESS);
  }

  void iMesh_getNumPrnt(iMesh_Instance instance,
                        /*in*/ const iBase_EntitySetHandle entity_set,
                        /*in*/ const int num_hops,
                        int *num_parent, int *err)
  {
    *num_parent = 0;
    MBErrorCode result = MBI->num_parent_meshsets
      (ENTITY_HANDLE(entity_set), num_parent, num_hops);

    if (result != MB_SUCCESS) {
      std::string msg("iMesh_getNumPrnt:\
           ERROR GetNumParents failed, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      RETURN(iBase_ERROR_MAP[result]);
    }

    RETURN(iBase_SUCCESS);
  }

  void iMesh_getChldn(iMesh_Instance instance,
                      /*in*/ const iBase_EntitySetHandle from_entity_set,
                      /*in*/ const int num_hops,
                      /*out*/ iBase_EntitySetHandle** entity_set_handles,
                      /*out*/ int* entity_set_handles_allocated,
                      /*out*/ int* entity_set_handles_size, int *err) 
  {
    std::vector<MBEntityHandle> children;

    MBErrorCode result = MBI->get_child_meshsets
      (ENTITY_HANDLE(from_entity_set), children, num_hops);

    if (result != MB_SUCCESS) {
      std::string msg("iMesh_getChldn:\
           ERROR getChildren failed, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      RETURN(iBase_ERROR_MAP[result]);
    }

    CHECK_SIZE(*entity_set_handles, *entity_set_handles_allocated,
               (int)children.size(), iBase_EntitySetHandle, iBase_MEMORY_ALLOCATION_FAILED);

    MBEntityHandle *ents = HANDLE_ARRAY_PTR(*entity_set_handles);
      // use a memcpy for efficiency
    memcpy(ents, &children[0], children.size()*sizeof(MBEntityHandle));

    *entity_set_handles_size = children.size();

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_getPrnts(iMesh_Instance instance,
                      /*in*/ const iBase_EntitySetHandle from_entity_set,
                      /*in*/ const int num_hops,
                      /*out*/ iBase_EntitySetHandle** entity_set_handles,
                      /*out*/ int* entity_set_handles_allocated,
                      /*out*/ int* entity_set_handles_size, int *err) 
  {
    std::vector<MBEntityHandle> parents;

    MBErrorCode result = MBI->get_parent_meshsets
      (ENTITY_HANDLE(from_entity_set), parents, num_hops);

    if (result != MB_SUCCESS) {
      std::string msg("iMesh_getPrnts:\
           ERROR getParents failed, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      RETURN(iBase_ERROR_MAP[result]);
    }

    CHECK_SIZE(*entity_set_handles, *entity_set_handles_allocated,
               (int)parents.size(), iBase_EntitySetHandle, iBase_MEMORY_ALLOCATION_FAILED);

    MBEntityHandle *ents = HANDLE_ARRAY_PTR(*entity_set_handles);
      // use a memcpy for efficiency
    memcpy(ents, &parents[0], parents.size()*sizeof(MBEntityHandle));

    *entity_set_handles_size = parents.size();

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_setVtxArrCoords (iMesh_Instance instance,
                              /*in*/ const iBase_EntityHandle* vertex_handles,
                              /*in*/ const int vertex_handles_size,
                              /*in*/ const int storage_order,
                              /*in*/ const double* new_coords,
                              /*in*/ const int new_coords_size, int *err) 
  {
    MBErrorCode result = MB_SUCCESS, tmp_result;
    if (storage_order == iBase_INTERLEAVED) {
      result = MBI->set_coords(CONST_HANDLE_ARRAY_PTR(vertex_handles),
                               vertex_handles_size, new_coords);
    }
    else {
      const MBEntityHandle *verts = CONST_HANDLE_ARRAY_PTR(vertex_handles);
      double dummy[3];
      for (int i = 0; i < vertex_handles_size; i++) {
        dummy[0] = new_coords[i]; dummy[1] = new_coords[vertex_handles_size+i]; 
        dummy[2] = new_coords[2*vertex_handles_size+i];
        tmp_result = MBI->set_coords(&verts[i], 1, dummy);
        if (MB_SUCCESS != tmp_result) result = tmp_result;
      }
    }
  
    if (MB_SUCCESS != result)
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_setVtxArrCoords: problem setting coordinates.");
  
    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_createVtxArr(iMesh_Instance instance,
                          /*in*/ const int num_verts,
                          /*in*/ const int storage_order,
                          /*in*/ const double* new_coords,
                          /*in*/ const int new_coords_size,
                          /*inout*/ iBase_EntityHandle** new_vertex_handles,
                          /*inout*/ int* new_vertex_handles_allocated,
                          /*inout*/ int* new_vertex_handles_size, int *err) 
  {
    if (new_coords_size != 3*num_verts) {
      iMesh_processError(iBase_INVALID_ARGUMENT, "iMesh_createVtxArr: Didn't get the right # coordinates.");
      RETURN(iBase_INVALID_ARGUMENT);
    }

      // if there aren't any elements in the array, allocate it
    CHECK_SIZE(*new_vertex_handles, *new_vertex_handles_allocated,
               num_verts, iBase_EntityHandle, iBase_MEMORY_ALLOCATION_FAILED);
  
      // make the entities
    MBEntityHandle *new_verts = HANDLE_ARRAY_PTR(*new_vertex_handles);
  
    for (int i = 0; i < num_verts; i++) {
      MBErrorCode result = MBI->create_vertex(&new_coords[3*i], new_verts[i]);
      if (MB_SUCCESS != result) {
        iMesh_processError(iBase_ERROR_MAP[result], "iMesh_createVtxArr: couldn't create vertex.");
        RETURN(iBase_ERROR_MAP[result]);
      }
    }  

    *new_vertex_handles_size = new_coords_size/3;

    RETURN(iBase_SUCCESS);
  }
                                                   
  void iMesh_createEntArr(iMesh_Instance instance,
                          /*in*/ const int new_entity_topology,
                          /*in*/ const iBase_EntityHandle* lower_order_entity_handles,
                          /*in*/ const int lower_order_entity_handles_size,
                          /*out*/ iBase_EntityHandle** new_entity_handles,
                          /*out*/ int* new_entity_handles_allocated,
                          /*out*/ int* new_entity_handles_size,
                          /*inout*/  int** status,
                          /*inout*/ int* status_allocated,
                          /*out*/ int* status_size, int *err) 
  {
      // for now, throw an error if lower order entity handles aren't vertices
    MBEntityType this_type = mb_topology_table[new_entity_topology];
    int num_ents = 0, num_verts;
    const MBEntityHandle *lower_ents;
    if (MBVERTEX != this_type) {
      num_verts = MBCN::VerticesPerEntity(this_type);
      num_ents = lower_order_entity_handles_size / num_verts;
      lower_ents = CONST_HANDLE_ARRAY_PTR(lower_order_entity_handles);
        // check that we have the right number of lower order entity handles
      if (lower_order_entity_handles_size % MBCN::VerticesPerEntity(this_type) != 0) {
        iMesh_processError(iBase_INVALID_ENTITY_COUNT, "iMesh_createEntArr: wrong # vertices for this entity type.");
        RETURN(iBase_INVALID_ENTITY_COUNT);
      }
    }
    else {
      iMesh_processError(iBase_INVALID_ARGUMENT, "iMesh_createEntArr: can't create vertices with this function, use createVtxArr instead.");
      RETURN(iBase_INVALID_ARGUMENT);
    }
  
    if (num_ents == 0) {
      iMesh_processError(iBase_INVALID_ENTITY_COUNT, 
                         "iMesh_createEntArr: called to create 0 entities.");
      RETURN(iBase_INVALID_ENTITY_COUNT);
    }

      // if there aren't any elements in the array, allocate it
    CHECK_SIZE(*new_entity_handles, *new_entity_handles_allocated,
               num_ents, iBase_EntityHandle, iBase_MEMORY_ALLOCATION_FAILED);
  
    CHECK_SIZE(*status, *status_allocated, num_ents, 
               int, iBase_MEMORY_ALLOCATION_FAILED);
  
      // make the entities
    MBEntityHandle *new_ents = HANDLE_ARRAY_PTR(*new_entity_handles);

    MBErrorCode tmp_result, result = MB_SUCCESS;
  
    for (int i = 0; i < num_ents; i++) {
      tmp_result = MBI->create_element(this_type, lower_ents, num_verts,
                                       new_ents[i]);
      if (MB_SUCCESS != tmp_result) {
        (*status)[i] = iBase_CREATION_FAILED;
        result = tmp_result;
      }
      else
        (*status)[i] = iBase_NEW;
    
      lower_ents += num_verts;
    }

    if (MB_SUCCESS != result)
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_createEntArr: couldn't create one of the entities.");

    if (MB_SUCCESS == result) {
      *new_entity_handles_size = num_ents;
      *status_size = num_ents;
    }

    if (MBimesh->AdjTable[5] || MBimesh->AdjTable[10]) {
      MBRange set_ents;
      std::copy(HANDLE_ARRAY_PTR(*new_entity_handles), 
                HANDLE_ARRAY_PTR(*new_entity_handles)+*new_entity_handles_size,
                mb_range_inserter(set_ents));
      result = create_int_ents(MBI, set_ents);
    }

    RETURN(iBase_ERROR_MAP[result]);
  }
                                                   
  void iMesh_deleteEntArr(iMesh_Instance instance,
                          /*in*/ const iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size, int *err) 
  {
    if (0 == entity_handles_size) {
      RETURN(iBase_SUCCESS);
    }

    MBErrorCode result = MBI->delete_entities(CONST_HANDLE_ARRAY_PTR(entity_handles),
                                              entity_handles_size);
    if (MB_SUCCESS != result)
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_deleteEntArr: trouble deleting entities.");

    RETURN(iBase_ERROR_MAP[result]);
  }
                                                 
  void iMesh_createTag(iMesh_Instance instance,
                       /*in*/ const char* tag_name,
                       /*in*/ const int tag_size,
                       /*in*/ const int tag_type,
                       /*out*/ iBase_TagHandle* tag_handle, 
                       int *err,
                       const int tag_name_size)
  {
    MBTag new_tag;
    int this_size = tag_size;

    std::string tmp_tagname(tag_name);
    eatwhitespace(tmp_tagname);

    switch (tag_type) {
      case iBase_INTEGER:
        this_size *= sizeof(int);
        break;
      case iBase_DOUBLE:
        this_size *= sizeof(double);
        break;
      case iBase_ENTITY_HANDLE:
        this_size *= sizeof(iBase_EntityHandle);
        break;
      case iBase_BYTES:
        break;
    }
      
    MBErrorCode result = MBI->tag_create(tmp_tagname.c_str(), this_size,
                                         MB_TAG_SPARSE, 
                                         mb_data_type_table[tag_type],
                                         new_tag,
                                         NULL);

    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) {
      std::string msg = std::string("iMesh_createTag: error creating tag with name '") +
        std::string(tag_name) + std::string("'.");
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    }
    else if (MB_ALREADY_ALLOCATED == result) {
      std::string msg = std::string("iMesh_createTag: tag with name '") +
        std::string(tag_name) + std::string("' already created.");
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    }

    *tag_handle = (iBase_TagHandle) new_tag;

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_destroyTag(iMesh_Instance instance,
                        /*in*/ iBase_TagHandle tag_handle,
                        /*in*/ const int forced, int *err)
  {
      // might need to check if it's used first
    if (false == forced) {
      MBRange ents;
      MBErrorCode result;
      MBTag this_tag = TAG_HANDLE(tag_handle);
      for (MBEntityType this_type = MBVERTEX; this_type != MBMAXTYPE; this_type++) {
        result = MBI->get_entities_by_type_and_tag(0, this_type, &this_tag, NULL, 1, 
                                                   ents, MBInterface::UNION);
        if (result != MB_SUCCESS) {
          std::string msg("iMesh_destroyTag: problem finding tag., with error type: ");
          msg += MBI->get_error_string(result);
          iMesh_processError(iBase_ERROR_MAP[result],
                             msg.c_str());
        }
        else if (!ents.empty()) {
          iMesh_processError(iBase_FAILURE, "iMesh_destroyTag: forced=false and entities"
                             " are still assigned this tag.");
          RETURN(iBase_FAILURE);
        }
      }
        // check if tag value is set on mesh
      const void* data_ptr;
      result = MBI->tag_get_data( this_tag, 0, 0, &data_ptr );
      if (MB_SUCCESS == result) {
        iMesh_processError(iBase_FAILURE, "iMesh_destroyTag: forced=false and mesh"
                           " is still assigned this tag.");
        RETURN(iBase_FAILURE);
      }
    }
  
      // ok, good to go - either forced or no entities with this tag
    MBErrorCode result = MBI->tag_delete(TAG_HANDLE(tag_handle));
    if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result)
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_destroyTag: problem deleting tag.");

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_getTagName(iMesh_Instance instance,
                        /*in*/ const iBase_TagHandle tag_handle,
                        char *out_data, int *err,
                        int out_data_len)
  {
    static ::std::string name;
    MBErrorCode result = MBI->tag_get_name(TAG_HANDLE(tag_handle), name);
    if (MB_SUCCESS != result)
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_getTagName: problem getting name.");
        else
          iMesh_LAST_ERROR.error_type = iBase_SUCCESS;

    strncpy(out_data, name.c_str(), out_data_len);

    RETURN(iMesh_LAST_ERROR.error_type);
  }

  void iMesh_getTagType (iMesh_Instance instance,
                         /*in*/ const iBase_TagHandle tag_handle,
                         int *value_type, int *err) 
  {
    MBDataType this_type;
    MBErrorCode result = MBI->tag_get_data_type(TAG_HANDLE(tag_handle),
                                                this_type);
    if (MB_SUCCESS != result) {
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_getTagType: problem getting type.");
      *value_type = iBase_BYTES;
      RETURN(iBase_ERROR_MAP[result]);
    }
  
    else
    {
      iMesh_LAST_ERROR.error_type = iBase_SUCCESS;
      *value_type = tstt_data_type_table[this_type];
      RETURN(iBase_ERROR_MAP[result]);
    }
  }

  void iMesh_getTagSizeValues(iMesh_Instance instance,
                              /*in*/ const iBase_TagHandle tag_handle,
                              int *tag_size_val, int *err)
  {
    MBErrorCode result = MBI->tag_get_size(TAG_HANDLE(tag_handle), *tag_size_val);
    if (MB_SUCCESS != result)
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_getTagSize: problem getting size.");
        else {
          int this_type;
          int success;
          iMesh_getTagType(instance, tag_handle, &this_type, &success);
          if (iBase_ERROR_MAP[success] != iBase_SUCCESS) 
            return;
          switch (this_type) {
            case iBase_INTEGER:
              *tag_size_val /= sizeof(int);
              break;
            case iBase_DOUBLE:
              *tag_size_val /= sizeof(double);
              break;
            case iBase_ENTITY_HANDLE:
              *tag_size_val /= sizeof(MBEntityHandle);
              break;
            case iBase_BYTES:
              break;
          }
        
          iMesh_LAST_ERROR.error_type = iBase_SUCCESS;
        }

  
    RETURN(iMesh_LAST_ERROR.error_type);
  }

  void iMesh_getTagSizeBytes(iMesh_Instance instance,
                             /*in*/ const iBase_TagHandle tag_handle,
                             int *tag_size_bytes, int *err)
  {
    MBErrorCode result = MBI->tag_get_size(TAG_HANDLE(tag_handle), *tag_size_bytes);
    if (MB_SUCCESS != result)
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_getTagSize: problem getting size.");
        else {
          iMesh_LAST_ERROR.error_type = iBase_SUCCESS;
        }
  
    RETURN(iMesh_LAST_ERROR.error_type);
  }

  void iMesh_getTagHandle(iMesh_Instance instance,
                          /*in*/ const char* tag_name,
                          iBase_TagHandle *tag_handle, int *err,
                          const int tag_name_len)
  {
    MBErrorCode result = MBI->tag_get_handle(tag_name, (MBTag&)*tag_handle);
    
    if (MB_SUCCESS != result) {
      std::string msg("iMesh_getTagHandle: problem getting handle for tag named '");
      msg += std::string(tag_name) + std::string("'");
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      *tag_handle = 0;
      RETURN(iBase_ERROR_MAP[result]);
    }

    iMesh_LAST_ERROR.error_type = iBase_SUCCESS;
    RETURN(iMesh_LAST_ERROR.error_type);
  }

  void iMesh_setEntSetData (iMesh_Instance instance,
                            /*in*/ iBase_EntitySetHandle entity_set_handle,
                            /*in*/ const iBase_TagHandle tag_handle,
                            /*in*/ const char* tag_value,
                            /*in*/ const int , int *err) 
  {
    MBErrorCode result;

    if (entity_set_handle == 0)
        // set the tag data on this entity set
      result = MBI->tag_set_data(TAG_HANDLE(tag_handle),
                                 NULL, 0, tag_value);
    else {
      MBEntityHandle set = ENTITY_HANDLE(entity_set_handle);
      result = MBI->tag_set_data(TAG_HANDLE(tag_handle), &set, 1, tag_value);
    }
    
    if (MB_SUCCESS != result)
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_setEntSetData: error");

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_setEntSetIntData (iMesh_Instance instance,
                               /*in*/ iBase_EntitySetHandle entity_set,
                               /*in*/ const iBase_TagHandle tag_handle,
                               /*in*/ const int tag_value, int *err) 
  {
    iMesh_setEntSetData(instance, entity_set, tag_handle, 
                        reinterpret_cast<const char*>(&tag_value), 
                        sizeof(int), err);
  }

  void iMesh_setEntSetDblData (iMesh_Instance instance,
                               /*in*/ iBase_EntitySetHandle entity_set,
                               /*in*/ const iBase_TagHandle tag_handle,
                               /*in*/ const double tag_value, int *err) 
  {
    iMesh_setEntSetData(instance, entity_set, tag_handle, 
                        reinterpret_cast<const char*>(&tag_value),
                        sizeof(double), err);
  }

  void iMesh_setEntSetBoolData (iMesh_Instance instance,
                                /*in*/ iBase_EntitySetHandle entity_set,
                                /*in*/ const iBase_TagHandle tag_handle,
                                /*in*/ const bool tag_value, int *err) 
  {
    iMesh_setEntSetData(instance, entity_set, tag_handle, 
                        reinterpret_cast<const char*>(&tag_value), 
                        sizeof(bool), err);
  }

  void iMesh_setEntSetEHData (iMesh_Instance instance,
                              /*in*/ iBase_EntitySetHandle entity_set,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*in*/ const iBase_EntityHandle tag_value, int *err) 
  {
    iMesh_setEntSetData(instance, entity_set, tag_handle, 
                        reinterpret_cast<const char*>(&tag_value), 
                        sizeof(iBase_EntityHandle), err);
  }

  void iMesh_getEntSetData (iMesh_Instance instance,
                            /*in*/ const iBase_EntitySetHandle entity_set_handle,
                            /*in*/ const iBase_TagHandle tag_handle,
                            /*inout*/ char** tag_value,
                            /*inout*/ int* tag_value_allocated,
                            /*inout*/ int* tag_value_size, int *err) 
  {
    MBEntityHandle eh = ENTITY_HANDLE(entity_set_handle);
    MBTag tag = TAG_HANDLE(tag_handle);

    int tag_size;
    MBErrorCode result = MBI->tag_get_size(tag, tag_size);
    if (MB_SUCCESS != result) {
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_getEntSetData: couldn't get tag size.");
    }

    else {
      TAG_CHECK_SIZE(*tag_value, *tag_value_allocated, tag_size);

      if (eh == 0)
        result = MBI->tag_get_data(tag, NULL, 0, *tag_value);
      else
        result = MBI->tag_get_data(tag, &eh, 1, *tag_value);

      if (MB_SUCCESS != result) {
        iMesh_processError(iBase_ERROR_MAP[result], "iMesh_getEntSetData didn't succeed.");
      }
      else
        *tag_value_size = tag_size;
    }
  
    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_getEntSetIntData (iMesh_Instance instance,
                               /*in*/ const iBase_EntitySetHandle entity_set,
                               /*in*/ const iBase_TagHandle tag_handle,
                               int *out_data, int *err) 
  {
    char *tag_ptr = reinterpret_cast<char*>(out_data);
    int dum_size = sizeof(int);
    iMesh_getEntSetData(instance, entity_set, tag_handle, &tag_ptr, 
                        &dum_size, &dum_size, err);
  }

  void iMesh_getEntSetDblData (iMesh_Instance instance,
                               /*in*/ const iBase_EntitySetHandle entity_set,
                               /*in*/ const iBase_TagHandle tag_handle,
                               double *out_data, int *err) 
  {
    char *tag_ptr = reinterpret_cast<char*>(out_data);
    int tag_size = sizeof(double);
    iMesh_getEntSetData(instance, entity_set, tag_handle, &tag_ptr, 
                        &tag_size, &tag_size, err);
  }

  void iMesh_getEntSetBoolData (iMesh_Instance instance,
                                /*in*/ const iBase_EntitySetHandle entity_set,
                                /*in*/ const iBase_TagHandle tag_handle,
                                int *out_data, int *err) 
  {
    char *tag_ptr = reinterpret_cast<char*>(out_data);
    int tag_size = sizeof(bool);
    iMesh_getEntSetData(instance, entity_set, tag_handle, &tag_ptr, 
                        &tag_size, &tag_size, err);
  }

  void iMesh_getEntSetEHData (iMesh_Instance instance,
                              /*in*/ const iBase_EntitySetHandle entity_set,
                              /*in*/ const iBase_TagHandle tag_handle,
                              iBase_EntityHandle *out_data, int *err) 
  {
    char *tag_ptr = reinterpret_cast<char*>(out_data);
    int tag_size = sizeof(MBEntityHandle);
    iMesh_getEntSetData(instance, entity_set, tag_handle, &tag_ptr, 
                        &tag_size, &tag_size, err);
  }

  void iMesh_getAllEntSetTags (iMesh_Instance instance,
                               /*in*/ const iBase_EntitySetHandle entity_set_handle,
                               /*out*/ iBase_TagHandle** tag_handles,
                               /*out*/ int* tag_handles_allocated,
                               /*out*/ int* tag_handles_size, int *err) 
  {
    MBEntityHandle eh = ENTITY_HANDLE(entity_set_handle);
    std::vector<MBTag> all_tags;
  
    MBErrorCode result = MBI->tag_get_tags_on_entity(eh, all_tags);
    if (MB_SUCCESS != result)
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_entitysetGetAllTagHandles failed.");
 
      // now put those tag handles into sidl array
    CHECK_SIZE(*tag_handles, *tag_handles_allocated, 
               (int)all_tags.size(), iBase_TagHandle, iBase_MEMORY_ALLOCATION_FAILED);
    memcpy(*tag_handles, &all_tags[0], all_tags.size()*sizeof(MBTag));

    if (MB_SUCCESS == result)
      *tag_handles_size = (int) all_tags.size();

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_rmvEntSetTag (iMesh_Instance instance,
                           /*in*/ iBase_EntitySetHandle entity_set_handle,
                           /*in*/ const iBase_TagHandle tag_handle, int *err) 
  {
    if (0 == entity_set_handle) {
      int success;
      iMesh_getRootSet(instance, &entity_set_handle, &success);
      if (iBase_ERROR_MAP[success] != iBase_SUCCESS) RETURN(iBase_ERROR_MAP[success]);
    }
    MBEntityHandle set = ENTITY_HANDLE(entity_set_handle);
    MBErrorCode result = MBI->tag_delete_data(TAG_HANDLE(tag_handle), &set, 1);
  
      // don't check return; this tag may have never been set on the entity set
    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_setVtxCoord (iMesh_Instance instance,
                          /*in*/ iBase_EntityHandle vertex_handle,
                          /*in*/ const double x, /*in*/ const double y, 
                          /*in*/ const double z, int *err)
                    
  {
    const double xyz[3] = {x, y, z};
  
    iMesh_setVtxArrCoords(instance, &vertex_handle, 1, iBase_BLOCKED,
                          xyz, 3, err);
  }

  void iMesh_createVtx(iMesh_Instance instance,
                       /*in*/ const double x, /*in*/ const double y, 
                       /*in*/ const double z,
                       /*out*/ iBase_EntityHandle* new_vertex_handle, int *err) 
  {
    int dum = 1;
    const double xyz[3] = {x, y, z};
    iMesh_createVtxArr(instance, 1, iBase_BLOCKED,
                       xyz, 3, &new_vertex_handle, &dum, &dum, err);
  }
                                                   
  void iMesh_createEnt(iMesh_Instance instance,
                       /*in*/ const int new_entity_topology,
                       /*in*/ const iBase_EntityHandle* lower_order_entity_handles,
                       /*in*/ const int lower_order_entity_handles_size,
                       /*out*/ iBase_EntityHandle* new_entity_handle,
                       /*out*/ int* status, int *err) 
  {
    if (0 == lower_order_entity_handles_size) {
      iMesh_processError(iBase_INVALID_ENTITY_COUNT, 
                         "iMesh_createEnt: need more than zero lower order entities.");
      RETURN(iBase_INVALID_ENTITY_COUNT);
    }

      // call MB directly to allow creation of higher-order entities
      // directly from connectivity
    MBEntityType this_type = mb_topology_table[new_entity_topology];
    MBEntityHandle tmp_ent;
    MBErrorCode result = MBI->create_element(this_type,
                                             CONST_HANDLE_ARRAY_PTR(lower_order_entity_handles), 
                                             lower_order_entity_handles_size, 
                                             tmp_ent);
    if (MB_SUCCESS != result)
      *status = iBase_CREATION_FAILED;
    else
      *status = iBase_SUCCESS;
    *new_entity_handle = reinterpret_cast<iBase_EntityHandle>(tmp_ent);

    *err = *status;
  }

  void iMesh_deleteEnt(iMesh_Instance instance,
                       /*in*/ iBase_EntityHandle entity_handle, int *err) 
  {
    iMesh_deleteEntArr(instance, &entity_handle, 1, err);
  }
                                                 
  void iMesh_getArrData (iMesh_Instance instance,
                         /*in*/ const iBase_EntityHandle* entity_handles,
                         /*in*/ const int entity_handles_size,
                         /*in*/ const iBase_TagHandle tag_handle,
                         /*inout*/ char** tag_values,
                         /*inout*/int* tag_values_allocated,
                         /*out*/ int* tag_values_size, int *err) 
  {
    const MBEntityHandle *ents = reinterpret_cast<const MBEntityHandle *>(entity_handles);
    MBTag tag = TAG_HANDLE(tag_handle);
    int tag_size;
    MBErrorCode result = MBI->tag_get_size(tag, tag_size);
    if (MB_SUCCESS != result) {
      int nerr=-1; char tagn[64], msg[256];
      iMesh_getTagName(instance, tag_handle, tagn, &nerr, sizeof(tagn));
      snprintf(msg, sizeof(msg), "iMesh_getArrData: couldn't get size for tag \"%s\"",
               nerr==0?tagn:"unknown");
      iMesh_processError(iBase_ERROR_MAP[result], msg);
      RETURN(iBase_ERROR_MAP[result]);
    }

    if (0 == entity_handles_size) {
      RETURN(iBase_ERROR_MAP[MB_SUCCESS]);
    }
  
    TAG_CHECK_SIZE(*tag_values, *tag_values_allocated, 
                   tag_size * entity_handles_size);

    result = MBI->tag_get_data(tag, ents, entity_handles_size,
                               *tag_values);

    if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result) {
      int nerr=-1; char tagn[64], msg[256];
      iMesh_getTagName(instance, tag_handle, tagn, &nerr, sizeof(tagn));
      snprintf(msg, sizeof(msg), "iMesh_getArrData: didn't succeed for tag \"%s\"",
               nerr==0?tagn:"unknown");
      iMesh_processError(iBase_ERROR_MAP[result], msg);
    }
    else if (MB_TAG_NOT_FOUND == result) {
      int nerr=-1; char tagn[64], msg[256];
      iMesh_getTagName(instance, tag_handle, tagn, &nerr, sizeof(tagn));
      snprintf(msg, sizeof(msg), "iMesh_getArrData: tag \"%s\" not found",
               nerr==0?tagn:"unknown");
      iMesh_processError(iBase_ERROR_MAP[result], msg);
    }

    if (MB_SUCCESS == result)
      *tag_values_size = tag_size * entity_handles_size;
  
    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_getIntArrData (iMesh_Instance instance,
                            /*in*/ const iBase_EntityHandle* entity_handles,
                            /*in*/ const int entity_handles_size,
                            /*in*/ const iBase_TagHandle tag_handle,
                            /*inout*/ int** tag_values,
                            /*inout*/ int* tag_values_allocated,
                            /*out*/ int* tag_values_size, int *err) 
  {
    *tag_values_allocated *= sizeof(int);
    *tag_values_size *= sizeof(int);
    iMesh_getArrData(instance, entity_handles, 
                     entity_handles_size, tag_handle,
                     reinterpret_cast<char**>(tag_values), 
                     tag_values_allocated, 
                     tag_values_size, err);
    *tag_values_allocated /= sizeof(int);
    *tag_values_size /= sizeof(int);
  }

  void iMesh_getDblArrData (iMesh_Instance instance,
                            /*in*/ const iBase_EntityHandle* entity_handles,
                            /*in*/ const int entity_handles_size,
                            /*in*/ const iBase_TagHandle tag_handle,
                            /*inout*/ double** tag_values,
                            /*inout*/ int* tag_values_allocated,
                            /*out*/ int* tag_values_size, int *err) 
  {
    *tag_values_allocated *= sizeof(double);
    *tag_values_size *= sizeof(double);
    iMesh_getArrData(instance, entity_handles, 
                     entity_handles_size, tag_handle,
                     reinterpret_cast<char**>(tag_values), 
                     tag_values_allocated, tag_values_size, err);
    *tag_values_allocated /= sizeof(double);
    *tag_values_size /= sizeof(double);
  }

  void iMesh_getBoolArrData (iMesh_Instance instance,
                             /*in*/ const iBase_EntityHandle* entity_handles,
                             /*in*/ const int entity_handles_size,
                             /*in*/ const iBase_TagHandle tag_handle,
                             /*inout*/ bool** tag_value,
                             /*inout*/ int* tag_value_allocated,
                             /*out*/ int* tag_value_size, int *err) 
  {
    *tag_value_allocated *= sizeof(bool);
    *tag_value_size *= sizeof(bool);
    iMesh_getArrData(instance, entity_handles, 
                     entity_handles_size, tag_handle,
                     reinterpret_cast<char**>(tag_value), 
                     tag_value_allocated, tag_value_size, err);
    *tag_value_allocated /= sizeof(bool);
    *tag_value_size /= sizeof(bool);
  }

  void iMesh_getEHArrData (iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*inout*/ iBase_EntityHandle** tag_value,
                           /*inout*/ int* tag_value_allocated,
                           /*out*/ int* tag_value_size, int *err) 
  {
    *tag_value_allocated *= sizeof(iBase_EntityHandle);
    *tag_value_size *= sizeof(iBase_EntityHandle);
    iMesh_getArrData(instance, entity_handles, 
                     entity_handles_size, tag_handle,
                     reinterpret_cast<char**>(tag_value), 
                     tag_value_allocated, 
                     tag_value_size, err);
    *tag_value_allocated /= sizeof(iBase_EntityHandle);
    *tag_value_size /= sizeof(iBase_EntityHandle);
  }

  void iMesh_setArrData (iMesh_Instance instance,
                         /*in*/ const iBase_EntityHandle* entity_handles,
                         /*in*/ const int entity_handles_size,
                         /*in*/ const iBase_TagHandle tag_handle,
                         /*in*/ const char* tag_values,
                         /*in*/ const int tag_values_size, int *err) 
  {
    if (0 == entity_handles_size) {
      RETURN(iBase_ERROR_MAP[MB_SUCCESS]);
    }

    MBErrorCode result = MBI->tag_set_data(TAG_HANDLE(tag_handle), 
                                           CONST_HANDLE_ARRAY_PTR(entity_handles),
                                           entity_handles_size,
                                           tag_values);
    if (MB_SUCCESS != result) {
      std::string msg("iMesh_setArrData didn't succeed, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    }
  
    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_setIntArrData (iMesh_Instance instance,
                            /*in*/ const iBase_EntityHandle* entity_handles,
                            /*in*/ const int entity_handles_size,
                            /*in*/ const iBase_TagHandle tag_handle,
                            /*in*/ const int* tag_values,
                            /*in*/ const int tag_values_size, int *err) 
  {
    iMesh_setArrData(instance, entity_handles, 
                     entity_handles_size, tag_handle, 
                     reinterpret_cast<const char*>(tag_values), 
                     sizeof(int)*tag_values_size, err);
  }

  void iMesh_setDblArrData (iMesh_Instance instance,
                            /*in*/ const iBase_EntityHandle* entity_handles,
                            /*in*/ const int entity_handles_size,
                            /*in*/ const iBase_TagHandle tag_handle,
                            /*in*/ const double* tag_values,
                            /*in*/ const int tag_values_size, int *err) 
  {
    iMesh_setArrData(instance, entity_handles, 
                     entity_handles_size, tag_handle, 
                     reinterpret_cast<const char*>(tag_values), 
                     sizeof(double)*tag_values_size, err);
  }

  void iMesh_setBoolArrData (iMesh_Instance instance,
                             /*in*/ const iBase_EntityHandle* entity_handles,
                             /*in*/ const int entity_handles_size,
                             /*in*/ const iBase_TagHandle tag_handle,
                             /*in*/ const bool* tag_values,
                             /*in*/ const int tag_values_size, int *err) 
  {
    iMesh_setArrData(instance, entity_handles, 
                     entity_handles_size, tag_handle, 
                     reinterpret_cast<const char*>(tag_values), 
                     sizeof(bool)*tag_values_size, err);
  }

  void iMesh_setEHArrData (iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*in*/ const iBase_EntityHandle* tag_values,
                           /*in*/ const int tag_values_size, int *err) 
  {
    iMesh_setArrData(instance, entity_handles, 
                     entity_handles_size, tag_handle, 
                     reinterpret_cast<const char*>(tag_values), 
                     sizeof(iBase_EntityHandle)*tag_values_size, err);
  }

  void iMesh_rmvArrTag (iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle* entity_handles,
                        /*in*/ const int entity_handles_size,
                        /*in*/ const iBase_TagHandle tag_handle, int *err) 
  {
    MBErrorCode result = MBI->tag_delete_data(TAG_HANDLE(tag_handle),
                                              CONST_HANDLE_ARRAY_PTR(entity_handles),
                                              entity_handles_size);
  
      // don't check return; this tag may have never been set on the entity
    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_getData (iMesh_Instance instance,
                      /*in*/ const iBase_EntityHandle entity_handle,
                      /*in*/ const iBase_TagHandle tag_handle,
                      /*out*/ char** tag_value,
                      /*inout*/ int *tag_value_allocated,
                      /*out*/ int *tag_value_size, int *err) 
  {
    iMesh_getArrData(instance, &entity_handle, 1,
                     tag_handle, tag_value, tag_value_allocated,
                     tag_value_size, err);
  }

  void iMesh_getIntData (iMesh_Instance instance,
                         /*in*/ const iBase_EntityHandle entity_handle,
                         /*in*/ const iBase_TagHandle tag_handle,
                         int *out_data, int *err) 
  {
    char *val_ptr = reinterpret_cast<char*>(out_data);
    int val_size = sizeof(int);
    iMesh_getArrData(instance, &entity_handle, 1,
                     tag_handle, &val_ptr, &val_size, &val_size, err);
  }

  void iMesh_getDblData (iMesh_Instance instance,
                         /*in*/ const iBase_EntityHandle entity_handle,
                         /*in*/ const iBase_TagHandle tag_handle,
                         double *out_data, int *err) 
  {
    char *val_ptr = reinterpret_cast<char*>(out_data);
    int val_size = sizeof(double);
    iMesh_getArrData(instance, &entity_handle, 1,
                     tag_handle, &val_ptr, &val_size, &val_size, err);
  }

  void iMesh_getBoolData (iMesh_Instance instance,
                          /*in*/ const iBase_EntityHandle entity_handle,
                          /*in*/ const iBase_TagHandle tag_handle,
                          int *out_data, int *err) 
  {
    char *val_ptr = reinterpret_cast<char*>(out_data);
      // make the data size a full word, because of sidl needing at least a full word
    int val_size = sizeof(int);
    iMesh_getArrData(instance, &entity_handle, 1,
                     tag_handle, &val_ptr, &val_size, &val_size, err);
  }

  void iMesh_getEHData (iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*in*/ const iBase_TagHandle tag_handle,
                        iBase_EntityHandle *out_data, int *err) 
  {
    char *val_ptr = reinterpret_cast<char*>(out_data);
    int dum = sizeof(iBase_EntityHandle);
    iMesh_getArrData(instance, &entity_handle, 1,
                     tag_handle, &val_ptr, &dum, &dum, err);
  }

  void iMesh_setData (iMesh_Instance instance,
                      /*in*/ iBase_EntityHandle entity_handle,
                      /*in*/ const iBase_TagHandle tag_handle,
                      /*in*/ const char* tag_value,
                      /*in*/ const int tag_value_size, int *err) 
  {
    iMesh_setArrData(instance, &entity_handle, 1,
                     tag_handle, tag_value, tag_value_size, err);
  }

  void iMesh_setIntData (iMesh_Instance instance,
                         /*in*/ iBase_EntityHandle entity_handle,
                         /*in*/ const iBase_TagHandle tag_handle,
                         /*in*/ const int tag_value, int *err) 
  {
    iMesh_setArrData(instance, &entity_handle, 1,
                     tag_handle, 
                     reinterpret_cast<const char*>(&tag_value), 
                     sizeof(int), err);
  }

  void iMesh_setDblData (iMesh_Instance instance,
                   
                         /*in*/ iBase_EntityHandle entity_handle,
                         /*in*/ const iBase_TagHandle tag_handle,
                         /*in*/ const double tag_value, int *err) 
  {
    iMesh_setArrData(instance, &entity_handle, 1,
                     tag_handle, 
                     reinterpret_cast<const char*>(&tag_value), 
                     sizeof(double), err);
  }

  void iMesh_setBoolData (iMesh_Instance instance,
                          /*in*/ iBase_EntityHandle entity_handle,
                          /*in*/ const iBase_TagHandle tag_handle,
                          /*in*/ const bool tag_value, int *err) 
  {
    iMesh_setArrData(instance, &entity_handle, 1,
                     tag_handle, 
                     reinterpret_cast<const char*>(&tag_value), 
                     sizeof(bool), err);
  }

  void iMesh_setEHData (iMesh_Instance instance,
                        /*in*/ iBase_EntityHandle entity_handle,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*in*/ const iBase_EntityHandle tag_value, int *err) 
  {
    iMesh_setArrData(instance, &entity_handle, 1,
                     tag_handle, 
                     reinterpret_cast<const char*>(&tag_value), 
                     sizeof(iBase_EntityHandle), err);
  }

  void iMesh_getAllTags (iMesh_Instance instance,
                         /*in*/ const iBase_EntityHandle entity_handle,
                         /*inout*/ iBase_TagHandle** tag_handles,
                         /*inout*/ int* tag_handles_allocated,
                         /*out*/ int* tag_handles_size, int *err) 
  {
    std::vector<MBTag> all_tags;
  
    MBErrorCode result = MBI->tag_get_tags_on_entity(ENTITY_HANDLE(entity_handle), all_tags);
    if (MB_SUCCESS != result)
      iMesh_processError(iBase_ERROR_MAP[result], "iMesh_getAllTags failed.");
    
      // now put those tag handles into sidl array
    CHECK_SIZE(*tag_handles, *tag_handles_allocated,
               (int)all_tags.size(), iBase_TagHandle, iBase_MEMORY_ALLOCATION_FAILED);
    memcpy(*tag_handles, &all_tags[0], all_tags.size()*sizeof(MBTag));

    if (MB_SUCCESS == result)
      *tag_handles_size = all_tags.size();

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_rmvTag (iMesh_Instance instance,
                     /*in*/ iBase_EntityHandle entity_handle,
                     /*in*/ const iBase_TagHandle tag_handle, int *err) 
  {
    iMesh_rmvArrTag(instance, &entity_handle, 1, tag_handle, err);
  }

  void iMesh_initEntIter (iMesh_Instance instance,
                          /*in*/ const iBase_EntitySetHandle entity_set_handle,
                          /*in*/ const int requested_entity_type,
                          /*in*/ const int requested_entity_topology,
                          /*out*/ iMesh_EntityIterator* entity_iterator,
                          int *err) 
  {
    iMesh_initEntArrIter(instance, entity_set_handle, requested_entity_type,
                         requested_entity_topology, 1, 
                         reinterpret_cast<iMesh_EntityArrIterator*>(entity_iterator),
                         err);
  }

  void iMesh_getNextEntIter (iMesh_Instance instance,
                             /*in*/ iMesh_EntityIterator entity_iterator,
                             /*out*/ iBase_EntityHandle* entity_handle, 
                             int *is_end, int *err) 
  {
    int eh_size = 1;
    iMesh_getNextEntArrIter(instance,
                            reinterpret_cast<iMesh_EntityArrIterator>(entity_iterator),
                            &entity_handle, &eh_size, &eh_size, is_end, err);
  
  }

  void iMesh_resetEntIter (iMesh_Instance instance,
                           /*in*/ iMesh_EntityIterator entity_iterator, int *err) 
  {
    iMesh_resetEntArrIter(instance,
                          reinterpret_cast<iMesh_EntityArrIterator>(entity_iterator),
                          err);  
  }

  void iMesh_endEntIter (iMesh_Instance instance,
                         /*in*/ iMesh_EntityIterator entity_iterator, int *err) 
  {
    iMesh_endEntArrIter(instance, 
                        reinterpret_cast<iMesh_EntityArrIterator>(entity_iterator),
                        err);
  }

  void iMesh_getEntTopo (iMesh_Instance instance,
                         /*in*/ const iBase_EntityHandle entity_handle,
                         int *out_topo, int *err) 
  {
    *out_topo = tstt_topology_table[MBI->type_from_handle(ENTITY_HANDLE(entity_handle))];
    RETURN(iBase_SUCCESS);
  }
  
  void iMesh_getEntType (iMesh_Instance instance,
                         /*in*/ const iBase_EntityHandle entity_handle,
                         int *out_type, int *err) 
  {
    *out_type = tstt_type_table[MBI->type_from_handle(ENTITY_HANDLE(entity_handle))];
    RETURN(iBase_SUCCESS);
  }

  void iMesh_getVtxCoord (iMesh_Instance instance,
                          /*in*/ const iBase_EntityHandle vertex_handle,
                          /*out*/ double *x, /*out*/ double *y, /*out*/ double *z, int *err)
  {
    int order = iBase_BLOCKED;
    double xyz[3], *tmp_xyz = xyz;
    int dum = 3;
  
    iMesh_getVtxArrCoords(instance,
                          &vertex_handle, 1, order,
                          &tmp_xyz, &dum, &dum, err);
    if (iBase_SUCCESS == *err) {
      *x = xyz[0]; *y = xyz[1]; *z = xyz[2];
    }
  }

  void iMesh_getEntAdj(iMesh_Instance instance,
                       /*in*/ const iBase_EntityHandle entity_handle,
                       /*in*/ const int entity_type_requested,
                       /*inout*/ iBase_EntityHandle** adj_entity_handles,
                       /*inout*/ int* adj_entity_handles_allocated,
                       /*out*/ int* adj_entity_handles_size, int *err)
  {
    int offsets[2];
    int *offsets_ptr = offsets;
    int offset_size, offset_allocated = 2;
  
    iMesh_getEntArrAdj(instance,
                       &entity_handle, 1, entity_type_requested,
                       adj_entity_handles, adj_entity_handles_allocated, 
                       adj_entity_handles_size, &offsets_ptr, &offset_allocated, 
                       &offset_size, err);
  }
 
  void iMesh_getEnt2ndAdj( iMesh_Instance instance,
                           iBase_EntityHandle entity_handle,
                           int order_adjacent_key,
                           int requested_entity_type,
                           iBase_EntityHandle** adj_entities,
                           int* adj_entities_allocated,
                           int* adj_entities_size,
                           int* err ) 
  {
    int offsets[2];
    int *offsets_ptr = offsets;
    int offset_size, offset_allocated = 2;
  
    iMesh_getEntArr2ndAdj(instance,
                          &entity_handle, 1, order_adjacent_key,
                          requested_entity_type,
                          adj_entities, adj_entities_allocated, 
                          adj_entities_size, &offsets_ptr, &offset_allocated, 
                          &offset_size, err);
  }

  void iMesh_subtract(iMesh_Instance instance,
                      /*in*/ const iBase_EntitySetHandle entity_set_1,
                      /*in*/ const iBase_EntitySetHandle entity_set_2,
                      /*out*/ iBase_EntitySetHandle* result_entity_set, int *err)
  {
    MBEntityHandle temp_set;
    MBEntityHandle set1 = ENTITY_HANDLE(entity_set_1), 
      set2 = ENTITY_HANDLE(entity_set_2);
    MBErrorCode result = MBI->create_meshset(MESHSET_SET, temp_set);
    if (MB_SUCCESS == result) result = MBI->unite_meshset(temp_set, set1);
    if (MB_SUCCESS == result) result = MBI->subtract_meshset(temp_set, set2);

    if (result != MB_SUCCESS) {
      std::string msg("iMesh_entitysetSubtract:\
           ERROR subtract failed, with error type: ");
      msg += MBI->get_error_string(result);
      MBI->delete_entities(&temp_set, 1);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    }

    *result_entity_set = (iBase_EntitySetHandle)temp_set;

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_intersect(iMesh_Instance instance,
                       /*in*/ const iBase_EntitySetHandle entity_set_1,
                       /*in*/ const iBase_EntitySetHandle entity_set_2,
                       /*out*/ iBase_EntitySetHandle* result_entity_set, int *err)
  {
    MBEntityHandle temp_set;
    MBEntityHandle set1 = ENTITY_HANDLE(entity_set_1), 
      set2 = ENTITY_HANDLE(entity_set_2);
    MBErrorCode result = MBI->create_meshset(MESHSET_SET, temp_set);
    if (MB_SUCCESS == result) result = MBI->unite_meshset(temp_set, set1);
    if (MB_SUCCESS == result) result = MBI->intersect_meshset(temp_set, set2);

    if (result != MB_SUCCESS) {
      std::string msg("iMesh_entitysetIntersect:\
           ERROR subtract failed, with error type: ");
      msg += MBI->get_error_string(result);
      MBI->delete_entities(&temp_set, 1);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    }

    *result_entity_set = (iBase_EntitySetHandle)temp_set;

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_unite(iMesh_Instance instance,
                   /*in*/ const iBase_EntitySetHandle entity_set_1,
                   /*in*/ const iBase_EntitySetHandle entity_set_2,
                   /*out*/ iBase_EntitySetHandle* result_entity_set, int *err)
  {  
    MBEntityHandle temp_set;
    MBEntityHandle set1 = ENTITY_HANDLE(entity_set_1), 
      set2 = ENTITY_HANDLE(entity_set_2);
    MBErrorCode result = MBI->create_meshset(MESHSET_SET, temp_set);
    if (MB_SUCCESS == result) result = MBI->unite_meshset(temp_set, set1);
    if (MB_SUCCESS == result) result = MBI->unite_meshset(temp_set, set2);

    if (result != MB_SUCCESS) {
      std::string msg("iMesh_entitysetIntersect:\
           ERROR subtract failed, with error type: ");
      msg += MBI->get_error_string(result);
      MBI->delete_entities(&temp_set, 1);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    }

    *result_entity_set = (iBase_EntitySetHandle)temp_set;

    RETURN(iBase_ERROR_MAP[result]);
  }

  void iMesh_getEntitiesRec(iMesh_Instance instance,
                            /*in*/ const iBase_EntitySetHandle entity_set_handle,
                            /*in*/ const int entity_type,
                            /*in*/ const int entity_topology,
                            /*in*/ const int recursive,
                            /*out*/ iBase_EntityHandle** entity_handles,
                            /*out*/ int* entity_handles_allocated,
                            /*out*/ int* entity_handles_size,
                            /*out*/ int *err) 
  {
    bool use_top = false;
    bool use_type = false;
      // initialize just to get rid of compiler warning
    MBEntityType type = mb_topology_table[iMesh_ALL_TOPOLOGIES];
    MBRange out_entities;
 
    if (entity_topology >= iMesh_POINT
        && entity_topology < iMesh_ALL_TOPOLOGIES) {
      type = mb_topology_table[entity_topology];
      use_top = true;
    }
    else if (entity_type >= iBase_VERTEX
             && entity_type <= iBase_ALL_TYPES)
      use_type = true;
    else {
      iMesh_processError(iBase_BAD_TYPE_AND_TOPO, 
                         "iMesh_getEntities:ERROR not valid entity type or topology");
      RETURN(iBase_BAD_TYPE_AND_TOPO);
    }

    MBEntityHandle handle = ENTITY_HANDLE(entity_set_handle);
    MBErrorCode result;

    if (use_top) {
      if (entity_topology == iMesh_SEPTAHEDRON)
        result = MB_SUCCESS;  // MOAB doesn't do septahedrons, so there are never any.
      else
        result = MBI->get_entities_by_type(handle, type, out_entities, recursive);
    }
    else if (use_type && entity_type != iBase_ALL_TYPES)
      result = MBI->get_entities_by_dimension(handle, entity_type, out_entities, recursive);
    else 
      result = MBI->get_entities_by_handle(handle, out_entities, recursive);

    if (result != MB_SUCCESS) {
      std::string msg("iMesh_GetEntities:ERROR getting entities, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      RETURN(iBase_ERROR_MAP[result]);
    }

    int num_ents = out_entities.size() - out_entities.num_of_type(MBENTITYSET);
    
    CHECK_SIZE(*entity_handles, *entity_handles_allocated, 
               num_ents, iBase_EntityHandle, iBase_MEMORY_ALLOCATION_FAILED);
  
    MBRange::iterator iter = out_entities.begin();
    MBRange::iterator end_iter = out_entities.end();
    int k = 0;

      // filter out entity sets here
    if (iBase_ALL_TYPES == entity_type && iMesh_ALL_TOPOLOGIES == entity_topology) {
      for (; iter != end_iter && MBI->type_from_handle(*iter) != MBENTITYSET; iter++)
        (*entity_handles)[k++] = (iBase_EntityHandle)*iter;
    }
    else {
      for (; iter != end_iter; iter++)
        (*entity_handles)[k++] = (iBase_EntityHandle)*iter;
    }

      // now it's safe to set the size; set it to k, not out_entities.size(), to
      // account for sets which might have been removed
    *entity_handles_size = k;

    RETURN(iBase_SUCCESS);
  }  

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
                             /*out*/ int *err) 
  {
    *num_type = 0;
    MBErrorCode result;
    if (entity_type == iBase_ALL_TYPES) {
      result = MBI->get_number_entities_by_handle
        (ENTITY_HANDLE(entity_set_handle), *num_type, recursive);
      if (MB_SUCCESS == result && !recursive) {
        int num_sets = 0;
        result = MBI->get_number_entities_by_type
          (ENTITY_HANDLE(entity_set_handle), MBENTITYSET, num_sets, recursive);
        *num_type -= num_sets;
      }
    } else {
      result = MBI->get_number_entities_by_dimension
        (ENTITY_HANDLE(entity_set_handle), entity_type, *num_type, recursive);
    }

    if (MB_SUCCESS != result) {
      std::string msg("iMesh_entitysetGetNumberEntityOfType: ERROR getting number of entities"
                      " by type, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      *num_type = 0;
      RETURN(iBase_ERROR_MAP[result]);
    }

    RETURN(iBase_SUCCESS);
  }


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
                             /*out*/ int *err) 
  {
    if (entity_topology == iMesh_SEPTAHEDRON) {
      *num_topo = 0;
      RETURN(iBase_SUCCESS);
    }

    *num_topo = 0;
    MBErrorCode result;
    if (iMesh_ALL_TOPOLOGIES == entity_topology) {
      result = MBI->get_number_entities_by_handle(ENTITY_HANDLE(entity_set_handle),
                                                  *num_topo, recursive);
                                       
      if (!recursive && MB_SUCCESS == result) { // remove entity sets from count
        int num_sets;
        result = MBI->get_number_entities_by_type
          (ENTITY_HANDLE(entity_set_handle), MBENTITYSET, num_sets, recursive);
        *num_topo -= num_sets;
      }
    }
    else {
      result = MBI->get_number_entities_by_type(ENTITY_HANDLE(entity_set_handle),
                                         mb_topology_table[entity_topology], 
                                         *num_topo, recursive);
    }

    if (MB_SUCCESS != result) {
      std::string msg("iMesh_entitysetGetNumberEntityOfTopology: ERROR getting "
                      "number of entities by topology., with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      *num_topo = -1;
      RETURN(iBase_ERROR_MAP[result]);
    }
    else {
      RETURN(iBase_SUCCESS);
    }
  }

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
                              /*out*/ int *err)
  {
    bool use_top = false;
    bool use_type = false;
      // initialize just to get rid of compiler warning
    MBEntityType type = mb_topology_table[iMesh_ALL_TOPOLOGIES];
    MBRange out_entities;
 
    if (entity_topology >= iMesh_POINT
        && entity_topology < iMesh_ALL_TOPOLOGIES) {
      type = mb_topology_table[entity_topology];
      use_top = true;
    }
    else if (entity_type >= iBase_VERTEX
             && entity_type <= iBase_ALL_TYPES)
      use_type = true;
    else {
      iMesh_processError(iBase_BAD_TYPE_AND_TOPO, 
                         "iMesh_getEntities:ERROR not valid entity type or topology");
      RETURN(iBase_BAD_TYPE_AND_TOPO);
    }

    MBEntityHandle handle = ENTITY_HANDLE(entity_set_handle);
    MBErrorCode result = MB_SUCCESS;

    if (use_top) {
      if (entity_topology == iMesh_SEPTAHEDRON)
        result = MB_SUCCESS;  // MOAB doesn't do septahedrons, so there are never any.
      else
        result = MBI->get_entities_by_type_and_tag(handle, type, (MBTag*)tag_handles, 
                                                   (const void* const *)tag_vals,
                                                   num_tags_vals, out_entities, 
                                                   MBInterface::INTERSECT, recursive);
    }
    else if (use_type && entity_type != iBase_ALL_TYPES) {
        // need to loop over all types of this dimension
      for (MBEntityType tp = MBCN::TypeDimensionMap[entity_type].first;
           tp <= MBCN::TypeDimensionMap[entity_type].second; tp++) {
        MBRange tmp_range;
        MBErrorCode tmp_result = MBI->get_entities_by_type_and_tag(handle, type, (MBTag*)tag_handles, 
                                                                   (const void* const *)tag_vals,
                                                                   num_tags_vals, tmp_range, 
                                                                   MBInterface::INTERSECT, recursive);
        if (MB_SUCCESS != tmp_result) result = tmp_result;
        else out_entities.merge(tmp_range);
      }
    }
    else 
      result = MBI->get_entities_by_type_and_tag(handle, type, (MBTag*)tag_handles, 
                                                 (const void* const *)tag_vals,
                                                 num_tags_vals, out_entities,
                                                 MBInterface::INTERSECT, recursive);

    if (result != MB_SUCCESS) {
      std::string msg("iMesh_GetEntities:ERROR getting entities, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      RETURN(iBase_ERROR_MAP[result]);
    }

    CHECK_SIZE(*entity_handles, *entity_handles_allocated, 
               (int)out_entities.size(), iBase_EntityHandle, iBase_MEMORY_ALLOCATION_FAILED);
  
    MBRange::iterator iter = out_entities.begin();
    MBRange::iterator end_iter = out_entities.end();
    int k = 0;

      // filter out entity sets here
    if (iBase_ALL_TYPES == entity_type && iMesh_ALL_TOPOLOGIES == entity_topology) {
      for (; iter != end_iter && MBI->type_from_handle(*iter) != MBENTITYSET; iter++)
        (*entity_handles)[k++] = (iBase_EntityHandle)*iter;
    }
    else {
      for (; iter != end_iter; iter++)
        (*entity_handles)[k++] = (iBase_EntityHandle)*iter;
    }

      // now it's safe to set the size; set it to k, not out_entities.size(), to
      // account for sets which might have been removed
    *entity_handles_size = k;

    RETURN(iBase_SUCCESS);
  }

  void iMesh_getEntSetsByTagsRec(iMesh_Instance instance,
                                 /*in*/ const iBase_EntitySetHandle entity_set_handle,
                                 /*in*/ const iBase_TagHandle *tag_handles,
                                 /*in*/ const char * const *tag_vals,
                                 /*in*/ const int num_tags_vals,
                                 /*in*/ const int recursive,
                                 /*out*/ iBase_EntitySetHandle** set_handles,
                                 /*out*/ int* set_handles_allocated,
                                 /*out*/ int* set_handles_size,
                                 /*out*/ int *err)
  {
    MBRange out_entities;
 
    MBEntityHandle handle = ENTITY_HANDLE(entity_set_handle);
    MBErrorCode result;

    result = MBI->get_entities_by_type_and_tag(handle, MBENTITYSET, (MBTag*)tag_handles, 
                                               (const void* const *)tag_vals,
                                               num_tags_vals, out_entities, 
                                               MBInterface::INTERSECT, recursive);
    if (result != MB_SUCCESS) {
      std::string msg("ERROR getting entities, with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      RETURN(iBase_ERROR_MAP[result]);
    }

    CHECK_SIZE(*set_handles, *set_handles_allocated, 
               (int)out_entities.size(), iBase_EntitySetHandle, iBase_MEMORY_ALLOCATION_FAILED);

    std::copy(out_entities.begin(), out_entities.end(), ((MBEntityHandle*) *set_handles));
  
      // now it's safe to set the size; set it to k, not out_entities.size(), to
      // account for sets which might have been removed
    *set_handles_size = out_entities.size();

    RETURN(iBase_SUCCESS);
  }

  void iMesh_MBCNType(/*in*/ const int imesh_entity_topology,
                      /*out*/ int *mbcn_type) 
  {
    if (iMesh_POINT > imesh_entity_topology ||
        iMesh_ALL_TOPOLOGIES <= imesh_entity_topology)
      *mbcn_type = -1;
    else
      *mbcn_type = mb_topology_table[imesh_entity_topology];
  }
    
#ifdef __cplusplus
} // extern "C"
#endif

MBErrorCode iMesh_tag_set_vertices(iMesh_Instance instance,
                                   MBEntityHandle in_set, 
                                   const int req_dimension, 
                                   const MBEntityType req_type,
                                   MBTag &tag, MBRange &req_entities, 
                                   int &num_verts, int *err) 
{
    // get all the entities then vertices
  MBRange vertices, entities;
  entities.clear();
  MBErrorCode result = MBI->get_entities_by_handle(in_set, entities, false);
  if (MB_SUCCESS != result) {
    std::string msg("MBMesh::tag_set_vertices: getting entities didn't succeed., with error type: ");
    msg += MBI->get_error_string(result);
    iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    *err = iBase_ERROR_MAP[result];
    return MB_FAILURE;
  }

    // remove any sets
  entities.erase(entities.lower_bound(MBENTITYSET),
                 entities.upper_bound(MBENTITYSET));

    // get vertices
  result = MBI->get_adjacencies(entities, 0, false, vertices,
                                MBInterface::UNION);
  if (MB_SUCCESS != result) {
    std::string msg("MBMesh::tag_set_vertices: getting vertices didn't succeed., with error type: ");
    msg += MBI->get_error_string(result);
    iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
    *err = iBase_ERROR_MAP[result];
    return result;
  }
  
    // tag each vertex with its index in this list
  int i = 0;
  if (tag == 0) {
    result = MBI->tag_create("__position_tag", 4, MB_TAG_DENSE, tag, &i);
    if (0 == tag) {
      std::string msg("MBMesh::tag_set_vertices: couldn't make tag., with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      *err = iBase_ERROR_MAP[result];
      return result;
    }
  }
  
  MBRange::iterator vit;
  for (vit = vertices.begin(), i = 0; vit != vertices.end(); vit++, i++) {
    result = MBI->tag_set_data(tag, &(*vit), 1, &i);
    if (MB_SUCCESS != result) {
      std::string msg("MBMesh::tag_set_vertices: couldn't set pos_tag., with error type: ");
      msg += MBI->get_error_string(result);
      iMesh_processError(iBase_ERROR_MAP[result], msg.c_str());
      *err = iBase_ERROR_MAP[result];
      return result;
    }
  }

    // winnow the list for entities of the desired type and/or topology
  num_verts = 0;
  MBRange::iterator ent_it;
  MBEntityType this_type;
  std::vector<MBEntityHandle> connect_v;
  const MBEntityHandle *connect;
  int num_connect;
  for (ent_it = entities.begin(); ent_it != entities.end(); ent_it++) {
    this_type = MBI->type_from_handle(*ent_it);
    if ((-1 == req_dimension || req_dimension == MBCN::Dimension(this_type)) &&
        (MBMAXTYPE == req_type || this_type == req_type)) {
      req_entities.insert(*ent_it);
      result = MBI->get_connectivity(*ent_it, connect, num_connect, false, 
                                     &connect_v);
      num_verts += num_connect;
    }
  }

  return result;
}

MBErrorCode create_int_ents(MBInterface *instance,
                            MBRange &from_ents,
                            MBEntityHandle in_set) 
{
  MBiMesh* mbimesh = dynamic_cast<MBiMesh*>(instance);
  assert(mbimesh->AdjTable[10] || mbimesh->AdjTable[5]);
  MBRange int_ents;
  MBErrorCode result;
  
  if (mbimesh->AdjTable[10]) {
    result = instance->get_adjacencies(from_ents, 2, true, int_ents,
                                  MBInterface::UNION);
    if (MB_SUCCESS != result) return result;
    unsigned int old_size = from_ents.size();
    from_ents.merge(int_ents);
    if (old_size != from_ents.size() && in_set) {
      result = instance->add_entities(in_set, int_ents);
      if (MB_SUCCESS != result) return result;
    }
  }
  
  if (mbimesh->AdjTable[5]) {
    int_ents.clear();
    result = instance->get_adjacencies(from_ents, 1, true, int_ents,
                                  MBInterface::UNION);
    if (MB_SUCCESS != result) return result;
    unsigned int old_size = from_ents.size();
    from_ents.merge(int_ents);
    if (old_size != from_ents.size() && in_set) {
      result = instance->add_entities(in_set, int_ents);
      if (MB_SUCCESS != result) return result;
    }
  }
    
  return MB_SUCCESS;
}

void eatwhitespace(std::string &this_string) 
{
  std::string::size_type len = this_string.find_last_not_of(" ");
  if (len != this_string.npos) this_string[len+1] = '\0';
}
  
void cfunc_(int arg3, char *mystr, char *mystr2, int arg2, 
            int strsz, int strsz2) 
{
  char tmpstr1[121], tmpstr2[121];
  strncpy(tmpstr1, mystr, strsz);
  tmpstr1[strsz] = '\0';
  strncpy(tmpstr2, mystr2, strsz2);
  tmpstr2[strsz2] = '\0';
  
  std::cout << "String1: " << tmpstr1 << ", string2: " << tmpstr2 << ", arg2 = " << arg2 
            << ", arg3 = " << arg3 << std::endl;
  return;
}

void cfptr_(void **instance) 
{
  *instance = malloc(sizeof(int));
  int *tmp_inst = (int*) instance;
  *tmp_inst = 6;
}

void cfptr2_(void *instance) 
{
  std::cout << "Instance ptr = " << (size_t) instance << std::endl;
}
