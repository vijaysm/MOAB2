#ifndef __iBase_INC__
#define __iBase_INC__

#ifdef __cplusplus

extern "C" 
{
#endif

  typedef int iBase_EntityHandle;
  typedef int iBase_EntitySetHandle;
  typedef int iBase_TagHandle;

  enum iBase_EntityType {       
    iBase_VERTEX, 
    iBase_EDGE, 
    iBase_FACE, 
    iBase_REGION,
    iBase_ALL_TYPES
  };

  enum iBase_CreationStatus {	
    iBase_NEW,               
    iBase_ALREADY_EXISTED,   
    iBase_CREATED_DUPLICATE,
    iBase_CREATION_FAILED
  };

  enum iBase_ErrorActions {
    iBase_SILENT, 
    iBase_WARN_ONLY, 
    iBase_ABORT_ON_ERROR, 
    iBase_PRINT_AND_THROW_ERROR, 
    iBase_THROW_ERROR
  };

   enum iBase_ErrorType {
      iBase_SUCCESS,	
      iBase_MESH_ALREADY_LOADED,	
      iBase_NO_MESH_DATA,	
      iBase_FILE_NOT_FOUND,	
      iBase_FILE_WRITE_ERROR,	
      iBase_NIL_ARRAY,	
      iBase_BAD_ARRAY_SIZE,	
      iBase_BAD_ARRAY_DIMENSION,	
      iBase_INVALID_ENTITY_HANDLE,	
      iBase_INVALID_ENTITY_COUNT,	
      iBase_INVALID_ENTITY_TYPE,	
      iBase_INVALID_ENTITY_TOPOLOGY,
      iBase_BAD_TYPE_AND_TOPO,	
      iBase_ENTITY_CREATION_ERROR,
      iBase_INVALID_TAG_HANDLE,
      iBase_TAG_NOT_FOUND,	
      iBase_TAG_ALREADY_EXISTS,
      iBase_TAG_IN_USE,
      iBase_INVALID_ENTITYSET_HANDLE,	
      iBase_INVALID_ITERATOR_HANDLE,	
      iBase_INVALID_ARGUMENT,	
      iBase_MEMORY_ALLOCATION_FAILED,
      iBase_NOT_SUPPORTED,	
      iBase_FAILURE	
   }; 

  struct iBase_Error 
  {
    char description[120];
    enum iBase_ErrorType error_type;
  };

  enum iBase_StorageOrder {
    iBase_BLOCKED,
    iBase_INTERLEAVED,
    iBase_UNDETERMINED
  };

  enum iBase_TagValueType {
    iBase_INTEGER,
    iBase_DOUBLE,
    iBase_ENTITY_HANDLE,
    iBase_BYTES
  };
#ifdef __cplusplus
}
#endif

#endif
