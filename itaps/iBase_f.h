#ifndef IBASE_F_H
#define IBASE_F_H

#define iBase_EntityHandle integer
#define iBase_EntitySetHandle integer
#define iBase_TagHandle integer

#endif

      integer iBase_VERTEX
      integer iBase_EDGE
      integer iBase_FACE
      integer iBase_REGION
      integer iBase_ALL_TYPES

      parameter (iBase_VERTEX = 0) 
      parameter (iBase_EDGE = 1) 
      parameter (iBase_FACE = 2) 
      parameter (iBase_REGION = 3)
      parameter (iBase_ALL_TYPES = 4)



      integer iBase_NEW
      integer iBase_ALREADY_EXISTED
      integer iBase_CREATED_DUPLICATE
      integer iBase_CREATION_FAILED

      parameter (iBase_NEW = 0)               
      parameter (iBase_ALREADY_EXISTED = 1)   
      parameter (iBase_CREATED_DUPLICATE = 2)
      parameter (iBase_CREATION_FAILED = 3)


      integer iBase_SILENT
      integer iBase_WARN_ONLY
      integer iBase_THROW_ERROR

      parameter (iBase_SILENT = 0) 
      parameter (iBase_WARN_ONLY = 1) 
      parameter (iBase_THROW_ERROR = 4)


      integer iBase_SUCCESS
      integer iBase_MESH_ALREADY_LOADED
      integer iBase_FILE_NOT_FOUND
      integer iBase_FILE_WRITE_ERROR
      integer iBase_NIL_ARRAY
      integer iBase_BAD_ARRAY_SIZE
      integer iBase_BAD_ARRAY_DIMENSION
      integer iBase_INVALID_ENTITY_HANDLE
      integer iBase_INVALID_ENTITY_COUNT
      integer iBase_INVALID_ENTITY_TYPE
      integer iBase_INVALID_ENTITY_TOPOLOGY
      integer iBase_BAD_TYPE_AND_TOPO
      integer iBase_ENTITY_CREATION_ERROR
      integer iBase_INVALID_TAG_HANDLE
      integer iBase_TAG_NOT_FOUND
      integer iBase_TAG_ALREADY_EXISTS
      integer iBase_TAG_IN_USE
      integer iBase_INVALID_ENTITYSET_HANDLE
      integer iBase_INVALID_ITERATOR_HANDLE
      integer iBase_INVALID_ARGUMENT
      integer iBase_MEMORY_ALLOCATION_FAILED
      integer iBase_NOT_SUPPORTED
      integer iBase_FAILURE

      parameter (iBase_SUCCESS = 0)
      parameter (iBase_MESH_ALREADY_LOADED = 1)
      parameter (iBase_FILE_NOT_FOUND = 2)
      parameter (iBase_FILE_WRITE_ERROR = 3)
      parameter (iBase_NIL_ARRAY = 4)
      parameter (iBase_BAD_ARRAY_SIZE = 5)
      parameter (iBase_BAD_ARRAY_DIMENSION = 6)
      parameter (iBase_INVALID_ENTITY_HANDLE = 7)
      parameter (iBase_INVALID_ENTITY_COUNT = 8)
      parameter (iBase_INVALID_ENTITY_TYPE = 9)
      parameter (iBase_INVALID_ENTITY_TOPOLOGY = 10)
      parameter (iBase_BAD_TYPE_AND_TOPO = 11)
      parameter (iBase_ENTITY_CREATION_ERROR = 12)
      parameter (iBase_INVALID_TAG_HANDLE = 13)
      parameter (iBase_TAG_NOT_FOUND = 14)
      parameter (iBase_TAG_ALREADY_EXISTS = 15)
      parameter (iBase_TAG_IN_USE = 16)
      parameter (iBase_INVALID_ENTITYSET_HANDLE = 17)
      parameter (iBase_INVALID_ITERATOR_HANDLE = 18)
      parameter (iBase_INVALID_ARGUMENT = 19)
      parameter (iBase_MEMORY_ALLOCATION_FAILED = 20)
      parameter (iBase_NOT_SUPPORTED = 21)
      parameter (iBase_FAILURE = 22)


      integer iBase_BLOCKED
      integer iBase_INTERLEAVED
      integer iBase_UNDETERMINED

      parameter (iBase_BLOCKED = 0)
      parameter (iBase_INTERLEAVED = 1)
      parameter (iBase_UNDETERMINED = 2)


      integer iBase_INTEGER
      integer iBase_DOUBLE
      integer iBase_ENTITY_HANDLE
      integer iBase_BYTES

      parameter (iBase_INTEGER = 0)
      parameter (iBase_DOUBLE = 1)
      parameter (iBase_ENTITY_HANDLE = 2)
      parameter (iBase_BYTES = 3)

