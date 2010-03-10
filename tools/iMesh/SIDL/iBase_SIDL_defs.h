/**
 * Copyright 2006 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */
#include "babel_config.h"
#ifdef __cplusplus
#include "sidl_cxx.hh"
#endif

#if defined(PACKAGE_VERSION)
#define ARRAY_PTR(array, type) reinterpret_cast<type*>(array._get_ior()->d_firstElement)

#define ARRAY_SIZE(array) (array._is_nil() ? 0 : (int)(array.upper(0) - array.lower(0) + 1))

#define ASSIGN_ARRAY(my_arr) \
          if (!my_arr ## _i_allocated && \
              ARRAY_SIZE(my_arr) < my_arr ## _allocated_size) {\
                PROCESS_ERROR_MSG(iBase_MEMORY_ALLOCATION_FAILED,\
                "Array passed in is non-zero length but not long enough.");} \
          else if (my_arr ## _i_allocated && my_arr ## _allocated_size > 0) {\
            my_arr = my_arr.create1d(0);\
            free(my_arr._get_ior()->d_firstElement);\
            my_arr._get_ior()->d_firstElement = my_arr ## _temp;\
            my_arr._get_ior()->d_lower[0] = 0;\
            my_arr._get_ior()->d_upper[0] = my_arr ## _allocated_size -1;\
            my_arr._get_ior()->d_stride[0] = 1;\
          }

#define CREATE_TEMP_ARRAY(this_type, my_arr) \
          bool my_arr ## _i_allocated = (my_arr._get_ior() == NULL); \
          this_type *my_arr ## _temp = (my_arr._get_ior() == NULL ? NULL : \
              my_arr._get_ior()->d_firstElement); \
          int my_arr ## _allocated_size = (my_arr._get_ior() == NULL ? 0 :\
          my_arr._get_ior()->d_upper[0] - my_arr._get_ior()->d_lower[0] + 1)

#define CREATE_TEMP_EH_ARRAY(my_arr) \
          bool my_arr ## _i_allocated = (my_arr._get_ior() == NULL); \
          iBase_EntityHandle *my_arr ## _temp = reinterpret_cast<iBase_EntityHandle*>((my_arr._get_ior() == NULL ? NULL : \
              my_arr._get_ior()->d_firstElement)); \
          int my_arr ## _allocated_size = (my_arr._get_ior() == NULL ? 0 :\
          my_arr._get_ior()->d_upper[0] - my_arr._get_ior()->d_lower[0] + 1)

#define CREATE_TEMP_ESH_ARRAY(my_arr) \
          bool my_arr ## _i_allocated = (my_arr._get_ior() == NULL); \
          iBase_EntitySetHandle *my_arr ## _temp = reinterpret_cast<iBase_EntitySetHandle*>((my_arr._get_ior() == NULL ? NULL : \
              my_arr._get_ior()->d_firstElement)); \
          int my_arr ## _allocated_size = (my_arr._get_ior() == NULL ? 0 :\
          my_arr._get_ior()->d_upper[0] - my_arr._get_ior()->d_lower[0] + 1)

#define CREATE_TEMP_TH_ARRAY(my_arr) \
          bool my_arr ## _i_allocated = (my_arr._get_ior() == NULL); \
          iBase_TagHandle *my_arr ## _temp = reinterpret_cast<iBase_TagHandle*>((my_arr._get_ior() == NULL ? NULL : \
              my_arr._get_ior()->d_firstElement)); \
          int my_arr ## _allocated_size = (my_arr._get_ior() == NULL ? 0 :\
          my_arr._get_ior()->d_upper[0] - my_arr._get_ior()->d_lower[0] + 1)

#define ASSIGN_TAG_ARRAY(my_arr) \
          if (!my_arr ## _i_allocated && \
              ARRAY_SIZE(my_arr) < my_arr ## _allocated_size) {\
                PROCESS_ERROR_MSG(iBase_MEMORY_ALLOCATION_FAILED,\
                "Array passed in is non-zero length but not long enough.")} \
          else if (my_arr ## _i_allocated && my_arr ## _allocated_size > 0) {\
            my_arr = my_arr.create1d(0);\
            free(my_arr._get_ior()->d_firstElement);\
            my_arr._get_ior()->d_firstElement = (char*) my_arr ## _temp;\
            my_arr._get_ior()->d_lower[0] = 0;\
            my_arr._get_ior()->d_upper[0] = (my_arr ## _allocated_size) -1;\
            my_arr._get_ior()->d_stride[0] = 1;\
          };

#define CREATE_TEMP_TAG_ARRAY(my_arr) \
          bool my_arr ## _i_allocated = (my_arr._get_ior() == NULL); \
          char *my_arr ## _temp = (my_arr._get_ior() == NULL ? NULL : \
              my_arr._get_ior()->d_firstElement); \
          int my_arr ## _allocated_size = (my_arr._get_ior() == NULL ? 0 :\
          (my_arr._get_ior()->d_upper[0] - my_arr._get_ior()->d_lower[0] + 1))

#define ASSIGN_ENUM_ARRAY(my_arr) \
          if (!my_arr ## _i_allocated && \
              ARRAY_SIZE(my_arr) < my_arr ## _allocated_size) {\
                PROCESS_ERROR_MSG(iBase_MEMORY_ALLOCATION_FAILED,\
                "Array passed in is non-zero length but not long enough.");} \
          else if (my_arr ## _i_allocated && my_arr ## _allocated_size > 0) {\
            my_arr = my_arr.create1d(0);\
            free(((sidl_int__array*)my_arr._get_ior())->d_firstElement);\
            ((sidl_int__array*)my_arr._get_ior())->d_firstElement = (int*)my_arr ## _temp;\
            ((sidl_int__array*)my_arr._get_ior())->d_lower[0] = 0;\
            ((sidl_int__array*)my_arr._get_ior())->d_upper[0] = my_arr ## _size -1;\
            ((sidl_int__array*)my_arr._get_ior())->d_stride[0] = 1; \
          }

#define CREATE_TEMP_ENUM_ARRAY(this_type, my_arr) \
          bool my_arr ## _i_allocated = (my_arr._get_ior() == NULL); \
          this_type *my_arr ## _temp = (this_type*)(my_arr._get_ior() == NULL ? NULL : \
                                                    ((sidl_int__array*)my_arr._get_ior())->d_firstElement); \
          int my_arr ## _allocated_size = (my_arr._get_ior() == NULL ? 0 : \
                                           ((sidl_int__array*)my_arr._get_ior())->d_upper[0] - \
                                           ((sidl_int__array*)my_arr._get_ior())->d_lower[0] + 1)

#define TEMP_ARRAY_INOUT(my_arr) &my_arr ## _temp, &my_arr ## _allocated_size, &my_arr ## _size

#define TEMP_EH_ARRAY_INOUT(my_arr) &my_arr ## _temp, &my_arr ## _allocated_size, &my_arr ## _size

#define TEMP_ARRAY_IN(my_arr) \
          (my_arr._get_ior() == NULL ? NULL : my_arr._get_ior()->d_firstElement), \
          my_arr ## _size

#define TEMP_TYPED_ARRAY_IN(type, my_arr) \
          reinterpret_cast<type*>(my_arr._get_ior() == NULL ? NULL : my_arr._get_ior()->d_firstElement), \
          my_arr ## _size

#define TEMP_TAG_ARRAY_IN(my_arr) \
          my_arr._get_ior()->d_firstElement, \
          (my_arr._get_ior()->d_upper[0] - my_arr._get_ior()->d_lower[0] + 1)

#else
#define ARRAY_PTR(array, type) reinterpret_cast<type*>(array._get_ior()->d_firstElement)

#define ARRAY_SIZE(array) (array._is_nil() ? 0 : (int)(array.upper(0) - array.lower(0) + 1))

#define ASSIGN_ARRAY(my_arr) ASSIGN_TYPED_ARRAY(void*, my_arr)
#define ASSIGN_TYPED_ARRAY(type, my_arr) \
          if (!my_arr ## _i_allocated && \
              ARRAY_SIZE(my_arr) < my_arr ## _allocated_size) {\
                PROCESS_ERROR_MSG(iBase_MEMORY_ALLOCATION_FAILED,\
                "Array passed in is non-zero length but not long enough.");} \
          else if (my_arr ## _i_allocated && my_arr ## _allocated_size > 0) {\
            my_arr = my_arr.create1d(0);\
            free(my_arr._get_ior()->d_firstElement);\
            my_arr._get_ior()->d_firstElement = reinterpret_cast<type*>(my_arr ## _temp);\
            my_arr._get_ior()->d_metadata.d_lower[0] = 0;\
            my_arr._get_ior()->d_metadata.d_upper[0] = my_arr ## _allocated_size -1;\
            my_arr._get_ior()->d_metadata.d_stride[0] = 1;\
          }

#define ASSIGN_ARRAY(my_arr) ASSIGN_TYPED_ARRAY(void*, my_arr)
#define ASSIGN_TYPED_ARRAY(type, my_arr) \
          if (!my_arr ## _i_allocated && \
              ARRAY_SIZE(my_arr) < my_arr ## _allocated_size) {\
                PROCESS_ERROR_MSG(iBase_MEMORY_ALLOCATION_FAILED,\
                "Array passed in is non-zero length but not long enough.");} \
          else if (my_arr ## _i_allocated && my_arr ## _allocated_size > 0) {\
            my_arr = my_arr.create1d(0);\
            free(my_arr._get_ior()->d_firstElement);\
            my_arr._get_ior()->d_firstElement = reinterpret_cast<type*>(my_arr ## _temp);\
            my_arr._get_ior()->d_metadata.d_lower[0] = 0;\
            my_arr._get_ior()->d_metadata.d_upper[0] = my_arr ## _allocated_size -1;\
            my_arr._get_ior()->d_metadata.d_stride[0] = 1;\
          }

#define CREATE_TEMP_ARRAY(this_type, my_arr) \
          bool my_arr ## _i_allocated = (my_arr._get_ior() == NULL); \
          this_type *my_arr ## _temp = (my_arr._get_ior() == NULL ? NULL : \
              my_arr._get_ior()->d_firstElement); \
          int my_arr ## _allocated_size = (my_arr._get_ior() == NULL ? 0 :\
          my_arr._get_ior()->d_metadata.d_upper[0] - my_arr._get_ior()->d_metadata.d_lower[0] + 1)

#define CREATE_TEMP_EH_ARRAY(my_arr) \
          bool my_arr ## _i_allocated = (my_arr._get_ior() == NULL); \
          iBase_EntityHandle *my_arr ## _temp = reinterpret_cast<iBase_EntityHandle*>((my_arr._get_ior() == NULL ? NULL : \
              my_arr._get_ior()->d_firstElement)); \
          int my_arr ## _allocated_size = (my_arr._get_ior() == NULL ? 0 :\
          my_arr._get_ior()->d_metadata.d_upper[0] - my_arr._get_ior()->d_metadata.d_lower[0] + 1)

#define CREATE_TEMP_ESH_ARRAY(my_arr) \
          bool my_arr ## _i_allocated = (my_arr._get_ior() == NULL); \
          iBase_EntitySetHandle *my_arr ## _temp = reinterpret_cast<iBase_EntitySetHandle*>((my_arr._get_ior() == NULL ? NULL : \
              my_arr._get_ior()->d_firstElement)); \
          int my_arr ## _allocated_size = (my_arr._get_ior() == NULL ? 0 :\
          my_arr._get_ior()->d_metadata.d_upper[0] - my_arr._get_ior()->d_metadata.d_lower[0] + 1)

#define CREATE_TEMP_TH_ARRAY(my_arr) \
          bool my_arr ## _i_allocated = (my_arr._get_ior() == NULL); \
          iBase_TagHandle *my_arr ## _temp = reinterpret_cast<iBase_TagHandle*>((my_arr._get_ior() == NULL ? NULL : \
              my_arr._get_ior()->d_firstElement)); \
          int my_arr ## _allocated_size = (my_arr._get_ior() == NULL ? 0 :\
          my_arr._get_ior()->d_metadata.d_upper[0] - my_arr._get_ior()->d_metadata.d_lower[0] + 1)

#define ASSIGN_TAG_ARRAY(my_arr) \
          if (!my_arr ## _i_allocated && \
              ARRAY_SIZE(my_arr) < my_arr ## _allocated_size) {\
                PROCESS_ERROR_MSG(iBase_MEMORY_ALLOCATION_FAILED,\
                "Array passed in is non-zero length but not long enough.");} \
          else if (my_arr ## _i_allocated && my_arr ## _allocated_size > 0) {\
            my_arr = my_arr.create1d(0);\
            free(my_arr._get_ior()->d_firstElement);\
            my_arr._get_ior()->d_firstElement = (char*) my_arr ## _temp;\
            my_arr._get_ior()->d_metadata.d_lower[0] = 0;\
            my_arr._get_ior()->d_metadata.d_upper[0] = (my_arr ## _allocated_size) -1;\
            my_arr._get_ior()->d_metadata.d_stride[0] = 1;\
          };

#define CREATE_TEMP_TAG_ARRAY(my_arr) \
          bool my_arr ## _i_allocated = (my_arr._get_ior() == NULL); \
          char *my_arr ## _temp = (my_arr._get_ior() == NULL ? NULL : \
              my_arr._get_ior()->d_firstElement); \
          int my_arr ## _allocated_size = (my_arr._get_ior() == NULL ? 0 :\
          (my_arr._get_ior()->d_metadata.d_upper[0] - my_arr._get_ior()->d_metadata.d_lower[0] + 1))

#define ASSIGN_ENUM_ARRAY(my_arr) \
          if (!my_arr ## _i_allocated && \
              ARRAY_SIZE(my_arr) < my_arr ## _allocated_size) {\
                PROCESS_ERROR_MSG(iBase_MEMORY_ALLOCATION_FAILED,\
                "Array passed in is non-zero length but not long enough.");} \
          else if (my_arr ## _i_allocated && my_arr ## _allocated_size > 0) {\
            my_arr = my_arr.create1d(0);\
            free(((sidl_int__array*)my_arr._get_ior())->d_firstElement);\
            ((sidl_int__array*)my_arr._get_ior())->d_firstElement = (int*)my_arr ## _temp;\
            ((sidl_int__array*)my_arr._get_ior())->d_metadata.d_lower[0] = 0;\
            ((sidl_int__array*)my_arr._get_ior())->d_metadata.d_upper[0] = my_arr ## _size -1;\
            ((sidl_int__array*)my_arr._get_ior())->d_metadata.d_stride[0] = 1; \
          }

#define CREATE_TEMP_ENUM_ARRAY(this_type, my_arr) \
          bool my_arr ## _i_allocated = (my_arr._get_ior() == NULL); \
          this_type *my_arr ## _temp = (this_type*)(my_arr._get_ior() == NULL ? NULL : \
                                                    ((sidl_int__array*)my_arr._get_ior())->d_firstElement); \
          int my_arr ## _allocated_size = (my_arr._get_ior() == NULL ? 0 : \
                                           ((sidl_int__array*)my_arr._get_ior())->d_metadata.d_upper[0] - \
                                           ((sidl_int__array*)my_arr._get_ior())->d_metadata.d_lower[0] + 1)

#define TEMP_ARRAY_INOUT(my_arr) &my_arr ## _temp, &my_arr ## _allocated_size, &my_arr ## _size

#define TEMP_ARRAY_IN(my_arr) \
          (my_arr._get_ior() == NULL ? NULL : my_arr._get_ior()->d_firstElement), \
          my_arr ## _size

#define TEMP_TYPED_ARRAY_IN(type, my_arr) \
          reinterpret_cast<type*>(my_arr._get_ior() == NULL ? NULL : my_arr._get_ior()->d_firstElement), \
          my_arr ## _size

#define TEMP_TAG_ARRAY_IN(my_arr) \
          my_arr._get_ior()->d_firstElement, \
          (my_arr._get_ior()->d_metadata.d_upper[0] - my_arr._get_ior()->d_metadata.d_lower[0] + 1)


#endif

#ifdef __cplusplus
template <class T> static inline sidl::array<T> alloc_sidl_vector( size_t size, T init )
{
  sidl::array<T> result = result.create1d(size);
  for (int32_t i = 0; i < (int32_t)size; ++i)
    result.set( i, init );
  return result;
}

template <class T> static inline 
sidl::array<T> convert_to_sidl_vector( T* array, size_t size )
{
  sidl::array<T> result;
  int32_t lower = 0, upper = size - 1, stride = 1;
  result.borrow( array, 1, &lower, &upper, &stride );
  return result;
}

template <class T> static inline 
sidl::array<T> convert_to_sidl_vector( std::vector<T>& vect, int offset = 0 )
{
  sidl::array<T> result;
  int32_t lower = 0, upper = vect.size() - offset - 1, stride = 1;
  result.borrow( &vect[offset], 1, &lower, &upper, &stride );
  return result;
}

template <class S, class T> static inline void copy_from_sidl( sidl::array<S>& source,
                                                        T* target )
{
  typename sidl::array<S>::iterator i = source.begin();
  for (; i != source.end(); ++i, ++target)
    *target = (T)*i;
}

template <class S, class T> static inline void append_from_sidl( sidl::array<S>& source,
                                                                 int source_size,
                                                                 std::vector<T>& target )
{
  for (int i = 0; i < source_size; ++i)
    target.push_back((T)source[i]);
}
#endif

#define CAST_IREL_INTERFACE(var_in, var_out, iface, ret) \
          iRel::iface var_out = var_in;\
          if (var_out._is_nil()) {\
            std::cerr << "Error: cast to iRel iface failed." << std::endl; \
            return ret;\
          }

#define CAST_IGEOM_INTERFACE(var_in, var_out, iface, ret) \
          iGeom::iface var_out = var_in;\
          if (var_out._is_nil()) {\
            std::cerr << "Error: cast to iGeom iface failed." << std::endl; \
            return ret;\
          }

#define CAST_IMESH_INTERFACE(var_in, var_out, iface, ret) \
          iMesh::iface var_out = var_in;\
          if (var_out._is_nil()) {\
            std::cerr << "Error: cast to iMesh iface failed." << std::endl; \
            return ret;\
          }

#define CAST_IBASE_INTERFACE(var_in, var_out, iface, ret) \
          iBase::iface var_out = var_in;\
          if (var_out._is_nil()) {\
            std::cerr << "Error: cast to iBase iface failed." << std::endl; \
            return ret;\
          }

#define CHECK_SIZE(array, size)  \
          if (array._is_nil() || ARRAY_SIZE(array) == 0) array = array.create1d(size); \
          else if (ARRAY_SIZE(array) < size) { \
             std::cerr << "   Array passed in is non-zero but too short." << std::cerr; \
           return false; }

#define CHECK_SIZE_VOID(type, array, allocated_size, size)  \
          if (NULL == *array || *allocated_size == 0) {\
            *array = (type *) malloc(sizeof(type) * size); \
            *allocated_size = size;} \
          else if (*allocated_size < size) { \
             std::cerr << "   Array passed in is non-zero but too short." << std::cerr; }

#ifdef TEST_DEFS
CREATE_TEMP_ARRAY(iBase_ErrorType, error_table);

TEMP_ARRAY_INOUT(error_table);

TEMP_ARRAY_IN(error_table);

ASSIGN_ARRAY(error_table);
#endif
