/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#include <stdlib.h>
#include <string.h>
#include <H5Tpublic.h>
#include <H5Gpublic.h>
#include <H5Dpublic.h>
#include "status.h"
#include "file-handle.h"
#include "mhdf.h"
#include "util.h"
#include "names-and-paths.h"

hid_t
mhdf_getNativeType( hid_t input_type,
                    int size,
                    mhdf_Status* status )
{
  hid_t result;
  H5T_sign_t sgn;
  H5T_class_t cls;
  hid_t tmp_id, type_id;
  
  mhdf_setOkay( status );
  
  cls = H5Tget_class( input_type );
  switch( cls )
  {
    case H5T_FLOAT:
      switch( size ) 
      {
        case  4: return H5T_NATIVE_FLOAT;
        case  8: return H5T_NATIVE_DOUBLE;
        case 16: return H5T_NATIVE_LDOUBLE;
        default:
          mhdf_setFail( status, "Invalid size for floating point type: %d", size );
          return -1;
      }

    case H5T_INTEGER:
      sgn = H5Tget_sign( input_type );
      if (H5T_SGN_ERROR == sgn)
      {
        mhdf_setFail( status, "Internall errror calling H5Tget_sign." );
        return -1;
      }
      if (     sizeof(      char ) == size)
        return sgn == H5T_SGN_NONE ? H5T_NATIVE_UCHAR  : H5T_NATIVE_SCHAR;
      else if (sizeof(     short ) == size)
        return sgn == H5T_SGN_NONE ? H5T_NATIVE_USHORT : H5T_NATIVE_SHORT;
      else if (sizeof(       int ) == size)
        return sgn == H5T_SGN_NONE ? H5T_NATIVE_UINT   : H5T_NATIVE_INT;
      else if (sizeof(      long ) == size)
        return sgn == H5T_SGN_NONE ? H5T_NATIVE_ULONG  : H5T_NATIVE_LONG;
      else if (sizeof( long long ) == size)
        return sgn == H5T_SGN_NONE ? H5T_NATIVE_ULLONG : H5T_NATIVE_LLONG;
      
      mhdf_setFail( status, "Invalid size for integer type: %d", size );
      return -1;
      
    case H5T_ENUM:
      tmp_id = H5Tget_super( input_type );
      if (tmp_id < 0)
      {
        mhdf_setFail( status, "Internal error calling H5Tget_super." );
        return -1;
      }
      type_id = mhdf_getNativeType( tmp_id, size, status );
      H5Tclose( tmp_id );
      return type_id;

    case H5T_TIME:
    case H5T_OPAQUE:
    case H5T_REFERENCE:
      mhdf_setFail( status, "Unsupported type class." );
      return -1;

    case H5T_COMPOUND:
    case H5T_VLEN:
    case H5T_ARRAY:
    case H5T_STRING:
      mhdf_setFail( status, "Only atomic types are supported." );
      return -1;

    default:
      mhdf_setFail( status, "Internal error calling H5Tget_class.  Bad handle?" );
      return -1;
  }
}

static hid_t get_tag( mhdf_FileHandle file_handle,
                      const char* tag_name,
                      mhdf_Status* status )
{
  hid_t group_id, tag_id;
  char* path;
  FileHandle* file_ptr;

  file_ptr = (FileHandle*)file_handle;
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;

  path = mhdf_name_to_path_copy( tag_name, status );
  if (NULL == path)
    return -1;
  
  group_id = H5Gopen( file_ptr->hdf_handle, TAG_GROUP );
  if (group_id < 0)
  {
    mhdf_setFail( status, "Failed to open tag group." );
    free( path );
    return -1;
  }
  
  tag_id = H5Gopen( group_id, path );
  H5Gclose( group_id );
  free( path );
  if (tag_id < 0)
  {
    mhdf_setFail( status, "Failed to open tag data for tag \"%s\".", tag_name );
    return -1;
  }
  
  mhdf_setOkay( status );
  return tag_id;
}
  

static hid_t get_tag_type( FileHandle* file_ptr,
                           const char* tag_path,
                           mhdf_Status* status )
{
  hid_t group_id, tag_id, type_id;
  
  group_id = H5Gopen( file_ptr->hdf_handle, TAG_GROUP );
  if (group_id < 0)
  {
    mhdf_setFail( status, "Failed to open tag group." );
    return -1;
  }
  
  tag_id = H5Gopen( group_id, tag_path );
  H5Gclose( group_id );
  if (tag_id < 0)
  {
    mhdf_setFail( status, "Failed to open group for tag \"%s\".", tag_path );
    return -1;
  }
  
  type_id = H5Topen( tag_id, TAG_TYPE_NAME );
  H5Gclose( tag_id );
  if (type_id < 0)
  {
    mhdf_setFail( status, "Failed to open type data for tag \"%s\".", tag_path );
    return -1;
  }
  
  return type_id;
}
  


void
mhdf_createOpaqueTag( mhdf_FileHandle file_handle,
                      const char* tag_name,
                      size_t tag_size,
                      void* default_value,
                      void* global_value,
                      int tstt_tag_type,
                      mhdf_Status* status )
{
  hid_t type_id;
  
  type_id = H5Tcreate( H5T_OPAQUE, tag_size );
  if (type_id < 0)
  {
    mhdf_setFail( status, "Internal failure in HT5create(H5T_OPAQUE,%d).", (int)tag_size);
    return;
  }
  
  if (0 > H5Tset_tag( type_id, tag_name ))
  {
    H5Tclose( type_id );
    mhdf_setFail( status, "H5Tset_tag failed.");
    return;
  }

  mhdf_createTypeTag( file_handle, tag_name, type_id,
                      default_value, global_value, 
                      tstt_tag_type, status );

  H5Tclose( type_id );
}


void
mhdf_createBitTag( mhdf_FileHandle file_handle,
                   const char* tag_name,
                   size_t tag_size,
                   void* default_value,
                   void* global_value,
                   int tstt_tag_type,
                   mhdf_Status* status )
{
  hid_t type_id, base_id;
  
  if (tag_size <= 8)
    base_id = H5T_NATIVE_B8;
  else if (tag_size <= 16)
    base_id = H5T_NATIVE_B16;
  else if (tag_size <= 32)
    base_id = H5T_NATIVE_B32;
  else if (tag_size <= 64)
    base_id = H5T_NATIVE_B64;
  else
  {
    mhdf_setFail( status, "Cannot createa a bit tag larger than 64-bits.  %d bits requested.\n", (int)tag_size);
    return;
  }
  
  type_id = H5Tcopy( base_id );
  if (type_id < 0)
  {
    mhdf_setFail( status, "Internal failure in HT5copy.", (int)tag_size);
    return;
  }
  
  if (0 > H5Tset_precision( type_id, tag_size ))
  {
    H5Tclose( type_id );
    mhdf_setFail( status, "H5Tset_precision failed.");
    return;
  }

  mhdf_createTypeTag( file_handle, tag_name, type_id,
                      default_value, global_value, 
                      tstt_tag_type, status );

  H5Tclose( type_id );
}


void
mhdf_createTypeTag( mhdf_FileHandle file_handle,
                    const char* tag_name,
                    hid_t hdf5_tag_type,
                    void* default_value,
                    void* global_value,
                    int tstt_tag_type,
                    mhdf_Status* status )
{
  hid_t type_id, group_id, tag_id;
  char* path;
  FileHandle* file_ptr;
  herr_t rval;
  API_BEGIN;

    /* Validate input */
  
  file_ptr = (FileHandle*)file_handle;
  if (!mhdf_check_valid_file( file_ptr, status ))
    return ;

  if (!tag_name || !*tag_name)
  {
    mhdf_setFail( status, "Invalid tag name" );
    return ;
  }
  
  
    /* Open the tag group */

  group_id = H5Gopen( file_ptr->hdf_handle, TAG_GROUP );
  if (group_id < 0)
  {
    mhdf_setFail( status, "H5Gopen(\"%s\") failed.", TAG_GROUP );
    return;
  }

    /* Create path string for tag object */

  path = mhdf_name_to_path_copy( tag_name, status );
  if (!path) 
  { 
    H5Gclose( group_id );
    return; 
  }
  
    /* Create group for this tag */

  tag_id = H5Gcreate( group_id, path, 3 );
  if (tag_id < 0)
  {
     mhdf_setFail( status, "H5Gcreate( \"%s\" ) failed.", path );
     free( path );
     H5Gclose( group_id );
     return;
  }
  
    /* Store the tag name as the comment on the group entry */
  
  rval = H5Gset_comment( group_id, path, tag_name );
  H5Gclose( group_id );
  free( path );
  if (rval < 0)
  {
    mhdf_setFail( status, "H5Gset_comment failed for tag \"%s\"", tag_name );
    H5Gclose( tag_id );
    return;
  }

    /* Store TSTT tag type as attribute */

  rval = mhdf_create_scalar_attrib( tag_id, 
                                   TAG_TYPE_ATTRIB,
                                   H5T_NATIVE_INT,
                                   &tstt_tag_type,
                                   status );
  if (!rval)
  {
    H5Gclose( tag_id );
    return;
  }


    /* Create tag type object, or write attribute if opaque */
 
  type_id = H5Tcopy( hdf5_tag_type );
  rval = H5Tcommit( tag_id, TAG_TYPE_NAME, type_id );
  if (rval < 0)
  {
    mhdf_setFail( status, "H5Tcommit failed for tag \"%s\"", tag_name );
    H5Gclose( tag_id );
    return;
  }

    /* Store the default value as a attribute of the tag group */

  if (default_value)
  {
    rval = mhdf_create_scalar_attrib( tag_id, 
                                     TAG_DEFAULT_ATTRIB, 
                                     type_id, 
                                     default_value, 
                                     status);
    if (!rval) 
    { 
      H5Gclose( tag_id );
      H5Tclose( type_id );
      return; 
    }
  }

    
    /* Store global tag value as attribute */
  
  if (global_value)
  {
    rval = mhdf_create_scalar_attrib( tag_id,  
                                     TAG_GLOBAL_ATTRIB,
                                     type_id,
                                     global_value, 
                                     status );
    if (!rval) 
    { 
      H5Gclose( tag_id );
      H5Tclose( type_id );
      return; 
    }
  }

  H5Gclose( tag_id );
  H5Tclose( type_id );
  mhdf_setOkay( status );
  API_END;
}

int
mhdf_getNumberTags( mhdf_FileHandle file_handle, mhdf_Status* status )
{
  hid_t group_id;
  hsize_t result;
  FileHandle* file_ptr;
  API_BEGIN;

    /* Validate input */
  
  file_ptr = (FileHandle*)file_handle;
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;
  
    /* Open the tags group */
  
  group_id = H5Gopen( file_ptr->hdf_handle, TAG_GROUP );
  if (group_id < 0)
  {
    mhdf_setFail( status, "H5Gopen(\"%s\") failed", TAG_GROUP);
    return -1;
  }
  
    /* Get number of objects in tags group */
  
  if (H5Gget_num_objs( group_id, &result ) < 0)
  {
    mhdf_setFail( status, "Internal failure calling H5Gget_num_objs.");
    H5Gclose(group_id);
    return -1;
  }
  
  H5Gclose( group_id );
  mhdf_setOkay( status );
  API_END;
  return (int)result;
}

char**
mhdf_getTagNames( mhdf_FileHandle file_handle,
                  int* num_names_out,
                  mhdf_Status* status )
{
  hid_t group_id;
  FileHandle* file_ptr;
  hsize_t count, index;
  char* name;
  char** result;
  ssize_t size;
  API_BEGIN;
  

    /* Validate input */
  
  file_ptr = (FileHandle*)file_handle;
  if (!mhdf_check_valid_file( file_ptr, status ))
    return NULL;
  
    /* Open the tags group */
  
  group_id = H5Gopen( file_ptr->hdf_handle, TAG_GROUP );
  if (group_id < 0)
  {
    mhdf_setFail( status, "H5Gopen(\"%s\") failed", TAG_GROUP);
    return NULL;
  }
  
    /* Get number of objects in tags group */
  
  if (H5Gget_num_objs( group_id, &count ) < 0)
  {
    mhdf_setFail( status, "Internal failure calling H5Gget_num_objs.");
    H5Gclose(group_id);
    return NULL;
  }
  
    /* No tags? */

  *num_names_out = (int)count;
  if (count == 0)
  {
    H5Gclose( group_id );
    mhdf_setOkay( status );
    return NULL;
  }
  
    /* Allocate string array */
  
  result = (char**)mhdf_malloc( sizeof(char*) * count, status );
  if (NULL == result)
  {
    H5Gclose( group_id );
    return NULL;
  }
  
    /* Get names */
  
  for (index = 0; index < count; ++index)
  {
    size = H5Gget_objname_by_idx( group_id, index, NULL, 0 );
    if (size < 1 || NULL == (name = (char*)mhdf_malloc( size+1, status )))
    {
      while ((--index) > 0)
        free( result[index] );
      free ( result );
      H5Gclose( group_id );
      mhdf_setFail( status, "Internal failure calling H5Gget_objname_by_idx.");
      return NULL;
    }
    
    H5Gget_objname_by_idx( group_id, index, name, size + 1 );
    if (!mhdf_path_to_name( name, name ))
    {
      mhdf_setFail( status, "Invalid character string in internal file path: \"%s\"\n",
        name );
      return NULL;
    }
    result[index] = name;
  }
  
  H5Gclose( group_id );
  mhdf_setOkay( status );
  API_END;
  return result;
}
  


void
mhdf_getTagInfo( mhdf_FileHandle file_handle,
                 const char* tag_name,
                 int* tag_data_len_out,
                 int* have_default_out,
                 int* have_global_out,
                 int* is_opaque_type_out,
                 int* have_sparse_data_out,
                 int* tstt_tag_class_out,
                 int* bit_tag_bits_out,
                 hid_t* hdf_type_out,
                 mhdf_Status* status )
{
  hid_t tag_id, type_id;
  int rval;
  unsigned int index;

  API_BEGIN;


    /* Validate input */
  if (NULL == tag_name             ||
      NULL == tag_data_len_out     ||
      NULL == have_default_out     ||
      NULL == have_global_out      ||
      NULL == is_opaque_type_out   ||
      NULL == have_sparse_data_out ||
      NULL == tstt_tag_class_out   ||
      NULL == bit_tag_bits_out     ||
      NULL == hdf_type_out          )
  {
    mhdf_setFail( status, "Invalid input." );
    return;
  }
  
    /* Get group for tag */
  tag_id = get_tag( file_handle, tag_name, status );
  if (tag_id < 0)
    return;
  
    /* Check for sparse data */
  rval = mhdf_is_in_group( tag_id, SPARSE_ENTITY_NAME, status );
  if (rval < 0)
  {
    H5Gclose( tag_id );
    return;
  }
  *have_sparse_data_out = rval ? 1 : 0;

    /* Check if have default value for tag */
  rval = mhdf_find_attribute( tag_id, TAG_DEFAULT_ATTRIB, &index, status );
  if (rval < 0)
  {
    H5Gclose( tag_id );
    return;
  }
  *have_default_out = rval ? 1 : 0;

    /* Check if have global value for tag */
  rval = mhdf_find_attribute( tag_id, TAG_GLOBAL_ATTRIB, &index, status );
  if (rval < 0)
  {
    H5Gclose( tag_id );
    return;
  }
  *have_global_out = rval ? 1 : 0;
  
    /* Get TSTT tag class */
  rval = mhdf_read_scalar_attrib( tag_id, TAG_TYPE_ATTRIB, 
                                  H5T_NATIVE_INT, tstt_tag_class_out,
                                  status );
  if (rval < 1)
  {
    H5Gclose( tag_id );
    return;
  }
  
    /* Get tag type */
  type_id = H5Topen( tag_id, TAG_TYPE_NAME );
  if (type_id < 0)
  {
    H5Gclose( tag_id );
    mhdf_setFail( status, "Failed to get type object for tag \"%s\".", tag_name );
    return ;
  }
  *hdf_type_out = H5Tcopy( type_id );
  H5Tclose( type_id );

  *tag_data_len_out = (int)H5Tget_size( *hdf_type_out );
  if (*tag_data_len_out < 1)
  {
    mhdf_setFail( status, "Invalid opaque tag size: %d.", *tag_data_len_out );
    H5Gclose( tag_id );
    H5Tclose( *hdf_type_out );
    return;
  }
  
    /* Check if tag is opaque or bitfield type */
  *is_opaque_type_out = H5Tget_class( *hdf_type_out ) == H5T_OPAQUE;
  *bit_tag_bits_out = 0;
  if (*is_opaque_type_out)
  {
    H5Tclose( *hdf_type_out );
    *hdf_type_out = 0;
  }
  else if (H5Tget_class( *hdf_type_out ) == H5T_BITFIELD)
  {
    *bit_tag_bits_out = H5Tget_precision( *hdf_type_out );
    H5Tclose( *hdf_type_out );
    *hdf_type_out = 0;
  }

  H5Gclose( tag_id );
  mhdf_setOkay( status );
  API_END;
}    

void
mhdf_getTagValues( mhdf_FileHandle file_handle,
                   const char* tag_name,
                   hid_t output_data_type,
                   void* default_value,
                   void* global_value,
                   mhdf_Status* status )
{
  hid_t tag_id;
  int rval;
  unsigned int junk;
  API_BEGIN;
  
    /* check args */
  if (NULL == tag_name || !*tag_name)
  {
    mhdf_setFail( status, "Invalid input." );
    return;
  }
  
    /* Get the tag group */
  tag_id = get_tag( file_handle, tag_name, status );
  if (tag_id < 0)
    return;
  
    /* Check if tag has default value */
  rval = mhdf_find_attribute( tag_id, TAG_DEFAULT_ATTRIB, &junk, status );
  if (rval < 0)
  {
    H5Gclose( tag_id );
    return;
  }
  
    /* Get default if there is one */
  if (rval)
  {
    if (NULL == default_value)
    {
      mhdf_setFail( status, "Invalid input." );
      return;
    }
    
    rval = mhdf_read_scalar_attrib( tag_id, TAG_DEFAULT_ATTRIB,
                                    output_data_type, default_value,
                                    status );
    if (!rval)
    {
      H5Gclose( tag_id );
      return;
    }
  }
  
    /* Check if tag has global value */
  rval = mhdf_find_attribute( tag_id, TAG_GLOBAL_ATTRIB, &junk, status );
  if (rval < 0)
  {
    H5Gclose( tag_id );
    return;
  }
  
    /* Get global value if there is one */
  if (rval)
  {
    if (NULL == global_value)
    {
      mhdf_setFail( status, "Invalid input." );
      return;
    }
    
    rval = mhdf_read_scalar_attrib( tag_id, TAG_GLOBAL_ATTRIB,
                                    output_data_type, global_value,
                                    status );
    if (!rval)
    {
      H5Gclose( tag_id );
      return;
    }
  }
  
  H5Gclose( tag_id );
  mhdf_setOkay( status );
  API_END;
}

int
mhdf_haveDenseTag( mhdf_FileHandle file_handle,
                   const char* tag_name,
                   const char* type_handle,
                   mhdf_Status* status )
{
  char* path;
  hid_t elem_id, group_id;
  FileHandle* file_ptr;
  int rval = 0;
  API_BEGIN;
  
  file_ptr = (FileHandle*)file_handle;
  if (!mhdf_check_valid_file( file_ptr, status )) return -1;
  
  if (type_handle == mhdf_node_type_handle())
  {
    elem_id = H5Gopen( file_ptr->hdf_handle, NODE_GROUP );
    if (elem_id < 0)
      mhdf_setFail( status, "Could not open node group." );
  }
  else if (type_handle == mhdf_set_type_handle())
  {
    elem_id = H5Gopen( file_ptr->hdf_handle, SET_GROUP );
    if (elem_id < 0)
      mhdf_setFail( status, "Could not open set group." );
  }
  else
  {
    elem_id = mhdf_elem_group_from_handle( file_ptr, type_handle, status );
  }
  if (elem_id < 0) return -1;
  
  rval = mhdf_is_in_group( elem_id, TAG_GROUP_NAME, status );
  if (rval < 0)
  {
    H5Gclose( elem_id );
    return -1;
  }
  else if (rval == 0)
  {
    H5Gclose( elem_id );
    mhdf_setOkay( status );
    return 0;
  }
  
  group_id = H5Gopen( elem_id, DENSE_TAG_SUBGROUP );
  H5Gclose( elem_id );
  if (group_id < 0)
  {
    mhdf_setFail( status, "Could not open tag subgroup." );
    return -1;
  }
  
  path = mhdf_name_to_path_copy( tag_name, status );
  if (NULL == path) { H5Gclose( group_id ); return -1; }
  
  rval = mhdf_is_in_group( group_id, path, status );
  H5Gclose( group_id );
  free( path );
  
  if (rval >= 0)
  {
    mhdf_setOkay( status );
  }
  
  API_END;
  return rval;
}

hid_t
mhdf_createDenseTagData( mhdf_FileHandle file_handle,
                         const char* tag_name,
                         const char* type_handle,
                         long num_values,
                         mhdf_Status* status )
{
  char* path;
  hid_t elem_id, data_id, type_id;
  FileHandle* file_ptr;
  size_t name_len, path_len, dir_len;
  hsize_t size;
  API_BEGIN;
  
  file_ptr = (FileHandle*)file_handle;
  if (!mhdf_check_valid_file( file_ptr, status )) return -1;
  
  if (type_handle == mhdf_node_type_handle())
  {
    elem_id = H5Gopen( file_ptr->hdf_handle, NODE_GROUP );
    if (elem_id < 0)
      mhdf_setFail( status, "Could not open node group." );
  }
  else if (type_handle == mhdf_set_type_handle())
  {
    elem_id = H5Gopen( file_ptr->hdf_handle, SET_GROUP );
    if (elem_id < 0)
      mhdf_setFail( status, "Could not open set group." );
  }
  else
  {
    elem_id = mhdf_elem_group_from_handle( file_ptr, type_handle, status );
  }
  if (elem_id < 0) return -1;
  
  dir_len = strlen( DENSE_TAG_SUBGROUP );
  name_len = mhdf_name_to_path( tag_name, NULL, 0 );
  path_len = dir_len + name_len + 1;
  path = (char*)mhdf_malloc( path_len, status );
  if (NULL == path) 
    { H5Gclose( elem_id ); return -1; }
  strcpy( path, DENSE_TAG_SUBGROUP );
  mhdf_name_to_path( tag_name, path + dir_len, name_len + 1 );

  type_id = get_tag_type( file_ptr, path + dir_len, status );
  if (type_id < 0) 
    { H5Gclose( elem_id ); return -1; }
  
  size = (hsize_t)num_values;
  data_id = mhdf_create_table( elem_id, path, type_id, 1, &size, status );
  free( path );
  H5Gclose( elem_id );
  H5Tclose( type_id );
  
  if (data_id > 0)
    mhdf_setOkay( status );
  
  API_END_H( 1 );
  return data_id;
}

hid_t
mhdf_openDenseTagData(  mhdf_FileHandle file_handle,
                        const char* tag_name,
                        const char* type_handle,
                        long* num_values_out,
                        mhdf_Status* status )
{
  char* path;
  hid_t elem_id, data_id;
  FileHandle* file_ptr;
  size_t name_len, path_len, dir_len;
  hsize_t size;
  API_BEGIN;
  
  file_ptr = (FileHandle*)file_handle;
  if (!mhdf_check_valid_file( file_ptr, status )) return -1;
  
  if (type_handle == mhdf_node_type_handle())
  {
    elem_id = H5Gopen( file_ptr->hdf_handle, NODE_GROUP );
    if (elem_id < 0)
      mhdf_setFail( status, "Could not open node group." );
  }
  else if (type_handle == mhdf_set_type_handle())
  {
    elem_id = H5Gopen( file_ptr->hdf_handle, SET_GROUP );
    if (elem_id < 0)
      mhdf_setFail( status, "Could not open set group." );
  }
  else
  {
    elem_id = mhdf_elem_group_from_handle( file_ptr, type_handle, status );
  }
  if (elem_id < 0) return -1;
  
  dir_len = strlen( DENSE_TAG_SUBGROUP );
  name_len = mhdf_name_to_path( tag_name, NULL, 0 );
  path_len = dir_len + name_len + 1;
  path = (char*)mhdf_malloc( path_len, status );
  if (NULL == path) 
    { H5Gclose( elem_id ); return -1; }
  strcpy( path, DENSE_TAG_SUBGROUP );
  mhdf_name_to_path( tag_name, path + dir_len, name_len + 1 );
  
  data_id = mhdf_open_table( elem_id, path, 1, &size, status );
  free( path );
  H5Gclose( elem_id );
  *num_values_out = (long)size;
  
  if (data_id > 0)
    mhdf_setOkay( status );
  
  API_END_H( 1 );
  return data_id;
}


void
mhdf_writeDenseTag( hid_t tag_table,
                    long offset,
                    long count,
                    hid_t type_id,
                    const void* tag_data,
                    mhdf_Status* status )
{
  hid_t my_type_id;
  API_BEGIN;
  
  if (type_id > 0)
  {
    my_type_id = type_id;
  }
  else
  {
    my_type_id = H5Dget_type( tag_table );
    if (my_type_id < 0)
    {
      mhdf_setFail( status, "Internal error calling H5Dget_type.  Bad handle?" );
      return;
    }
  }

  mhdf_write_data( tag_table, offset, count, my_type_id, tag_data, status );

  if (type_id < 1)
    H5Tclose( my_type_id );
  API_END;
}

void
mhdf_readDenseTag( hid_t tag_table,
                   long offset,
                   long count,
                   hid_t type_id,
                   void* tag_data,
                   mhdf_Status* status )

{
  hid_t my_type_id;
  API_BEGIN;
  
  if (type_id > 0)
  {
    my_type_id = type_id;
  }
  else
  {
    my_type_id = H5Dget_type( tag_table );
    if (my_type_id < 0)
    {
      mhdf_setFail( status, "Internal error calling H5Dget_type.  Bad handle?" );
      return;
    }
  }

  mhdf_read_data( tag_table, offset, count, my_type_id, tag_data, status );

  if (type_id < 1)
    H5Tclose( my_type_id );
  API_END;
}

void
mhdf_createSparseTagData( mhdf_FileHandle file_handle,
                          const char* tag_name,
                          long num_values,
                          hid_t handles_out[2],
                          mhdf_Status* status )
{
  hid_t tag_id, index_id, data_id, type_id;
  hsize_t count = (hsize_t)num_values;
  API_BEGIN;
  
  tag_id = get_tag( file_handle, tag_name, status );
  if (tag_id < 0) return ;
  
  type_id = H5Topen( tag_id, TAG_TYPE_NAME );
  if (type_id < 0)
  {
    H5Gclose( tag_id );
    mhdf_setFail( status, "Failed to get type object for tag \"%s\".", tag_name );
    return ;
  }
  
  index_id = mhdf_create_table( tag_id, SPARSE_ENTITY_NAME,
                                H5T_NATIVE_LONG, 1, &count,
                                status );
  if (index_id < 0) 
  { 
    H5Gclose( tag_id ); 
    H5Tclose( type_id );
    return ; 
  }
  
  data_id = mhdf_create_table( tag_id, SPARSE_VALUES_NAME,
                               type_id, 1, &count, status );
  H5Tclose( type_id );
  H5Gclose( tag_id ); 
  if (data_id < 0) 
  { 
    H5Dclose( index_id );
    return ; 
  }
  
  handles_out[0] = index_id;
  handles_out[1] = data_id;
  mhdf_setOkay( status );
  API_END_H(2);
}


void
mhdf_openSparseTagData( mhdf_FileHandle file_handle,
                        const char* tag_name,
                        long* num_values_out,
                        hid_t handles_out[2],
                        mhdf_Status* status )
{
  hid_t tag_id, index_id, data_id;
  hsize_t size1, size2;
  API_BEGIN;
  
  tag_id = get_tag( file_handle, tag_name, status );
  if (tag_id < 0) return ;
 
  index_id = mhdf_open_table( tag_id, SPARSE_ENTITY_NAME, 1, &size1, status );
  if (index_id < 0) 
  { 
    H5Gclose( tag_id ); 
    return ; 
  }
  
  data_id = mhdf_open_table( tag_id, SPARSE_VALUES_NAME, 1, &size2, status );
  H5Gclose( tag_id ); 
  if (data_id < 0) 
  { 
    H5Dclose( index_id );
    return ; 
  }
  
  if (size1 != size2)
  {
    mhdf_setFail( status, "Data length mismatch for sparse tag data -- invalid file.");
    H5Dclose( index_id );
    H5Dclose( data_id );
    return ;
  }
  *num_values_out = (long)size1;
  
  handles_out[0] = index_id;
  handles_out[1] = data_id;
  mhdf_setOkay( status );
  API_END_H(2);
}

void
mhdf_writeSparseTagEntities( hid_t table_id,
                             long offset,
                             long count,
                             hid_t int_type,
                             const void* id_list,
                             mhdf_Status* status )
{
  API_BEGIN;
  mhdf_write_data( table_id, offset, count, int_type, id_list, status );
  API_END;
}
                        
void
mhdf_writeSparseTagValues( hid_t table_id,
                           long offset,
                           long count,
                           hid_t tag_type,
                           const void* tag_data,
                           mhdf_Status* status )
{
  hid_t type_id;
  API_BEGIN;
  
  if (tag_type > 0)
  {
    type_id = tag_type;
  }
  else
  {
    type_id = H5Dget_type( table_id );
    if (type_id < 0)
    {
      mhdf_setFail( status, "Internal error calling H5Dget_type.  Bad handle?" );
      return;
    }
  }
  
  mhdf_write_data( table_id, offset, count, type_id, tag_data, status );

  if (tag_type < 1)
    H5Tclose( type_id );
  API_END;
}

void
mhdf_readSparseTagEntities( hid_t table_id,
                            long offset,
                            long count,
                            hid_t int_type,
                            void* id_list,
                            mhdf_Status* status )
{
  API_BEGIN;
  mhdf_read_data( table_id, offset, count, int_type, id_list, status );
  API_END;
}
                        
void
mhdf_readSparseTagValues( hid_t table_id,
                          long offset,
                          long count,
                          hid_t tag_type,
                          void* tag_data,
                          mhdf_Status* status )
{
  hid_t type_id;
  API_BEGIN;
  
  if (tag_type > 0)
  {
    type_id = tag_type;
  }
  else
  {
    type_id = H5Dget_type( table_id );
    if (type_id < 0)
    {
      mhdf_setFail( status, "Internal error calling H5Dget_type.  Bad handle?" );
      return;
    }
  }
  
  mhdf_read_data( table_id, offset, count, type_id, tag_data, status );

  if (tag_type < 1)
    H5Tclose( type_id );
  API_END;
}
