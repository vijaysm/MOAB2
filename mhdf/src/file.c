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
#include <assert.h>
#include <H5Fpublic.h>
#include <H5Ppublic.h>
#include <H5Gpublic.h>
#include <H5Spublic.h>
#include <H5Tpublic.h>
#include <H5Apublic.h>
#include "mhdf.h"
#include "status.h"
#include "names-and-paths.h"
#include "util.h"
#include "file-handle.h"

static int
make_hdf_group( const char* path, hid_t file, size_t size, mhdf_Status* status );

mhdf_FileHandle
mhdf_createFile( const char* filename, 
                 int overwrite, 
                 const char** elem_type_list,
                 size_t elem_list_len,
                 mhdf_Status* status )
{
  FileHandle* file_ptr;
  unsigned int flags;
  unsigned char index;
  size_t i;
  hid_t enum_id, root_id, group_id;
  int rval;
  API_BEGIN;
  
  if (elem_list_len > 255)
  {
    mhdf_setFail( status, "Element type list too long." );
    return NULL;
  }
  mhdf_setOkay( status );
  
    /* Create struct to hold working data */
  file_ptr = mhdf_alloc_FileHandle( 0, status );
  if (!file_ptr) return NULL;

    /* Create the file */
  flags = overwrite ? H5F_ACC_TRUNC : H5F_ACC_EXCL;
  file_ptr->hdf_handle = H5Fcreate( filename, flags, H5P_DEFAULT, H5P_DEFAULT );
  if (file_ptr->hdf_handle < 0)
  {
    mhdf_setFail( status, "Failed to create file \"%s\"", filename );
    free( file_ptr );
    return NULL;
  }
  
  
    /* Create file structure */
  if (!make_hdf_group(     ROOT_GROUP, file_ptr->hdf_handle, 6, status )
   || !make_hdf_group(      TAG_GROUP, file_ptr->hdf_handle, 0, status )
   || !make_hdf_group(  ELEMENT_GROUP, file_ptr->hdf_handle, 8, status )
   || !make_hdf_group(     NODE_GROUP, file_ptr->hdf_handle, 3, status )
   || !make_hdf_group(      SET_GROUP, file_ptr->hdf_handle, 5, status )
   || !make_hdf_group( NODE_TAG_GROUP, file_ptr->hdf_handle, 0, status )
   || !make_hdf_group(  SET_TAG_GROUP, file_ptr->hdf_handle, 0, status ))
  {
    H5Fclose( file_ptr->hdf_handle );
    free( file_ptr );
    return NULL;
  }
  
    /* Store the max ID as an attribite on the /tstt/ group */
  group_id = H5Gopen( file_ptr->hdf_handle, ROOT_GROUP );
  rval = mhdf_create_scalar_attrib( group_id, 
                                    MAX_ID_ATTRIB, 
                                    H5T_NATIVE_ULONG, 
                                    &file_ptr->max_id,
                                    status );
  H5Gclose( group_id );
  if (!rval)
  {
    H5Fclose( file_ptr->hdf_handle );
    free( file_ptr );
    return NULL;
  }
  
    /* Create the type name list in file */
  enum_id = H5Tenum_create( H5T_NATIVE_UCHAR );
  if (enum_id < 0)
  {
    mhdf_setFail( status, "Failed to store elem type list." );
    H5Fclose( file_ptr->hdf_handle );
    free( file_ptr );
    return NULL;
  }
  for (i = 0; i < elem_list_len; ++i)
  {
    if (!elem_type_list[i] || !*elem_type_list[i])
      continue;
      
    index = (unsigned char)i;
    if ( H5Tenum_insert( enum_id, elem_type_list[i], &index ) < 0)
    {
      mhdf_setFail( status, "Failed to store elem type list." );
      H5Fclose( file_ptr->hdf_handle );
      free( file_ptr );
      return NULL;
    }
  }
  if (H5Tcommit( file_ptr->hdf_handle, TYPE_ENUM_PATH, enum_id ) < 0)  
  {
    mhdf_setFail( status, "Failed to store elem type list." );
    H5Fclose( file_ptr->hdf_handle );
    free( file_ptr );
    return NULL;
  }
  H5Tclose( enum_id );
  
  API_END_H( 1 );
  return file_ptr;
}


mhdf_FileHandle
mhdf_openFile( const char* filename, 
               int writeable, 
               unsigned long* max_id_out,
               mhdf_Status* status )
{
  return mhdf_openFileWithOpt( filename, 
                               writeable, 
                               max_id_out,
                               H5P_DEFAULT, 
                               status );
}


static herr_t get_max_id( hid_t group_id, 
                          const char* subgroup, 
                          const char* datatable,
                          unsigned long* data )
{
  unsigned long id;
  hid_t elem_id, conn_id, attr_id, space_id, poly_id = -1;
  herr_t rval;
  int rank;
  hsize_t dims[2];
  mhdf_Status status;
  
  elem_id = H5Gopen( group_id, subgroup );
  if (elem_id < 0) return (herr_t)-1;
  
  conn_id = H5Dopen( elem_id, datatable );
  H5Gclose( elem_id );
  if (conn_id < 0) return (herr_t)-1;
  
  space_id = H5Dget_space( conn_id );
  if (space_id < 0) { H5Dclose( conn_id ); return -1; }
  
  rank = H5Sget_simple_extent_ndims( space_id );
  if (rank <= 0 || rank > 2) { H5Dclose(conn_id); H5Sclose(space_id); return -1; }
  
  rval = H5Sget_simple_extent_dims( space_id, dims, NULL );
  H5Sclose( space_id );
  if (rval < 0) { H5Dclose( conn_id ); return -1; }
  
  attr_id = H5Aopen_name( conn_id, START_ID_ATTRIB );
  H5Dclose( conn_id );
  if (attr_id < 0) return (herr_t)-1;
  
  rval = H5Aread( attr_id, H5T_NATIVE_ULONG, &id );
  H5Aclose( attr_id );
  if (rval < 0) return rval;
  
  id += dims[0];
  if (id > *data)
    *data = id;
  return 0;
}

static herr_t max_id_iter( hid_t group_id, const char* name, void* data )
{
  return get_max_id( group_id, name, CONNECTIVITY_NAME, (unsigned long*)data );
}

static int
scan_for_max_id( FileHandle* file_ptr, mhdf_Status* status )
{
  hid_t group_id;
  herr_t rval;
  unsigned long tmp;
  
    /* Check for new format, with max_id as attrib of root group */
  group_id = H5Gopen( file_ptr->hdf_handle, ROOT_GROUP );
  if (group_id < 0)
  {
    mhdf_setFail( status, "Internal error - invalid file.");
    return 0;
  }
  if (mhdf_read_scalar_attrib( group_id, MAX_ID_ATTRIB,
                               H5T_NATIVE_ULONG, &file_ptr->max_id,
                               status ))
  {
    H5Gclose( group_id );
    return 1;
  }
  
    /* Didn't find it, scan the elements group */
  rval = H5Giterate( group_id, ELEMENT_GROUP_NAME, 0, &max_id_iter, &file_ptr->max_id );
  if (rval)
  {
    H5Gclose( group_id );
    mhdf_setFail( status, "Internal error -- invalid file." );
    return 0;
  }
  
    /* Check node table too */
  rval = get_max_id( group_id, NODE_GROUP_NAME, "coordinates", &file_ptr->max_id );
  if (rval)
  {
    H5Gclose( group_id );
    mhdf_setFail( status, "Internal error -- invalid file." );
    return 0;
  }
  
    /* Check set table, if it exists */
  rval = mhdf_is_in_group( group_id, SET_GROUP_NAME, status );
  if (rval < 1)
  {
    H5Gclose( group_id );
    return !rval;
  }
  rval = get_max_id( group_id, SET_GROUP_NAME, SET_META_NAME, &file_ptr->max_id );
  H5Gclose( group_id );
  if (rval)
  {
    mhdf_setFail( status, "Internal error -- invalid file." );
    return 0;
  }

  return 1;
}    
  
   
   


mhdf_FileHandle
mhdf_openFileWithOpt( const char* filename, 
                      int writable, 
                      unsigned long* max_id_out,
                      hid_t access_prop,
                      mhdf_Status* status )
{
  FileHandle* file_ptr;
  unsigned int flags;
  size_t size;
  ssize_t tmp;
  hsize_t count, index;
  hid_t group_id, elem_id;
  char* name;
  unsigned long long max_id = 0;
  unsigned long long gp_max_id;
  API_BEGIN;
  
    /* Check if file is HDF5 */
  if (H5Fis_hdf5( filename ) <= 0)
    return NULL;
  
    /* Create struct to hold working data */
  file_ptr = mhdf_alloc_FileHandle( 0, status );
  if (!file_ptr) return NULL;

    /* Create the file */
  flags = writable ? H5F_ACC_RDWR : H5F_ACC_RDONLY;
  file_ptr->hdf_handle = H5Fopen( filename, flags, access_prop );
  if (file_ptr->hdf_handle < 0)
  {
    mhdf_setFail( status, "Failed to open file \"%s\"", filename );
    free( file_ptr );
    return NULL;
  }
  
    /* Check for TSTT data in file */
  group_id = H5Gopen( file_ptr->hdf_handle, ROOT_GROUP );
  if (group_id < 0)
  {
    mhdf_setFail( status, "Invalid file \"%s\"\n", filename );
    H5Fclose( file_ptr->hdf_handle );
    free( file_ptr );
    return NULL;
  }
  H5Gclose( group_id );
  
    /* Get max id */
  if (!scan_for_max_id( file_ptr, status ))
  {
    H5Fclose( file_ptr->hdf_handle );
    free( file_ptr );
    return NULL;
  }
  
  if (max_id_out)
    *max_id_out = max_id;
    
  mhdf_setOkay( status );
  API_END_H(1);
  return file_ptr;
}


void
mhdf_getElemName( mhdf_FileHandle file_handle,
                  unsigned int type_index,
                  char* buffer,
                  size_t buf_size,
                  mhdf_Status* status )
{
  FileHandle* file_ptr;
  herr_t rval;
  hid_t enum_id;
  API_BEGIN;
  
  if (type_index > 255)
  {
    mhdf_setFail( status, "Type index out of bounds." );
    return;
  }
  
  file_ptr = (FileHandle*)(file_handle);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return;

  enum_id = get_elem_type_enum( file_ptr, status );
  if (enum_id < 0)
    return;
  
  rval = H5Tconvert( H5T_NATIVE_UINT, H5Tget_super(enum_id), 1, &index, NULL, H5P_DEFAULT );
  if (rval < 0)
  {
    H5Tclose( enum_id );
    mhdf_setFail( status, "Internal error converting to enum type." );
    return;
  }
  
  rval = H5Tenum_nameof( enum_id, &index, buffer, buf_size );
  H5Tclose( enum_id );
  if (rval < 0)
    mhdf_setFail( status, "H5Tenum_nameof failed.  Invalid type index?" );
  else
    mhdf_setOkay( status );
    
  API_END;
}

void 
mhdf_closeFile( mhdf_FileHandle handle,
                mhdf_Status* status )
{
  FileHandle* file_ptr;
  API_BEGIN;
  
  file_ptr = (FileHandle*)(handle);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return;
/* 
  if (file_ptr->open_handle_count)
  {
    mhdf_setError( status, "Cannot close file with %d open data handles.", 
      file_ptr->open_handle_count );
    return;
  }
*/  
  
  /* Check for open handles.  HDF5 will not actually close the
     file until all handles are closed. */
  if (H5Fget_obj_count( file_ptr->hdf_handle, H5F_OBJ_ALL ) != 1)
  {
    mhdf_setFail( status, "Cannot close file with open handles: "
                 "%d file, %d data, %d group, %d type, %d attr\n",
                  H5Fget_obj_count( file_ptr->hdf_handle, H5F_OBJ_FILE ) - 1,
                  H5Fget_obj_count( file_ptr->hdf_handle, H5F_OBJ_DATASET ),
                  H5Fget_obj_count( file_ptr->hdf_handle, H5F_OBJ_GROUP ),
                  H5Fget_obj_count( file_ptr->hdf_handle, H5F_OBJ_DATATYPE ),
                  H5Fget_obj_count( file_ptr->hdf_handle, H5F_OBJ_ATTR ) );
    return;
  }
 
  if (0 > H5Fclose( file_ptr->hdf_handle ))
  {
    mhdf_setFail( status, "H5FClose failed.  Invalid handle?" );
    return;
  } 
  
  bzero( file_ptr, sizeof(FileHandle) );
  free( file_ptr );
  mhdf_setOkay( status );
  API_END_H( -1 );
}

void
mhdf_closeData( mhdf_FileHandle file, hid_t handle, mhdf_Status* status )
{
  FileHandle* file_ptr;
  herr_t rval = -1;
  
  file_ptr = (FileHandle*)(file);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return;
  
  switch ( H5Iget_type( handle ) )
  {
    case H5I_GROUP    :  rval = H5Gclose( handle );  break;
    case H5I_DATATYPE :  rval = H5Tclose( handle );  break;
    case H5I_DATASPACE:  rval = H5Sclose( handle );  break;
    case H5I_DATASET  :  rval = H5Dclose( handle );  break;
    default           :  rval = -1;
  }
  
  if (rval < 0)
  {
    mhdf_setFail( status, "H5Xclose failed.  Invalid handle?\n");
  }
  else
  {
    file_ptr->open_handle_count--;
    mhdf_setOkay( status );
  }
}


void
mhdf_addElement( mhdf_FileHandle file_handle, 
                 const char* name, 
                 unsigned int elem_type,
                 mhdf_Status* status )
{
  FileHandle* file_ptr = (FileHandle*)file_handle;
  hid_t group_id, tag_id, enum_id;
  char* path, *ptr;
  size_t name_len;
  herr_t rval;
  API_BEGIN;
  
  if (!mhdf_check_valid_file( file_ptr, status ))
    return;
  
  name_len = mhdf_name_to_path( name, NULL, 0 );
  name_len += strlen(ELEMENT_GROUP) + 1;
  path = (char*)mhdf_malloc( name_len, status );
  if (!path)
    return;
  
  strcpy( path, ELEMENT_GROUP );
  ptr = path + strlen(ELEMENT_GROUP);
  if (!mhdf_path_to_name( name, ptr ))
  {
    mhdf_setFail( status, "Invalid character string in internal file path: \"%s\"\n",
      name );
    return;
  }
  
  group_id = H5Gcreate( file_ptr->hdf_handle, path, 3 );
  if (group_id < 0)
  {
    mhdf_setFail( status, "Creation of \"%s\" group failed.\n", path );
    free( path );
    return;
  }
  free( path );
  
  tag_id = H5Gcreate( group_id, DENSE_TAG_SUBGROUP, 0 );
  if (tag_id < 0)
  {
    H5Gclose( group_id );
    mhdf_setFail( status, "Creation of tag subgroup failed.\n" );
    return;
  }
  H5Gclose( tag_id );
  
  enum_id = get_elem_type_enum( file_ptr, status );
  if (enum_id < 0)
  {
    H5Gclose( group_id );
    return;
  }
  
  rval = H5Tconvert( H5T_NATIVE_UINT, H5Tget_super(enum_id), 1, &elem_type, NULL, H5P_DEFAULT );
  if (rval < 0)
  {
    H5Gclose( group_id );
    H5Tclose( enum_id );
    mhdf_setFail( status, "Internal error converting to enum type." );
    return;
  }
  
  rval = mhdf_create_scalar_attrib( group_id, ELEM_TYPE_ATTRIB, enum_id,
                                 &elem_type, status );
  H5Tclose( enum_id );
  if (rval < 0)
  {
    H5Gclose( group_id );
    return;
  }
  
  H5Gclose( group_id );
  mhdf_setOkay( status );
  API_END;
}


char**
mhdf_getElemHandles( mhdf_FileHandle file_handle,
                     unsigned int* count_out,
                     mhdf_Status* status )
{
  hsize_t idx, count, length, i;
  char** buffer;
  char* current;
  hid_t group_id;
  herr_t rval;
  ssize_t rlen;
  size_t remaining;
  FileHandle* file_ptr = (FileHandle*)file_handle;
  if (!mhdf_check_valid_file( file_ptr, status ))
    return NULL;
  
  group_id = H5Gopen( file_ptr->hdf_handle, ELEMENT_GROUP );
  if (group_id < 0) 
  {
    mhdf_setFail( status, "Invalid file -- element group does not exist." );
    return NULL;
  }
  
  rval = H5Gget_num_objs( group_id, &count );
  if (rval < 0) 
  {
    H5Gclose( group_id );
    mhdf_setFail( status, "Internal error calling H5Gget_num_objs." );
    return NULL;
  }
  *count_out = count;
  
  for (i = 0; i < count; ++i)
  {
    rlen += H5Gget_objname_by_idx( group_id, i, NULL, 0 ) + 1;
  }
  
  length = count * sizeof(char*) + rlen;
  buffer = (char**)mhdf_malloc( length, status );
  if (!buffer) { H5Gclose( group_id ); return NULL; }
  current = (char*)(buffer + count);
  remaining = rlen;
  
  for (i = 0; i < count; ++i)
  {
    buffer[i] = current;
    rlen = H5Gget_objname_by_idx( group_id, i, current, remaining ) + 1;
    if (rlen < 0)
    {
      H5Gclose( group_id );
      free( buffer );
      mhdf_setFail( status, "Internal error calling H5Gget_objname_by_idx." );
      return NULL;
    }
    
    mhdf_path_to_name( current, current );
    remaining -= rlen;
    current += rlen;
  }
  
  H5Gclose( group_id );
  mhdf_setOkay( status );
  return buffer;
}


void 
mhdf_getElemTypeName( mhdf_FileHandle file_handle,
                      const char* elem_handle,
                      char* buffer, size_t buf_len,
                      mhdf_Status* status )
{
  FileHandle* file_ptr;
  hid_t elem_id, type_id, attr_id;
  char bytes[16];
  herr_t rval;
  API_BEGIN;
 
  if (NULL == buffer || buf_len < 2)
  {
    mhdf_setFail( status, "invalid input" );
    return;
  }
  buffer[0] = '\0';
  
  file_ptr = (FileHandle*)(file_handle);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return;
  
  elem_id = mhdf_elem_group_from_handle( file_ptr, elem_handle, status );
  if (elem_id < 0) return;
  
  attr_id = H5Aopen_name( elem_id, ELEM_TYPE_ATTRIB );
  H5Gclose( elem_id );
  if (attr_id < 0)
  {
    mhdf_setFail( status, "Missing element type attribute.  Invalid file." );
    return;
  }
  
  type_id = H5Aget_type( attr_id );
  assert( type_id > 0 );
  
  rval = H5Aread( attr_id, type_id, bytes );
  H5Aclose( attr_id );
  if (rval < 0)
  {
    H5Tclose( type_id );
    mhdf_setFail( status, "Failed to read element type attribute.  Invalid file." );
    return;
  }
  
  rval = H5Tenum_nameof( type_id, bytes, buffer, buf_len );
  H5Tclose( type_id );
  if (rval < 0)
  {
    mhdf_setFail( status, "Invalid datatype for element type attribute.  Invalid file." );
    return;
  }
  
  mhdf_setOkay( status );  
  API_END;
  return ;
}


static int
make_hdf_group( const char* path, hid_t file, size_t size, mhdf_Status* status )
{
  hid_t handle = H5Gcreate( file, path, size );
  if (handle < 0)
  {
    mhdf_setFail( status, "Failed to create \"%s\" group.", path );
    return 0;
  }
  else
  {
    H5Gclose( handle );
    return 1;
  }
}

const char* 
mhdf_node_type_handle(void)
{
  static const char rval[] = "nodes";
  return rval;
}

const char*
mhdf_set_type_handle(void)
{
  static const char rval[] = "sets";
  return rval;
}

int
mhdf_isPolyElement( mhdf_FileHandle file_handle,
                    const char* elem_handle,
                    mhdf_Status* status )
{
  FileHandle* file_ptr;
  hid_t elem_id;
  int rval;
  API_BEGIN;
  
  file_ptr = (FileHandle*)(file_handle);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;
  
  elem_id = mhdf_elem_group_from_handle( file_ptr, elem_handle, status );
  if (elem_id < 0) return -1;
  
  mhdf_setOkay( status );
  rval = mhdf_is_in_group( elem_id, POLY_INDEX_NAME, status );
  H5Gclose( elem_id );
  API_END;
  return rval;
}

void
mhdf_writeHistory( mhdf_FileHandle file_handle, 
                   const char** strings, 
                   int num_strings,
                   mhdf_Status* status )
{
  FileHandle* file_ptr;
  hid_t data_id, type_id, space_id;
  hsize_t dim = (hsize_t)num_strings;
  herr_t rval;
  API_BEGIN;
  
  file_ptr = (FileHandle*)(file_handle);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return;
    
  type_id = H5Tcopy( H5T_C_S1 );
  if (type_id < 0 || H5Tset_size( type_id, H5T_VARIABLE ) < 0)
  {
    if (type_id >= 0) H5Tclose(type_id);
    mhdf_setFail( status, "Could not create variable length string type." );
    return;
  }
  
  space_id = H5Screate_simple( 1, &dim, NULL );
  if (space_id < 0)
  {
    H5Tclose( type_id );
    mhdf_setFail( status, "H5Screate_simple failed." );
    return;
  }
  
  data_id = H5Dcreate( file_ptr->hdf_handle, HISTORY_PATH, type_id, space_id, H5P_DEFAULT );
  H5Sclose( space_id );
  if (data_id < 0)
  {
    H5Tclose( type_id );
    mhdf_setFail( status, "Failed to create \"%s\".", HISTORY_PATH );
    return;
  }
    
  rval = H5Dwrite( data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, strings );
  H5Dclose( data_id );
  H5Tclose( type_id );
  if (rval < 0)
  {
    H5Gunlink( file_ptr->hdf_handle, HISTORY_PATH );
    mhdf_setFail( status, "Failed to write history data." );
    return;
  }
  
  mhdf_setOkay( status );
  API_END;
}

char**
mhdf_readHistory( mhdf_FileHandle file_handle, 
                  int* num_strings,
                  mhdf_Status* status )
{
  FileHandle* file_ptr;
  hid_t data_id, type_id, space_id, group_id;
  hsize_t dim;
  herr_t rval;
  char** array;
  API_BEGIN;
  
  file_ptr = (FileHandle*)(file_handle);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return NULL;
  
    /* check if file contains history data */
  group_id = H5Gopen( file_ptr->hdf_handle, ROOT_GROUP );
  if (group_id < 0)
  {
    mhdf_setFail( status, "Could not open root group.  Invalid file." );
    return NULL;
  }
  
  rval = mhdf_is_in_group( group_id, HISTORY_NAME, status );
  if (rval < 1)
  {
    H5Gclose( group_id );
    *num_strings = 0;
    if (0 == rval)
      mhdf_setOkay( status );
    return NULL;
  }
  
  data_id = H5Dopen( group_id, HISTORY_NAME );
  H5Gclose( group_id );
  if (data_id < 0)
  {
    mhdf_setFail( status, "Failed to open \"%s\".", HISTORY_PATH );
    return NULL;
  }
  
  space_id = H5Dget_space( data_id );
  if (space_id < 0)
  {
    H5Dclose( data_id );
    mhdf_setFail( status, "Internal error calling H5Dget_space.");
    return NULL;
  }
  
  if (1 != H5Sget_simple_extent_ndims( space_id ) ||
      1 != H5Sget_simple_extent_dims( space_id, &dim, NULL ))
  {
    H5Dclose( data_id );
    mhdf_setFail( status, "Invalid dimension for \"%s\".", HISTORY_PATH );
    return NULL;
  }
  H5Sclose( space_id );
  
  if (0 == dim)
  {
    H5Dclose( data_id );
    *num_strings = 0;
    mhdf_setOkay( status );
    return NULL;
  }
  
  array = (char**)mhdf_malloc( dim * sizeof(char*), status );
  if (!array)
  {
    H5Dclose( data_id );
    return NULL;
  }
    
  type_id = H5Tcopy( H5T_C_S1 );
  if (type_id < 0 || H5Tset_size( type_id, H5T_VARIABLE ) < 0)
  {
    H5Dclose( data_id );
    if (type_id >= 0) H5Tclose(type_id);
    mhdf_setFail( status, "Could not create variable length string type." );
    free( array );
    return NULL;
  }
  
  rval = H5Dread( data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, array );
  H5Tclose( type_id );
  H5Dclose( data_id );
  if (rval < 0)
  {
    free( array );
    mhdf_setFail( status, "H5Dread failed." );
    return NULL;
  }
   
  *num_strings = dim;
  mhdf_setOkay( status );
  API_END;
  return array;
}
