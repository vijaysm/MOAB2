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
  hid_t enum_id;
  API_BEGIN;
  
  if (elem_list_len > 255)
  {
    mhdf_setFail( status, "Element type list too long." );
    return NULL;
  }
  
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


static long
mhdf_get_max_elem_id( hid_t elem_id, mhdf_Status* status )
{
  hid_t table_id, space_id;
  hsize_t dims[2];
  int rank;
  long result;
  int ispoly;

    /* Get connectivity table for element */
  table_id = H5Dopen( elem_id, CONNECTIVITY_NAME );
  if (table_id < 0)
  {
    mhdf_setFail( status, "Missing connectivity data for element." );
    return 0L;
  }
  
  space_id = H5Dget_space( table_id );
  if (space_id < 0)
  {
    mhdf_setFail( status, "Internal error in H5Dget_space.");
    H5Dclose( table_id );
    return 0L;
  }
  
  ispoly = mhdf_is_in_group( elem_id, POLY_INDEX_NAME, status );
  if (ispoly < 0)
  {
    H5Sclose( space_id );
    H5Dclose( table_id );
    return 0L;
  }
  
  rank = H5Sget_simple_extent_ndims( space_id );
  if (rank != (ispoly ? 1 : 2))
  {
    if (rank < 0)
      mhdf_setFail( status, "Internal error in H5Sget_simple_extent_ndims." );
    else if (ispoly)
      mhdf_setFail( status, "Invalid file: Poly(gon|hedron) connectivity with rank %d.\n", rank );
    else
      mhdf_setFail( status, "Invalid file: Element connectivity data has rank %d.\n", rank );
    H5Sclose( space_id );
    H5Dclose( table_id );
    return 0L;
  }
  
  rank = H5Sget_simple_extent_dims( space_id, dims, NULL );
  H5Sclose( space_id );
  if (rank < 0)
  {
    mhdf_setFail( status, "Internal error calling H5Sget_simple_extent_dims.");
    H5Dclose( table_id );
    return 0L;
  }
  
  rank = mhdf_read_scalar_attrib( table_id, START_ID_ATTRIB, 
                                  H5T_NATIVE_LONG, &result, 
                                  status );
  H5Dclose( table_id );
  if (!rank)
  {
    return 0L;
  }
  if (result < 1)
  {
    mhdf_setFail( status, "Invalid start ID for element data: %ld\n", result );
    return 0L;
  }
  
  result += dims[0];
  return result;
}


mhdf_FileHandle
mhdf_openFile( const char* filename, int writable, mhdf_Status* status )
{
  FileHandle* file_ptr;
  unsigned int flags;
  size_t size;
  ssize_t tmp;
  hsize_t count, index;
  hid_t group_id, elem_id;
  char* name;
  long max_id;
  API_BEGIN;
  
    /* Check if file is HDF5 */
  if (H5Fis_hdf5( filename ) <= 0)
    return NULL;
  
    /* Create struct to hold working data */
  file_ptr = mhdf_alloc_FileHandle( 0, status );
  if (!file_ptr) return NULL;

    /* Create the file */
  flags = writable ? H5F_ACC_RDWR : H5F_ACC_RDONLY;
  file_ptr->hdf_handle = H5Fopen( filename, flags, H5P_DEFAULT );
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
  
    /* Check for subgroups */
  mhdf_setFail( status, "Invalid file." );
  if (1 > mhdf_is_in_group( group_id, NODE_GROUP_NAME, status ) ||
      1 > mhdf_is_in_group( group_id, "elements", status ) ||
      1 > mhdf_is_in_group( group_id, TAG_GROUP_NAME, status) ||
      1 > mhdf_is_in_group( group_id, SET_GROUP_NAME, status))
  {
    H5Gclose( group_id );
    H5Fclose( file_ptr->hdf_handle );
    free( file_ptr );
    return NULL;
  }
  H5Gclose( group_id );
  
    /* Get element group */
  group_id = H5Gopen( file_ptr->hdf_handle, ELEMENT_GROUP );
  if (group_id < 0)
  {
    H5Fclose( file_ptr->hdf_handle );
    free( file_ptr );
    return NULL;
  }
  
    /* Get number of entries in element group */
  if (H5Gget_num_objs( group_id, &count ) < 0)
  {
    H5Gclose( group_id );
    H5Fclose( file_ptr->hdf_handle );
    free( file_ptr );
    return NULL;
  }
  
    /* Pre-allocate array for element handles, as the size is known */
  file_ptr->type_array_len = (int)count;
  file_ptr->num_type_handles = 0;
  file_ptr->elem_type_handles = (hid_t*)mhdf_malloc( count * sizeof(hid_t), status );
  if (NULL == file_ptr->elem_type_handles)
  {
    H5Gclose( group_id );
    H5Fclose( file_ptr->hdf_handle );
    free( file_ptr );
    return NULL;
  }
  
    /* Determine the maximum name length for element types,
       and make sure each entry is a group */
  size = 0;
  for (index = 0; index < count; ++index)
  {
    if (H5Gget_objtype_by_idx( group_id, index ) != H5G_GROUP)
    {
      H5Gclose( group_id );
      H5Fclose( file_ptr->hdf_handle );
      free( file_ptr );
      return NULL;
    }
    
    tmp = H5Gget_objname_by_idx( group_id, index, NULL, 0 );
    if (tmp < 1)
    {
      H5Gclose( group_id );
      H5Fclose( file_ptr->hdf_handle );
      free( file_ptr );
      return NULL;
    }
    
    if ((unsigned)tmp > size)
      size = tmp;
  }
  
    /* Allocate buffer to hold element names */
  ++size;
  name = (char*)mhdf_malloc( size, status );
  if (!name)
  {
    H5Gclose( group_id );
    H5Fclose( file_ptr->hdf_handle );
    free( file_ptr );
    return NULL;
  }
  
    /* Open each entry in the elements group by name, check for 
       the existance of a connectivity table and if one exists,
       aquire the ID range for the elements from said table. */
  mhdf_setFail( status, "Internal error traversing elem group contents." );
  for (index = 0; index < count; ++index)
  {
      /* Open element group */
    tmp = H5Gget_objname_by_idx( group_id, index, name, size );
    if (tmp < 1 || (unsigned)tmp > size  || 
        (elem_id = H5Gopen( group_id, name )) < 0)
    {
      free( name );
      H5Gclose( group_id );
      mhdf_closeFile( file_ptr, status );
      mhdf_setFail( status, "Invalid file." );
      return NULL;
    }
    
      /* Add element type to FileHandle */
    if (0 > mhdf_add_elem_type( file_ptr, elem_id, status ))
    {
      free( name );
      H5Gclose( group_id );
      mhdf_closeFile( file_ptr, status );
      mhdf_setFail( status, "Invalid file." );
      return NULL;
    }
    
    max_id = mhdf_get_max_elem_id( elem_id, status );
    if (max_id <= 0)
    {
      free( name );
      H5Gclose( group_id );
      mhdf_closeFile( file_ptr, NULL );
      return NULL;
    }
    
    if (max_id > file_ptr->max_id)
      file_ptr->max_id = max_id;
  }
  
  free( name );
  H5Gclose( group_id );
  mhdf_setOkay( status );
  API_END_H(file_ptr->num_type_handles+1);
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
  int i;
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
  
  for (i = 0; i < file_ptr->num_type_handles; i++)
  {
    H5Gclose( file_ptr->elem_type_handles[i] );
  }
  
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
  
  free( file_ptr->elem_type_handles );
  bzero( file_ptr, sizeof(FileHandle) );
  free( file_ptr );
  mhdf_setOkay( status );
  API_END_H( -i - 1 );
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


mhdf_ElemHandle
mhdf_addElement( mhdf_FileHandle file_handle, 
                 const char* name, 
                 unsigned int elem_type,
                 mhdf_Status* status )
{
  FileHandle* file_ptr = (FileHandle*)file_handle;
  hid_t group_id, tag_id, enum_id;
  char* path, *ptr;
  size_t name_len;
  mhdf_ElemHandle result;
  herr_t rval;
  API_BEGIN;
  
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;
  
  name_len = mhdf_name_to_path( name, NULL, 0 );
  name_len += strlen(ELEMENT_GROUP) + 1;
  path = (char*)mhdf_malloc( name_len, status );
  if (!path)
    return -1;
  
  strcpy( path, ELEMENT_GROUP );
  ptr = path + strlen(ELEMENT_GROUP);
  if (!mhdf_path_to_name( name, ptr ))
  {
    mhdf_setFail( status, "Invalid character string in internal file path: \"%s\"\n",
      name );
    return -1;
  }
  
  group_id = H5Gcreate( file_ptr->hdf_handle, path, 3 );
  if (group_id < 0)
  {
    mhdf_setFail( status, "Creation of \"%s\" group failed.\n", path );
    free( path );
    return -1;
  }
  free( path );
  
  tag_id = H5Gcreate( group_id, DENSE_TAG_SUBGROUP, 0 );
  if (tag_id < 0)
  {
    H5Gclose( group_id );
    mhdf_setFail( status, "Creation of tag subgroup failed.\n" );
    return -1;
  }
  H5Gclose( tag_id );
  
  enum_id = get_elem_type_enum( file_ptr, status );
  if (enum_id < 0)
  {
    H5Gclose( group_id );
    return -1;
  }
  
  rval = H5Tconvert( H5T_NATIVE_UINT, H5Tget_super(enum_id), 1, &elem_type, NULL, H5P_DEFAULT );
  if (rval < 0)
  {
    H5Gclose( group_id );
    H5Tclose( enum_id );
    mhdf_setFail( status, "Internal error converting to enum type." );
    return -1;
  }
  
  rval = mhdf_write_scalar_attrib( group_id, ELEM_TYPE_ATTRIB, enum_id,
                                 &elem_type, status );
  H5Tclose( enum_id );
  if (rval < 0)
  {
    H5Gclose( group_id );
    return -1;
  }
  
  result = mhdf_add_elem_type( file_ptr, group_id, status );
  if (result < 0)
  {
    H5Gclose( group_id );
    return -1;
  }
  
  mhdf_setOkay( status );
  API_END_H(1);
  return result;
}


int
mhdf_numElemGroups( mhdf_FileHandle file_handle, mhdf_Status* status )
{
  FileHandle* file_ptr = (FileHandle*)file_handle;
  
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;
  
  return file_ptr->num_type_handles;
}


void
mhdf_getElemGroups( mhdf_FileHandle file_handle, 
                    mhdf_ElemHandle* out, 
                    mhdf_Status* status )
{
  mhdf_ElemHandle i;
  FileHandle* file_ptr = (FileHandle*)file_handle;
  
  if (!mhdf_check_valid_file( file_ptr, status ))
    return;
  
  for (i = 0; i < file_ptr->num_type_handles; ++i)
    out[i] = i;
}


void 
mhdf_getElemTypeName( mhdf_FileHandle file_handle,
                      mhdf_ElemHandle elem_handle,
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
  
  elem_id = mhdf_handle_from_type_index( file_ptr, elem_handle, status );
  if (elem_id < 0) return;
  
  attr_id = H5Aopen_name( elem_id, ELEM_TYPE_ATTRIB );
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

mhdf_ElemHandle 
mhdf_node_type_handle(void)
{
  return 1 << 15;
}

mhdf_ElemHandle
mhdf_set_type_handle(void)
{
  return (1 << 15) - 1;
}

int
mhdf_isPolyElement( mhdf_FileHandle file_handle,
                    mhdf_ElemHandle elem_handle,
                    mhdf_Status* status )
{
  FileHandle* file_ptr;
  hid_t elem_id;
  int rval;
  API_BEGIN;
  
  file_ptr = (FileHandle*)(file_handle);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;
  
  elem_id = mhdf_handle_from_type_index( file_ptr, elem_handle, status );
  if (elem_id < 0) return -1;
  
  mhdf_setOkay( status );
  rval = mhdf_is_in_group( elem_id, POLY_INDEX_NAME, status );
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
