#include <H5Tpublic.h>
#include <H5Gpublic.h>
#include "mhdf.h"
#include "util.h"
#include "file-handle.h"
#include "status.h"
#include "names-and-paths.h"

int
mhdf_haveAdjacency( mhdf_FileHandle file,
                    mhdf_ElemHandle elem_group,
                    mhdf_Status* status )
{
  FileHandle* file_ptr;
  hid_t elem_id;
  int result;
  API_BEGIN;
  
  file_ptr = (FileHandle*)(file);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;
  
  if (elem_group == mhdf_node_type_handle())
  {
    elem_id = H5Gopen( file_ptr->hdf_handle, NODE_GROUP );
    if (elem_id < 0)
    {
      mhdf_setFail( status, "H5Gopen( \"%s\" ) failed.\n", NODE_GROUP );
      return -1;
    }
  }
  else
  {
    elem_id = mhdf_handle_from_type_index( file_ptr, elem_group, status );
    if (elem_id < 0)
      return -1;
  }
  
  result = mhdf_is_in_group( elem_id, ADJACENCY_NAME, status );
  if (elem_group == mhdf_node_type_handle())
    H5Gclose( elem_id );
  
  mhdf_setOkay( status );
  API_END;
  return result;
}


hid_t
mhdf_createAdjacency( mhdf_FileHandle file,
                      mhdf_ElemHandle element_handle,
                      long adj_list_size,
                      mhdf_Status* status )
{
  FileHandle* file_ptr;
  hid_t elem_id, table_id;
  hsize_t dim = (hsize_t)adj_list_size;
  API_BEGIN;
  
  file_ptr = (FileHandle*)(file);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;

  if (adj_list_size < 1)
  {
    mhdf_setFail( status, "Invalid argument.\n" );
    return -1;
  }
  
  if (element_handle == mhdf_node_type_handle())
  {
    table_id = mhdf_create_table( file_ptr->hdf_handle,
                                  NODE_ADJCY_PATH,
                                  H5T_NATIVE_LONG,
                                  1, &dim,
                                  status );
  }
  else
  {
    elem_id = mhdf_handle_from_type_index( file_ptr, element_handle, status );
    if (elem_id < 0)
      return -1;
    
    table_id = mhdf_create_table( elem_id,
                                  ADJACENCY_NAME,
                                  H5T_NATIVE_LONG,
                                  1, &dim,
                                  status );
  }
  
  API_END_H(1);
  return table_id;
}


  
hid_t
mhdf_openAdjacency( mhdf_FileHandle file,
                    mhdf_ElemHandle element_handle,
                    long* adj_list_size_out,
                    mhdf_Status* status )

{
  FileHandle* file_ptr;
  hid_t elem_id, table_id;
  hsize_t dim;
  API_BEGIN;
  
  file_ptr = (FileHandle*)(file);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;

  if (!adj_list_size_out)
  {
    mhdf_setFail( status, "Invalid argument.\n" );
    return -1;
  }
  
  if (element_handle == mhdf_node_type_handle())
  {
    table_id = mhdf_open_table( file_ptr->hdf_handle,
                                NODE_ADJCY_PATH,
                                1, &dim,
                                status );
  }
  else
  {
    elem_id = mhdf_handle_from_type_index( file_ptr, element_handle, status );
    if (elem_id < 0)
      return -1;
    
    table_id = mhdf_open_table( elem_id,
                                ADJACENCY_NAME,
                                1, &dim,
                                status );
  }
  
  *adj_list_size_out = (long)dim;
  API_END_H(1);
  return table_id;
}

void
mhdf_writeAdjacency( hid_t table_id,
                     long offset,
                     long count,
                     hid_t type,
                     const void* data,
                     mhdf_Status* status )
{
  API_BEGIN;
  mhdf_write_data( table_id, offset, count, type, data, status );
  API_END;
}

void
mhdf_readAdjacency( hid_t table_id,
                    long offset,
                    long count,
                    hid_t type,
                    void* data,
                    mhdf_Status* status )
{
  API_BEGIN;
  mhdf_read_data( table_id, offset, count, type, data, status );
  API_END;
}



