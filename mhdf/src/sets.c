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

#include <H5Tpublic.h>
#include <H5Dpublic.h>
#include <H5Gpublic.h>
#include "mhdf.h"
#include "util.h"
#include "file-handle.h"
#include "status.h"
#include "names-and-paths.h"

int
mhdf_haveSets( mhdf_FileHandle file,
               int* have_data,
               int* have_child,
               mhdf_Status* status )
{
  FileHandle* file_ptr = (FileHandle*)file;
  hid_t root_id, set_id;
  int result;
  API_BEGIN;
  
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;
  
  root_id = H5Gopen( file_ptr->hdf_handle, ROOT_GROUP );
  if (root_id < 0)
  {
    mhdf_setFail( status, "H5Gopen( \"%s\" ) failed.", ROOT_GROUP );
    return -1;
  }
  
  result = mhdf_is_in_group( root_id, SET_GROUP_NAME, status );
  if (result < 1)
  {
    H5Gclose(root_id);
    return result;
  }
  
  set_id = H5Gopen( root_id, SET_GROUP_NAME );
  H5Gclose( root_id );
  if (set_id < 0)
  {
    mhdf_setFail( status, "H5Gopen( \"%s\" ) failed.", SET_GROUP );
    return -1;
  }
  
  result = mhdf_is_in_group( set_id, SET_META_NAME, status );
  if (result < 0)
  {
    H5Gclose(set_id);
    return result;
  }
  
  *have_data = mhdf_is_in_group( set_id, SET_DATA_NAME, status );
  if (*have_data < 0)
  {
    H5Gclose(set_id);
    return *have_data;
  }
  
  *have_child = mhdf_is_in_group( set_id, SET_CHILD_NAME, status );
  if (*have_child < 0)
  {
    H5Gclose(set_id);
    return *have_child;
  }
  
  mhdf_setOkay( status );
  H5Gclose( set_id );
  API_END;
  return result;
}
  

hid_t
mhdf_createSetMeta( mhdf_FileHandle file,
                    long num_sets,
                    long* first_id_out,
                    mhdf_Status* status )
{
  FileHandle* file_ptr = (FileHandle*)file;
  hid_t table_id;
  hsize_t dims[2];
  long first_id;
  API_BEGIN;
  
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;
  
  dims[0] = (hsize_t)num_sets;
  dims[1] = (hsize_t)3;
  table_id = mhdf_create_table( file_ptr->hdf_handle,
                                SET_META_PATH,
                                H5T_NATIVE_LONG,
                                2, dims,
                                status );
  if (table_id < 0)
    return -1;
  
  first_id = file_ptr->max_id + 1;
  if (!mhdf_write_scalar_attrib( table_id, 
                                 START_ID_ATTRIB, 
                                 H5T_NATIVE_LONG,
                                 &first_id,
                                 status ))
  {
    H5Dclose( table_id );
    return -1;
  }
  
  *first_id_out = first_id;
  file_ptr->max_id += num_sets;
  file_ptr->open_handle_count++;
  mhdf_setOkay( status );
 
  API_END_H(1);
  return table_id;
}
  

hid_t
mhdf_openSetMeta( mhdf_FileHandle file,
                  long* num_sets,
                  long* first_id_out,
                  mhdf_Status* status )
{
  FileHandle* file_ptr = (FileHandle*)file;
  hid_t table_id;
  hsize_t dims[2];
  API_BEGIN;
  
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;
  
  table_id = mhdf_open_table2( file_ptr->hdf_handle,
                               SET_META_PATH, 2,
                               dims, first_id_out, status );
  if (table_id < 0)
    return -1;
  
  if (dims[1] != 3)
  { 
    mhdf_setFail( status, "Invalid format for meshset table.\n" );
    H5Dclose( table_id );
    return -1;
  }

 
  *num_sets = dims[0];
  file_ptr->open_handle_count++;
  mhdf_setOkay( status );
  API_END_H(1);
  return table_id;
}




void
mhdf_readSetMeta( hid_t table_id,
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

void
mhdf_writeSetMeta( hid_t table_id,
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

                  

hid_t
mhdf_createSetData( mhdf_FileHandle file_handle,
                    long data_list_size,
                    mhdf_Status* status )
{
  FileHandle* file_ptr;
  hid_t table_id;
  hsize_t dim = (hsize_t)data_list_size;
  API_BEGIN;
  
  file_ptr = (FileHandle*)(file_handle);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;

  if (data_list_size < 1)
  {
    mhdf_setFail( status, "Invalid argument.\n" );
    return -1;
  }
  
  table_id = mhdf_create_table( file_ptr->hdf_handle,
                                SET_DATA_PATH,
                                H5T_NATIVE_LONG,
                                1, &dim,
                                status );
  
  API_END_H(1);
  return table_id;
}

hid_t
mhdf_openSetData( mhdf_FileHandle file_handle,
                  long* data_list_size_out,
                  mhdf_Status* status )
{
  FileHandle* file_ptr;
  hid_t table_id;
  hsize_t dim;
  API_BEGIN;
  
  file_ptr = (FileHandle*)(file_handle);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;

  if (!data_list_size_out)
  {
    mhdf_setFail( status, "Invalid argument.\n" );
    return -1;
  }
  
  table_id = mhdf_open_table( file_ptr->hdf_handle,
                              SET_DATA_PATH,
                              1, &dim,
                              status );
 
  *data_list_size_out = (long)dim;
  API_END_H(1);
  return table_id;
}


void
mhdf_writeSetData( hid_t table_id,
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
mhdf_readSetData( hid_t table_id,
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

hid_t
mhdf_createSetChildren( mhdf_FileHandle file_handle,
                        long child_list_size,
                        mhdf_Status* status )
{
  FileHandle* file_ptr;
  hid_t table_id;
  hsize_t dim = (hsize_t)child_list_size;
  API_BEGIN;
  
  file_ptr = (FileHandle*)(file_handle);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;

  if (child_list_size < 1)
  {
    mhdf_setFail( status, "Invalid argument.\n" );
    return -1;
  }
  
  table_id = mhdf_create_table( file_ptr->hdf_handle,
                                SET_CHILD_PATH,
                                H5T_NATIVE_LONG,
                                1, &dim,
                                status );
  
  API_END_H(1);
  return table_id;
}

hid_t
mhdf_openSetChildren( mhdf_FileHandle file_handle,
                      long* child_list_size,
                      mhdf_Status* status )
{
  FileHandle* file_ptr;
  hid_t table_id;
  hsize_t dim;
  API_BEGIN;
  
  file_ptr = (FileHandle*)(file_handle);
  if (!mhdf_check_valid_file( file_ptr, status ))
    return -1;

  if (!child_list_size)
  {
    mhdf_setFail( status, "Invalid argument.\n" );
    return -1;
  }
  
  table_id = mhdf_open_table( file_ptr->hdf_handle,
                              SET_CHILD_PATH,
                              1, &dim,
                              status );
 
  *child_list_size = (long)dim;
  API_END_H(1);
  return table_id;
}

void
mhdf_writeSetChildren( hid_t table_id,
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
mhdf_readSetChildren( hid_t table_id,
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
