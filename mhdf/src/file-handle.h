#ifndef MHDF_FILE_HANDLE_H
#define MHDF_FILE_HANDLE_H

#include "mhdf.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct struct_FileHandle {
  uint32_t magic;
  hid_t hdf_handle;
  int open_handle_count;

  hid_t* elem_type_handles;
  int num_type_handles;
  int type_array_len;
  
  long max_id;
} FileHandle;

FileHandle* mhdf_alloc_FileHandle( hid_t hdf_handle, mhdf_Status* status );

int mhdf_check_valid_file( FileHandle* handle, mhdf_Status* status );

int mhdf_add_elem_type( FileHandle* handle, 
                        hid_t hdf_handle, 
                        mhdf_Status* status );

hid_t mhdf_handle_from_type_index( FileHandle* handle, 
                                   int index,
                                   mhdf_Status* status );

#ifdef __cplusplus
} // extern "C"
#endif

#endif
