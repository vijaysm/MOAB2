#ifndef mhdf_FILE_UTIL_H
#define mhdf_FILE_UTIL_H

#include <sys/types.h>
#include <H5Ipublic.h>
#include "status.h"
#include "file-handle.h"

#ifdef __cplusplus
extern "C" {
#endif


void* mhdf_malloc( size_t size, mhdf_Status* status );

size_t mhdf_name_to_path( const char* name, char* path, size_t path_len );

int mhdf_path_to_name( const char* path, char* name );

char* mhdf_name_to_path_copy( const char* name, mhdf_Status* status );

int mhdf_write_scalar_attrib( hid_t object,
                              const char* name,
                              hid_t type,
                              const void* value,
                              mhdf_Status* status );

/* If type is zero, assumes opaque type.
   On error, sets status and returns zero.
   On success, returns non-zero and does not modify status */
int mhdf_read_scalar_attrib( hid_t object,
                             const char* name,
                             hid_t type,
                             void* value,
                             mhdf_Status* status );

/* Search the specified object to see if it contains an 
   an attribute with passed name.  Returns -1 on error, 1 
   if attribute was found, and zero if attribute was not 
   found.
*/
int mhdf_find_attribute( hid_t object, 
                         const char* attrib_name,
                         unsigned int* index_out,
                         mhdf_Status* status );

int mhdf_is_in_group( hid_t group, const char* name, mhdf_Status* status );

int
mhdf_read_data( hid_t data_table,
                long offset,
                long count,
                hid_t type,
                void* array,
                mhdf_Status* status );


int
mhdf_write_data( hid_t data_table,
                 long offset,
                 long count,
                 hid_t type,
                 const void* array,
                 mhdf_Status* status );

int
mhdf_read_column( hid_t data_table,
                  int column,
                  long offset,
                  long count,
                  hid_t type,
                  void* array,
                  mhdf_Status* status );


int
mhdf_write_column( hid_t data_table,
                   int column,
                   long offset,
                   long count,
                   hid_t type,
                   const void* array,
                   mhdf_Status* status );

hid_t 
mhdf_create_table( hid_t group,
                   const char* path,
                   hid_t type,
                   int rank,
                   hsize_t* dims,
                   mhdf_Status* status );

hid_t
mhdf_open_table( hid_t group,
                 const char* path,
                 int columns,
                 hsize_t* rows_out,
                 mhdf_Status* status );

hid_t
mhdf_open_table2( hid_t group,
                  const char* path,
                  int rank,
                  hsize_t* dims_out,
                  long* start_id_out,
                  mhdf_Status* status );

int
mhdf_compact_to_ranges( int* length_in_out, int* ids_in, int ordered );

hid_t 
get_elem_type_enum( FileHandle* file_ptr, mhdf_Status* status );


#ifdef __cplusplus
} // extern "C"
#endif

#endif
 
