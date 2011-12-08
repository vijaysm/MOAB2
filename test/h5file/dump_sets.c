#include "mhdf.h"
#include <stdio.h>
#include <stdlib.h>
#include <H5Tpublic.h>

#define CHECK \
  do { if (mhdf_isError(&status)) { \
    fprintf(stderr,"mhdf error at %s:%d : \"%s\"\n", __FILE__, __LINE__, \
      mhdf_message(&status)); \
    exit(2); \
  } } while(0)

void dump_sets( mhdf_FileHandle file );
void dump_set_contents( long first,
                        long num_sets,
                        const int* flags,
                        const long* end_indices,
                        const long* contents );
void dump_set_list( long id, const long* array, long count );
void dump_set_ranges( long id, const long* array, long count );

int WIDTH;

int main( int argc, char* argv[] )
{
  mhdf_Status status;
  mhdf_FileHandle file;
  unsigned long max_id;

  if (argc != 2) {
    fprintf(stderr,"Usage: %s <input_file>\n", argv[0]);
    return 1;
  }
  
  file = mhdf_openFile( argv[1], 0, &max_id, -1, &status ); CHECK;
  
  WIDTH = 1;
  while (max_id >= 10) {
    max_id /= 10;
    ++WIDTH;
  }
  dump_sets( file );
  mhdf_closeFile( file, &status ); CHECK;
  return 0;
}

void dump_sets( mhdf_FileHandle file )
{
  long *end_indices, *data;
  int* flags;
  long num_sets, first, num_data;
  hid_t meta, handle;
  mhdf_Status status;
  
  meta = mhdf_openSetMeta( file, &num_sets, &first, &status ); CHECK;
  end_indices = malloc( num_sets * sizeof(long) );

  printf("\nSet IDs: %ld - %ld\n", first, first+num_sets - 1);

  puts("\nSets:\n");
  flags = malloc( num_sets * sizeof(int) );
  mhdf_readSetFlags( meta, 0, num_sets, H5T_NATIVE_INT, flags, &status ); CHECK;
  mhdf_readSetContentEndIndices( meta, 0, num_sets, H5T_NATIVE_LONG, end_indices, &status ); CHECK;
  handle = mhdf_openSetData( file, &num_data, &status ); CHECK;
  data = malloc( num_data * sizeof(long) ); 
  mhdf_readSetData( handle, 0, num_data, H5T_NATIVE_LONG, data, &status ); CHECK;
  mhdf_closeData( file, handle, &status ); CHECK;
  dump_set_contents( first, num_sets, flags, end_indices, data );
  free(flags);
  
  puts("\nSet Children:\n");
  mhdf_readSetChildEndIndices( meta, 0, num_sets, H5T_NATIVE_LONG, end_indices, &status ); CHECK;
  handle = mhdf_openSetChildren( file, &num_data, &status ); CHECK;
  data = realloc( data, num_data * sizeof(long) );
  mhdf_readSetParentsChildren( handle, 0, num_data, H5T_NATIVE_LONG, data, &status); CHECK;
  mhdf_closeData( file, handle, &status ); CHECK;
  dump_set_contents( first, num_sets, NULL, end_indices, data );
  
  puts("\nSet Parents:\n");
  mhdf_readSetParentEndIndices( meta, 0, num_sets, H5T_NATIVE_LONG, end_indices, &status ); CHECK;
  handle = mhdf_openSetParents( file, &num_data, &status ); CHECK;
  data = realloc( data, num_data * sizeof(long) );
  mhdf_readSetParentsChildren( handle, 0, num_data, H5T_NATIVE_LONG, data, &status); CHECK;
  mhdf_closeData( file, handle, &status ); CHECK;
  dump_set_contents( first, num_sets, NULL, end_indices, data );

  free(end_indices);
  free(data);
  mhdf_closeData( file, meta, &status ); CHECK;
}

void dump_set_contents( long first,
                        long num_sets,
                        const int* flags,
                        const long* end_indices,
                        const long* contents )
{
  long prev = -1;
  long i, end;
  
  for (i = 0; i < num_sets; ++i) {
    end = end_indices[i];
    if (end < prev) {
      fprintf(stderr,"Invalid end index for set %ld (ID %ld): %ld following %ld\n",
        i, first+i, end, prev );
      exit(2);
    }
    if (flags && (flags[i] & mhdf_SET_RANGE_BIT)) {
      dump_set_ranges( i+first, contents + prev + 1, end - prev );
    }
    else {
      dump_set_list( i+first, contents + prev + 1, end - prev );
    }
    prev = end;
  }
}

void dump_set_list( long id, const long* array, long count )
{
  long i;

  if (!count)
    return; /* skip empty sets */
  
  printf("%*ld: %ld", WIDTH, id, array[0] );
  for (i = 1; i < count; ++i)
    printf(", %ld", array[i] );
  printf(" (%ld)\n", count );
}

void dump_set_ranges( long id, const long* array, long count )
{
  long i, n = 0;
  
  if (!count)
    return; /* skip empty sets */
  printf("%*ld:", WIDTH, id );
  if (count % 2) {
    puts(" <INVALID DATA LENGTH FOR RANGED LIST>\n" );
    return;
  }
  
  for (i = 0; i < count; i += 2) {
    if (array[i+1] > 1) 
      printf(" %ld-%ld", array[i], array[i]+array[i+1]-1);
    else
      printf(" %ld", array[i]);
    if (i+2 < count)
      putchar( ',' );
    n += array[i+1];
  }
  printf(" (%ld)\n", n );
}


