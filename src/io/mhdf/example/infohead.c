#include "mhdf.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <H5Tpublic.h>

int main( int argc, char* argv[] )
{
  mhdf_FileHandle file;
  mhdf_Status status;
  unsigned long max_id;
  struct mhdf_FileDesc* data;
  /* find PARALLEL_PARTITION tag index */
  const char * pname = "PARALLEL_PARTITION";
  int i, maxc, ne;
  long int nval, junk;
  hid_t table[3];
  
  if (argc != 2) {
    fprintf( stderr,"Usage: %s <filename>\n", argv[0] );
    return 1;
  }
  
  file = mhdf_openFile( argv[1], 0, &max_id, -1, &status );
  if (mhdf_isError( &status )) {
    fprintf( stderr,"%s: %s\n", argv[1], mhdf_message( &status ) );
    return 1;
  }
  
  data = mhdf_getFileSummary( file, H5T_NATIVE_ULONG, &status );
  if (mhdf_isError( &status )) {
    fprintf( stderr,"%s: %s\n", argv[1], mhdf_message( &status ) );
    return 1;
  }
  
  printf(" num dimensions: %d\n", data->nodes.vals_per_ent);
  printf(" num vertices:   %d\n", (int)data->nodes.count);

  /* max connectivity number of elem */
  maxc = -1;
  ne = 0;
  for (i=0; i<data->num_elem_desc; i++)
  {
    struct mhdf_ElemDesc * el_desc = &(data->elems[i]);
    struct mhdf_EntDesc * ent_d = &(el_desc->desc);
    if (ent_d->vals_per_ent > maxc)
    {
      ne = ent_d->count;
      maxc = ent_d->vals_per_ent;
    }
  }
  printf(" num elements:   %d\n", ne);

  for (i=0; i<data->num_tag_desc; i++)
  {
    struct mhdf_TagDesc * tag_desc = &(data->tags[i]);
    if (strcmp(pname,tag_desc->name)==0)
    {
      /*printf(" tag index %d is parallel partition tag\n", i);*/
      if (tag_desc->have_sparse) {
        mhdf_openSparseTagData(file, pname, &nval, &junk, table, &status);
        if (mhdf_isError( &status )) {
          fprintf( stderr,"%s: %s\n", argv[1], mhdf_message( &status ) );
          return 1;
        }
      }
      else
      {
        /* could be dense tags on sets */
        mhdf_openDenseTagData(file, pname, mhdf_set_type_handle(), &nval, &status);
        if (mhdf_isError( &status )) {
          fprintf( stderr,"%s: %s\n", argv[1], mhdf_message( &status ) );
          return 1;
        }
      }
      printf(" num parallel partitions: %d\n", (int)nval);
    }
  }

  mhdf_closeFile( file, &status );
  
  free( data );
  return 0;
}
