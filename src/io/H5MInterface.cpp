/*
 * H5MInterface.cpp
 */

#include "moab/H5MInterface.hpp"

#include "mhdf.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <H5Tpublic.h>

namespace moab {

H5MInterface::H5MInterface(const char * filena) : filename(filena) {

}

H5MInterface::~H5MInterface() {

}

int H5MInterface::get_mesh_info( int* verts, int *edges, int*faces, int* regions, int *numdim, int* parts)
{
  *verts = 0;
  *edges = 0;
  *faces = 0;
  *regions = 0;
  *numdim = 0;
  *parts = 0;

  mhdf_FileHandle file;
  mhdf_Status status;
  unsigned long max_id;
  struct mhdf_FileDesc* data;
  /* find PARALLEL_PARTITION tag index */
  const char * pname = "PARALLEL_PARTITION";

  long int nval, junk;
  hid_t table[3];


  file = mhdf_openFile( filename, 0, &max_id, -1, &status );
  if (mhdf_isError( &status )) {
    fprintf( stderr,"%s: %s\n", filename, mhdf_message( &status ) );
    return 1;
  }

  data = mhdf_getFileSummary( file, H5T_NATIVE_ULONG, &status );
  if (mhdf_isError( &status )) {
    fprintf( stderr,"%s: %s\n", filename, mhdf_message( &status ) );
    return 1;
  }
  *numdim = data->nodes.vals_per_ent;
  *verts = (int)data->nodes.count;

  for (int i=0; i<data->num_elem_desc; i++)
  {
    struct mhdf_ElemDesc * el_desc = &(data->elems[i]);
    struct mhdf_EntDesc * ent_d = &(el_desc->desc);
    if (0==strcmp(el_desc->type, mhdf_EDGE_TYPE_NAME)) *edges += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_TRI_TYPE_NAME))  *faces += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_QUAD_TYPE_NAME)) *faces += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_POLYGON_TYPE_NAME)) *faces += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_TET_TYPE_NAME)) *regions += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_PYRAMID_TYPE_NAME)) *regions += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_PRISM_TYPE_NAME)) *regions += ent_d->count;
    if (0==strcmp(el_desc->type, mdhf_KNIFE_TYPE_NAME)) *regions += ent_d->count;
    if (0==strcmp(el_desc->type, mdhf_HEX_TYPE_NAME)) *regions += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_POLYHEDRON_TYPE_NAME)) *regions += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_SEPTAHEDRON_TYPE_NAME)) *regions += ent_d->count;
  }
  for (int i=0; i<data->num_tag_desc; i++)
  {
    struct mhdf_TagDesc * tag_desc = &(data->tags[i]);
    if (strcmp(pname,tag_desc->name)==0)
    {
      /*printf(" tag index %d is parallel partition tag\n", i);*/
      if (tag_desc->have_sparse) {
        mhdf_openSparseTagData(file, pname, &nval, &junk, table, &status);
        if (mhdf_isError( &status )) {
          fprintf( stderr,"%s: %s\n", filename, mhdf_message( &status ) );
          return 1;
        }
      }
      else
      {
        /* could be dense tags on sets */
        mhdf_openDenseTagData(file, pname, mhdf_set_type_handle(), &nval, &status);
        if (mhdf_isError( &status )) {
          fprintf( stderr,"%s: %s\n", filename, mhdf_message( &status ) );
          return 1;
        }
      }

      *parts = (int)nval;
    }
  }

  if (*faces >0 ) *numdim = 2;
  if (*regions>0) *numdim = 3;
  mhdf_closeFile( file, &status );

  free( data );
  return 0;
}
} /* namespace moab */
