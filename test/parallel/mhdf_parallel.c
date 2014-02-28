/* Test mhdf library on top of Parallel HDF5.
 *
 * Call mhdf API similar to how WriteHDF5Parallel does,
 * but avoid any of our own parallel communiation.
 *
 * Output will be composed of:
 *  a 1-D array of hexes, one per processor, where the first
 *    process writes all the nodes for its hex and every other
 *    process writes the right 4 nodes of its hex.
 *  a set containing all of the hexes
 *  a set containing all of the nodes
 *  a set per processor containing all of the entities written by that proc
 *  a tag specifying the process ID that wrote each entity.
 */

#include "mhdf.h"

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <H5Ppublic.h>
#include <H5Tpublic.h>
#include <H5Epublic.h>
#include <H5FDmpi.h>
#include <H5FDmpio.h>

#ifdef H5_HAVE_PARALLEL

const char filename[] = "mhdf_ll.h5m";
const char proc_tag_name[] = "proc_id";
const char elem_handle[] = "Hex";

int RANK;
int NUM_PROC;

#define CHECK( A ) check( &(A), __LINE__ )
void check( mhdf_Status* status, int line )
{
  if (mhdf_isError( status )) {
    fprintf(stderr,"%d: ERROR at line %d: \"%s\"\n", RANK, line, mhdf_message( status ) );
    abort();
  }
}

typedef int handle_id_t;
#define id_type H5T_NATIVE_INT

/* create file layout in serial */
void create_file()
{
  const char* elem_types[1] = {elem_handle};
  const int num_elem_types = sizeof(elem_types)/sizeof(elem_types[0]);
  const int num_nodes = 4+4*NUM_PROC;
  const int num_hexes = NUM_PROC;
  const int num_sets = 2+NUM_PROC;
  const int set_data_size = num_hexes + num_nodes + num_hexes + num_nodes;
  const int num_tag_vals = num_hexes + num_nodes + num_sets - 2;
  const char* history[] = { __FILE__, NULL };
  const int history_size = sizeof(history) / sizeof(history[0]);
  int default_tag_val = NUM_PROC;
  hid_t data, tagdata[2];
  long junk;
  
  mhdf_Status status;
  mhdf_FileHandle file;
  
  time_t t;
  t = time(NULL);
  history[1] = ctime(&t);
  
  file = mhdf_createFile( filename, 1, elem_types, num_elem_types, id_type, &status );
  CHECK(status);
  
  mhdf_writeHistory( file, history, history_size, &status );
  CHECK(status);
  
  data = mhdf_createNodeCoords( file, 3, num_nodes, &junk, &status );
  CHECK(status);
  mhdf_closeData( file, data, &status );
  CHECK(status);
  
  mhdf_addElement( file, elem_types[0], 0, &status );
  CHECK(status);
  data = mhdf_createConnectivity( file, elem_handle, 8, num_hexes, &junk, &status );
  CHECK(status);
  mhdf_closeData( file, data, &status );
  CHECK(status);
  
  data = mhdf_createSetMeta( file, num_sets, &junk, &status );
  CHECK(status);
  mhdf_closeData( file, data, &status );
  CHECK(status);

  data = mhdf_createSetData( file, set_data_size, &status );
  CHECK(status);
  mhdf_closeData( file, data, &status );
  CHECK(status);
  
  mhdf_createTag( file, proc_tag_name, mhdf_INTEGER, 1, mhdf_DENSE_TYPE, 
                  &default_tag_val, &default_tag_val, 0, 0, &status );
  CHECK(status);
  
  mhdf_createSparseTagData( file, proc_tag_name, num_tag_vals, tagdata, &status );
  CHECK(status);
  mhdf_closeData( file, tagdata[0], &status );
  CHECK(status);
  mhdf_closeData( file, tagdata[1], &status );
  CHECK(status);
 
  mhdf_closeFile( file, &status );
  CHECK(status);
}

void write_file_data()
{
  const int total_num_nodes = 4+4*NUM_PROC;
  const int total_num_hexes = NUM_PROC;
  long first_node, first_elem, first_set, count, ntag;
  unsigned long ucount;
  mhdf_index_t set_desc[4] = { 0, -1, -1, 0 };
  hid_t handle, handles[2];
  mhdf_Status status;
  mhdf_FileHandle file;
  int num_node, offset, dim;
  double coords[3*8] = { 0.0, 0.0, 0.0, 
                         0.0, 1.0, 0.0, 
                         0.0, 1.0, 1.0, 
                         0.0, 0.0, 1.0,
                         0.0, 0.0, 0.0, 
                         0.0, 1.0, 0.0, 
                         0.0, 1.0, 1.0, 
                         0.0, 0.0, 1.0 };
  handle_id_t list[10];
  int i, tagdata[10];
  for (i = 0; i < 10; i++) tagdata[i] = RANK;
  
  
    /* open file */
  handle = H5Pcreate( H5P_FILE_ACCESS );
  H5Pset_fapl_mpio( handle, MPI_COMM_WORLD, MPI_INFO_NULL );
  file = mhdf_openFileWithOpt( filename, 1, &ucount, id_type, handle, &status );
  CHECK(status);
  H5Pclose( handle );
  
    /* write node coordinates */
  if (0 == RANK) {
    num_node = 8;
    offset = 0;
    coords[12] = coords[15] = coords[18] = coords[21] = 1.0;
  }
  else {
    num_node = 4;
    offset = 4 + 4*RANK;
    coords[0] = coords[3] = coords[6] = coords[9] = 1.0 + RANK;
  }
  handle = mhdf_openNodeCoords( file, &count, &dim, &first_node, &status );
  CHECK(status);
  assert( count == total_num_nodes );
  assert( dim == 3 );
  mhdf_writeNodeCoords( handle, offset, num_node, coords, &status );
  CHECK(status);
  mhdf_closeData( file, handle, &status );
  CHECK(status);
  
    /* write hex connectivity */
  for (i = 0; i < 8; ++i)
    list[i] = 4*RANK + i + first_node;
  handle = mhdf_openConnectivity( file, elem_handle, &dim, &count, &first_elem, &status );
  CHECK(status);
  assert( count == total_num_hexes );
  assert( 8 == dim );
  mhdf_writeConnectivity( handle, RANK, 1, id_type, list, &status );
  CHECK(status);
  mhdf_closeData( file, handle, &status );
  CHECK(status);
  
    /* write set descriptions */
  handle = mhdf_openSetMeta( file, &count, &first_set, &status );
  CHECK(status);
  assert( count == 2+NUM_PROC );
  
    /* write set descriptions for per-proc sets */
  set_desc[0] = total_num_nodes + total_num_hexes + 9 + 5*RANK - 1;
  mhdf_writeSetMeta( handle, 2+RANK, 1, MHDF_INDEX_TYPE, set_desc, &status );
  CHECK(status);
  
    /* write set descriptions for multi-proc sets */
  if (0 == RANK) {
    set_desc[0] = total_num_nodes - 1;
    mhdf_writeSetMeta( handle, 0, 1, MHDF_INDEX_TYPE, set_desc, &status );
    CHECK(status);
    set_desc[0] += total_num_hexes;
    mhdf_writeSetMeta( handle, 1, 1, MHDF_INDEX_TYPE, set_desc, &status );
    CHECK(status);
  }

  mhdf_closeData( file, handle, &status );
  CHECK(status);
  
    /* write set contents */
  
  handle = mhdf_openSetData( file, &count, &status );
  CHECK(status);
  assert( 2*total_num_nodes + 2*total_num_hexes == count );
  
    /* write per-proc sets */
  if (0 == RANK) {
    count = 9;
    offset = total_num_nodes + total_num_hexes;
    for (i = 0; i < 8; ++i)
      list[i] = first_node + i;
    list[8] = first_elem;
  }
  else {
    count = 5;
    offset = total_num_nodes + total_num_hexes + 4 + 5*RANK;
    for (i = 0; i < 4; ++i)
      list[i] = 4 + 4*RANK + first_node + i;
    list[4] = RANK + first_elem;
  }
  mhdf_writeSetData( handle, offset, count, id_type, list, &status );
  CHECK(status);
  
    /* write multi-proc sets */
    /* write nodes */
  offset = (0 == RANK) ? 0 : 4 + 4*RANK;
  mhdf_writeSetData( handle, offset, count-1, id_type, list, &status );
  CHECK(status);
    /* write hex */
  offset = total_num_nodes + RANK;
  mhdf_writeSetData( handle, offset, 1, id_type, list + count - 1, &status );
  CHECK(status);
  
  mhdf_closeData( file, handle, &status );
  CHECK(status);
  
    /* write tag data */
  offset = (0 == RANK) ? 0 : 4 + 4*RANK;
  offset += 2*RANK;
  list[count++] = first_set + 2 + RANK;
  mhdf_openSparseTagData( file, proc_tag_name, &ntag, &ntag, handles, &status );
  CHECK(status);
  assert( ntag == total_num_nodes + total_num_hexes + NUM_PROC );
  mhdf_writeSparseTagEntities( handles[0], offset, count, id_type, list, &status );
  CHECK(status);
  mhdf_writeTagValues( handles[1], offset, count, H5T_NATIVE_INT, tagdata, &status );
  CHECK(status);
  mhdf_closeData( file, handles[0], &status );
  CHECK(status);
  mhdf_closeData( file, handles[1], &status );
  CHECK(status);
  
    /* done */
  mhdf_closeFile( file, &status );
  CHECK(status);
}

/* Define a dummy error handler to register with HDF5. 
 * This function doesn't do anything except pass the
 * error on to the default handler that would have been
 * called anyway.  It's only purpose is to provide a 
 * spot to set a break point so we can figure out where
 * (in our code) that we made an invalid call into HDF5
 */
#if defined(H5E_auto_t_vers) && H5E_auto_t_vers > 1
herr_t (*default_handler)( hid_t, void* );
static herr_t handle_hdf5_error( hid_t stack, void* data )
#else
herr_t (*default_handler)( void* );
static herr_t handle_hdf5_error( void* data )
#endif
{
  herr_t result = 0;
  if (default_handler)
#if defined(H5E_auto_t_vers) && H5E_auto_t_vers > 1
    result = (default_handler)(stack,data);
#else
    result = (default_handler)(data);
#endif
  assert(0);
  return result;
}

#endif /* #ifdef H5_HAVE_PARALLEL */
 
int main( int argc, char* argv[] )
{
#ifdef H5_HAVE_PARALLEL
  int rval;
  void* data;
  herr_t err;

#if defined(H5Eget_auto_vers) && H5Eget_auto_vers > 1
  err = H5Eget_auto( H5E_DEFAULT, &default_handler, &data );
#else
  err = H5Eget_auto( &default_handler, &data );
#endif
  if (err >= 0) {
#if defined(H5Eset_auto_vers) && H5Eset_auto_vers > 1
    H5Eset_auto( H5E_DEFAULT, &handle_hdf5_error, data );
#else
    H5Eset_auto( &handle_hdf5_error, data );
#endif
  }
  
  rval = MPI_Init( &argc, &argv );
  if (rval)
    return rval;
  rval = MPI_Comm_rank( MPI_COMM_WORLD, &RANK );
  if (rval)
    return rval;
  rval = MPI_Comm_size( MPI_COMM_WORLD, &NUM_PROC );
  if (rval)
    return rval;
  
  if (RANK == 0)
    create_file();

  /* Wait for rank 0 to finish creating the file, otherwise rank 1 may find it to be invalid */
  rval = MPI_Barrier(MPI_COMM_WORLD);
  if (rval)
    return rval;

  write_file_data();
  
  MPI_Finalize();
#endif
  return 0;
}
