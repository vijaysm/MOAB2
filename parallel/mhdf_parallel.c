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
#include "MBTypes.h"

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <H5Ppublic.h>
#include <H5Tpublic.h>

const char* filename = "mhdf_ll.h5m";
const char* proc_tag_name = "proc_id";
const char* elem_handle = "Hex";

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

/* create file layout in serial */
void create_file()
{
  const char* elem_types[] = { elem_handle };
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
  
  time_t t = time(NULL);
  history[1] = ctime(&t);
  
  mhdf_Status status;
  mhdf_FileHandle file;
  
  file = mhdf_createFile( filename, 1, elem_types, num_elem_types, &status );
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
  long set_desc[4] = { 0, -1, -1, MESHSET_SET };
  hid_t handle, handles[2];
  mhdf_Status status;
  mhdf_FileHandle file;
  double coords[3*8] = { 0.0, 0.0, 0.0, 
                         0.0, 1.0, 0.0, 
                         0.0, 1.0, 1.0, 
                         0.0, 0.0, 1.0,
                         0.0, 0.0, 0.0, 
                         0.0, 1.0, 0.0, 
                         0.0, 1.0, 1.0, 
                         0.0, 0.0, 1.0 };
  int i, list[10], tagdata[10] = {RANK,RANK,RANK,RANK,RANK,RANK,RANK,RANK,RANK,RANK};
  int num_node, offset, dim;
  
  
    /* open file */
  handle = H5Pcreate( H5P_FILE_ACCESS );
  H5Pset_fapl_mpio( handle, MPI_COMM_WORLD, MPI_INFO_NULL );
  file = mhdf_openFileWithOpt( filename, 1, &count, handle, &status );
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
    list[i] = 4*RANK + i + 1;
  handle = mhdf_openConnectivity( file, elem_handle, &dim, &count, &first_elem, &status );
  CHECK(status);
  assert( count == total_num_hexes );
  assert( 8 == dim );
  mhdf_writeConnectivity( handle, RANK, 1, H5T_NATIVE_INT, list, &status );
  CHECK(status);
  mhdf_closeData( file, handle, &status );
  CHECK(status);
  
    /* write set descriptions */
  handle = mhdf_openSetMeta( file, &count, &first_set, &status );
  CHECK(status);
  assert( count == 2+NUM_PROC );
  
    /* write set descriptions for per-proc sets */
  set_desc[0] = total_num_nodes + total_num_hexes + 9 + 5*RANK - 1;
  mhdf_writeSetMeta( handle, 2+RANK, 1, H5T_NATIVE_LONG, set_desc, &status );
  CHECK(status);
  
    /* write set descriptions for multi-proc sets */
  if (0 == RANK) {
    set_desc[0] = total_num_nodes - 1;
    mhdf_writeSetMeta( handle, 0, 1, H5T_NATIVE_LONG, set_desc, &status );
    CHECK(status);
    set_desc[0] += total_num_hexes;
    mhdf_writeSetMeta( handle, 1, 1, H5T_NATIVE_LONG, set_desc, &status );
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
  mhdf_writeSetData( handle, offset, count, H5T_NATIVE_INT, list, &status );
  CHECK(status);
  
    /* write multi-proc sets */
    /* write nodes */
  offset = (0 == RANK) ? 0 : 4 + 4*RANK;
  mhdf_writeSetData( handle, offset, count-1, H5T_NATIVE_INT, list, &status );
  CHECK(status);
    /* write hex */
  offset = total_num_nodes + RANK;
  mhdf_writeSetData( handle, offset, 1, H5T_NATIVE_INT, list + count - 1, &status );
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
  mhdf_writeSparseTagEntities( handles[0], offset, count, H5T_NATIVE_INT, list, &status );
  CHECK(status);
  mhdf_writeSparseTagValues( handles[1], offset, count, H5T_NATIVE_INT, tagdata, &status );
  CHECK(status);
  mhdf_closeData( file, handles[0], &status );
  CHECK(status);
  mhdf_closeData( file, handles[1], &status );
  CHECK(status);
  
    /* done */
  mhdf_closeFile( file, &status );
  CHECK(status);
}
 
 
int main( int argc, char* argv[] )
{
  mhdf_Status status;
  int rval;
  
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
  
  write_file_data();
  
  MPI_Finalize();
  return 0;
}
