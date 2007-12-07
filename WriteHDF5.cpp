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

//-------------------------------------------------------------------------
// Filename      : WriteHDF5.cpp
//
// Purpose       : TSTT HDF5 Writer 
//
// Special Notes : WriteSLAC used as template for this
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/01/04
//-------------------------------------------------------------------------

#include "MBEntityHandle.h"
#ifndef HDF5_FILE
#  error Attempt to compile WriteHDF5 with HDF5 support disabled
#endif

#include <assert.h>
#ifndef _MSC_VER
#include <sys/time.h>
#endif
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <limits>
#include <H5Tpublic.h>
#include "MBInterface.hpp"
#include "MBInternals.hpp"
#include "MBTagConventions.hpp"
#include "MBCN.hpp"
#include "WriteHDF5.hpp"
#include "MBWriteUtilIface.hpp"
#include "mhdf.h"
/* Access HDF5 file handle for debugging
#include <H5Fpublic.h>
struct file { uint32_t magic; hid_t handle; };
*/
#undef DEBUG

#ifdef DEBUG
#  define DEBUGOUT(A) fputs( A, stderr )
#  include <stdio.h>
#else
#  define DEBUGOUT(A)
#endif

#ifdef DEBUG
/*
# include <H5Epublic.h>
  extern "C" herr_t hdf_error_handler( void*  )
  {
    H5Eprint( stderr );
    assert( 0 );
  }
*/
# define myassert(A) assert(A)
#else
# define myassert(A)
#endif

#define WRITE_HDF5_BUFFER_SIZE (40*1024*1024)

  // This is the tag used to store file offsets (file ids)
  // during export.  It is temporary data.
#define HDF5_ID_TAG_NAME "hdf5_id"

  // This is the HDF5 type for the data of the "hdf5_id" tag.
const hid_t WriteHDF5::id_type = H5T_NATIVE_INT;

  // Some macros to handle error checking.  The
  // CHK_MHDF__ERR* macros check the value of an mhdf_Status 
  // object.  The CHK_MB_ERR_* check the value of an MBStatus.
  // The *_0 macros accept no other arguments. The *_1
  // macros accept a single hdf5 handle to close on error.
  // The *_2 macros accept an array of two hdf5 handles to
  // close on error.  The _*2C macros accept one hdf5 handle
  // to close on error and a bool and an hdf5 handle where
  // the latter handle is conditionally closed depending on
  // the value of the bool.  All macros contain a "return"
  // statement.
#define CHK_MHDF_ERR_0( A )                                 \
do if ( mhdf_isError( &(A) )) {                             \
    writeUtil->report_error( "%s\n", mhdf_message( &(A) ) );\
    myassert(0);                                            \
    return MB_FAILURE;                                      \
} while(false)                                               

#define CHK_MHDF_ERR_1( A, B )                              \
do if ( mhdf_isError( &(A) )) {                             \
    writeUtil->report_error( "%s\n", mhdf_message( &(A) ) );\
    myassert(0);                                            \
    mhdf_closeData( filePtr, (B), &(A) );                   \
    return MB_FAILURE;                                      \
} while(false)                                               

#define CHK_MHDF_ERR_2( A, B )                              \
do if ( mhdf_isError( &(A) )) {                             \
    writeUtil->report_error( "%s\n", mhdf_message( &(A) ) );\
    myassert(0);                                            \
    mhdf_closeData( filePtr, (B)[0], &(A) );                \
    mhdf_closeData( filePtr, (B)[1], &(A) );                \
    return MB_FAILURE;                                      \
} while(false)                                               

#define CHK_MHDF_ERR_2C( A, B, C, D )                       \
do if ( mhdf_isError( &(A) )) {                             \
    writeUtil->report_error( "%s\n", mhdf_message( &(A) ) );\
    myassert(0);                                            \
    mhdf_closeData( filePtr, (B), &(A) );                   \
    if (C) mhdf_closeData( filePtr, (D), &(A) );            \
    return MB_FAILURE;                                      \
} while(false)                                               


#define CHK_MB_ERR_0( A ) \
do if (MB_SUCCESS != (A)) return (A); while(false)

#define CHK_MB_ERR_1( A, B, C )         \
do if (MB_SUCCESS != (A)) {             \
  mhdf_closeData( filePtr, (B), &(C) ); \
  myassert(0);                          \
  return (A);                           \
} while(false)

#define CHK_MB_ERR_2( A, B, C )            \
do if (MB_SUCCESS != (A)) {                \
  mhdf_closeData( filePtr, (B)[0], &(C) ); \
  mhdf_closeData( filePtr, (B)[1], &(C) ); \
  write_finished();                        \
  myassert(0);                             \
  return (A);                              \
} while(false)

#define CHK_MB_ERR_2C( A, B, C, D, E )          \
do if (MB_SUCCESS != (A)) {                     \
  mhdf_closeData( filePtr, (B), &(E) );         \
  if (C) mhdf_closeData( filePtr, (D), &(E) );  \
  write_finished();                             \
  myassert(0);                                  \
  return (A);                                   \
} while(false)

/** When writing tags containing MBEntityHandles to file, need to convert tag
 *  data from MBEntityHandles to file IDs.  This function does that. 
 *
 * If the handle is not valid or does not correspond to an entity that will
 * be written to the file, the file ID is set to zero.
 *\param iface MOAB instance
 *\param tag   The TAG containing file IDs for entities
 *\param data  The data buffer.  As input, an array of MBEntityHandles.  As
 *             output an array of file IDS, where the size of each integral
 *             file ID is the same as the size of MBEntityHandle.
 *\param count The number of handles in the buffer.
 *\return true if at least one of the handles is valid and will be written to
 *             the file or at least one of the handles is NULL (zero). false
 *             otherwise
 */
static bool convert_handle_tag( MBInterface* iface, MBTag tag, void* data, size_t count )
{
  bool some_valid = false;
  MBErrorCode rval;
  
  // if same saize
  MBEntityHandle *buffer, *end;
  WriteHDF5::id_t* witer;
  int step;
  if (sizeof(MBEntityHandle) >= sizeof(WriteHDF5::id_t)) {
    buffer = (MBEntityHandle*)data;
    end = buffer + count;
    witer = (WriteHDF5::id_t*)data;
    step = 1;
  }
  else {
    // iterate in reverse order if sizeof(id_t) > sizeof(MBEntityHandle)
    buffer = (MBEntityHandle*)data + count - 1;
    end = (MBEntityHandle*)data - 1;
    witer = (WriteHDF5::id_t*)data + count - 1;
    step = -1;
  }
  for ( ; buffer != end; buffer += step, witer += step) {
    if (!*buffer) {
      some_valid = true;
      *witer = 0;
    }
    else {
      int id;
      rval = iface->tag_get_data( tag, buffer, 1, &id );
      if (MB_SUCCESS == rval && id > 0) {
        some_valid = true;
        *witer = id;
      }
      else {
        *buffer = 0;
        *witer = 0;
      }
    }
  }

  return some_valid;
}

const char* WriteHDF5::ExportSet::name() const
{
  static char buffer[32];
  switch (type) {
    case MBVERTEX:
      return mhdf_node_type_handle();
    case MBENTITYSET:
      return mhdf_set_type_handle();
    default:
      sprintf( buffer, "%s%d", MBCN::EntityTypeName( type ), num_nodes );
      return buffer;
  }
}
  

MBWriterIface* WriteHDF5::factory( MBInterface* iface )
  { return new WriteHDF5( iface ); }

WriteHDF5::WriteHDF5( MBInterface* iface )
  : bufferSize( WRITE_HDF5_BUFFER_SIZE ),
    dataBuffer( 0 ),
    iFace( iface ), 
    writeUtil( 0 ), 
    filePtr( 0 ), 
    createdIdTag(false),
    idTag( 0 ),
    setContentsOffset( 0 ),
    setChildrenOffset( 0 ),
    setParentsOffset( 0 ),
    writeSets(false),
    writeSetContents(false),
    writeSetChildren(false),
    writeSetParents(false)
{
}

MBErrorCode WriteHDF5::init()
{
  MBErrorCode rval;
  id_t zero_int = -1;

  if (writeUtil) // init has already been called
    return MB_SUCCESS;
/* 
#ifdef DEBUG
  H5Eset_auto( &hdf_error_handler, writeUtil );  // HDF5 callback for errors
#endif
*/ 
    // For known tag types, store the corresponding HDF5 in which
    // the tag data is to be written in the file.
  //register_known_tag_types( iFace ); 
 
    // Get the util interface
  rval = iFace->query_interface( "MBWriteUtilIface", (void**)&writeUtil );
  CHK_MB_ERR_0(rval);
  
    // Get the file id tag if it exists or create it otherwise.
  rval = iFace->tag_get_handle( HDF5_ID_TAG_NAME, idTag );
  if (MB_TAG_NOT_FOUND == rval)
  {
    rval = iFace->tag_create( HDF5_ID_TAG_NAME, sizeof(id_t), 
                              MB_TAG_DENSE, idTag, &zero_int );
    if (MB_SUCCESS == rval)
      createdIdTag = true;
  }
  
  if (MB_SUCCESS != rval)
  {
    iFace->release_interface( "MBWriteUtilIFace", writeUtil );
    writeUtil = 0;
    return rval;
  }

  return MB_SUCCESS;
}
  
MBErrorCode WriteHDF5::write_finished()
{
    // release memory allocated in lists
  exportList.clear();
  nodeSet.range.clear();
  setSet.range.clear();
  rangeSets.clear();
  tagList.clear();
  return MB_SUCCESS;
}
  

WriteHDF5::~WriteHDF5()
{
  if (!writeUtil) // init() failed.
    return;

  iFace->release_interface( "MBWriteUtilIface", writeUtil );
  if (createdIdTag)
    iFace->tag_delete( idTag );
}


MBErrorCode WriteHDF5::write_file( const char* filename,
                                   bool overwrite,
                                   const FileOptions& ,
                                   const MBEntityHandle* set_array,
                                   const int num_sets,
                                   std::vector<std::string>& qa_records,
                                   int user_dimension )
{
  MBErrorCode result;
  mhdf_Status rval;
  std::list<SparseTag>::const_iterator t_itor;
  std::list<ExportSet>::iterator ex_itor;
  MBEntityHandle elem_count, max_id;
  
  if (MB_SUCCESS != init())
    return MB_FAILURE;

DEBUGOUT("Gathering Mesh\n");
  
    // Gather mesh to export
  exportList.clear();
  if (0 == num_sets || 1 == num_sets && set_array[0] == 0)
  {
    result = gather_all_mesh( );
    CHK_MB_ERR_0(result);
  }
  else
  {
    std::vector<MBEntityHandle> passed_export_list(num_sets);
    memcpy( &passed_export_list[0], set_array, sizeof(MBEntityHandle)*num_sets );
    result = gather_mesh_info( passed_export_list );
    if (MB_SUCCESS != result) goto write_fail;
    
      // Mark all entities invalid.  Later the ones we are
      // exporting will be marked valid.  This way we can
      // distinguish adjacent entities and such that aren't
      // being exported.
      // Don't to this, just set the default value to -1 when
      // the tag is created in the init() function.
    //result = clear_all_id_tags();
    //if (MB_SUCCESS != result) goto write_fail;
  }
  
  //if (nodeSet.range.size() == 0)
  //  goto write_fail;
  
DEBUGOUT("Checking ID space\n");

    // Make sure ID space is sufficient
  elem_count = nodeSet.range.size() + setSet.range.size();
  for (ex_itor = exportList.begin(); ex_itor != exportList.end(); ++ex_itor)
    elem_count += ex_itor->range.size();
  max_id = (MBEntityHandle)1 << (8*sizeof(id_t)-1);
  if (elem_count > max_id)
  {
    writeUtil->report_error("ID space insufficient for mesh size.\n");
    goto write_fail;
  }

DEBUGOUT( "Creating File\n" );  

    // Figure out the dimension in which to write the mesh.  
  int mesh_dim;
  result = iFace->get_dimension( mesh_dim );
  CHK_MB_ERR_0(result);
  
  if (user_dimension < 1) 
    user_dimension = mesh_dim;
  user_dimension = user_dimension > mesh_dim ? mesh_dim : user_dimension;
  
    // Allocate internal buffer to use when gathering data to write.
  dataBuffer = (char*)malloc( bufferSize );
  if (!dataBuffer)
    goto write_fail;
  
    // Create the file layout, including all tables (zero-ed) and
    // all structure and meta information.
  result = create_file( filename, overwrite, qa_records, user_dimension );
  if (MB_SUCCESS != result)
    goto write_fail;

DEBUGOUT("Writing Nodes.\n");
  
    // Write node coordinates
  if (!nodeSet.range.empty() && write_nodes() != MB_SUCCESS)
    goto write_fail;

DEBUGOUT("Writing connectivity.\n");
  
    // Write element connectivity
  for (ex_itor = exportList.begin(); ex_itor != exportList.end(); ++ex_itor)
    if (MB_SUCCESS != write_elems( *ex_itor ))
      goto write_fail;

DEBUGOUT("Writing sets.\n");
  
    // Write meshsets
  if (write_sets() != MB_SUCCESS)
    goto write_fail;

DEBUGOUT("Writing adjacencies.\n");
  
    // Write adjacencies
  // Tim says don't save node adjacencies!
#ifdef WRITE_NODE_ADJACENCIES
  if (write_adjacencies( nodeSet ) != MB_SUCCESS)
    goto write_fail;
#endif
  for (ex_itor = exportList.begin(); ex_itor != exportList.end(); ++ex_itor)
    if (write_adjacencies( *ex_itor ) != MB_SUCCESS)
      goto write_fail;

DEBUGOUT("Writing tags.\n");
  

    // Write tags
  for (t_itor = tagList.begin(); t_itor != tagList.end(); ++t_itor)
    if (t_itor->write)
      if (write_sparse_tag( *t_itor ) != MB_SUCCESS)
        goto write_fail;

DEBUGOUT("Closing file.\n");

    // Clean up and exit.
  free( dataBuffer );
  dataBuffer = 0;
  mhdf_closeFile( filePtr, &rval );
  filePtr = 0;
  result = write_finished();
  CHK_MHDF_ERR_0( rval );
  return result;
  
write_fail:
  
  if (dataBuffer)
  {
    free( dataBuffer );
    dataBuffer = 0;
  }
  mhdf_closeFile( filePtr, &rval );
  filePtr = 0;
  write_finished();
  return MB_FAILURE;
}

  // Initialize all file ids to -1.  We do this so that
  // if we are writing only a part of the mesh, we can
  // easily identify adjacencies, set contents, etc. which
  // are not in the mesh we are writing.
MBErrorCode WriteHDF5::clear_all_id_tags()
{
  id_t cleared_value = -1;
  MBRange range;
  for (MBEntityType type = MBVERTEX; type < MBMAXTYPE; ++type)
  {
    range.clear();
    MBErrorCode rval = iFace->get_entities_by_type( 0, type, range, false );
    CHK_MB_ERR_0(rval);
    
    for (MBRange::iterator itor = range.begin(); itor != range.end(); ++itor)
    {
      MBEntityHandle handle = *itor;
      rval = iFace->tag_set_data( idTag, &handle, 1, &cleared_value );
      CHK_MB_ERR_0(rval);
    }
  }
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5::initialize_mesh( const MBRange ranges[5] )
{
  MBErrorCode rval;
  
  if (!ranges[0].all_of_type(MBVERTEX))
    return MB_FAILURE;
  nodeSet.range = ranges[0];
  nodeSet.type = MBVERTEX;
  nodeSet.num_nodes = 1;
  
  if (!ranges[4].all_of_type(MBENTITYSET))
    return MB_FAILURE;
  setSet.range = ranges[4];
  setSet.type = MBENTITYSET;
  setSet.num_nodes = 0;

  exportList.clear();
  std::vector<MBRange> bins(1024); // sort entities by connectivity length
                                   // resize is expensive due to MBRange copy, so start big
  for (MBEntityType type = MBEDGE; type < MBENTITYSET; ++type)
  {
    ExportSet set;
    const int dim = MBCN::Dimension(type);

      // Group entities by connectivity length
    bins.clear();
    std::pair<MBRange::const_iterator,MBRange::const_iterator> p = ranges[dim].equal_range(type);
    MBRange::const_iterator i = p.first;
    while (i != p.second) {
      MBRange::const_iterator first = i;
      MBEntityHandle const* conn;
      int len, firstlen;
      rval = iFace->get_connectivity( *i, conn, firstlen );
      if (MB_SUCCESS != rval)
        return rval;
      
      for (++i; i != p.second; ++i) {
        rval = iFace->get_connectivity( *i, conn, len );
        if (MB_SUCCESS != rval)
          return rval;
        
        if (len != firstlen)
          break;
      }
      
      if (firstlen >= (int)bins.size())
        bins.resize(firstlen+1);
      bins[firstlen].merge( first, i );
    }

      // Create ExportSet for each group
    for (std::vector<MBRange>::iterator j = bins.begin(); j != bins.end(); ++j) {
      if (j->empty())
        continue;
        
      set.range.clear();
      set.type = type;
      set.num_nodes = j - bins.begin();
      exportList.push_back( set );
      exportList.back().range.swap( *j );
    }
  }
    
  return MB_SUCCESS;  
}

                                         
  // Gather the mesh to be written from a list of owning meshsets.
MBErrorCode WriteHDF5::gather_mesh_info( 
                           const std::vector<MBEntityHandle>& export_sets )
{
  MBErrorCode rval;
  
  int dim;
  MBRange range;      // temporary storage
  MBRange ranges[5];  // lists of entities to export, grouped by dimension
  
    // Gather list of all related sets
  std::vector<MBEntityHandle> stack(export_sets);
  std::copy( export_sets.begin(), export_sets.end(), stack.begin() );
  std::vector<MBEntityHandle> set_children;
  while( !stack.empty() )
  {
    MBEntityHandle meshset = stack.back(); stack.pop_back();
    ranges[4].insert( meshset );
  
      // Get contained sets
    range.clear();
    rval = iFace->get_entities_by_type( meshset, MBENTITYSET, range );
    CHK_MB_ERR_0(rval);
    for (MBRange::iterator ritor = range.begin(); ritor != range.end(); ++ritor)
      if (ranges[4].find( *ritor ) == ranges[4].end())
        stack.push_back( *ritor );
    
      // Get child sets
    set_children.clear();
    rval = iFace->get_child_meshsets( meshset, set_children, 1 );
    CHK_MB_ERR_0(rval);
    for (std::vector<MBEntityHandle>::iterator vitor = set_children.begin();
         vitor != set_children.end(); ++vitor )
      if (ranges[4].find( *vitor ) == ranges[4].end())
        stack.push_back( *vitor );
  }
  
    // Gather list of all mesh entities from list of sets,
    // grouped by dimension.
  for (MBRange::iterator setitor = ranges[4].begin();
       setitor != ranges[4].end(); ++setitor)
  {
    for (dim = 0; dim < 4; ++dim)
    {
      range.clear();
      rval = iFace->get_entities_by_dimension( *setitor, dim, range, false );
      CHK_MB_ERR_0(rval);

      ranges[dim].merge(range);
    }
  }
  
    // For each list of elements, append adjacent children and
    // nodes to lists.
  for (dim = 3; dim > 0; --dim)
  {
    for (int cdim = 1; cdim < dim; ++cdim)
    {
      range.clear();
      rval = iFace->get_adjacencies( ranges[dim], cdim, false, range );
      CHK_MB_ERR_0(rval);
      ranges[cdim].merge( range );
    }  
    range.clear();
    rval = writeUtil->gather_nodes_from_elements( ranges[dim], 0, range );
    CHK_MB_ERR_0(rval);
    ranges[0].merge( range );      
  }
  
  return initialize_mesh( ranges );
}

  // Gather all the mesh and related information to be written.
MBErrorCode WriteHDF5::gather_all_mesh( )
{
  MBErrorCode rval;
  MBRange ranges[5];

  rval = iFace->get_entities_by_type( 0, MBVERTEX, ranges[0] );
  if (MB_SUCCESS != rval)
    return rval;

  rval = iFace->get_entities_by_dimension( 0, 1, ranges[1] );
  if (MB_SUCCESS != rval)
    return rval;

  rval = iFace->get_entities_by_dimension( 0, 2, ranges[2] );
  if (MB_SUCCESS != rval)
    return rval;

  rval = iFace->get_entities_by_dimension( 0, 3, ranges[3] );
  if (MB_SUCCESS != rval)
    return rval;

  rval = iFace->get_entities_by_type( 0, MBENTITYSET, ranges[4] );
  if (MB_SUCCESS != rval)
    return rval;

  return initialize_mesh( ranges );
}
  
MBErrorCode WriteHDF5::write_nodes( )
{
  mhdf_Status status;
  int dim, mesh_dim;
  MBErrorCode rval;
  hid_t node_table;
  long first_id, num_nodes;
  
  rval = iFace->get_dimension( mesh_dim );
  CHK_MB_ERR_0(rval);
  
  node_table = mhdf_openNodeCoords( filePtr, &num_nodes, &dim, &first_id, &status );
  CHK_MHDF_ERR_0(status);
  
  double* buffer = (double*)dataBuffer;
  int chunk_size = bufferSize / sizeof(double);
  
  long remaining = nodeSet.range.size();
  long offset = nodeSet.offset;
  MBRange::const_iterator iter = nodeSet.range.begin();
  while (remaining)
  {
    long count = chunk_size < remaining ? chunk_size : remaining;
    remaining -= count;
    MBRange::const_iterator end = iter;
    end += count;
    
    for (int d = 0; d < dim; d++)
    {
      if (d < mesh_dim)
      {
        rval = writeUtil->get_node_array( d, iter, end, count, buffer );
        CHK_MB_ERR_1(rval, node_table, status);
      }
      else
      {
        memset( buffer, 0, count * sizeof(double) );
      }
    
      mhdf_writeNodeCoord( node_table, offset, count, d, buffer, &status );
      CHK_MHDF_ERR_1(status, node_table);
    }
    
    iter = end;
    offset += count;
  }
  
  mhdf_closeData( filePtr, node_table, &status );
  CHK_MHDF_ERR_0(status);
 
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5::write_elems( ExportSet& elems )
{
  mhdf_Status status;
  MBErrorCode rval;
  long first_id;
  int nodes_per_elem;
  long table_size;

  hid_t elem_table = mhdf_openConnectivity( filePtr, 
                                            elems.name(), 
                                            &nodes_per_elem,
                                            &table_size,
                                            &first_id,
                                            &status );
                                            
  CHK_MHDF_ERR_0(status);
  assert (first_id <= elems.first_id);
  assert ((unsigned long)table_size >= elems.offset + elems.range.size());
  
  
  id_t* buffer = (id_t*)dataBuffer;
  int chunk_size = bufferSize / (elems.num_nodes * sizeof(id_t));
  long offset = elems.offset;
  long remaining = elems.range.size();
  MBRange::iterator iter = elems.range.begin();
  
  while (remaining)
  {
    long count = chunk_size < remaining ? chunk_size : remaining;
    remaining -= count;
  
    MBRange::iterator next = iter;
    next += count;
    rval = writeUtil->get_element_array( iter, next, elems.num_nodes, 
                                         idTag, count * elems.num_nodes, buffer );
    CHK_MB_ERR_1(rval, elem_table, status);
    iter = next;
    
    mhdf_writeConnectivity( elem_table, offset, count, 
                            id_type, buffer, &status );
    CHK_MHDF_ERR_1(status, elem_table);
    
    offset += count;
  }

  mhdf_closeData( filePtr, elem_table, &status );
  CHK_MHDF_ERR_0(status);
 
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5::get_set_info( MBEntityHandle set,
                                     long& num_entities,
                                     long& num_children,
                                     long& num_parents,
                                     unsigned long& flags )
{
  MBErrorCode rval;
  int i;
  unsigned int u;
  
  rval = iFace->get_number_entities_by_handle( set, i, false );
  CHK_MB_ERR_0(rval);
  num_entities = i;

  rval = iFace->num_child_meshsets( set, &i );
  CHK_MB_ERR_0(rval);
  num_children = i;

  rval = iFace->num_parent_meshsets( set, &i );
  CHK_MB_ERR_0(rval);
  num_parents = i;

  rval = iFace->get_meshset_options( set, u );
  CHK_MB_ERR_0(rval);
  flags = u;
  
  return MB_SUCCESS;
}


MBErrorCode WriteHDF5::write_sets( )
{
  mhdf_Status status;
  MBRange& sets = setSet.range;
  MBErrorCode rval;
  long first_id, meta_size, data_size, child_size, parent_size;
  hid_t set_table = 0, content_table = 0, child_table = 0, parent_table = -1;
  
  /* If no sets, just return success */
  if (!writeSets)
    return MB_SUCCESS;
  
  /* Write set description table and set contents table */
  
  /* Create the table */
  set_table = mhdf_openSetMeta( filePtr, &meta_size, &first_id, &status );
  CHK_MHDF_ERR_0(status);
  
  if (writeSetContents)
  {
    content_table = mhdf_openSetData( filePtr, &data_size, &status );
    CHK_MHDF_ERR_1(status, set_table);
  }
  
    
  MBRange set_contents;
  MBRange::const_iterator iter = sets.begin();
  MBRange::const_iterator comp = rangeSets.begin();
  const MBRange::const_iterator end = sets.end();
  long set_data[4];
  long set_offset = setSet.offset;
  long content_offset = setContentsOffset;
  long child_offset = setChildrenOffset;
  long parent_offset = setParentsOffset;
  unsigned long flags;
  std::vector<id_t> id_list;
  std::vector<MBEntityHandle> handle_list;
  for ( ; iter != end; ++iter )
  {
    rval = get_set_info( *iter, data_size, child_size, parent_size, flags );
    CHK_MB_ERR_2C(rval, set_table, writeSetContents, content_table, status);
    
    id_list.clear();
    if (*iter == *comp)
    {
      set_contents.clear();
      
      rval = iFace->get_entities_by_handle( *iter, set_contents, false );
      CHK_MB_ERR_2C(rval, set_table, writeSetContents, content_table, status);

      rval = range_to_id_list( set_contents, id_list );
      CHK_MB_ERR_2C(rval, set_table, writeSetContents, content_table, status);

      assert (id_list.size() < (unsigned long)data_size);
      flags |= mhdf_SET_RANGE_BIT;
      data_size = id_list.size();
      ++comp;
    }
    else
    {
      handle_list.clear();
      
      rval = iFace->get_entities_by_handle( *iter, handle_list, false );
      CHK_MB_ERR_2C(rval, set_table, writeSetContents, content_table, status);
      
      rval = vector_to_id_list( handle_list, id_list );
      CHK_MB_ERR_2C(rval, set_table, writeSetContents, content_table, status);
    }
    
    child_offset += child_size;
    parent_offset += parent_size;
    set_data[0] = content_offset + data_size - 1;
    set_data[1] = child_offset - 1;
    set_data[2] = parent_offset - 1;
    set_data[3] = flags;

    mhdf_writeSetMeta( set_table, set_offset++, 1L, H5T_NATIVE_LONG, set_data, &status );
    CHK_MHDF_ERR_2C(status, set_table, writeSetContents, content_table );
    
    if (id_list.size())
    {
      mhdf_writeSetData( content_table, 
                         content_offset,
                         id_list.size(),
                         id_type,
                         &id_list[0],
                         &status );
      CHK_MHDF_ERR_2C(status, set_table, writeSetContents, content_table );
      content_offset += data_size;
    }
  }
  
  rval = write_shared_set_descriptions( set_table );
  CHK_MB_ERR_2C(rval, set_table, writeSetContents, content_table, status);
  mhdf_closeData( filePtr, set_table, &status );
  
  if (writeSetContents)
  {
    rval = write_shared_set_contents( content_table );
    mhdf_closeData( filePtr, content_table, &status );
    CHK_MB_ERR_0( rval );
  }
  
    /* Write set children */
  if (writeSetChildren)
  {

    child_offset = setChildrenOffset;
    child_table = mhdf_openSetChildren( filePtr, &child_size, &status );
    CHK_MHDF_ERR_0(status);

    for (iter = sets.begin(); iter != end; ++iter)
    {
      handle_list.clear();
      id_list.clear();
      rval = iFace->get_child_meshsets( *iter, handle_list, 1 );
      CHK_MB_ERR_1(rval, child_table, status);

      if (handle_list.size() == 0)
        continue;

      rval = vector_to_id_list( handle_list, id_list );
      CHK_MB_ERR_1(rval, child_table, status);


      mhdf_writeSetParentsChildren( child_table, 
                                    child_offset, 
                                    id_list.size(), 
                                    id_type, 
                                    &id_list[0], 
                                    &status );
      CHK_MHDF_ERR_1(status, child_table);
      child_offset += id_list.size();
    }

    rval = write_shared_set_children( child_table );
    mhdf_closeData( filePtr, child_table, &status );
    CHK_MB_ERR_0(rval);
  }
  
    /* Write set parents */
  if (writeSetParents)
  {

    parent_offset = setParentsOffset;
    parent_table = mhdf_openSetParents( filePtr, &parent_size, &status );
    CHK_MHDF_ERR_0(status);

    for (iter = sets.begin(); iter != end; ++iter)
    {
      handle_list.clear();
      id_list.clear();
      rval = iFace->get_parent_meshsets( *iter, handle_list, 1 );
      CHK_MB_ERR_1(rval, parent_table, status);

      if (handle_list.size() == 0)
        continue;

      rval = vector_to_id_list( handle_list, id_list );
      CHK_MB_ERR_1(rval, parent_table, status);


      mhdf_writeSetParentsChildren( parent_table, 
                                    parent_offset, 
                                    id_list.size(), 
                                    id_type, 
                                    &id_list[0], 
                                    &status );
      CHK_MHDF_ERR_1(status, parent_table);
      parent_offset += id_list.size();
    }

    rval = write_shared_set_parents( parent_table );
    mhdf_closeData( filePtr, parent_table, &status );
    CHK_MB_ERR_0(rval);
  }

  return MB_SUCCESS;
}


MBErrorCode WriteHDF5::range_to_id_list( const MBRange& input_range,
                                         std::vector<id_t>& output_id_list )
{
  MBRange::const_iterator r_iter;
  MBRange::const_iterator const r_end = input_range.end();
  std::vector<id_t>::iterator i_iter, w_iter;
  MBErrorCode rval;
  
    // Get file IDs from handles
  output_id_list.resize( input_range.size() );
  rval = iFace->tag_get_data( idTag, input_range, &output_id_list[0] );
  CHK_MB_ERR_0(rval);
  std::sort( output_id_list.begin(), output_id_list.end() );
  
    // Count the number of ranges in the id list
  unsigned long count = 0;
  bool need_to_copy = false;
  std::vector<id_t>::iterator const i_end = output_id_list.end();
  i_iter = output_id_list.begin();
  while (i_iter != i_end)
  {
    ++count;
    id_t prev = *i_iter;
    for (++i_iter; (i_iter != i_end) && (++prev == *i_iter); ++i_iter);
    if (i_iter - output_id_list.begin() < (long)(2*count))
      need_to_copy = true;
  }
  
    // If the range format is larger than half the size of the
    // the simple list format, just keep the list format
  if (4*count >= output_id_list.size())
    return MB_SUCCESS;
  
    // Convert to ranged format
  std::vector<id_t>* range_list = &output_id_list;
  if (need_to_copy)
    range_list = new std::vector<id_t>( 2*count );

  w_iter = range_list->begin();
  i_iter = output_id_list.begin();
  while (i_iter != i_end)
  {
    unsigned long range_size = 1;
    id_t prev = *w_iter = *i_iter;
    w_iter++;
    for (++i_iter; (i_iter != i_end) && (++prev == *i_iter); ++i_iter)
      ++range_size;
    *w_iter = range_size;
    ++w_iter;
  }

  if (need_to_copy)
  {
    std::swap( *range_list, output_id_list );
    delete range_list;
  }
  else
  {
    assert( w_iter - output_id_list.begin() == (long)(2*count) );
    output_id_list.resize( 2*count );
  }
  
  return MB_SUCCESS;
}
 
MBErrorCode WriteHDF5::vector_to_id_list( 
                                 const std::vector<MBEntityHandle>& input,
                                 std::vector<id_t>& output )
{
  id_t id;
  MBErrorCode rval;
  MBEntityHandle handle;
  std::vector<MBEntityHandle>::const_iterator i_iter = input.begin();
  const std::vector<MBEntityHandle>::const_iterator i_end = input.end();
  output.resize(input.size());
  std::vector<id_t>::iterator o_iter = output.begin();
  
  while (i_iter != i_end)
  {
    handle = *i_iter;
    rval = iFace->tag_get_data( idTag, &handle, 1, &id );
    CHK_MB_ERR_0(rval);
    *o_iter = id;
    ++o_iter;
    ++i_iter;
  }

  return MB_SUCCESS;
}


template <class T> static void 
erase_from_vector( std::vector<T>& vector, T value )
{
  typename std::vector<T>::iterator r_itor, w_itor;
  const typename std::vector<T>::iterator end = vector.end();
  unsigned count = 0;
  
  for (r_itor = vector.begin(); 
       r_itor != end && *r_itor != value; 
       ++r_itor)
    count++;
  
  w_itor = r_itor;
  for ( ; r_itor != end; ++r_itor)
  {
    if (*r_itor != value)
    {
      *w_itor = *r_itor;
      ++w_itor;
      count++;
    }
  }
  
  vector.resize( count );
}
 

  

inline MBErrorCode WriteHDF5::get_adjacencies( MBEntityHandle entity,
                                        std::vector<id_t>& adj )
{
  MBErrorCode rval = writeUtil->get_adjacencies( entity, idTag, adj );
  //erase_from_vector( adj, (id_t)0 );
  erase_from_vector( adj, (id_t)-1 );
  return rval;
}


MBErrorCode WriteHDF5::write_adjacencies( const ExportSet& elements )
{
  MBErrorCode rval;
  mhdf_Status status;
  MBRange::const_iterator iter;
  const MBRange::const_iterator end = elements.range.end();
  std::vector<int> adj_list;
  
  /* Count Adjacencies */
  long count = 0;
  //for (iter = elements.range.begin(); iter != end; ++iter)
  //{
  //  adj_list.clear();
  //  rval = get_adjacencies( *iter, adj_list);
  //  CHK_MB_ERR_0(rval);
  //
  //  if (adj_list.size() > 0)
  //    count += adj_list.size() + 2;
  //}
  
  //if (count == 0)
  //  return MB_SUCCESS;

  long offset = elements.adj_offset;
  if (offset < 0)
    return MB_SUCCESS;
  
  /* Create data list */
  hid_t table = mhdf_openAdjacency( filePtr, elements.name(), &count, &status );
  CHK_MHDF_ERR_0(status);
  
  /* Write data */
  id_t* buffer = (id_t*)dataBuffer;
  long chunk_size = bufferSize / sizeof(id_t); 
  count = 0;
  for (iter = elements.range.begin(); iter != end; ++iter)
  {
    adj_list.clear();
    rval = get_adjacencies( *iter, adj_list );
    CHK_MB_ERR_1(rval, table, status);
    if (adj_list.size() == 0)
      continue;
    
    if (count + adj_list.size() + 2 > (unsigned long)chunk_size)
    {
      mhdf_writeAdjacency( table, offset, count, id_type, buffer, &status );
      CHK_MHDF_ERR_1(status, table);
      
      offset += count;
      count = 0;
    }
    
    rval = iFace->tag_get_data( idTag, &*iter, 1, buffer + count );
    if (MB_SUCCESS != rval || buffer[count] < 1)
    {
      mhdf_closeData( filePtr, table, &status );
      return MB_FAILURE;
    }
    buffer[++count] = adj_list.size();
    ++count;
    
    assert (adj_list.size()+2 < (unsigned long)chunk_size);
    memcpy( buffer + count, &adj_list[0], adj_list.size() * sizeof(id_t) );
    count += adj_list.size();
  }
  
  if (count)
  {
    mhdf_writeAdjacency( table, offset, count, id_type, buffer, &status );
    CHK_MHDF_ERR_1(status, table);

    offset += count;
    count = 0;
  }
  
  mhdf_closeData( filePtr, table, &status );
  CHK_MHDF_ERR_0(status);
  
  return MB_SUCCESS;
}

/*

MBErrorCode WriteHDF5::write_tag( MBTag tag_handle )
{
  MBErrorCode rval;
  MBTagType tag_type;
  MBTag type_handle;
  int tag_size, mem_size;
  mhdf_Status status;
  hid_t hdf_tag_type;
  
  rval = iFace->tag_get_type( tag_handle, tag_type ); CHK_MB_ERR_0(rval);
  rval = iFace->tag_get_size( tag_handle, tag_size ); CHK_MB_ERR_0(rval);
  
  bool sparse = true;
  bool have_type = false;
  std::string tag_type_name = "__hdf5_tag_type_";
  std::string tag_name;
  rval = iFace->tag_get_name( tag_handle, tag_name ); CHK_MB_ERR_0(rval);
  
  tag_type_name += tag_name;
  rval = iFace->tag_get_handle( tag_type_name.c_str(), type_handle );
  if (MB_SUCCESS == rval)
  {
    rval = iFace->tag_get_data( type_handle, 0, 0, &hdf_tag_type );
    if (rval != MB_SUCCESS)
      return rval;
    have_type = true;
  }
  else if (MB_TAG_NOT_FOUND != rval)
    return rval;
  
  mem_size = tag_size;
  switch ( tag_type )
  {
    case MB_TAG_BIT:
      sparse = true;
      assert( tag_size < 9 );
      mem_size = 1;
      break;
    case MB_TAG_SPARSE:
      sparse = true;
      break;
    case MB_TAG_DENSE:
      sparse = false;
      break;
    case MB_TAG_MESH:
      sparse = true;
      break;
    default:
      return MB_FAILURE;
  }
  
  assert( 2*tag_size + sizeof(long) < (unsigned long)bufferSize );
  bool have_default = true;
  rval = iFace->tag_get_default_value( tag_handle, dataBuffer );
  if (MB_ENTITY_NOT_FOUND == rval)
    have_default = false;
  else if (MB_SUCCESS != rval)
    return rval;
  rval = iFace->tag_get_data( tag_handle, 0, 0, dataBuffer + mem_size );
  bool have_global = true;
  if (MB_TAG_NOT_FOUND == rval)
    have_global = false;
  else if (MB_SUCCESS != rval)
    return rval;
  
  if (have_type)
  {
    mhdf_createTypeTag( filePtr, tag_name.c_str(),
                        hdf_tag_type, have_default ? dataBuffer : 0, 
                        have_global ? dataBuffer + mem_size : 0,
                        tag_type, &status );
  }
  else if (MB_TAG_BIT == tag_type)
  {
    mhdf_createBitTag( filePtr, tag_name.c_str(), 
                       tag_size, have_default ? dataBuffer : 0,
                       have_global ? dataBuffer + mem_size : 0,
                       tag_type, &status );
    hdf_tag_type = H5T_NATIVE_B8;
  }  
  else 
  {
    mhdf_createOpaqueTag( filePtr, tag_name.c_str(),
                          tag_size, have_default ? dataBuffer : 0,
                          have_global ? dataBuffer + mem_size : 0,
                          tag_type, &status );
    hdf_tag_type = 0;
  }

  CHK_MHDF_ERR_0(status);
  
  // FIX ME
  // Always write tags as sparse to work around MOAB issues
  //   with dense tags.  (Can't determine which entities tag
  //   is actually set for.)
  //if (sparse)
  //  rval = write_sparse_tag( tag_handle, hdf_tag_type );
  //else
  //  rval = write_dense_tag( tag_handle, hdf_tag_type );
  //
  rval = write_sparse_tag( tag_handle hdf_tag_type );
    
  return rval;
}


MBErrorCode WriteHDF5::write_dense_tag( MBTag handle,
                                        hid_t type )
{
  MBErrorCode rval = MB_SUCCESS;
  
  if (!nodeSet.range.empty())
    rval = write_dense_tag( nodeSet, handle, type );
  CHK_MB_ERR_0(rval);
  
  std::list<ExportSet>::iterator iter, end = exportList.end();
  for (iter = exportList.begin(); iter != end; ++iter)
  {
    MBErrorCode rval = write_dense_tag( *iter, handle, type );
    CHK_MB_ERR_0(rval);
  }
  
  if (!setSet.range.empty())
    rval = write_dense_tag( setSet, handle, type );
  return rval;
}

MBErrorCode WriteHDF5::write_dense_tag( ExportSet& set,
                                        MBTag handle,
                                        hid_t type )
{
  MBRange sub_range;
  MBErrorCode rval;
  mhdf_Status status;
  hid_t data_handle;
  std::string name;
  int tag_size;
  MBTagType mb_type;
  
    //get tag properties
  if (MB_SUCCESS != iFace->tag_get_name( handle, name )    ||
      MB_SUCCESS != iFace->tag_get_type( handle, mb_type ) ||
      MB_SUCCESS != iFace->tag_get_size( handle, tag_size ))
    return MB_FAILURE;
  
  if (mb_type == MB_TAG_BIT)
    tag_size = 1;
  assert( type == 0 || H5Tget_size(type) == (unsigned)tag_size );

  data_handle = mhdf_createDenseTagData( filePtr, name.c_str(), set.type2, set.range.size(), &status );
  CHK_MHDF_ERR_0(status);
  
  long chunk_size = bufferSize / tag_size;
  long offset = 0;
  
  MBRange::const_iterator iter = set.range.begin();
  long remaining = set.range.size();
  while (remaining)
  {
    long count = remaining > chunk_size ? chunk_size : remaining;
    MBRange::const_iterator next = iter;
    next += count;
    sub_range.clear();
    sub_range.merge( iter, next );
    iter = next;
    remaining -= count;
    
    rval = iFace->tag_get_data( handle, sub_range, dataBuffer );
    if (MB_TAG_NOT_FOUND == rval)
    {
        // Dense tag that doesn't have a default value -- use zero.
      memset( dataBuffer, 0, bufferSize );
    }
    else CHK_MB_ERR_1( rval, data_handle, status );
    
    mhdf_writeDenseTag( data_handle, offset, count, type, dataBuffer, &status );
    CHK_MHDF_ERR_1( status, data_handle );
    
    offset += count;
  }
  
  mhdf_closeData( filePtr, data_handle, &status );
  CHK_MHDF_ERR_0(status);
  
  return MB_SUCCESS;
}

*/

MBErrorCode WriteHDF5::write_sparse_tag( const SparseTag& tag_data )
{
  MBErrorCode rval;
  mhdf_Status status;
  hid_t tables[2];
  std::string name;
  int mb_size;
  MBTagType mb_type;
  MBDataType mb_data_type;
  long table_size;
  
    //get tag properties from moab
  if (MB_SUCCESS != iFace->tag_get_name( tag_data.tag_id, name )    ||
      MB_SUCCESS != iFace->tag_get_type( tag_data.tag_id, mb_type ) ||
      MB_SUCCESS != iFace->tag_get_size( tag_data.tag_id, mb_size ) ||
      MB_SUCCESS != iFace->tag_get_data_type( tag_data.tag_id, mb_data_type ))
    return MB_FAILURE;
  if (mb_type == MB_TAG_BIT)
    mb_size = 1;

DEBUGOUT((std::string("Tag: ") + name + "\n").c_str());
  
    //open tables to write info
  mhdf_openSparseTagData( filePtr,
                          name.c_str(),
                          &table_size,
                          tables,
                          &status);
  CHK_MHDF_ERR_0(status);
  assert( tag_data.range.size() + tag_data.offset <= (unsigned long)table_size );

    // Set up data buffer for writing IDs
  size_t chunk_size = bufferSize / sizeof(id_t);
  id_t* id_buffer = (id_t*)dataBuffer;
  
    // Write IDs of tagged entities.
  MBRange range;
  long remaining = tag_data.range.size();
  long offset = tag_data.offset;
  MBRange::const_iterator iter = tag_data.range.begin();
  while (remaining)
  {
      // write "chunk_size" blocks of data
    long count = (unsigned long)remaining > chunk_size ? chunk_size : remaining;
    remaining -= count;
    MBRange::const_iterator stop = iter;
    stop += count;
    range.clear();
    range.merge( iter, stop );
    iter = stop;
    assert(range.size() == (unsigned)count);
    
    rval = iFace->tag_get_data( idTag, range, id_buffer );
    CHK_MB_ERR_2( rval, tables, status );
    
      // write the data
    mhdf_writeSparseTagEntities( tables[0], offset, count, id_type, 
                                 id_buffer, &status );
    CHK_MHDF_ERR_2( status, tables );
   
    offset += count;
  } // while (remaining)
  mhdf_closeData( filePtr, tables[0], &status );
  CHK_MHDF_ERR_0(status);
  
    // Set up data buffer for writing tag values
  chunk_size = bufferSize / mb_size;
  assert( chunk_size > 0 );
  char* tag_buffer = (char*)dataBuffer;
  
    // Write the tag values
  remaining = tag_data.range.size();
  offset = tag_data.offset;
  iter = tag_data.range.begin();
  while (remaining)
  {
      // write "chunk_size" blocks of data
    long count = (unsigned long)remaining > chunk_size ? chunk_size : remaining;
    remaining -= count;
    memset( tag_buffer, 0, count * mb_size );
    MBRange::const_iterator stop = iter;
    stop += count;
    range.clear();
    range.merge( iter, stop );
    iter = stop;
    assert(range.size() == (unsigned)count);
 
 /** Fix me - stupid API requires we get these one at a time for BIT tags */
    if (mb_type == MB_TAG_BIT)
    {
      rval = MB_SUCCESS;
      char* buf_iter = tag_buffer;
      for (MBRange::const_iterator it = range.begin(); 
           MB_SUCCESS == rval && it != range.end(); 
           ++it, buf_iter += mb_size)
        rval = iFace->tag_get_data( tag_data.tag_id, &*it, 1, buf_iter );
    }
    else
    {
      rval = iFace->tag_get_data( tag_data.tag_id, range, tag_buffer );
    }
    CHK_MB_ERR_1(rval, tables[1], status);
    
      // Convert MBEntityHandles to file ids
    if (mb_data_type == MB_TYPE_HANDLE)
      convert_handle_tag( iFace, idTag, tag_buffer, count * mb_size / sizeof(MBEntityHandle) );
    
      // write the data
    mhdf_writeSparseTagValues( tables[1], offset, count,
                               0, tag_buffer, &status );
    CHK_MHDF_ERR_1(status, tables[1]);
   
    offset += count;
  } // while (remaining)
  
  mhdf_closeData( filePtr, tables[1], &status );
  CHK_MHDF_ERR_0(status);
  
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5::write_qa( std::vector<std::string>& list )
{
  const char* app = "MOAB";
  const char* vers = "0";
  char date_str[64];
  char time_str[64];
  
  std::vector<const char*> strs(list.size() ? list.size() : 4);
  if (list.size() == 0)
  {
    time_t t = time(NULL);
    tm* lt = localtime( &t );
    strftime( date_str, sizeof(date_str), "%D", lt );
    strftime( time_str, sizeof(time_str), "%T", lt );
    
    strs[0] = app;
    strs[1] = vers;
    strs[2] = date_str;
    strs[3] = time_str;
  }
  else
  {
    for (unsigned int i = 0; i < list.size(); ++i)
      strs[i] = list[i].c_str();
  }
  
  mhdf_Status status;
  mhdf_writeHistory( filePtr, &strs[0], strs.size(), &status );
  CHK_MHDF_ERR_0(status);
  
  return MB_SUCCESS;
}

/*
MBErrorCode WriteHDF5::register_known_tag_types( MBInterface* iface )
{
  hid_t int4, double16;
  hsize_t dim[1];
  int error = 0;
  MBErrorCode rval;
  
  dim[0] = 4;
  int4 = H5Tarray_create( H5T_NATIVE_INT, 1, dim, NULL );
  
  dim[0] = 16;
  double16 = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, dim, NULL );
  
  if (int4 < 0 || double16 < 0)
    error = 1;
  
  struct { const char* name; hid_t type; } list[] = {
    { GLOBAL_ID_TAG_NAME, H5T_NATIVE_INT } ,
    { MATERIAL_SET_TAG_NAME, H5T_NATIVE_INT },
    { DIRICHLET_SET_TAG_NAME, H5T_NATIVE_INT },
    { NEUMANN_SET_TAG_NAME, H5T_NATIVE_INT },
    { HAS_MID_NODES_TAG_NAME, int4 },
    { GEOM_DIMENSION_TAG_NAME, H5T_NATIVE_INT },
    { MESH_TRANSFORM_TAG_NAME, double16 },
    { 0, 0 } };
  
  for (int i = 0; list[i].name; ++i)
  {
    if (list[i].type < 1)
      { ++error; continue; }
    
    MBTag handle;
    
    std::string name("__hdf5_tag_type_");
    name += list[i].name;
    
    rval = iface->tag_get_handle( name.c_str(), handle );
    if (MB_TAG_NOT_FOUND == rval)
    {
      rval = iface->tag_create( name.c_str(), sizeof(hid_t), MB_TAG_SPARSE, handle, NULL );
      if (MB_SUCCESS != rval)
        { ++error; continue; }
      
      hid_t copy_id = H5Tcopy( list[i].type );
      rval = iface->tag_set_data( handle, 0, 0, &copy_id );
      if (MB_SUCCESS != rval)
        { ++error; continue; }
    }
  }
  
  H5Tclose( int4 );
  H5Tclose( double16 );
  return error ? MB_FAILURE : MB_SUCCESS;
}
*/

MBErrorCode WriteHDF5::gather_tags()
{
  MBErrorCode result;
  std::string tagname;
  std::vector<MBTag> tag_list;
  std::vector<MBTag>::iterator t_itor;
  MBRange range;
    
    // Get list of Tags to write
  result = iFace->tag_get_tags( tag_list );
  CHK_MB_ERR_0(result);

    // Get list of tags
  for (t_itor = tag_list.begin(); t_itor != tag_list.end(); ++t_itor)
  {
      // Don't write global ID tag
    if (*t_itor == idTag)
      continue;
    
      // Don't write tags that have name beginning with "__"
    result = iFace->tag_get_name( *t_itor, tagname );
    if (MB_SUCCESS != result)
      return result;
      // Skip anonymous tags
    if (tagname.empty())
      continue;
      // skip tags for which the name begins with two underscores
    if (tagname.size() >= 2 && tagname[0] == '_' && tagname[1] == '_')
      continue;
  
      // Add tag to export list
    SparseTag tag_data;
    tag_data.tag_id = *t_itor;
    tag_data.offset = 0;
    tagList.push_back( tag_data );
  }
  
    // Get entities for each tag
  std::list<SparseTag>::iterator td_iter = tagList.begin();
  const std::list<SparseTag>::iterator td_end = tagList.end();
  for ( ; td_iter != td_end; ++td_iter)
  {  
    MBTag handle = td_iter->tag_id;
      // Get list of entities for which tag is set
    std::list<ExportSet>::iterator e_iter, e_end = exportList.end();
    for (e_iter = exportList.begin(); e_iter != e_end; ++e_iter)
    {
      range.clear();
      result = iFace->get_entities_by_type_and_tag( 0, e_iter->type, &handle, NULL, 1, range );
      CHK_MB_ERR_0(result);
      td_iter->range.merge( range.intersect( e_iter->range ) );
    }
    
    range.clear();
    result = iFace->get_entities_by_type_and_tag( 0, MBVERTEX, &handle, NULL, 1, range );
    CHK_MB_ERR_0(result);
    td_iter->range.merge( range.intersect( nodeSet.range ) );


    range.clear();
    result = iFace->get_entities_by_type_and_tag( 0, MBENTITYSET, &handle, NULL, 1, range );
    CHK_MB_ERR_0(result);
    td_iter->range.merge( range.intersect( setSet.range ) );
    
      // For tags containing entity handles, skip values if
      // handle doesn't reference something being written to the file.
      // If the tag contains multiple handle values, write it if any one
      // of those handles is valid.  Consider special case of 0 handle as
      // valid.
    MBDataType data_type;
    iFace->tag_get_data_type( handle, data_type );
    if (MB_TYPE_HANDLE == data_type)
    {
      int tag_size;
      result = iFace->tag_get_size( handle, tag_size );
      CHK_MB_ERR_0(result);
      if (tag_size % sizeof(MBEntityHandle)) // not an even multiple?
        td_iter->range.clear(); // don't write any values
      
      std::vector<MBEntityHandle> values(tag_size / sizeof(MBEntityHandle));
      MBRange::iterator i = td_iter->range.begin();
      while (i != td_iter->range.end())
      {
        result = iFace->tag_get_data( handle, &*i, 1, &values[0] );
        CHK_MB_ERR_0(result);

        if (convert_handle_tag( iFace, idTag, &values[0], values.size() ))
          ++i;
        else
          i = td_iter->range.erase( i );
      }
    }
  
    td_iter->write = !td_iter->range.empty();
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5::create_file( const char* filename,
                                    bool overwrite,
                                    std::vector<std::string>& qa_records,
                                    int dimension )
{
  long first_id;
  mhdf_Status status;
  hid_t handle;
  std::list<ExportSet>::iterator ex_itor;
  MBErrorCode rval;
  
  const char* type_names[MBMAXTYPE];
  memset( type_names, 0, MBMAXTYPE * sizeof(char*) );
  for (MBEntityType i = MBEDGE; i < MBENTITYSET; ++i)
    type_names[i] = MBCN::EntityTypeName( i );
 
    // Create the file
  filePtr = mhdf_createFile( filename, overwrite, type_names, MBMAXTYPE, &status );
  CHK_MHDF_ERR_0(status);
  assert(!!filePtr);

  rval = write_qa( qa_records );
  CHK_MB_ERR_0(rval);
  
    // Create node table
  if (nodeSet.range.size()) {
    handle = mhdf_createNodeCoords( filePtr, dimension, nodeSet.range.size(), &first_id, &status );
    CHK_MHDF_ERR_0(status);
    mhdf_closeData( filePtr, handle, &status );
    CHK_MHDF_ERR_0(status);
    writeUtil->assign_ids( nodeSet.range, idTag, (id_t)first_id );
    nodeSet.first_id = (id_t)first_id;
  }
  else {
    nodeSet.first_id = std::numeric_limits<id_t>::max();
  } 
  nodeSet.offset = 0;

    // Create element tables
  for (ex_itor = exportList.begin(); ex_itor != exportList.end(); ++ex_itor)
  {
    rval = create_elem_tables( ex_itor->type,
                               ex_itor->num_nodes,
                               ex_itor->range.size(),
                               first_id );
    CHK_MB_ERR_0(rval);
      
    writeUtil->assign_ids( ex_itor->range, idTag, (id_t)first_id );
    ex_itor->first_id = (id_t)first_id;
    ex_itor->offset = 0;
  }

    // create node adjacency table
  id_t num_adjacencies;
#ifdef WRITE_NODE_ADJACENCIES  
  rval = count_adjacencies( nodeSet.range, num_adjacencies );
  CHK_MB_ERR_0(rval);
  if (num_adjacencies > 0)
  {
    handle = mhdf_createAdjacency( filePtr,
                                   mhdf_node_type_handle(),
                                   num_adjacencies,
                                   &status );
    CHK_MHDF_ERR_0(status);
    mhdf_closeData( filePtr, handle, &status );
    nodeSet.adj_offset = 0;
  }
  else
    nodeSet.adj_offset = -1;
#endif
  
    // create element adjacency tables
  for (ex_itor = exportList.begin(); ex_itor != exportList.end(); ++ex_itor)
  {
    rval = count_adjacencies( ex_itor->range, num_adjacencies );
    CHK_MB_ERR_0(rval);
    
    if (num_adjacencies > 0)
    {
      handle = mhdf_createAdjacency( filePtr,
                                     ex_itor->name(),
                                     num_adjacencies,
                                     &status );
      CHK_MHDF_ERR_0(status);
      mhdf_closeData( filePtr, handle, &status );
      ex_itor->adj_offset = 0;
    }
    else
      ex_itor->adj_offset = -1;
  }
  
    // create set tables
  writeSets = !setSet.range.empty();
  if (writeSets)
  {
    long contents_len, children_len, parents_len;
    writeSets = true;
    
    rval = create_set_meta( setSet.range.size(), first_id );
    writeUtil->assign_ids( setSet.range, idTag, (id_t)first_id );
    
    rval = count_set_size( setSet.range, rangeSets, contents_len, children_len, parents_len );
    CHK_MB_ERR_0(rval);
    
    rval = create_set_tables( contents_len, children_len, parents_len );
    CHK_MB_ERR_0(rval);
   
    setSet.first_id = (id_t)first_id;
    setSet.offset = 0;
    setContentsOffset = 0;
    setChildrenOffset = 0;
    setParentsOffset = 0;
    writeSetContents = !!contents_len;
    writeSetChildren = !!children_len;
    writeSetParents = !!parents_len;
  } // if(!setSet.range.empty())
  
  
DEBUGOUT( "Gathering Tags\n" );
  
  rval = gather_tags();
  CHK_MB_ERR_0(rval);

    // Create the tags and tag data tables
  std::list<SparseTag>::iterator tag_iter = tagList.begin();
  const std::list<SparseTag>::iterator tag_end = tagList.end();
  for ( ; tag_iter != tag_end; ++tag_iter)
  {
    rval = create_tag( tag_iter->tag_id, tag_iter->range.size() );
    CHK_MB_ERR_0(rval);
  } // for(tags)
  
  return MB_SUCCESS;
}


MBErrorCode WriteHDF5::count_adjacencies( const MBRange& set, id_t& result )
{
  MBErrorCode rval;
  std::vector<int> adj_list;
  MBRange::const_iterator iter = set.begin();
  const MBRange::const_iterator end = set.end();
  result = 0;
  for ( ; iter != end; ++iter )
  {
    adj_list.clear();
    rval = get_adjacencies( *iter, adj_list );
    CHK_MB_ERR_0(rval);
    
    if (adj_list.size() > 0)
      result += 2 + adj_list.size();
  }
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5::create_elem_tables( MBEntityType mb_type,
                                           int nodes_per_elem,
                                           id_t num_elements,
                                           long& first_id_out )
{
  char name[64];
  mhdf_Status status;
  hid_t handle;
  
  sprintf( name, "%s%d", MBCN::EntityTypeName(mb_type), nodes_per_elem );
  mhdf_addElement( filePtr, name, mb_type, &status );
  CHK_MHDF_ERR_0(status);
  
  handle = mhdf_createConnectivity( filePtr, 
                                    name,
                                    nodes_per_elem,
                                    num_elements,
                                    &first_id_out,
                                    &status );
  CHK_MHDF_ERR_0(status);
  mhdf_closeData( filePtr, handle, &status );
  CHK_MHDF_ERR_0(status);
  
  return MB_SUCCESS;
}


MBErrorCode WriteHDF5::count_set_size( const MBRange& sets, 
                                       MBRange& compressed_sets,
                                       long& contents_length_out,
                                       long& children_length_out,
                                       long& parents_length_out )
{
  MBErrorCode rval;
  MBRange set_contents;
  MBRange::const_iterator iter = sets.begin();
  const MBRange::const_iterator end = sets.end();
  long contents_length_set, children_length_set, parents_length_set;
  unsigned long flags;
  std::vector<id_t> set_contents_ids;
  
  contents_length_out = 0;
  children_length_out = 0;
  parents_length_out = 0;
  
  for (; iter != end; ++iter)
  {
    rval = get_set_info( *iter, contents_length_set, children_length_set,
                         parents_length_set, flags );
    CHK_MB_ERR_0(rval);
    
      // check if can and should compress as ranges
    if ((flags&MESHSET_SET) && !(flags&MESHSET_ORDERED))
    {
      set_contents.clear();
      rval = iFace->get_entities_by_handle( *iter, set_contents, false );
      CHK_MB_ERR_0(rval);
      
      rval = range_to_id_list( set_contents, set_contents_ids );
      CHK_MB_ERR_0(rval);
      
      if (set_contents_ids.size() < (unsigned long)contents_length_set)
      {
        contents_length_set = set_contents_ids.size();
        compressed_sets.insert( *iter );
      }
    }
    
    contents_length_out += contents_length_set;
    children_length_out += children_length_set;
    parents_length_out += parents_length_set;
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5::create_set_meta( id_t num_sets, long& first_id_out )
{
  hid_t handle;
  mhdf_Status status;
  
  handle = mhdf_createSetMeta( filePtr, num_sets, &first_id_out, &status );
  CHK_MHDF_ERR_0(status);
  mhdf_closeData( filePtr, handle, &status );
  
  return MB_SUCCESS;
}


MBErrorCode WriteHDF5::create_set_tables( long num_set_contents,
                                          long num_set_children,
                                          long num_set_parents )
{
  hid_t handle;
  mhdf_Status status;
  
  if (num_set_contents > 0)
  {
    handle = mhdf_createSetData( filePtr, num_set_contents, &status );
    CHK_MHDF_ERR_0(status);
    mhdf_closeData( filePtr, handle, &status );
  }
  
  if (num_set_children > 0)
  {
    handle = mhdf_createSetChildren( filePtr, num_set_children, &status );
    CHK_MHDF_ERR_0(status);
    mhdf_closeData( filePtr, handle, &status );
  }
  
  if (num_set_parents > 0)
  {
    handle = mhdf_createSetParents( filePtr, num_set_parents, &status );
    CHK_MHDF_ERR_0(status);
    mhdf_closeData( filePtr, handle, &status );
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5::create_tag( MBTag tag_id, id_t num_sparse_entities )
{
  MBTagType storage;
  MBDataType mb_type;
  MBTag type_handle;
  mhdf_TagDataType mhdf_type;
  int tag_size, elem_size = 0, mhdf_size;
  hid_t hdf_type = (hid_t)0;
  hid_t handles[2];
  std::string tag_name, tag_type_name = "__hdf5_tag_type_";
  MBErrorCode rval;
  mhdf_Status status;
  

    // get tag properties
  rval = iFace->tag_get_type( tag_id, storage  ); CHK_MB_ERR_0(rval);
  rval = iFace->tag_get_size( tag_id, tag_size ); CHK_MB_ERR_0(rval);
  rval = iFace->tag_get_name( tag_id, tag_name ); CHK_MB_ERR_0(rval);
  rval = iFace->tag_get_data_type( tag_id, mb_type ); CHK_MB_ERR_0(rval);
  
  
    // get type-specific parameters
  if (MB_TAG_BIT == storage)
  {
    mhdf_type = mhdf_BITFIELD;
  }
  else 
  {
    switch (mb_type)
    {
    case MB_TYPE_INTEGER:
      elem_size = sizeof(int);
      mhdf_type = mhdf_INTEGER;
      break;
    case MB_TYPE_DOUBLE:
      elem_size = sizeof(double);
      mhdf_type = mhdf_FLOAT;
      break;
    case MB_TYPE_BIT:
      elem_size = sizeof(bool);
      mhdf_type = mhdf_BOOLEAN;
      break;
    case MB_TYPE_HANDLE:
      elem_size = sizeof(MBEntityHandle);
      mhdf_type = mhdf_ENTITY_ID;
      break;
    case MB_TYPE_OPAQUE:
    default:
      mhdf_type = mhdf_OPAQUE;
      
      tag_type_name = "__hdf5_tag_type_";
      tag_type_name += tag_name;
      rval = iFace->tag_get_handle( tag_type_name.c_str(), type_handle );
      if (MB_SUCCESS == rval)
      {
        rval = iFace->tag_get_data( type_handle, 0, 0, &hdf_type );
        if (rval != MB_SUCCESS || H5Tget_size(hdf_type) != (unsigned)tag_size) 
          return MB_FAILURE;
       }
      else if(MB_TAG_NOT_FOUND != rval)
        return rval;
    }
  }
    
    // if a basic type, check if it is an array of them
  if (elem_size)
  {
    if (tag_size % elem_size)  // tag_size must be a multiple of elem_size
      return MB_FAILURE;
    mhdf_size = tag_size / elem_size;
  }
  else
  {
    mhdf_size = tag_size;
  }
  
  
    // check for default and global/mesh values
  assert( 2*tag_size + sizeof(long) < (unsigned long)bufferSize );
 
  bool have_default = false;
  rval = iFace->tag_get_default_value( tag_id, dataBuffer );
  if (MB_SUCCESS == rval) {
    have_default = true;
    if (mb_type == MB_TYPE_HANDLE) 
      have_default = convert_handle_tag( iFace, idTag, dataBuffer, mhdf_size );
  }
  else if(MB_ENTITY_NOT_FOUND != rval)
    return rval;

  bool have_global = false;
  rval = iFace->tag_get_data( tag_id, 0, 0, dataBuffer + tag_size );
  if (MB_SUCCESS == rval) {
    have_global = true;
    if (mb_type == MB_TYPE_HANDLE) 
      have_default = convert_handle_tag( iFace, idTag, dataBuffer+tag_size, mhdf_size );
  }
  else if(MB_TAG_NOT_FOUND != rval)
    return rval;


    // write the tag description to the file
  mhdf_createTag( filePtr,
                  tag_name.c_str(),
                  mhdf_type,
                  mhdf_size,
                  storage,
                  have_default ? dataBuffer : 0,
                  have_global ? dataBuffer + tag_size : 0,
                  hdf_type,
                  &status );
  CHK_MHDF_ERR_0(status);

  
    // create empty table for tag data
  if (num_sparse_entities)
  {
    mhdf_createSparseTagData( filePtr, 
                              tag_name.c_str(), 
                              num_sparse_entities,
                              handles,
                              &status );
    CHK_MHDF_ERR_0(status);
    mhdf_closeData( filePtr, handles[0], &status );
    mhdf_closeData( filePtr, handles[1], &status );
  }
  
  return MB_SUCCESS;
}
