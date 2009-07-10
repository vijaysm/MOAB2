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
// Filename      : ReadHDF5.cpp
//
// Purpose       : TSTT HDF5 Writer 
//
// Special Notes : WriteSLAC used as template for this
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 04/18/04
//-------------------------------------------------------------------------

#include <assert.h>
#include <H5Tpublic.h>
#include <H5Ppublic.h>
#include "MBInterface.hpp"
#include "MBInternals.hpp"
#include "MBTagConventions.hpp"
#include "ReadHDF5.hpp"
#include "MBCN.hpp"
#include "FileOptions.hpp"
#ifdef HDF5_PARALLEL
#include "ReadParallel.hpp"
#include <H5FDmpi.h>
#include <H5FDmpio.h>
#endif
//#include "WriteHDF5.hpp"

#include <stdlib.h>
#include <string.h>
#include <limits>
#include <functional>

#undef DEBUG

#ifdef DEBUG
#  define DEBUGOUT(A) fputs( A, stderr )
#  include <stdio.h>
#else
#  define DEBUGOUT(A)
#endif

#define READ_HDF5_BUFFER_SIZE (40*1024*1024)

// This function doesn't do anything useful.  It's just a nice
// place to set a break point to determine why the reader fails.
static inline MBErrorCode error( MBErrorCode rval )
  { return rval; }

static void copy_sorted_file_ids( const MBEntityHandle* sorted_ids, 
                                  long num_ids,
                                  MBRange& results )
{
  MBRange::iterator hint = results.begin();
  long i = 0;
  while (i < num_ids) {
    MBEntityHandle start = sorted_ids[i];
    for (++i; i < num_ids && sorted_ids[i] == 1+sorted_ids[i-1]; ++i);
    hint = results.insert( hint, start, sorted_ids[i-1] );
  }
}

static void intersect( const mhdf_EntDesc& group, const MBRange& range, MBRange& result )
{
  MBRange::const_iterator s, e;
  s = MBRange::lower_bound( range.begin(), range.end(), group.start_id );
  e = MBRange::lower_bound( s, range.end(), group.start_id + group.count );
  result.merge( s, e );
}

MBReaderIface* ReadHDF5::factory( MBInterface* iface )
  { return new ReadHDF5( iface ); }

ReadHDF5::ReadHDF5( MBInterface* iface )
  : bufferSize( READ_HDF5_BUFFER_SIZE ),
    dataBuffer( 0 ),
    iFace( iface ), 
    filePtr( 0 ), 
    fileInfo( 0 ), 
    readUtil( 0 ),
    handleType( 0 ),
    ioProp( H5P_DEFAULT )
{
}

MBErrorCode ReadHDF5::init()
{
  MBErrorCode rval;

  if (readUtil) 
    return MB_SUCCESS;
  
  ioProp = H5P_DEFAULT;
  //WriteHDF5::register_known_tag_types( iFace );
  
  handleType = H5Tcopy( H5T_NATIVE_ULONG );
  if (handleType < 0)
    return error(MB_FAILURE);
  
  if (H5Tset_size( handleType, sizeof(MBEntityHandle)) < 0)
  {
    H5Tclose( handleType );
    return error(MB_FAILURE);
  }
  
  void* ptr = 0;
  rval = iFace->query_interface( "MBReadUtilIface", &ptr );
  if (MB_SUCCESS != rval)
  {
    H5Tclose( handleType );
    return error(rval);
  }
  readUtil = reinterpret_cast<MBReadUtilIface*>(ptr);
  
  idMap.clear();
  fileInfo = 0;
  
  return MB_SUCCESS;
}
  

ReadHDF5::~ReadHDF5()
{
  if (!readUtil) // init() failed.
    return;

  iFace->release_interface( "MBReadUtilIface", readUtil );
  H5Tclose( handleType );
}

MBErrorCode ReadHDF5::load_file( const char* filename, 
                                 MBEntityHandle& file_set, 
                                 const FileOptions& opts,
                                 const char* name,
                                 const int* id_list, 
                                 const int num_ids )
{
  MBErrorCode rval;
  mhdf_Status status;
  ioProp = H5P_DEFAULT;

  if (MB_SUCCESS != init())
    return error(MB_FAILURE);

  bool use_mpio = (MB_SUCCESS == opts.get_null_option("USE_MPIO"));
  if (use_mpio) {
#ifndef HDF5_PARALLEL
    return error(MB_NOT_IMPLEMENTED);
#else
    int parallel_mode;
    rval = opts.match_option( "PARALLEL", 
                              ReadParallel::parallelOptsNames, 
                              parallel_mode );
    if (MB_FAILURE == rval) {
      readUtil->report_error("Unexpected value for 'PARALLEL' option\n");
      return error(MB_FAILURE);
    }
    else if (MB_SUCCESS != rval ||
             parallel_mode != ReadParallel::POPT_READ_DELETE) {
      use_mpio = false;
    }
#endif
  }
  
  // This option is intended for testing purposes only, and thus
  // is not documented anywhere.  Decreasing the buffer size can
  // expose bugs that would otherwise only be seen when reading
  // very large files.
  rval = opts.get_int_option( "BUFFER_SIZE", bufferSize );
  if (MB_SUCCESS != rval) {
    bufferSize = READ_HDF5_BUFFER_SIZE;
  }
  else if (bufferSize < (int)std::max( sizeof(MBEntityHandle), sizeof(void*) )) {
    return error(MB_INVALID_SIZE);
  }
  
  dataBuffer = (char*)malloc( bufferSize );
  if (!dataBuffer)
    return error(MB_MEMORY_ALLOCATION_FAILED);

  rval = iFace->create_meshset( MESHSET_SET, file_set );
  if (MB_SUCCESS != rval) {
    free(dataBuffer);
    return error(rval);
  }
  
    // Open the file
  hid_t file_prop = H5P_DEFAULT;
#ifdef HDF5_PARALLEL
  if (use_mpio) {
    file_prop = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(file_prop, MPI_COMM_WORLD, MPI_INFO_NULL);
    ioProp = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(ioProp, H5FD_MPIO_COLLECTIVE);
  }
#endif
    
  
  filePtr = mhdf_openFileWithOpt( filename, 0, NULL, file_prop, &status );
  if (file_prop != H5P_DEFAULT)
    H5Pclose( file_prop );
  if (!filePtr)
  {
    readUtil->report_error( mhdf_message( &status ));
    free( dataBuffer );
    if (ioProp != H5P_DEFAULT)
      H5Pclose( ioProp );
    iFace->delete_entities( &file_set, 1 );
    return error(MB_FAILURE);
  }

    // get file info
  fileInfo = mhdf_getFileSummary( filePtr, handleType, &status );
  if (is_error(status)) {
    free( dataBuffer );
    if (ioProp != H5P_DEFAULT)
      H5Pclose( ioProp );
    mhdf_closeFile( filePtr, &status );
    iFace->delete_entities( &file_set, 1 );
    return error(MB_FAILURE);
  }
 
  if (name) 
    rval = load_file_partial( file_set, name, id_list, num_ids, opts );
  else
    rval = load_file_impl( file_set, opts );
  mhdf_closeFile( filePtr, &status );
  filePtr = 0;
  if (ioProp != H5P_DEFAULT)
    H5Pclose( ioProp );
  if (MB_SUCCESS == rval && is_error(status))
    rval = MB_FAILURE;

  if (file_set && MB_SUCCESS == rval) {
    DEBUGOUT("Creating entity set for file contents\n")
    MBRange range;
    MBRange::iterator hint = range.begin();
    for (IDMap::iterator j = idMap.begin(); j != idMap.end(); ++j)
      hint = range.insert( hint, j->value, j->value + j->count - 1);
    rval = iFace->add_entities( file_set, range );
  }
        
      // delete everything that was read in if read failed part-way through
  if (MB_SUCCESS != rval) {
    MBRange range;
    range.insert( file_set );
    for (IDMap::iterator i = idMap.begin(); i != idMap.end(); ++i) 
      range.insert( i->value, i->value + i->count - 1 );
    iFace->delete_entities( range );
  }
  
  free( dataBuffer );
  free( fileInfo );
  return rval;
}
  


MBErrorCode ReadHDF5::load_file_impl( MBEntityHandle file_set, const FileOptions& opts )
{
  MBErrorCode rval;
  mhdf_Status status;
  std::string tagname;
  int i;

DEBUGOUT("Reading Nodes.\n");
  
  MBRange ids;
  if (fileInfo->nodes.count) {
    ids.insert( fileInfo->nodes.start_id,
                fileInfo->nodes.start_id + fileInfo->nodes.count - 1);
    rval = read_nodes( ids );
    if (MB_SUCCESS != rval)
      return error(rval);
  }

DEBUGOUT("Reading element connectivity.\n");

  std::vector<int> polyhedra; // need to do these last so that faces are loaded
  for (i = 0; i < fileInfo->num_elem_desc; ++i) {
    if (MBCN::EntityTypeFromName(fileInfo->elems[i].type) == MBPOLYHEDRON) {
      polyhedra.push_back(i);
      continue;
    }
    
    rval = read_elems( i );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  for (std::vector<int>::iterator it = polyhedra.begin();
       it != polyhedra.end(); ++it) {
    rval = read_elems( *it );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
DEBUGOUT("Reading sets.\n");
  
  ids.clear();
  if (fileInfo->sets.count) {
    ids.insert( fileInfo->sets.start_id,
                fileInfo->sets.start_id + fileInfo->sets.count - 1);
    rval = read_sets( ids );
    if (rval != MB_SUCCESS) {
      return error(rval);
    }
  }
  
DEBUGOUT("Reading adjacencies.\n");

  for (i = 0; i < fileInfo->num_elem_desc; ++i) {
    if (!fileInfo->elems[i].have_adj)
      continue;
    
    long table_len;
    hid_t table = mhdf_openAdjacency( filePtr, 
                                      fileInfo->elems[i].handle,
                                      &table_len,
                                      &status );
    if (is_error(status))
      return error(MB_FAILURE);
      
    rval = read_adjacencies( table, table_len );
    mhdf_closeData( filePtr, table, &status );
    if (MB_SUCCESS != rval)
      return error(rval);
    if (is_error(status))
      return error(MB_FAILURE);
  }

DEBUGOUT("Reading tags.\n");
  
  for (i = 0; i < fileInfo->num_tag_desc; ++i) {
    rval = read_tag( i );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
DEBUGOUT("Finishing read.\n");
  rval = read_qa( file_set );
  return rval;
}

MBErrorCode ReadHDF5::load_file_partial( MBEntityHandle file_set, 
                                         const char* name,
                                         const int* id_list,
                                         int num_ids,
                                         const FileOptions& opts )
{
  mhdf_Status status;
  
DEBUGOUT( "RETREIVING TAGGED ENTITIES" );
    
  int tag_index;
  for (tag_index = 0; tag_index < fileInfo->num_tag_desc; ++tag_index) 
    if (!strcmp( name, fileInfo->tags[tag_index].name))
      break;
  if (tag_index == fileInfo->num_tag_desc) {
    readUtil->report_error( "File does not contain subset tag '%s'", name );
    return error(MB_TAG_NOT_FOUND);
  }
  if (fileInfo->tags[tag_index].type != mhdf_INTEGER ||
      fileInfo->tags[tag_index].size != 1) {
    readUtil->report_error( "Tag '%s' does not containa single integer value", name );
    return error(MB_TYPE_OUT_OF_RANGE);
  }
  
    
  MBErrorCode rval;
  MBRange file_ids;
  std::vector<int> ids( id_list, id_list+num_ids );
  std::sort( ids.begin(), ids.end() );
  rval = search_tag_values( tag_index, ids, file_ids );
  if (MB_SUCCESS != rval) {
    return error(rval);
  }
  
DEBUGOUT( "GATHERING ADDITIONAL ENTITIES" );
  
  const char* const set_opts[] = { "NONE", "SETS", "CONTENTS" };
  int child_mode;
  rval = opts.match_option( "CHILDREN", set_opts, child_mode );
  if (MB_ENTITY_NOT_FOUND == rval)
    child_mode = 2;
  else if (MB_SUCCESS != rval) {
    readUtil->report_error( "Invalid value for 'CHILDREN' option" );
    return error(rval);
  }
  int content_mode;
  rval = opts.match_option( "SETS", set_opts, content_mode );
  if (MB_ENTITY_NOT_FOUND == rval)
    content_mode = 2;
  else if (MB_SUCCESS != rval) {
    readUtil->report_error( "Invalid value for 'SETS' option" );
    return error(rval);
  }
  
    // If we want the contents of contained/child sets, 
    // search for them now (before gathering the non-set contents
    // of the sets.)
  MBRange sets;
  intersect( fileInfo->sets, file_ids, sets );
  if (content_mode == 2 || child_mode == 2) {
    rval = read_set_ids_recursive( sets, content_mode == 2, child_mode == 2 );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
    // get elements and vertices contained in sets
  rval = get_set_contents( sets, file_ids );
  if (MB_SUCCESS != rval)
    return error(rval);

    
DEBUGOUT( "READING NODES" );
  
    // if input contained any polyhedra, need to get faces
  for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
    MBEntityType type = MBCN::EntityTypeFromName( fileInfo->elems[i].type );
    if (type != MBPOLYHEDRON)
      continue;
    
    MBRange polyhedra;
    intersect( fileInfo->elems[i].desc, file_ids, polyhedra );
    rval = read_elems( i, polyhedra, file_ids );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
    // get node file ids for all elements
  MBRange nodes;
  intersect( fileInfo->nodes, file_ids, nodes );
  for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
    MBEntityType type = MBCN::EntityTypeFromName( fileInfo->elems[i].type );
    if (type <= MBVERTEX || type >= MBENTITYSET) {
      assert( false ); // for debug code die for unknown element tyoes
      continue; // for release code, skip unknown element types
    }
    if (MBPOLYHEDRON == type)
      continue;
    
    MBRange subset;
    intersect( fileInfo->elems[i].desc, file_ids, subset );
    rval = read_elems( i, subset, nodes );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
    // Read node coordinates and create vertices in MOAB
    // NOTE:  This populates the RangeMap with node file ids,
    //        which is expected by read_node_adj_elems.
  rval = read_nodes( nodes );
  if (MB_SUCCESS != rval)
    return error(rval);

 
DEBUGOUT( "READING ELEMENTS" );
 
    // decide if we need to read additional elements
  int side_mode;
  const char* const options[] = { "EXPLICIT", "NODES", "SIDES", 0 };
  rval = opts.match_option( "ELEMENTS", options, side_mode );
  if (MB_ENTITY_NOT_FOUND == rval) {
      // Chose default based on whether or not any elements have been
      // specified.
    MBRange tmp;
    intersect( fileInfo->nodes, file_ids, tmp );
      // If only nodes were specified, then default to "NODES", otherwise
      // default to "SIDES".
    if (!tmp.empty() && (tmp.size() + sets.size()) == file_ids.size())
      side_mode = 1;
    else
      side_mode = 2;
  }
  else if (MB_SUCCESS != rval) {
    readUtil->report_error( "Invalid value for 'ELEMENTS' option" );
    return error(rval);
  }
  
  switch (side_mode) {
    case 0: // ELEMENTS=EXPLICIT : read only specified element IDS
    for (int dim = 1; dim <= 3; ++dim) {
      for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
        MBEntityType type = MBCN::EntityTypeFromName( fileInfo->elems[i].type );
        if (MBCN::Dimension(type) == dim) {
          MBRange subset;
          intersect( fileInfo->elems[i].desc, file_ids, subset );
          rval = read_elems( fileInfo->elems[i],  subset );
          if (MB_SUCCESS != rval)
            return error(rval);
        }
      }
    }
    break;
    
    case 1: // ELEMENTS=NODES : read all elements for which all nodes have been read
    for (int dim = 1; dim <= 3; ++dim) {
      for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
        MBEntityType type = MBCN::EntityTypeFromName( fileInfo->elems[i].type );
        if (MBCN::Dimension(type) == dim) {
          rval = read_node_adj_elems( fileInfo->elems[i] );
          if (MB_SUCCESS != rval)
            return error(rval);
        }
      }
    }
    break;
    
    case 2: // ELEMENTS=SIDES : read explicitly specified elems and any sides of those elems
    rval = read_elements_and_sides( file_ids );
    if (MB_SUCCESS != rval)
      return error(rval);
    break;
  }
  
  
DEBUGOUT( "READING SETS" );
    
    // If reading contained/child sets but not their contents then find
    // them now. If we were also reading their contents we would
    // have found them already.
  if (content_mode == 1 || child_mode == 1) {
    rval = read_set_ids_recursive( sets, content_mode != 0, child_mode != 0 );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
    // Append file IDs of sets containing any of the nodes or elements
    // we've read up to this point.
  rval = find_sets_containing( sets );
  if (MB_SUCCESS != rval)
    return error(rval);
    // Now actually read all set data and instantiate sets in MOAB.
    // Get any contained sets out of file_ids.
  MBEntityHandle first_set = fileInfo->sets.start_id;
  sets.merge( file_ids.lower_bound( first_set ),
              file_ids.lower_bound( first_set + fileInfo->sets.count ) );
  rval = read_sets( sets );
  if (MB_SUCCESS != rval)
    return error(rval);

  
DEBUGOUT( "READING ADJACENCIES" );
    
  for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
    if (fileInfo->elems[i].have_adj &&
        idMap.intersects( fileInfo->elems[i].desc.start_id, fileInfo->elems[i].desc.count )) {
      long len;
      hid_t th = mhdf_openAdjacency( filePtr, fileInfo->elems[i].handle, &len, &status );
      if (is_error(status))
        return error(MB_FAILURE);

      rval = read_adjacencies( th, len );
      mhdf_closeData( filePtr, th, &status );
      if (MB_SUCCESS != rval)
        return error(rval);
    }
  }
  
DEBUGOUT( "READING TAGS" );
  
  for (int i = 0; i < fileInfo->num_tag_desc; ++i) {
    rval = read_tag( i );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
  return MB_SUCCESS;
}

MBErrorCode ReadHDF5::search_tag_values( int tag_index,
                                         const std::vector<int>& sorted_values,
                                         MBRange& file_ids )
{
  MBErrorCode rval;
  mhdf_Status status;
  std::vector<MBEntityHandle>::iterator iter;
  const mhdf_TagDesc& tag = fileInfo->tags[tag_index];
  long size;
  long start_id;
   
    // do dense data
    
  hid_t table;
  const char* name;
  std::vector<MBEntityHandle> indices;
    // These are probably in order of dimension, so iterate
    // in reverse order to make MBRange insertions more efficient.
  std::vector<int> grp_indices( tag.dense_elem_indices, tag.dense_elem_indices+tag.num_dense_indices );
  for (std::vector<int>::reverse_iterator i = grp_indices.rbegin(); i != grp_indices.rend(); ++i)
  {
    int idx = *i;
    if (idx == -2) {
      name = mhdf_set_type_handle();
      start_id = fileInfo->sets.start_id;
    }
    else if (idx == -1) {
      name = mhdf_node_type_handle();
     start_id = fileInfo->nodes.start_id;
    }
    else {
      if (idx < 0 || idx >= fileInfo->num_elem_desc) 
        return error(MB_FAILURE);
      name = fileInfo->elems[idx].handle;
      start_id = fileInfo->elems[idx].desc.start_id;
    }
    table = mhdf_openDenseTagData( filePtr, tag.name, name, &size, &status );
    if (is_error(status))
      return error(MB_FAILURE);
    rval = search_tag_values( table, size, sorted_values, indices );
    mhdf_closeData( filePtr, table, &status );
    if (MB_SUCCESS != rval || is_error(status))
      return error(MB_FAILURE);
      // Convert from table indices to file IDs and add to result list
    std::sort( indices.begin(), indices.end(), std::greater<MBEntityHandle>() );
    std::transform( indices.begin(), indices.end(), mb_range_inserter(file_ids),
                    std::bind1st( std::plus<long>(), start_id ) );
    indices.clear();
  }
  
  if (!tag.have_sparse)
    return MB_SUCCESS;
  
    // do sparse data
    
  hid_t tables[2]; 
  long junk; // redundant value for non-variable-length tags
  mhdf_openSparseTagData( filePtr, tag.name, &size, &junk, tables, &status );
  if (is_error(status))
    return error(MB_FAILURE);
  rval = search_tag_values( tables[1], size, sorted_values, indices );
  mhdf_closeData( filePtr, tables[1], &status );
  if (MB_SUCCESS != rval || is_error(status)) {
    mhdf_closeData( filePtr, tables[0], &status );
    return error(MB_FAILURE);
  }
    // convert to ranges
  std::sort( indices.begin(), indices.end() );
  std::vector<MBEntityHandle> ranges;
  iter = indices.begin();
  while (iter != indices.end()) {
    ranges.push_back( *iter );
    MBEntityHandle last = *iter;
    for (++iter; iter != indices.end() && (last + 1) == *iter; ++iter, ++last);
    ranges.push_back( last );
  }
    // read file ids
  iter = ranges.begin();
  unsigned long offset = 0;
  while (iter != ranges.end()) {
    long begin = *iter; ++iter;
    long end   = *iter; ++iter;
    mhdf_readSparseTagEntities( tables[0], begin, end - begin + 1, 
                                handleType, &indices[offset], &status );
    if (is_error(status)) {
      mhdf_closeData( filePtr, tables[0], &status );
      return error(MB_FAILURE);
    }
    offset += end - begin + 1;
  }
  mhdf_closeData( filePtr, tables[0], &status );
  if (is_error(status))
    return error(MB_FAILURE);
  assert( offset == indices.size() );
  std::sort( indices.begin(), indices.end() );
  copy_sorted_file_ids( &indices[0], indices.size(), file_ids );
  
  return MB_SUCCESS;  
}

MBErrorCode ReadHDF5::search_tag_values( hid_t tag_table, 
                                         unsigned long table_size,
                                         const std::vector<int>& sorted_values,
                                         std::vector<MBEntityHandle>& value_indices )
{
  mhdf_Status status;
  size_t chunk_size = bufferSize / sizeof(unsigned);
  unsigned * buffer = reinterpret_cast<unsigned*>(dataBuffer);
  size_t remaining = table_size, offset = 0;
  while (remaining) {
      // Get a block of tag values
    size_t count = std::min( chunk_size, remaining );
    mhdf_readDenseTag( tag_table, offset, count, H5T_NATIVE_UINT, buffer, &status );
    if (is_error(status))
      return error(MB_FAILURE);
    
      // search tag values
    for (size_t i = 0; i < count; ++i)
      if (std::binary_search( sorted_values.begin(), sorted_values.end(), buffer[i] ))
        value_indices.push_back( i + offset );
    
    offset += count;
    remaining -= count;
  }

  return MB_SUCCESS;
}

MBErrorCode ReadHDF5::read_nodes( const MBRange& node_file_ids )
{
  MBErrorCode rval;
  mhdf_Status status;
  const int dim = fileInfo->nodes.vals_per_ent;
  MBRange range;
  
  if (node_file_ids.empty())
    return MB_SUCCESS;
  
  int cdim;
  rval = iFace->get_dimension( cdim );
  if (MB_SUCCESS != rval)
    return error(rval);
  
  if (cdim < dim)
  {
    rval = iFace->set_dimension( dim );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
  hid_t data_id = mhdf_openNodeCoordsSimple( filePtr, &status );
  if (is_error(status))
    return error(MB_FAILURE);

  
  MBEntityHandle handle;
  std::vector<double*> arrays(dim);
  rval = readUtil->get_node_arrays( dim, (int)node_file_ids.size(), 0, handle, arrays );
  if (MB_SUCCESS != rval)
  {
    mhdf_closeData( filePtr, data_id, &status );
    return error(rval);
  }
  
    // read blocks of coordinates
  MBRange::const_pair_iterator p;
  for (p = node_file_ids.const_pair_begin(); p != node_file_ids.const_pair_end(); ++p)
  {
    long count = p->second - p->first + 1;
    long offset = p->first - fileInfo->nodes.start_id;
    for (int i = 0; i < dim; ++i) {
      mhdf_readNodeCoord( data_id, offset, count, i, arrays[i], &status );
      if (is_error(status)) {
        mhdf_closeData( filePtr, data_id, &status );
        return error(MB_FAILURE);
      }
      arrays[i] += count;
    }
    for (int i = dim; i < cdim; ++i) {
      memset( arrays[i], 0, count*sizeof(double) );
      arrays[i] += count;
    }
    IDMap::iterator ins = idMap.insert( p->first, handle, count );
    if (ins == idMap.end()) {
      mhdf_closeData( filePtr, data_id, &status );
      return error(MB_FAILURE);
    }
    handle += count;
  }
  
  mhdf_closeData( filePtr, data_id, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message(&status) );
    return error(MB_FAILURE);
  }
  
  return MB_SUCCESS;
}

MBErrorCode ReadHDF5::read_elems( int i )
{
  MBRange ids;
  ids.insert( fileInfo->elems[i].desc.start_id,
              fileInfo->elems[i].desc.start_id + fileInfo->elems[i].desc.count - 1);
  return read_elems( i, ids );
}

MBErrorCode ReadHDF5::read_elems( int i, const MBRange& file_ids )
{
  if (fileInfo->elems[i].desc.vals_per_ent < 0)
    return read_poly( fileInfo->elems[i], file_ids );
  else
    return read_elems( fileInfo->elems[i], file_ids );
}

MBErrorCode ReadHDF5::read_elems( const mhdf_ElemDesc& elems, const MBRange& file_ids )
{
  MBErrorCode rval = MB_SUCCESS;
  mhdf_Status status;
  
  MBEntityType type = MBCN::EntityTypeFromName( elems.type );
  if (type == MBMAXTYPE)
  {
    readUtil->report_error( "Unknown element type: \"%s\".\n", elems.type );
    return error(MB_FAILURE);
  }
  
  const int nodes_per_elem = elems.desc.vals_per_ent;
  const long first_id = elems.desc.start_id;
  hid_t data_id = mhdf_openConnectivitySimple( filePtr, elems.handle, &status );
  if (is_error(status))
    return error(MB_FAILURE);

  MBRange::const_pair_iterator p;
  for (p = file_ids.const_pair_begin(); p != file_ids.const_pair_end(); ++p) {
    const long count = p->second - p->first + 1;
    
    MBEntityHandle handle;
    MBEntityHandle* array;
    rval = readUtil->get_element_array( (int)count,
                                         nodes_per_elem,
                                         type,
                                         0,
                                         handle, 
                                         array );
    if (MB_SUCCESS != rval) 
      break;
    IDMap::iterator ins = idMap.insert( p->first, handle, count );
    if (ins == idMap.end()) 
      { rval = MB_FAILURE; break; }

    mhdf_readConnectivityWithOpt( data_id, p->first - first_id, count, handleType, array, ioProp, &status );
    if (is_error(status)) 
      break;

    rval = convert_id_to_handle( array, (size_t)(nodes_per_elem*count) );
    if (MB_SUCCESS != rval) 
      break;

      // notify MOAB of the new elements
    rval = readUtil->update_adjacencies(handle, count, nodes_per_elem, array);
    if (MB_SUCCESS != rval) 
      break;
  }


  mhdf_closeData( filePtr, data_id, &status );
  if (is_error(status) && MB_SUCCESS == rval)
    rval = error(MB_FAILURE);
  return rval;
}

MBErrorCode ReadHDF5::read_node_adj_elems( const mhdf_ElemDesc& group )
{
  mhdf_Status status;
  MBErrorCode rval;
  
  hid_t table = mhdf_openConnectivitySimple( filePtr, group.handle, &status );
  if (is_error(status))
    return error(MB_FAILURE);
    
  rval = read_node_adj_elems( group, table );
  
  mhdf_closeData( filePtr, table, &status );
  if (MB_SUCCESS == rval && is_error(status))
    return error(rval = MB_FAILURE);
  return rval;
}

MBErrorCode ReadHDF5::read_node_adj_elems( const mhdf_ElemDesc& group, 
                                           hid_t table_handle )
{
  mhdf_Status status;
  MBErrorCode rval;

    // copy data to local variables (makes other code clearer)
  const int node_per_elem = group.desc.vals_per_ent;
  long start_id = group.desc.start_id;
  long remaining = group.desc.count;
  const MBEntityType type = MBCN::EntityTypeFromName( group.type );
  
    // figure out how many elements we can read in each pass
  long* const buffer = reinterpret_cast<long*>( dataBuffer );
  const long buffer_size = bufferSize / (node_per_elem * sizeof(buffer[0]));
    // read all element connectivity in buffer_size blocks
  long offset = 0;
  while (remaining) {
      // read a block of connectivity data
    const long count = std::min( remaining, buffer_size );
    mhdf_readConnectivity( table_handle, offset, count, H5T_NATIVE_LONG, buffer, &status );
    if (is_error(status))
      return error(MB_FAILURE);
    offset += count;
    remaining -= count;
    
      // count the number of elements in the block that we want,
      // zero connectivity for other elements
    long num_elem = 0;
    long* iter = buffer;
    for (long i = 0; i < count; ++i) {
      for (int j = 0; j < node_per_elem; ++j) {
        iter[j] = (long)idMap.find( iter[j] );
        if (!iter[j]) {
          iter[0] = 0;
          break;
        }
      }
      if (iter[0])
        ++num_elem;
      iter += node_per_elem;
    }
    
    if (!num_elem) {
      start_id += count;
      continue;
    }
    
      // create elements
    MBEntityHandle handle;
    MBEntityHandle* array;
    rval = readUtil->get_element_array( (int)num_elem,
                                         node_per_elem,
                                         type,
                                         0,
                                         handle, 
                                         array );
    if (MB_SUCCESS != rval)
      return error(rval);
    
      // copy all non-zero connectivity values
    iter = buffer;
    MBEntityHandle* iter2 = array;
    MBEntityHandle h = handle;
    for (long i = 0; i < count; ++i) {
      if (!*iter) {
        iter += node_per_elem;
        continue;
      }
      IDMap::iterator ins = idMap.insert( start_id + i, h++, 1 );
      if (ins == idMap.end()) 
        return error(MB_FAILURE);
        
      long* const end = iter + node_per_elem;
      for (; iter != end; ++iter, ++iter2)
        *iter2 = (MBEntityHandle)*iter;
    }
    assert( iter2 - array == num_elem * node_per_elem );
    start_id += count;
  }
  
  return MB_SUCCESS;
}
  

MBErrorCode ReadHDF5::read_elems( int i, const MBRange& elems_in, MBRange& nodes )
{
  MBEntityHandle* const buffer = reinterpret_cast<MBEntityHandle*>(dataBuffer);
  const size_t buffer_size = bufferSize / sizeof(MBEntityHandle);
  int node_per_elem = fileInfo->elems[i].desc.vals_per_ent;
  
  if (elems_in.empty())
    return MB_SUCCESS;
    
  assert( (long)elems_in.front() >= fileInfo->elems[i].desc.start_id );
  assert( (long)elems_in.back() - fileInfo->elems[i].desc.start_id < fileInfo->elems[i].desc.count );
  
    // we don't support version 3 style poly element data
  if (fileInfo->elems[i].desc.vals_per_ent <= 0)
    return error(MB_TYPE_OUT_OF_RANGE);
  
  mhdf_Status status;
  hid_t table = mhdf_openConnectivitySimple( filePtr, fileInfo->elems[i].handle, &status );
  if (is_error(status))
    return error(MB_FAILURE);
    
  MBRange elements( elems_in );
  while (!elements.empty()) {
    MBEntityHandle file_id = elements.front();
    long count = elements.const_pair_begin()->second - file_id + 1;
    long offset = file_id - fileInfo->elems[i].desc.start_id;
    if (count*node_per_elem > (long)buffer_size)
      count = buffer_size/node_per_elem;
    elements.erase( elements.begin(), elements.begin()+count );
    
      // read element connectivity
    mhdf_readConnectivity( table, offset, count, handleType, buffer, &status );
    if (is_error(status)) {
      mhdf_closeData( filePtr, table, &status );
      return error(MB_FAILURE);
    }
    
    count *= node_per_elem;
    std::sort( buffer, buffer + count );
    count = std::unique( buffer, buffer + count ) - buffer;
    copy_sorted_file_ids( buffer, count, nodes );
  }
  
  mhdf_closeData( filePtr, table, &status );
  return is_error(status) ? error(MB_FAILURE) : MB_SUCCESS;
}

MBErrorCode ReadHDF5::read_poly( const mhdf_ElemDesc& elems, const MBRange& file_ids )
{
  class PolyReader : public ContentReader {
    private:
      const MBEntityType type;
      hid_t indexHandle, connHandle, handleType;
      MBInterface *const mb;
      IDMap& idMap;
    public:
    PolyReader( MBEntityType elem_type, hid_t idx, hid_t conn, 
                hid_t handle_type, MBInterface* iface, IDMap& id_map )
               : type(elem_type), indexHandle(idx), connHandle(conn), 
                 handleType(handle_type), mb(iface), idMap(id_map) 
               {}
    void read_indices( long offset, long count, long* buffer, mhdf_Status& status ) 
      { mhdf_readPolyConnIndices( indexHandle, offset, count, H5T_NATIVE_LONG, buffer, &status ); }
    void read_contents( long offset, long count, MBEntityHandle* buffer, mhdf_Status& status )
      { mhdf_readPolyConnIDs( connHandle, offset, count, handleType, buffer, &status ); }
    MBErrorCode store_data( MBEntityHandle, long file_id, MBEntityHandle* conn, long len, bool )
    {
      size_t valid;
      convert_id_to_handle( conn, len, valid, idMap );
      if (valid != (size_t)len)
        return error(MB_ENTITY_NOT_FOUND);
      MBEntityHandle handle;
      MBErrorCode rval = mb->create_element( type, conn, len, handle );
      if (MB_SUCCESS != rval)
        return error(rval);
      IDMap::iterator ins = idMap.insert( file_id, handle, 1 );
      if (ins == idMap.end()) 
        return error(MB_FAILURE);
      return MB_SUCCESS;
    }
  };
  
  MBEntityType type = MBCN::EntityTypeFromName( elems.type );
  if (type == MBMAXTYPE)
  {
    readUtil->report_error( "Unknown element type: \"%s\".\n", elems.type );
    return error(MB_FAILURE);
  }
  
  hid_t handles[2];
  mhdf_Status status;
  long num_poly, num_conn, first_id;
  mhdf_openPolyConnectivity( filePtr, elems.handle, &num_poly, &num_conn, &first_id, 
                             handles, &status );
  if (is_error(status))
    return error(MB_FAILURE);

  PolyReader tool( type, handles[0], handles[1], handleType, iFace, idMap );
  MBRange empty;
  MBErrorCode rval = read_contents( tool, file_ids, first_id, 0, num_poly, num_conn, empty );
  mhdf_closeData( filePtr, handles[0], &status );
  if (MB_SUCCESS == rval && is_error(status))
    rval = error(MB_FAILURE);
  mhdf_closeData( filePtr, handles[1], &status );
  if (MB_SUCCESS == rval && is_error(status))
    rval = error(MB_FAILURE);

  return rval;
}
/*
MBErrorCode ReadHDF5::read_poly( const char* elem_group )
{
  MBErrorCode rval;
  mhdf_Status status;
  char name[64];
   
  mhdf_getElemTypeName( filePtr, elem_group, name, sizeof(name), &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return error(MB_FAILURE);
  }
  MBEntityType type = MBCN::EntityTypeFromName( name );

  long count, first_id, data_len;
  hid_t handles[2];
  mhdf_openPolyConnectivity( filePtr, elem_group, &count, &data_len,
                             &first_id, handles, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return error(MB_FAILURE);
  }

  ElemSet empty_set;
  empty_set.type = MBCN::EntityTypeFromName( name );
  empty_set.type2 = elem_group;
  
  MBEntityHandle h;
  bool first = true;
  long connend = -1;
  std::vector<MBEntityHandle> connectivity; 
  for (long i = 0; i < count; ++i) {
    long prevend = connend;
    mhdf_readPolyConnIndicesWithOpt( handles[0], i, 1, H5T_NATIVE_LONG, &connend, ioProp, &status );
    if (mhdf_isError( &status ))
    {
      readUtil->report_error( mhdf_message( &status ) );
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      return error(MB_FAILURE);
    }
    
    connectivity.resize( connend - prevend );
    mhdf_readPolyConnIDsWithOpt( handles[1], prevend+1, connectivity.size(), handleType,
                          &connectivity[0], ioProp, &status );
    if (mhdf_isError( &status ))
    {
      readUtil->report_error( mhdf_message( &status ) );
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      return error(MB_FAILURE);
    }
    
    rval= convert_id_to_handle( &connectivity[0], connectivity.size() );
    if (MB_SUCCESS != rval)     {
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      return error(rval);
    }
    
    rval = iFace->create_element( type, &connectivity[0], connectivity.size(), h );
    if (MB_SUCCESS != rval) 
    {
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      return error(rval);
    }
    
    if (first || elemList.back().range.back() + 1 >= h) {
      elemList.push_back( empty_set );
      elemList.back().first_id = first_id + i;
      first = false;
    }
    elemList.back().range.insert( h );
  }
 
  MBErrorCode result = MB_SUCCESS;
  mhdf_closeData( filePtr, handles[0], &status );
  if (mhdf_isError( &status )) {
    readUtil->report_error( mhdf_message( &status ));
    result = error(MB_FAILURE);
  }
  mhdf_closeData( filePtr, handles[1], &status );
  if (mhdf_isError( &status )) {
    readUtil->report_error( mhdf_message( &status ));
    result = error(MB_FAILURE);
  }
  return result;
}
*/

MBErrorCode ReadHDF5::read_elements_and_sides( const MBRange& file_ids )
{
  MBErrorCode rval;
  MBEntityType type;
    
    // determine largest element dimension
  int max_dim = 0;
  for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
    type = MBCN::EntityTypeFromName( fileInfo->elems[i].type );
    int dim = MBCN::Dimension(type);
    if (dim > max_dim) {
      MBEntityHandle start = (MBEntityHandle)fileInfo->elems[i].desc.start_id;
      MBRange::iterator it = file_ids.lower_bound( start );
      if (it != file_ids.end() && (long)(*it - start ) < fileInfo->elems[i].desc.count)
        max_dim = dim;
    }
  }
  
    // Get MBRange of element IDs only
  MBRange elem_ids( file_ids );
  MBRange::iterator s, e;
  if (fileInfo->nodes.count) {
    MBEntityHandle first = (MBEntityHandle)fileInfo->nodes.start_id;
    s = elem_ids.lower_bound( first );
    e = MBRange::lower_bound( s, elem_ids.end(), first + fileInfo->nodes.count );
    elem_ids.erase( s, e );
  }
  if (fileInfo->sets.count) {
    MBEntityHandle first = (MBEntityHandle)fileInfo->sets.start_id;
    s = elem_ids.lower_bound( first );
    e = MBRange::lower_bound( s, elem_ids.end(), first + fileInfo->sets.count );
    elem_ids.erase( s, e );
  }
  
    // read all node-adjacent elements of smaller dimensions
  for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
    MBEntityType type = MBCN::EntityTypeFromName( fileInfo->elems[i].type );
    if (MBCN::Dimension(type) < max_dim) {
      rval = read_node_adj_elems( fileInfo->elems[i] );
      if (MB_SUCCESS != rval)
        return error(rval);
    }
  }
  
    // Read the subset of the explicitly sepecified elements that are of the
    // largest dimension.   We could read all explicitly specified elements,
    // but then the lower-dimension elements will be read a second when reading
    // node-adjacent elements.
  for (int i = 0; i < fileInfo->num_elem_desc; ++i) {
    type = MBCN::EntityTypeFromName( fileInfo->elems[i].type );
    if (MBCN::Dimension(type) == max_dim) {
      MBRange subset;
      intersect( fileInfo->elems[i].desc, elem_ids, subset );
      if (!subset.empty()) {
        elem_ids.subtract( subset );
        rval = read_elems( fileInfo->elems[i],  subset );
        if (MB_SUCCESS != rval)
          return error(rval); 
      } 
    }
  } 
      
    // delete anything we read in that we should not have
    // (e.g. an edge spanning two disjoint blocks of elements)
    
    // get MBRange of all explicitly specified elements, grouped by dimension
  MBRange explicit_elems[4];
  MBRange::iterator hints[4] = { explicit_elems[0].begin(),
                                explicit_elems[1].begin(),
                                explicit_elems[2].begin(),
                                explicit_elems[3].begin() };
  RangeMap<long, MBEntityHandle>::iterator rit;
  while (!elem_ids.empty()) {
    long start = elem_ids.front();
    long count = elem_ids.const_pair_begin()->second - start + 1;
    rit = idMap.lower_bound( start );
    assert( rit->begin <= start && rit->begin + rit->count > start );
    long offset = start - rit->begin;
    long avail = rit->count - offset;
    if (avail > count)
      count = avail;
    MBEntityHandle start_h = rit->value + offset;
    int d = MBCN::Dimension( TYPE_FROM_HANDLE( start_h ) );
    if (MBCN::Dimension( TYPE_FROM_HANDLE( start_h + count - 1 ) ) != d) 
      count = start - LAST_HANDLE( TYPE_FROM_HANDLE(start_h) ) + 1;
    
    elem_ids.erase( elem_ids.begin(), elem_ids.begin() + count );
    hints[d] = explicit_elems[d].insert( hints[d], rit->value + offset, rit->value + offset + count - 1 );
  }
  
    // get MBRange of all read elements, and remove any that were explicity specified
    // get handles for everything we've read
  MBRange all_elems;
  MBRange::iterator hint = all_elems.begin();
  for (rit = idMap.begin(); rit != idMap.end(); ++rit)
    hint = all_elems.insert( hint, rit->value, rit->value + rit->count - 1 );
    // remove any vertex handles
  all_elems.erase( all_elems.begin(), all_elems.upper_bound( MBVERTEX ) );
    // remove any set handles and any handles for elements >= max_dim
  all_elems.erase( all_elems.lower_bound( MBCN::TypeDimensionMap[max_dim].first ), all_elems.end() );
    // remove explicit elements < max_dim
  if (!explicit_elems[1].empty() && max_dim > 1)
    all_elems = all_elems.subtract( explicit_elems[1] );
  if (!explicit_elems[2].empty() && max_dim > 2)
    all_elems = all_elems.subtract( explicit_elems[2] );
 
    // remove any elements that are adjacent to some explicity specified element.
  if (max_dim > 1 && !explicit_elems[2].empty()) {
    MBRange adj;
    rval = iFace->get_adjacencies( explicit_elems[2], 1, false, adj, MBInterface::UNION );
    if (MB_SUCCESS != rval)
      return error(rval);
    if (!adj.empty())
      all_elems = all_elems.subtract( adj );
  }
  if (max_dim == 3) {
    MBRange adj;
    rval = iFace->get_adjacencies( explicit_elems[3], 1, false, adj, MBInterface::UNION );
    if (MB_SUCCESS != rval)
      return error(rval);
    if (!adj.empty())
      all_elems = all_elems.subtract( adj );
    adj.clear();
    rval = iFace->get_adjacencies( explicit_elems[3], 2, false, adj, MBInterface::UNION );
    if (MB_SUCCESS != rval)
      return error(rval);
    if (!adj.empty())
      all_elems = all_elems.subtract( adj );
  }
  
    // now delete anything remaining in all_elems
  rval = iFace->delete_entities( all_elems );
  if (MB_SUCCESS != rval)
    return error(rval);
  
    // remove dead entities from ID map
  while (!all_elems.empty()) {
    MBEntityHandle start = all_elems.front();
    MBEntityID count = all_elems.const_pair_begin()->second - start + 1;
    for (rit = idMap.begin(); rit != idMap.end(); ++rit) 
      if (rit->value <= start && (long)(start - rit->value) < rit->count)
        break;
    if (rit == idMap.end())
      return error(MB_FAILURE);
  
    MBEntityID offset = start - rit->value;
    MBEntityID avail = rit->count - offset;
    if (avail < count)
      count = avail;
    
    all_elems.erase( all_elems.begin(), all_elems.begin() + count );
    idMap.erase( rit->begin + offset, count );
  }
  
  return MB_SUCCESS;
}

MBErrorCode ReadHDF5::read_sets( const MBRange& file_ids )
{
  mhdf_Status status;
  MBErrorCode rval;
  if (fileInfo->sets.count == 0 || file_ids.empty()) {
    assert(file_ids.empty());
    return MB_SUCCESS;
  }
  
  hid_t meta_handle = mhdf_openSetMetaSimple( filePtr, &status );
  if (is_error(status))
    return error(MB_FAILURE);


    // create sets 
  MBRange ranged_set_ids;
  MBEntityHandle start_handle;
  rval = read_sets( file_ids, meta_handle, ranged_set_ids, start_handle );
  if (MB_SUCCESS != rval) {
    mhdf_closeData( filePtr, meta_handle, &status );
    return error(rval);
  }
    
    // read contents
  if (fileInfo->have_set_contents) {
    long len = 0;
    hid_t handle = mhdf_openSetData( filePtr, &len, &status );
    if (is_error(status)) {
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(MB_FAILURE);
    }
    
    rval = read_contents( file_ids, start_handle, meta_handle, handle, len,
                          ranged_set_ids );
    mhdf_closeData( filePtr, handle, &status );
    if (MB_SUCCESS == rval && is_error(status))
      rval = MB_FAILURE;
    if (MB_SUCCESS != rval) {
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(rval);
    }
  }
  
    // read set child lists
  if (fileInfo->have_set_children) {
    long len = 0;
    hid_t handle = mhdf_openSetChildren( filePtr, &len, &status );
    if (is_error(status)) {
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(MB_FAILURE);
    }
    
    rval = read_children( file_ids, start_handle, meta_handle, handle, len );
    mhdf_closeData( filePtr, handle, &status );
    if (MB_SUCCESS == rval && is_error(status))
      rval = MB_FAILURE;
    if (MB_SUCCESS != rval) {
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(rval);
    }
  }
  
    // read set parent lists
  if (fileInfo->have_set_parents) {
    long len = 0;
    hid_t handle = mhdf_openSetParents( filePtr, &len, &status );
    if (is_error(status)) {
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(MB_FAILURE);
    }
    
    rval = read_parents( file_ids, start_handle, meta_handle, handle, len );
    mhdf_closeData( filePtr, handle, &status );
    if (MB_SUCCESS == rval && is_error(status))
      rval = MB_FAILURE;
    if (MB_SUCCESS != rval) {
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(rval);
    }
  }
    
  mhdf_closeData( filePtr, meta_handle, &status );
  return is_error(status) ? error(MB_FAILURE) : MB_SUCCESS;
}

MBErrorCode ReadHDF5::read_set_ids_recursive( MBRange& sets_in_out,
                                              bool contained_sets,
                                              bool child_sets )
{
  if (!fileInfo->have_set_children)
    child_sets = false;
  if (!fileInfo->have_set_contents)
    contained_sets = false;
  if (!child_sets && !contained_sets)
    return MB_SUCCESS;

    // open data tables
  if (fileInfo->sets.count == 0) {
    assert( sets_in_out.empty() );
    return MB_SUCCESS;
  }
  
  if (!contained_sets && !child_sets)
    return MB_SUCCESS;
  
  hid_t meta_handle, content_handle = 0, child_handle = 0;
  
  mhdf_Status status;
  meta_handle = mhdf_openSetMetaSimple( filePtr, &status );
  if (is_error(status))
    return error(MB_FAILURE);
  
  if (contained_sets) {
    long content_len = 0;
    content_handle = mhdf_openSetData( filePtr, &content_len, &status );
    if (is_error(status)) {
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(MB_FAILURE);
    }
  }
  
  if (child_sets) {
    long child_len = 0;
    child_handle = mhdf_openSetChildren( filePtr, &child_len, &status );
    if (is_error(status)) {
      if (contained_sets)
        mhdf_closeData( filePtr, content_handle, &status );
      mhdf_closeData( filePtr, meta_handle, &status );
      return error(MB_FAILURE);
    }
  }
  
  MBErrorCode rval = MB_SUCCESS;
  MBRange children, new_children(sets_in_out);
  do {
    children.clear();
    if (child_sets) {
      rval = read_child_ids( new_children, meta_handle, child_handle, children );
      if (MB_SUCCESS != rval)
        break;
    }
    if (contained_sets) {
      rval = read_contained_set_ids( new_children, meta_handle, content_handle, children );
      if (MB_SUCCESS != rval)
        break;
    }
    new_children = children.subtract( sets_in_out );
    sets_in_out.merge( new_children );
  } while (!new_children.empty());
  
  if (child_sets) {
    mhdf_closeData( filePtr, child_handle, &status );
    if (MB_SUCCESS == rval && is_error(status))
      rval = error(MB_FAILURE);
  }
  if (contained_sets) {
    mhdf_closeData( filePtr, content_handle, &status );
    if (MB_SUCCESS == rval && is_error(status))
      rval = error(MB_FAILURE);
  }
  mhdf_closeData( filePtr, meta_handle, &status );
  if (MB_SUCCESS == rval && is_error(status))
    rval = error(MB_FAILURE);
  
  return rval;
}

MBErrorCode ReadHDF5::find_sets_containing( MBRange& sets_out )
{
  MBErrorCode rval;
  mhdf_Status status;

  if (!fileInfo->have_set_contents)
    return MB_SUCCESS;
  assert( fileInfo->sets.count );

    // open data tables
  hid_t meta_handle = mhdf_openSetMetaSimple( filePtr, &status );
  if (is_error(status))
    return error(MB_FAILURE);
  long content_len = 0;
  hid_t content_handle = mhdf_openSetData( filePtr, &content_len, &status );
  if (is_error(status)) {
    mhdf_closeData( filePtr, meta_handle, &status );
    return error(MB_FAILURE);
  }

  rval = find_sets_containing( meta_handle, content_handle, content_len, sets_out );

  mhdf_closeData( filePtr, content_handle, &status );
  if(MB_SUCCESS == rval && is_error(status))
    rval = error(MB_FAILURE);
  mhdf_closeData( filePtr, meta_handle, &status );
  if(MB_SUCCESS == rval && is_error(status))
    rval = error(MB_FAILURE);
    
  return rval;
}

static bool set_map_intersect( unsigned short flags,
                               const long* contents,
                               int content_len,
                               const RangeMap<long,MBEntityHandle>& id_map  )
{
  if (flags & mhdf_SET_RANGE_BIT) {
    if (!content_len || id_map.empty())
      return false;
      
    const long* j = contents;
    const long* const end = contents + content_len;
    assert(content_len % 2 == 0);
    while (j != end) {
      long start = *(j++);
      long count = *(j++);
      if (id_map.intersects( start, count ))
        return true;
    }
  }
  else {
    const long* const end = contents + content_len;
    for (const long* i = contents; i != end; ++i)
      if (id_map.exists( *i ))
        return true;
  }
  return false;
}

MBErrorCode ReadHDF5::find_sets_containing( hid_t meta_handle,
                                            hid_t contents_handle, 
                                            long contents_len,
                                            MBRange& file_ids )
{
  const long avg_set_len = contents_len / fileInfo->sets.count;
  long sets_per_buffer = bufferSize / (sizeof(short) + sizeof(long) * (2+avg_set_len));
    // round to down multiple of 8 to avoid alignment issues
  sets_per_buffer = 8 * (sets_per_buffer / 8);
  if (sets_per_buffer < 10) // just in case there's one huge set
    sets_per_buffer = 10;  
  unsigned short* flag_buffer = (unsigned short*)dataBuffer;
  long* offset_buffer = (long*)(flag_buffer + sets_per_buffer);
  long* content_buffer = offset_buffer + sets_per_buffer;
  assert(bufferSize % sizeof(long) == 0);
  long content_len = (long*)(dataBuffer + bufferSize) - content_buffer;
  assert(dataBuffer + bufferSize >= (char*)(content_buffer + content_len));
    // scan set table  
  mhdf_Status status;
  MBRange::iterator hint = file_ids.begin();
  long remaining = fileInfo->sets.count;
  long offset = 0;
  long prev_idx = -1;
  while (remaining) {
    long count = std::min( remaining, sets_per_buffer );
    mhdf_readSetFlags( meta_handle, offset, count, H5T_NATIVE_USHORT, flag_buffer, &status );
    if (is_error(status)) 
      return error(MB_FAILURE);
    mhdf_readSetContentEndIndices( meta_handle, offset, count, H5T_NATIVE_LONG, offset_buffer, &status );
    if (is_error(status))
      return error(MB_FAILURE);
    
    long sets_remaining = count;
    long sets_offset = 0;
    while (sets_remaining) {
        // figure how many of the remaining sets are required to 
        // fill the set contents buffer.
      long sets_count = std::lower_bound( offset_buffer + sets_offset, 
                          offset_buffer + count, content_len + prev_idx )
                          - offset_buffer - sets_offset;
      if (!sets_count) { // contents of single set don't fit in buffer
        long content_remaining = offset_buffer[sets_offset] - prev_idx;
        long content_offset = prev_idx+1;
        while (content_remaining) {
          long content_count = content_len < content_remaining ?
                               2*(content_len/2) : content_remaining;
          mhdf_readSetData( contents_handle, content_offset,
                            content_count, H5T_NATIVE_LONG, 
                            content_buffer, &status );
          if (is_error(status))
            return error(MB_FAILURE);
          if (set_map_intersect( flag_buffer[sets_offset],
                                 content_buffer, content_count, idMap )) {
            long id = fileInfo->sets.start_id + offset + sets_offset;
            hint = file_ids.insert( hint, id, id );
            break;
          }
          content_remaining -= content_count;
          content_offset += content_count;
        }
        prev_idx = offset_buffer[sets_offset];
        sets_count = 1;
      }
      else if (long read_num = offset_buffer[sets_offset + sets_count - 1] - prev_idx) {
        assert(sets_count > 0);
        mhdf_readSetData( contents_handle, prev_idx+1, read_num, 
                          H5T_NATIVE_LONG, content_buffer, &status );
        if (is_error(status))
          return error(MB_FAILURE);
        
        long* buff_iter = content_buffer;
        for (long i = 0; i < sets_count; ++i) {
          long set_size = offset_buffer[i+sets_offset] - prev_idx;
          prev_idx += set_size;
          if (set_map_intersect( flag_buffer[sets_offset+i],
                                 buff_iter, set_size, idMap )) {
            long id = fileInfo->sets.start_id + offset + sets_offset + i;
            hint = file_ids.insert( hint, id, id );
          }
          buff_iter += set_size;
        }
      }
    
      sets_offset += sets_count;
      sets_remaining -= sets_count;
    }
    
    offset += count;
    remaining -= count;
  }
  
  return MB_SUCCESS;
}

MBErrorCode ReadHDF5::read_child_ids( const MBRange& input_file_ids,
                                      hid_t meta_handle,
                                      hid_t child_handle,
                                      MBRange& child_file_ids )
{
  mhdf_Status status;
  long* buffer = reinterpret_cast<long*>(dataBuffer);
  long buffer_size = bufferSize / sizeof(long);
  long first, range[2], count, remaining;
  MBRange sets(input_file_ids);
  MBRange::iterator hint;
  while (!sets.empty()) {
    count = (long)sets.const_pair_begin()->second - sets.front() + 1;
    first = (long)sets.front() - fileInfo->sets.start_id;
    sets.erase( sets.begin(), sets.begin() + count );
    
    if (!first) {
      range[0] = -1;
      mhdf_readSetChildEndIndices( meta_handle, first+count-1, 1, H5T_NATIVE_LONG, range+1, &status );
      if (is_error(status))
        return error(MB_FAILURE);
    }
    else if (count == 1) {
      mhdf_readSetChildEndIndices( meta_handle, first-1, 2, H5T_NATIVE_LONG, range, &status );
      if (is_error(status))
        return error(MB_FAILURE);
    }
    else {
      mhdf_readSetChildEndIndices( meta_handle, first-1, 1, H5T_NATIVE_LONG, range+1, &status );
      if (is_error(status))
        return error(MB_FAILURE);
      mhdf_readSetChildEndIndices( meta_handle, first+count-1, 1, H5T_NATIVE_LONG, range+1, &status );
      if (is_error(status))
        return error(MB_FAILURE);
    }
    
    remaining = range[1] - range[0];
    long offset = range[0] + 1;
    while (remaining) {
      count = std::min( buffer_size, remaining );
      remaining -= count;
      mhdf_readSetParentsChildren( child_handle, offset, count, H5T_NATIVE_LONG, buffer, &status );
  
      std::sort( buffer, buffer + count );
      count = std::unique( buffer, buffer + count ) - buffer;
      hint = child_file_ids.begin();
      for (long i = 0; i < count; ++i) {
        MBEntityHandle h = (MBEntityHandle)buffer[i];
        hint = child_file_ids.insert( hint, h, h );
      }
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode ReadHDF5::read_contained_set_ids( const MBRange& input_file_ids,
                                              hid_t meta_handle,
                                              hid_t content_handle,
                                              MBRange& contained_set_file_ids )
{
  mhdf_Status status;
  long buffer_size = bufferSize / (sizeof(long) + sizeof(short));
    // don't want to worry about reading half of a range pair later
  if (buffer_size % 2) --buffer_size;
  long* content_buffer = reinterpret_cast<long*>(dataBuffer);
  unsigned short* flag_buffer = reinterpret_cast<unsigned short*>(content_buffer + buffer_size);
  long first, range[2], count, remaining, sets_offset;

  MBRange sets(input_file_ids);
  MBRange::iterator hint;
  while (!sets.empty()) {
    count = (long)sets.const_pair_begin()->second - sets.front() + 1;
    first = (long)sets.front() - fileInfo->sets.start_id;
    if (count > buffer_size)
      count = buffer_size;
    sets.erase( sets.begin(), sets.begin() + count );
    
    mhdf_readSetFlags( meta_handle, first, count, H5T_NATIVE_USHORT, flag_buffer, &status );
    if (is_error(status))
      return MB_FAILURE;
    
    sets_offset = 0;
    while (sets_offset < count) {
        // Find block of sets with same value for ranged flag
      long start_idx = sets_offset;
      unsigned short ranged = flag_buffer[start_idx] & mhdf_SET_RANGE_BIT;
      for (++sets_offset; sets_offset < count; ++sets_offset)
        if ((flag_buffer[sets_offset] & mhdf_SET_RANGE_BIT) != ranged)
          break;
          
      if (!first && !start_idx) { // first set
        range[0] = -1;
        mhdf_readSetContentEndIndices( meta_handle, first+sets_offset-1, 1, 
                                       H5T_NATIVE_LONG, range+1, &status );
        if (is_error(status))
          return error(MB_FAILURE);
      }
      else if (count == 1) {
        mhdf_readSetContentEndIndices( meta_handle, first+start_idx-1, 2, 
                                       H5T_NATIVE_LONG, range, &status );
        if (is_error(status))
          return error(MB_FAILURE);
      }
      else {
        mhdf_readSetContentEndIndices( meta_handle, first+start_idx-1, 1, 
                                       H5T_NATIVE_LONG, range, &status );
        if (is_error(status))
          return error(MB_FAILURE);
        mhdf_readSetContentEndIndices( meta_handle, first+sets_offset-1, 1, 
                                       H5T_NATIVE_LONG, range+1, &status );
        if (is_error(status))
          return error(MB_FAILURE);
      }
    
      remaining = range[1] - range[0];
      long offset = range[0] + 1;
      while (remaining) {
        assert( !ranged || !(remaining % 2) );
        long content_count = std::min( buffer_size, remaining );
        remaining -= content_count;
        mhdf_readSetData( content_handle, offset, content_count, H5T_NATIVE_LONG, content_buffer, &status );
  
        if (ranged) {
          hint = contained_set_file_ids.begin();
          for (long i = 0; i < content_count; i += 2) {
            MBEntityHandle s = (MBEntityHandle)content_buffer[i];
            MBEntityHandle e = s + content_buffer[i+1];
            if ((long)s < fileInfo->sets.start_id)
              s = fileInfo->sets.start_id;
            if ((long)e > fileInfo->sets.start_id + fileInfo->sets.count)
              e = fileInfo->sets.start_id + fileInfo->sets.count;
            if (s < e) 
              hint = contained_set_file_ids.insert( hint, s, e - 1 );
          }
        }
        else {
          std::sort( content_buffer, content_buffer + content_count );
          long* s = std::lower_bound( content_buffer, content_buffer + content_count,
                                      fileInfo->sets.start_id );
          long* e = std::lower_bound( s, content_buffer + content_count, 
                                      fileInfo->sets.start_id + fileInfo->sets.count );
          e = std::unique( s, e );
          hint = contained_set_file_ids.begin();
          for ( ; s != e; ++s) {
            MBEntityHandle h = *s;
            hint = contained_set_file_ids.insert( hint, h, h );
          }
        }
      }
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode ReadHDF5::read_sets( const MBRange& file_ids,
                                 hid_t meta_handle, 
                                 MBRange& ranged_file_ids,
                                 MBEntityHandle& start_handle,
                                 bool create )
{
  MBErrorCode rval;
  mhdf_Status status;

  size_t count = file_ids.size();
  if (!count)
    return MB_SUCCESS;

  std::vector<unsigned> tmp_buffer;
  unsigned* buffer;
  if (count > (bufferSize/sizeof(unsigned))) {
    tmp_buffer.resize( count );
    buffer = &tmp_buffer[0];
  }
  else {
    buffer = reinterpret_cast<unsigned*>(dataBuffer);
  }
  
  unsigned* buff_iter = buffer;
  MBRange::const_pair_iterator p;
  for (p = file_ids.const_pair_begin(); p != file_ids.const_pair_end(); ++p) {
    long offset = p->first - fileInfo->sets.start_id;
    long count = p->second - p->first + 1;
    mhdf_readSetFlags( meta_handle, offset, count, H5T_NATIVE_UINT, buff_iter, &status );
    if (is_error(status))
      return error(MB_FAILURE);
    offset += count;
    buff_iter += count;
  }
  
  buff_iter = buffer;
  MBRange::iterator hint = ranged_file_ids.begin();
  for (MBRange::iterator i = file_ids.begin(); i != file_ids.end(); ++i, ++buff_iter) {
    if ((*buff_iter) & mhdf_SET_RANGE_BIT) {
      *buff_iter &= ~(unsigned)mhdf_SET_RANGE_BIT;
      hint = ranged_file_ids.insert( hint, *i, *i );
    }
  }
  
  if (create) {
    rval = readUtil->create_entity_sets( count, buffer, 0, start_handle );
    if (MB_SUCCESS != rval)
      return error(rval);
      
    MBEntityHandle h = start_handle;
    for (p = file_ids.const_pair_begin(); p != file_ids.const_pair_end(); ++p) {
      long count = p->second - p->first + 1;
      IDMap::iterator ins = idMap.insert( p->first, h, count );
      if (ins == idMap.end()) 
        return error(MB_FAILURE);
      h += count;
    }
  }
  
  return MB_SUCCESS;
}


MBErrorCode ReadHDF5::read_contents( ContentReader& tool,
                                     const MBRange& file_ids,
                                     const long start_id,
                                     const MBEntityHandle start_handle,
                                     const long entity_count,
                                     const long content_len,
                                     const MBRange& ranged_ids_in )
{
  MBErrorCode rval;
  mhdf_Status status;
  if (file_ids.empty())
    return MB_SUCCESS;
    
    // things will get messed up if this isn't true
  assert( ranged_ids_in.subtract( file_ids ).empty() );
  assert( file_ids.front() >= (MBEntityHandle)start_id );
  assert( file_ids.back() - start_id < (MBEntityHandle)entity_count );

  const long avg_set_len = content_len / entity_count;
  long sets_per_buffer = bufferSize / (sizeof(long) + (avg_set_len+1)*sizeof(MBEntityHandle));
    // round to down multiple of 8 to avoid alignment issues
  sets_per_buffer = 8 * (sets_per_buffer / 8);
  if (sets_per_buffer < 10) // just in case there's one huge set
    sets_per_buffer = 10;  
  long* offset_buffer = (long*)dataBuffer;
  MBEntityHandle* content_buffer = (MBEntityHandle*)(offset_buffer + sets_per_buffer);
  assert(bufferSize % sizeof(long) == 0);
  long content_size = (MBEntityHandle*)(dataBuffer + bufferSize) - content_buffer;
  assert(dataBuffer + bufferSize >= (char*)(content_buffer + content_size));
 
  MBRange ranged_ids(ranged_ids_in);
  MBEntityHandle h = start_handle;
  MBRange sets(file_ids), range;
  long prev_start = start_id;
  long prev_count = 0;
  while (!sets.empty()) {
    long start = sets.front();
    long count = sets.const_pair_begin()->second - start + 1;
    long file_id = start;
    if (count >= sets_per_buffer)
      count = sets_per_buffer - 1; // less one because we need the previous offset
    sets.erase( sets.begin(), sets.begin() + count );

    if (start == start_id) {
      offset_buffer[0] = -1;
      tool.read_indices( 0, count, offset_buffer + 1, status );
    }
    else if (prev_start + prev_count == start) {
      offset_buffer[0] = offset_buffer[prev_count];
      tool.read_indices( start - start_id, count, offset_buffer + 1, status );
    }
    else {
      tool.read_indices( start - start_id - 1, count+1, offset_buffer, status );
    }
    if (is_error(status))
      return error(MB_FAILURE);
    prev_start = start;
    prev_count = count;
     
      // read set contents
    long remaining = count; //offset_buffer[count] - offset_buffer[0];
    long offset = 0;
    while (remaining) {
        // Figure out how many set contents list we can fit in the
        // buffer.
      long read_count;
      for (read_count = 1; read_count < remaining; ++read_count) 
        if (offset_buffer[read_count+offset+1] - offset_buffer[offset] >= content_size)
          break;
        // If no set contents list fit in the buffer, special case:
        // iteratively read subsets of the content list.
      if (read_count == 1) { 
        assert( 0 == content_size % 2 ); // otherwise problems with ranged data
        long rem_contents = offset_buffer[offset+1] - offset_buffer[offset];
        long set_offset = 0;
        const bool ranged = !ranged_ids.empty() && ((long)ranged_ids.front() == file_id);
        if (ranged) 
          ranged_ids.pop_front();
        while (rem_contents) {
          long content_count = std::min( content_size, rem_contents );
          tool.read_contents( offset_buffer[offset]+1 + set_offset, content_count, content_buffer, status );
          if (is_error(status)) 
            return error(MB_FAILURE);
         
          rval = tool.store_data( h, file_id, content_buffer, content_count, ranged );
          if (MB_SUCCESS != rval)
            return error(rval);
         
          set_offset += content_count;
          rem_contents -= content_count;
        }
        ++h;
        ++file_id;
      }
        // Read contents for 'read_count' sets.
      else if (long read_num = offset_buffer[read_count+offset] - offset_buffer[offset]){
        tool.read_contents( offset_buffer[offset]+1, read_num, 
                            content_buffer, status );
        if (is_error(status))
          return error(MB_FAILURE);
        MBEntityHandle* content_iter = content_buffer;
        for (long i = 0; i < read_count; ++i) {
          long content_count = offset_buffer[i+1] - offset_buffer[i];
          bool ranged = !ranged_ids.empty() && ((long)ranged_ids.front() == file_id);
          rval = tool.store_data( h, file_id, content_iter, content_count, ranged );
          if (MB_SUCCESS != rval)
            return error(rval);
         
          content_iter += content_count;
          ++file_id;
          ++h;
          
          if (ranged)
            ranged_ids.pop_front();
        }
      }
      
      remaining -= read_count;
      offset += read_count;
    }
  }
  
  return MB_SUCCESS;
}



MBErrorCode ReadHDF5::read_contents( const MBRange& set_file_ids,
                                     MBEntityHandle start_handle,
                                     hid_t set_meta_data_table,
                                     hid_t set_contents_table,
                                     long set_contents_length,
                                     const MBRange& ranged_set_file_ids )
{

  class ReadSetContents : public ReadHDF5::ContentReader {
    const hid_t metaHandle, contentHandle, handleType;
    MBInterface *const mb;
    const IDMap& idMap;
  public:
    ReadSetContents( hid_t meta, hid_t content, hid_t handle_type, MBInterface* iface, 
                    const IDMap& id_map )
                    : metaHandle(meta), contentHandle(content), 
                      handleType(handle_type), mb(iface), idMap(id_map) 
                    {}
    void read_indices( long offset, long count, long* buffer, mhdf_Status& status )
      { mhdf_readSetContentEndIndices( metaHandle, offset, count, H5T_NATIVE_LONG, buffer, &status ); }
    void read_contents( long offset, long count, MBEntityHandle* buffer, mhdf_Status& status )
      { mhdf_readSetData( contentHandle, offset, count, handleType, buffer, &status ); }
    MBErrorCode store_data( MBEntityHandle set, long, MBEntityHandle* array, long len, bool ranged ) 
    {
      if (ranged) {
        if (len % 2) 
          return error(MB_INDEX_OUT_OF_RANGE);
        MBRange range;
        convert_range_to_handle( array, len/2, idMap, range );
        return mb->add_entities( set, range );
      }
      else {
        if (!len)
          return MB_SUCCESS;
        size_t valid;
        convert_id_to_handle( array, len, valid, idMap );
        return mb->add_entities( set, array, valid );
      }
    }
  };

  ReadSetContents tool( set_meta_data_table, set_contents_table, 
                        handleType, iFace, idMap );
  return read_contents( tool, set_file_ids, fileInfo->sets.start_id, start_handle, 
                        fileInfo->sets.count, set_contents_length, ranged_set_file_ids );
}


MBErrorCode ReadHDF5::read_children( const MBRange& set_file_ids,
                                     MBEntityHandle start_handle,
                                     hid_t set_meta_data_table,
                                     hid_t set_contents_table,
                                     long set_contents_length )
{
  class ReadSetChildren : public ReadHDF5::ContentReader {
    const hid_t metaHandle, contentHandle, handleType;
    MBInterface *const mb;
    const IDMap& idMap;
  public:
    ReadSetChildren( hid_t meta, hid_t content, hid_t handle_type, MBInterface* iface, 
                    const IDMap& id_map )
                    : metaHandle(meta), contentHandle(content), 
                      handleType(handle_type), mb(iface), idMap(id_map) 
                    {}
    void read_indices( long offset, long count, long* buffer, mhdf_Status& status )
      { mhdf_readSetChildEndIndices( metaHandle, offset, count, H5T_NATIVE_LONG, buffer, &status ); }
    void read_contents( long offset, long count, MBEntityHandle* buffer, mhdf_Status& status )
      { mhdf_readSetParentsChildren( contentHandle, offset, count, handleType, buffer, &status ); }
    MBErrorCode store_data( MBEntityHandle set, long, MBEntityHandle* array, long len, bool ranged ) 
    {
      assert(!ranged);
      size_t valid;
      convert_id_to_handle( array, len, valid, idMap );
      return mb->add_child_meshsets( set, array, valid );
    }
  };

  MBRange empty;
  ReadSetChildren tool( set_meta_data_table, set_contents_table, 
                        handleType, iFace, idMap );
  return read_contents( tool, set_file_ids, fileInfo->sets.start_id, start_handle, 
                        fileInfo->sets.count, set_contents_length, empty );
}

MBErrorCode ReadHDF5::read_parents( const MBRange& set_file_ids,
                                    MBEntityHandle start_handle,
                                    hid_t set_meta_data_table,
                                    hid_t set_contents_table,
                                    long set_contents_length )
{
  class ReadSetParents : public ReadHDF5::ContentReader {
    const hid_t metaHandle, contentHandle, handleType;
    MBInterface *const mb;
    const IDMap& idMap;
  public:
    ReadSetParents( hid_t meta, hid_t content, hid_t handle_type, MBInterface* iface, 
                    const IDMap& id_map )
                    : metaHandle(meta), contentHandle(content), 
                      handleType(handle_type), mb(iface), idMap(id_map) 
                    {}
    void read_indices( long offset, long count, long* buffer, mhdf_Status& status )
      { mhdf_readSetParentEndIndices( metaHandle, offset, count, H5T_NATIVE_LONG, buffer, &status ); }
    void read_contents( long offset, long count, MBEntityHandle* buffer, mhdf_Status& status )
      { mhdf_readSetParentsChildren( contentHandle, offset, count, handleType, buffer, &status ); }
    MBErrorCode store_data( MBEntityHandle set, long, MBEntityHandle* array, long len, bool ranged ) 
    {
      assert(!ranged);
      size_t valid;
      convert_id_to_handle( array, len, valid, idMap );
      return mb->add_parent_meshsets( set, array, valid );
    }
  };

  MBRange empty;
  ReadSetParents tool( set_meta_data_table, set_contents_table, 
                        handleType, iFace, idMap );
  return read_contents( tool, set_file_ids, fileInfo->sets.start_id, start_handle, 
                        fileInfo->sets.count, set_contents_length, empty );
}

static void copy_set_contents( int ranged,
                               const MBEntityHandle* contents,
                               long length,
                               MBRange& results )
{
  if (ranged) {
    assert( length%2 == 0 );
    MBRange::iterator hint = results.begin();
    for (long i = 0; i < length; i += 2)
      hint = results.insert( hint, contents[i], contents[i] + contents[i+1] - 1 );
  }
  else {
    for (long i = 0; i < length; ++i)
      results.insert( contents[i] );
  }
}


MBErrorCode ReadHDF5::get_set_contents( const MBRange& sets, MBRange& file_ids )
{
  class GetContentList : public ReadHDF5::ContentReader {
    const hid_t metaHandle, contentHandle, handleType;
    MBRange *const resultList;
  public:
    GetContentList( hid_t meta, hid_t content, hid_t handle_type, MBRange* result_set )
                    : metaHandle(meta), contentHandle(content), 
                      handleType(handle_type), resultList(result_set) 
                    {}

    void read_indices( long offset, long count, long* buffer, mhdf_Status& status )
      { mhdf_readSetContentEndIndices( metaHandle, offset, count, H5T_NATIVE_LONG, buffer, &status ); }
    void read_contents( long offset, long count, MBEntityHandle* buffer, mhdf_Status& status )
      { mhdf_readSetData( contentHandle, offset, count, handleType, buffer, &status ); }
    MBErrorCode store_data( MBEntityHandle set, long, MBEntityHandle* array, long len, bool ranged ) 
    {
      if (ranged) {
        if (len % 2)
          return error(MB_INDEX_OUT_OF_RANGE);
        copy_set_contents( 1, array, len, *resultList );
      }
      else {
        std::sort( array, array+len );
        copy_sorted_file_ids( array, len, *resultList );
      }
      return MB_SUCCESS;
    }
  };

  MBErrorCode rval;
  if (sets.empty())
    return MB_SUCCESS;

  mhdf_Status status;
  long content_len;
  hid_t meta, contents;
  meta = mhdf_openSetMetaSimple( filePtr, &status );
  if (is_error(status))
    return error(MB_FAILURE);
  contents = mhdf_openSetData( filePtr, &content_len, &status );
  if (is_error(status)) {
    mhdf_closeData( filePtr, meta, &status );
    return error(MB_FAILURE);
  }

  MBEntityHandle junk;
  MBRange ranged;
  rval = read_sets( sets, meta, ranged, junk, false );
  if (MB_SUCCESS != rval) {
    mhdf_closeData( filePtr, meta, &status );
    mhdf_closeData( filePtr, contents, &status );
    return error(rval);
  }
  
  GetContentList tool( meta, contents, handleType, &file_ids );
  rval = read_contents( tool, sets, fileInfo->sets.start_id, junk, 
                        fileInfo->sets.count, content_len, ranged );
  mhdf_closeData( filePtr, meta, &status );
  mhdf_closeData( filePtr, contents, &status );
  return rval;
}


MBErrorCode ReadHDF5::read_adjacencies( hid_t table, long table_len )
{
  MBErrorCode rval;
  mhdf_Status status;
  
  MBEntityHandle* buffer = (MBEntityHandle*)dataBuffer;
  size_t chunk_size = bufferSize / sizeof(MBEntityHandle);
  size_t remaining = table_len;
  size_t left_over = 0;
  size_t offset = 0;
  while (remaining)
  {
    size_t count = std::min( chunk_size, remaining );
    count -= left_over;
    remaining -= count;
    
    mhdf_readAdjacencyWithOpt( table, offset, count, handleType, buffer + left_over,
                               ioProp, &status );
    if (is_error(status))
      return error(MB_FAILURE);
    
    MBEntityHandle* iter = buffer;
    MBEntityHandle* end = buffer + count + left_over;
    while (end - iter >= 3)
    {
      MBEntityHandle h = idMap.find( *iter++ );
      MBEntityHandle count2 = *iter++;
      if (!h) {
        iter += count2;
        continue;
      }

      if (count2 < 1)
        return error(MB_FAILURE);

      if (end < count2 + iter)
      {
        iter -= 2;
        break;
      }
      
      size_t valid;
      convert_id_to_handle( iter, count2, valid, idMap );
      rval = iFace->add_adjacencies( h, iter, valid, false );
      if (MB_SUCCESS != rval)
        return error(rval);
     
      iter += count2;
    }
    
    left_over = end - iter;
    memmove( buffer, iter, left_over );
  }
  
  assert(!left_over);  // unexpected truncation of data
  
  return MB_SUCCESS;  
}


MBErrorCode ReadHDF5::read_tag( int tag_index )
{
  MBErrorCode rval;
  mhdf_Status status;
  MBTag tag = 0;
  hid_t read_type;
  rval = create_tag( fileInfo->tags[tag_index], tag, read_type ); 
  if (MB_SUCCESS != rval)
    return error(rval);

  if (fileInfo->tags[tag_index].have_sparse) {
    hid_t handles[3];
    long num_ent, num_val;
    mhdf_openSparseTagData( filePtr, 
                            fileInfo->tags[tag_index].name,
                            &num_ent, &num_val,
                            handles, &status );
    if (is_error(status)) {
      if (read_type) H5Tclose( read_type );
      return error(MB_FAILURE);
    }

    if (fileInfo->tags[tag_index].size > 0)
      rval = read_sparse_tag( tag, read_type, handles[0], handles[1], num_ent );
    else
      rval = read_var_len_tag( tag, read_type, handles[0], handles[1], handles[2], num_ent, num_val );

    mhdf_closeData( filePtr, handles[0], &status );
    if (MB_SUCCESS == rval && is_error(status))
      rval = MB_FAILURE;
    mhdf_closeData( filePtr, handles[1], &status );
    if (MB_SUCCESS == rval && is_error(status))
      rval = MB_FAILURE;
    if (fileInfo->tags[tag_index].size <= 0) {
      mhdf_closeData( filePtr, handles[2], &status );
      if (MB_SUCCESS == rval && is_error(status))
        rval = MB_FAILURE;
    }
    if (MB_SUCCESS != rval) {
      if (read_type) H5Tclose( read_type );
      return error(rval);
    }
  }
  
  for (int j = 0; j < fileInfo->tags[tag_index].num_dense_indices; ++j) {
    long count;
    const char* name = 0;
    mhdf_EntDesc* desc;
    int elem_idx = fileInfo->tags[tag_index].dense_elem_indices[j];
    if (elem_idx == -2) {
      desc = &fileInfo->sets;
      name = mhdf_set_type_handle();
    }
    else if (elem_idx == -1) {
      desc = &fileInfo->nodes;
      name = mhdf_node_type_handle();
    }
    else if (elem_idx >= 0 && elem_idx < fileInfo->num_elem_desc) {
      desc = &fileInfo->elems[elem_idx].desc;
      name = fileInfo->elems[elem_idx].handle;
    }
    else {
      return error(MB_FAILURE);
    }
    
    hid_t handle = mhdf_openDenseTagData( filePtr, 
                                          fileInfo->tags[tag_index].name,
                                          name,
                                          &count, &status );
    if (is_error(status)) {
      rval = error(MB_FAILURE);
      break;
    }
    
    if (count > desc->count) {
      readUtil->report_error( "Invalid data length for dense tag data: %s/%s\n",
                              name, fileInfo->tags[tag_index].name );
      mhdf_closeData( filePtr, handle, &status );
      rval = error(MB_FAILURE);
      break;
    }

    rval = read_dense_tag( tag, read_type, handle, desc->start_id, count );
    mhdf_closeData( filePtr, handle, &status );
    if (MB_SUCCESS != rval)
      break;
    if (is_error(status)) {
      rval = error(MB_FAILURE);
      break;
    }
  }
  
  if (read_type) 
    H5Tclose( read_type );
  return rval;
}
                              
                                            
    


MBErrorCode ReadHDF5::create_tag( const mhdf_TagDesc& info,
                                  MBTag& handle,
                                  hid_t& hdf_type )
{
  MBErrorCode rval;
  mhdf_Status status;
  MBTagType storage;
  MBDataType mb_type;
  bool re_read_default = false;

  switch (info.storage) {
    case mhdf_DENSE_TYPE : storage = MB_TAG_DENSE ; break;
    case mhdf_SPARSE_TYPE: storage = MB_TAG_SPARSE; break;
    case mhdf_BIT_TYPE   : storage = MB_TAG_BIT;    break;
    case mhdf_MESH_TYPE  : storage = MB_TAG_MESH;   break;
    default:
      readUtil->report_error( "Invalid storage type for tag '%s': %d\n", info.name, info.storage );
      return error(MB_FAILURE);
  }

    // Type-specific stuff
  if (info.type == mhdf_BITFIELD) {
    if (info.size < 1 || info.size > 8)
    {
      readUtil->report_error( "Invalid bit tag:  class is MB_TAG_BIT, num bits = %d\n", info.size );
      return error(MB_FAILURE);
    }
    hdf_type = H5Tcopy(H5T_NATIVE_B8);
    mb_type = MB_TYPE_BIT;
    if (hdf_type < 0)
      return error(MB_FAILURE);
  }
  else if (info.type == mhdf_OPAQUE) {
    mb_type = MB_TYPE_OPAQUE;

      // Check for user-provided type
    MBTag type_handle;
    std::string tag_type_name = "__hdf5_tag_type_";
    tag_type_name += info.name;
    rval = iFace->tag_get_handle( tag_type_name.c_str(), type_handle );
    if (MB_SUCCESS == rval) {
      rval = iFace->tag_get_data( type_handle, 0, 0, &hdf_type );
      if (MB_SUCCESS != rval)
        return error(rval);
      hdf_type = H5Tcopy( hdf_type );
      re_read_default = true;
    }
    else if (MB_TAG_NOT_FOUND == rval) {
      hdf_type = 0;
    }
    else
      return error(rval);
      
    if (hdf_type < 0)
      return error(MB_FAILURE);
  }
  else {
    switch (info.type)
    {
      case mhdf_INTEGER:
        hdf_type = H5T_NATIVE_INT;
        mb_type = MB_TYPE_INTEGER;
        break;

      case mhdf_FLOAT:
        hdf_type = H5T_NATIVE_DOUBLE;
        mb_type = MB_TYPE_DOUBLE;
        break;

      case mhdf_BOOLEAN:
        hdf_type = H5T_NATIVE_UINT;
        mb_type = MB_TYPE_INTEGER;
        break;

      case mhdf_ENTITY_ID:
        hdf_type = handleType;
        mb_type = MB_TYPE_HANDLE;
        break;

      default:
        return error(MB_FAILURE);
    }
    
    if (info.size > 1) { // array
        hsize_t tmpsize = info.size;
#if defined(H5Tarray_create_vers) && H5Tarray_create_vers > 1  
        hdf_type = H5Tarray_create2( hdf_type, 1, &tmpsize );
#else
        hdf_type = H5Tarray_create( hdf_type, 1, &tmpsize, NULL );
#endif
    }
    else {
      hdf_type = H5Tcopy( hdf_type );
    }
    if (hdf_type < 0)
      return error(MB_FAILURE);
  }

  
    // If default or global/mesh value in file, read it.
  if (info.default_value || info.global_value)
  {
    if (re_read_default) {
      mhdf_getTagValues( filePtr, info.name, hdf_type, info.default_value, info.global_value, &status );
      if (mhdf_isError( &status ))
      {
        readUtil->report_error( mhdf_message( &status ) );
        if (hdf_type) H5Tclose( hdf_type );
        return error(MB_FAILURE);
      }
    }
    
    if (MB_TYPE_HANDLE == mb_type) {
      if (info.default_value) {
        rval = convert_id_to_handle( (MBEntityHandle*)info.default_value, info.default_value_size );
        if (MB_SUCCESS != rval) {
          if (hdf_type) H5Tclose( hdf_type );
          return error(rval);
        }
      }
      if (info.global_value) {
        rval = convert_id_to_handle( (MBEntityHandle*)info.global_value, info.global_value_size );
        if (MB_SUCCESS != rval) {
          if (hdf_type) H5Tclose( hdf_type );
          return error(rval);
        }
      }
    }
  }
  
  
    // Check if tag already exists
  rval = iFace->tag_get_handle( info.name, handle );
  if (MB_SUCCESS == rval) {
    // If tag exists, make sure it is consistant with the type in the file
    int curr_size;
    MBDataType curr_type;
    MBTagType curr_store;
    
    rval = iFace->tag_get_size( handle, curr_size );
    if (MB_VARIABLE_DATA_LENGTH == rval)
      curr_size = -1;
    else if (MB_SUCCESS != rval) {
      if (hdf_type) H5Tclose( hdf_type );
      return error(rval);
    }
    
    rval = iFace->tag_get_data_type( handle, curr_type );
    if (MB_SUCCESS != rval) {
      if (hdf_type) H5Tclose( hdf_type );
      return error(rval);
    }
    
    rval = iFace->tag_get_type( handle, curr_store );
    if (MB_SUCCESS != rval) {
      if (hdf_type) H5Tclose( hdf_type );
      return error(rval);
    }
    
    if ((curr_store != MB_TAG_BIT && curr_size != info.bytes) || curr_type != mb_type ||
        ((curr_store == MB_TAG_BIT || storage == MB_TAG_BIT) && 
          curr_store != storage))
    {
      readUtil->report_error( "Tag type in file does not match type in "
                              "database for \"%s\"\n", info.name );
      if (hdf_type) H5Tclose( hdf_type );
      return error(MB_FAILURE);
    }
  }
    // Create the tag if it doesn't exist
  else if (MB_TAG_NOT_FOUND == rval)
  {
    if (info.size < 0) {
      size_t size = hdf_type ? H5Tget_size(hdf_type) : 1;
      rval = iFace->tag_create_variable_length( info.name, storage, mb_type,
                                                handle, info.default_value, 
                                                info.default_value_size * size );
    }
    else
      rval = iFace->tag_create( info.name, info.bytes, storage, mb_type,
                                handle, info.default_value );
    if (MB_SUCCESS != rval) {
      if (hdf_type) H5Tclose( hdf_type );
      return error(rval);
    }
  }
    // error
  else {
    if (hdf_type) H5Tclose( hdf_type );
    return error(rval);
  }
    
  if (info.global_value) {
    int type_size = hdf_type ? H5Tget_size(hdf_type) : 1;
    int tag_size = info.global_value_size * type_size;
    rval = iFace->tag_set_data( handle, 0, 0, &info.global_value, &tag_size );
    if (MB_SUCCESS != rval) {
      if (hdf_type) H5Tclose( hdf_type );
      return error(rval);
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode ReadHDF5::read_dense_tag( MBTag tag_handle,
                                      hid_t hdf_read_type,
                                      hid_t data,
                                      long start_id,
                                      long num_values )
{
  mhdf_Status status;
  MBErrorCode rval;
  MBDataType mb_type;
  
  rval = iFace->tag_get_data_type( tag_handle, mb_type );
  if (MB_SUCCESS != rval) 
    return error(rval);
  
  hsize_t read_size = 0;
  if (mb_type == MB_TYPE_HANDLE) {
    read_size = H5Tget_size( hdf_read_type );
    if (!read_size)
      return error(MB_FAILURE);
  }
  
    // get actual entities read from file
  MBRange file_ids;
  IDMap::iterator l, u;
  l = idMap.lower_bound( start_id );
  u = idMap.lower_bound( start_id + num_values - 1 );
  if (l != idMap.end() && start_id + num_values > l->begin) {
    if (l == u) {
      file_ids.insert( std::max(start_id, l->begin), 
                       std::min(start_id + num_values, u->begin + u->count) - 1 );
    }
    else {
      file_ids.insert( std::max(start_id, l->begin), l->begin + l->count - 1 );
      if (u != idMap.end() && u->begin < start_id + num_values)
        file_ids.insert( u->begin, std::min( start_id + num_values, u->begin + u->count) - 1 );
      for (++l; l != u; ++l) 
        file_ids.insert( l->begin, l->begin + l->count - 1 );
    }
  }

  long buffer_size = bufferSize / read_size;
  while (!file_ids.empty()) {
    long first = file_ids.front();
    long count = file_ids.const_pair_begin()->second - first + 1;
    if (buffer_size < count)
      count = buffer_size;
    file_ids.erase( file_ids.begin(), file_ids.begin() + (count - 1) );
    
    mhdf_readDenseTagWithOpt( data, first - start_id, count, hdf_read_type, 
                              dataBuffer, ioProp, &status );
    if (is_error(status))
      return error(MB_FAILURE);
  
    if (MB_TYPE_HANDLE == mb_type) {
      rval = convert_id_to_handle( (MBEntityHandle*)dataBuffer, count * read_size / sizeof(MBEntityHandle) );
      if (MB_SUCCESS != rval)
        return error(rval);
    }
    
    MBRange ents;
    ents.insert( idMap.find( first ), idMap.find( first + count - 1 ) );
    rval = iFace->tag_set_data( tag_handle, ents, dataBuffer );
    if (MB_SUCCESS != rval)
      return error(MB_FAILURE);
  }
    
  return MB_SUCCESS;
}


MBErrorCode ReadHDF5::read_sparse_tag( MBTag tag_handle,
                                       hid_t hdf_read_type,
                                       hid_t id_table,
                                       hid_t value_table,
                                       long num_values )
{
  mhdf_Status status;
  MBErrorCode rval;
  MBDataType mbtype;
  
  rval = iFace->tag_get_data_type( tag_handle, mbtype );
  if (MB_SUCCESS != rval) 
    return error(rval);
  
  int read_size;
  rval = iFace->tag_get_size( tag_handle, read_size );
  if (MB_SUCCESS != rval) // wrong function for variable-length tags
    return error(rval);
  if (MB_TYPE_BIT == mbtype) 
    read_size = (read_size + 7)/8; // convert bits to bytes, plus 7 for ceiling
    
  if (hdf_read_type) { // if not opaque
    hsize_t hdf_size = H5Tget_size( hdf_read_type );
    if (hdf_size != (hsize_t)read_size) 
      return error(MB_FAILURE);
  }

  const int handles_per_tag = read_size/sizeof(MBEntityHandle);
  
  
    // Split buffer into two portions: one for handles and one for data.
    
    // We want to read the same number of handles as data values for each
    // iteration.  Calculate the total number of entries to read in each
    // pass as the size of the buffer over the sum of the size of a handle
    // and value.  Subtract off the size of one value so we reserve space
    // for adjusting for data alignment.
  long chunk_size = (bufferSize - read_size) / (sizeof(MBEntityHandle) + read_size);
  
    // Use the first half of the buffer for the handles.
  MBEntityHandle* idbuf = (MBEntityHandle*)dataBuffer;
    // Use the latter portion of the buffer for data
  char* databuf = dataBuffer + (chunk_size * sizeof(MBEntityHandle));
    // To be safe, align tag data to the size of an entire tag value
  if ((size_t)databuf % read_size)
    databuf += read_size - ((size_t)databuf % read_size);
      // Make sure the above calculations are correct
  assert( databuf + chunk_size*read_size < dataBuffer + bufferSize );
  
  long remaining = num_values;
  long offset = 0;
  while (remaining)
  {
    long count = std::min( chunk_size, remaining );
    long next_offset = offset + count;
    remaining -= count;
    
    mhdf_readSparseTagEntitiesWithOpt( id_table, offset, count, handleType, idbuf,
                                       ioProp, &status );
    if (is_error(status))
      return error(MB_FAILURE);
    
    rval = convert_id_to_handle( idbuf, count );
    if (MB_SUCCESS != rval)
     return error(rval);
     
      // Find contiguous blocks of values to read.
      // (convert_id_to_handle will return zero-valued handles for entities
      //  not read from the file--i.e. values to skip).
    long i = 0, j;
    for (;;) {
        // skip initial zeros
      while (i < count && 0 == idbuf[i])
        ++i;
        // if no more, stop
      if (i == count)
        break;
        // advance to end of range of non-zero handles
      j = i +1;
      while (j < count && 0 != idbuf[j])
        ++j;
        
        // if we've reached the end and the chunk is fairly small,
        // shift to the front of the buffer and read more to avoid
        // small data reads when blocks get split due to the buffer
        // size
      if (j == count && 2*i >= count && remaining > count) {
        long keep = j - i;
        memmove( idbuf, idbuf + i, sizeof(MBEntityHandle)*keep );
        long count2 = std::min( remaining - keep, count - keep );
        mhdf_readSparseTagEntitiesWithOpt( id_table, next_offset, count2, handleType, 
                                           idbuf+keep, ioProp, &status );
        if (is_error(status))
          return error(MB_FAILURE);
        offset += i;
        next_offset += count2;
        remaining -= count2;
        i = 0;

        rval = convert_id_to_handle( idbuf + keep, count2 );
        if (MB_SUCCESS != rval)
         return error(rval);
      
        j = keep;
        while (j < count && 0 != idbuf[j])
          ++j;
      }
      
        // read tag values 
      mhdf_readSparseTagValuesWithOpt( value_table, offset + i, j - i,
                                       hdf_read_type, databuf + i*read_size, 
                                       ioProp, &status );
      if (is_error(status))
        return error(MB_FAILURE);
      
      if (MB_TYPE_HANDLE == mbtype) {
        rval = convert_id_to_handle( (MBEntityHandle*)databuf + i*handles_per_tag, 
                                     (j-i) * handles_per_tag );
        if (MB_SUCCESS != rval)
          return error(rval);
      }
      
/*** FIX ME - need to do one at a time for BIT tags!  This is stupid. ***/
      if (mbtype == MB_TYPE_BIT)
      {
        rval = MB_SUCCESS;
        for ( ; MB_SUCCESS == rval && i < j; ++i)
          rval = iFace->tag_set_data( tag_handle, idbuf + i, 1, databuf + i );
      }
      else
      {
        rval = iFace->tag_set_data( tag_handle, idbuf + i, j - i, databuf + i*read_size );
        i = j;
      }
      if (MB_SUCCESS != rval)
        return error(rval);
    } // for(ever)
    
    
    offset = next_offset;
  } // while(remaining)
  
  return MB_SUCCESS;
}

#define assert_range( ARRAY, BYTES ) \
  assert( (char*)(ARRAY) >= dataBuffer && ((char*)(ARRAY)) + (BYTES) <= dataBuffer + bufferSize )

MBErrorCode ReadHDF5::read_var_len_tag( MBTag tag_handle,
                                        hid_t hdf_read_type,
                                        hid_t ent_table,
                                        hid_t val_table,
                                        hid_t off_table,
                                        long num_entities,
                                        long num_values )
{
  std::vector<char> tmp_buffer;
  mhdf_Status status;
  MBErrorCode rval;
  MBDataType mbtype;
  
  rval = iFace->tag_get_data_type( tag_handle, mbtype );
  if (MB_SUCCESS != rval) 
    return error(rval);
    
    // can't do variable-length bit tags
  if (MB_TYPE_BIT == mbtype)
    return error(MB_VARIABLE_DATA_LENGTH);

    // if here, MOAB tag must be variable-length
  int mbsize;
  if (MB_VARIABLE_DATA_LENGTH != iFace->tag_get_size( tag_handle, mbsize )) {
    assert(false);
    return error(MB_VARIABLE_DATA_LENGTH);
  }
  
  int read_size;
  if (hdf_read_type) {
    hsize_t hdf_size = H5Tget_size( hdf_read_type );
    if (hdf_size < 1)
      return error(MB_FAILURE);
    read_size = hdf_size;
  }
  else {
    // opaque
    read_size = 1;
  }
  
    // Subdivide buffer into 4 chunks:
    // Be careful of order so alignment is valid
    // 1) pointer array (input for MOAB) / offset array (read from file)
    // 2) entity handles (read as file IDs and converted to handles)
    // 3) tag sizes (calculated from file data)
    // 4) tag data
  const long avg_data_size = (num_values * read_size) / num_entities;
  const long per_ent_size = sizeof(long) + sizeof(MBEntityHandle) + sizeof(int);
  long num_ent = bufferSize / (per_ent_size + avg_data_size + avg_data_size/2);
  if (num_ent == 0) 
    num_ent = bufferSize / 4 / sizeof(void*);
  const long data_buffer_size = (bufferSize - num_ent * per_ent_size) / read_size;
  long* const end_idx_buffer = reinterpret_cast<long*>(dataBuffer);
  char* handle_buffer_start = reinterpret_cast<char*>(end_idx_buffer + num_ent+1);
  if (((size_t)handle_buffer_start) % sizeof(MBEntityHandle))
    handle_buffer_start += sizeof(MBEntityHandle) - ((size_t)handle_buffer_start);
  MBEntityHandle* const handle_buffer = reinterpret_cast<MBEntityHandle*>(handle_buffer_start);
  int* const size_buffer = reinterpret_cast<int*>(handle_buffer + num_ent);
  char* const data_buffer = reinterpret_cast<char*>(size_buffer + num_ent);
  
  
    // do num_ent blocks of entities
  long remaining = num_entities;
  long offset = 0;
  long prev_end_idx = -1;
  long prev_end = 0;
  while (remaining) {
    const long count = std::min( remaining, num_ent );
    long next_offset = offset + count;
    remaining -= count;
 
      // read entity IDs
    assert_range( handle_buffer, count * sizeof(MBEntityHandle) );
    mhdf_readSparseTagEntitiesWithOpt( ent_table, offset, count, handleType, 
                                       handle_buffer, ioProp, &status );
    if (is_error(status))
      return error(MB_FAILURE);
    
      // convert entity ID to MBEntityHandle
    rval = convert_id_to_handle( handle_buffer, count );
    if (MB_SUCCESS != rval)
      return error(rval);
     
      // Find contiguous blocks of values to read.
      // (convert_id_to_handle will return zero-valued handles for entities
      //  not read from the file--i.e. tag values to skip).
    long i = 0, j;
    for (;;) {
        // skip initial zeros
      while (i < count && 0 == handle_buffer[i])
        ++i;
        // if no more, stop
      if (i == count)
        break;
        // advance to end of range of non-zero handles
      j = i +1;
      while (j < count && 0 != handle_buffer[j])
        ++j;
        
        // if we've reached the end and the chunk is fairly small,
        // shift to the front of the buffer and read more to avoid
        // small data reads when blocks get split due to the buffer
        // size
      if (j == count && 2*i >= count && remaining > count) {
        long keep = j - i;
        memmove( handle_buffer, handle_buffer + i, sizeof(MBEntityHandle)*keep );
          // read more entity IDs
        long count2 = std::min( remaining - keep, count - keep );
        mhdf_readSparseTagEntitiesWithOpt( ent_table, next_offset, count2, handleType, 
                                           handle_buffer+keep, ioProp, &status );
        if (is_error(status))
          return error(MB_FAILURE);
        offset += i;
        next_offset += count2;
        remaining -= count2;
        i = 0;

        rval = convert_id_to_handle( handle_buffer + keep, count2 );
        if (MB_SUCCESS != rval)
         return error(rval);
      
        j = keep;
        while (j < count && 0 != handle_buffer[j])
          ++j;
      }
 
        // Read end-indices for each tagged entity.  Also need the
        // end of the previous entity.  If no previous entity, then -1.
      long remaining_ent = j - i;     
      if (offset + i == prev_end) {
        end_idx_buffer[0] = prev_end_idx;
        mhdf_readSparseTagIndicesWithOpt( off_table, offset + i, 
                                          remaining_ent, H5T_NATIVE_LONG, 
                                          end_idx_buffer + 1, ioProp, &status );
      }
      else {
        mhdf_readSparseTagIndicesWithOpt( off_table, offset + i - 1, 
                                          remaining_ent + 1, H5T_NATIVE_LONG, 
                                          end_idx_buffer, ioProp, &status );
      }
      if (is_error(status))
        return error(MB_FAILURE);
      prev_end = offset + i + remaining_ent;
        
        
      prev_end_idx = end_idx_buffer[0];
      long* idx_buffer = end_idx_buffer + 1;
      while (remaining_ent) {
          // Count the subset of remaining_ent entities for which the
          // concatentated tag data will fit in the data buffer.
        long tag_count = 0;
        while (tag_count < remaining_ent && 
               idx_buffer[tag_count] - prev_end_idx < data_buffer_size)
          ++tag_count;
        
          // If the tag value for a single entity doesn't fit in the buffer,
          // then we need to allocate a temporary buffer for it.
        long val_count;
        char* memptr = data_buffer;
        if (tag_count == 0) {
          tag_count = 1;
          val_count = idx_buffer[0] - prev_end_idx;
          tmp_buffer.resize(val_count * read_size);
          memptr = &tmp_buffer[0];
        }
        else {
          val_count = idx_buffer[tag_count-1] - prev_end_idx;
        }

          // Read the tag data
        mhdf_readSparseTagValues( val_table, prev_end_idx + 1,
                                  val_count, hdf_read_type, memptr, &status );
        if (is_error(status))
          return error(MB_FAILURE);
      
        if (MB_TYPE_HANDLE == mbtype) {
          rval = convert_id_to_handle( (MBEntityHandle*)memptr, val_count );
          if (MB_SUCCESS != rval)
            return error(rval);
        }
        
          // Calculate sizes from end indices
        size_buffer[0] = read_size * (idx_buffer[0] - prev_end_idx);
        for (long k = 1; k < tag_count; ++k)
          size_buffer[k] = read_size * (idx_buffer[k] - idx_buffer[k-1]);
          // save prev_end_index for next iteration because it will
          // get overwritten in the buffer next.
        prev_end_idx = idx_buffer[tag_count-1];
          // store pointers to individual tag values in end_idx_buffer
        assert(sizeof(long) == sizeof(void*));
        const void** ptr_buffer = reinterpret_cast<const void**>(idx_buffer);
        const char* ptr_iter = memptr;
        ptr_buffer[0] = ptr_iter;
        ptr_iter += size_buffer[0];
        for (long k = 1; k < tag_count; ++k) {
          ptr_buffer[k] = ptr_iter;
          ptr_iter += size_buffer[k];
        }
      
          // Pass data to MOAB to store for each tag
        rval = iFace->tag_set_data( tag_handle, handle_buffer+i, tag_count, 
                                    ptr_buffer, size_buffer );
        if (MB_SUCCESS != rval)
          return error(rval);
        remaining_ent -= tag_count;
        i +=  tag_count;
        idx_buffer += tag_count;
      }
        
    } // for(ever)
    
    
    offset = next_offset;
  } // while(remaining)
  
  return MB_SUCCESS;
}

MBErrorCode ReadHDF5::convert_id_to_handle( MBEntityHandle* array, 
                                            size_t size )
{
  convert_id_to_handle( array, size, idMap );
  return MB_SUCCESS;
}

void ReadHDF5::convert_id_to_handle( MBEntityHandle* array, 
                                     size_t size,
                                     const RangeMap<long,MBEntityHandle>& id_map )
{
  for (MBEntityHandle* const end = array + size; array != end; ++array)
    *array = id_map.find( *array );
}

void ReadHDF5::convert_id_to_handle( MBEntityHandle* array, 
                                     size_t size, size_t& new_size,
                                     const RangeMap<long,MBEntityHandle>& id_map )
{
  RangeMap<long,MBEntityHandle>::const_iterator it;
  new_size = 0;
  for (size_t i = 0; i < size; ++i) {
    it = id_map.lower_bound( array[i] );
    if (it != id_map.end() && it->begin <= (long)array[i])
      array[new_size++] = it->value + (array[i] - it->begin);
  }
}

void ReadHDF5::convert_range_to_handle( const MBEntityHandle* ranges,
                                        size_t num_ranges,
                                        const RangeMap<long,MBEntityHandle>& id_map,
                                        MBRange& merge )
{
  RangeMap<long,MBEntityHandle>::iterator it = id_map.begin();
  for (size_t i = 0; i < num_ranges; ++i) {
    long id = ranges[2*i];
    const long end = id + ranges[2*i+1];
    it = id_map.lower_bound( it, id_map.end(), id );
    if (id < it->begin)
      id = it->begin;
    const long off = id - it->begin;
    while (id < end) {
      long count = std::min( it->count - off,  end - id );
      merge.insert( it->value + off, it->value + off + count - 1 );
      id += count;
      if (id < end)
        ++it;
    }
  }
}

MBErrorCode ReadHDF5::convert_range_to_handle( const MBEntityHandle* array,
                                               size_t num_ranges,
                                               MBRange& range )
{
  convert_range_to_handle( array, num_ranges, idMap, range );
  return MB_SUCCESS;
}
  

MBErrorCode ReadHDF5::read_qa( MBEntityHandle import_set )
{
  mhdf_Status status;
  std::vector<std::string> qa_list;
  
  int qa_len;
  char** qa = mhdf_readHistory( filePtr, &qa_len, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( "%s", mhdf_message( &status ) );
    return error(MB_FAILURE);
  }
  qa_list.resize(qa_len);
  for (int i = 0; i < qa_len; i++)
  {
    qa_list[i] = qa[i];
    free( qa[i] );
  }
  free( qa );
  
  /** FIX ME - how to put QA list on set?? */

  return MB_SUCCESS;
}

  
    
