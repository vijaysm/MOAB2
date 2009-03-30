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

#undef DEBUG

#ifdef DEBUG
#  define DEBUGOUT(A) fputs( A, stderr )
#  include <stdio.h>
#else
#  define DEBUGOUT(A)
#endif

#define READ_HDF5_BUFFER_SIZE (40*1024*1024)

MBReaderIface* ReadHDF5::factory( MBInterface* iface )
  { return new ReadHDF5( iface ); }

ReadHDF5::ReadHDF5( MBInterface* iface )
  : bufferSize( READ_HDF5_BUFFER_SIZE ),
    dataBuffer( 0 ),
    iFace( iface ), 
    filePtr( 0 ), 
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
    return MB_FAILURE;
  
  if (H5Tset_size( handleType, sizeof(MBEntityHandle)) < 0)
  {
    H5Tclose( handleType );
    return MB_FAILURE;
  }
  
  void* ptr = 0;
  rval = iFace->query_interface( "MBReadUtilIface", &ptr );
  if (MB_SUCCESS != rval)
  {
    H5Tclose( handleType );
    return rval;
  }
  readUtil = reinterpret_cast<MBReadUtilIface*>(ptr);
  
  setSet.first_id = 0;
  setSet.type2 = mhdf_set_type_handle();
  setSet.type = MBENTITYSET;
  nodeSet.first_id = 0;
  nodeSet.type2 = mhdf_node_type_handle();
  nodeSet.type = MBVERTEX;
  
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
                                 const int*, const int )
{
  MBErrorCode rval;
  mhdf_Status status;
  ioProp = H5P_DEFAULT;

  if (name) {
    readUtil->report_error( "Reading subset of files not supported for HDF5." );
    return MB_UNSUPPORTED_OPERATION;
  }

  if (MB_SUCCESS != init())
    return MB_FAILURE;

  bool use_mpio = (MB_SUCCESS == opts.get_null_option("USE_MPIO"));
  if (use_mpio) {
#ifndef HDF5_PARALLEL
    return MB_NOT_IMPLEMENTED;
#else
    int parallel_mode;
    rval = opts.match_option( "PARALLEL", 
                              ReadParallel::parallelOptsNames, 
                              parallel_mode );
    if (MB_FAILURE == rval) {
      readUtil->report_error("Unexpected value for 'PARALLEL' option\n");
      return MB_FAILURE;
    }
    else if (MB_SUCCESS != rval ||
             parallel_mode != ReadParallel::POPT_READ_DELETE) {
      use_mpio = false;
    }
#endif
  }
  
  dataBuffer = (char*)malloc( bufferSize );
  if (!dataBuffer)
    return MB_MEMORY_ALLOCATION_FAILED;

  rval = iFace->create_meshset( MESHSET_SET, file_set );
  if (MB_SUCCESS != rval) {
    free(dataBuffer);
    return rval;
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
    return MB_FAILURE;
  }
    
  
  rval = load_file_impl( file_set, use_mpio );
  mhdf_closeFile( filePtr, &status );
  filePtr = 0;
  if (ioProp != H5P_DEFAULT)
    H5Pclose( ioProp );
  if (mhdf_isError( &status )) {
    if (MB_SUCCESS == rval)
      rval = MB_FAILURE;
    readUtil->report_error( mhdf_message( &status ));
  }
        
      // delete everything that was read in if read failed part-way through
  if (MB_SUCCESS != rval) {
    iFace->delete_entities( &file_set, 1 );
    file_set = 0;
    iFace->delete_entities( setSet.range );
    for (std::list<ElemSet>::reverse_iterator rel_itor = elemList.rbegin();
         rel_itor != elemList.rend(); ++rel_itor)
      iFace->delete_entities( rel_itor->range );
    iFace->delete_entities( nodeSet.range );
  }
  
  elemList.clear();
  nodeSet.range.clear();
  setSet.range.clear();
  free( dataBuffer );
  return rval;
}
  


MBErrorCode ReadHDF5::load_file_impl( 
                                 MBEntityHandle file_set, 
                                 bool use_mpio )
{
  MBErrorCode rval;
  mhdf_Status status;
  std::string tagname;
  int num_tags = 0;
  char** tag_names = NULL;
  char** groups = NULL;
  std::list<ElemSet>::iterator el_itor;
  unsigned int i, num_groups;
  bool have_nodes = true;

DEBUGOUT("Reading Nodes.\n");
  
  rval = read_nodes();
  if (MB_FILE_WRITE_ERROR == rval) {
    DEBUGOUT("No nodes in file.!\n");
    have_nodes = false;
  }
  else if (MB_SUCCESS != rval)
    return rval;

DEBUGOUT("Reading element connectivity.\n");

  groups = mhdf_getElemHandles( filePtr, &num_groups, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ));
    return MB_FAILURE;
  }

  for (i = 0; i < num_groups; ++i)
  {
    int poly = mhdf_isPolyElement( filePtr, groups[i], &status );
    if (mhdf_isError( &status ))
    {
      readUtil->report_error( mhdf_message( &status ));
      free(groups);
      return MB_FAILURE;
    }

    if (poly)
      rval = read_poly( groups[i] );
    else
      rval = read_elems( groups[i] );
      
    if (MB_SUCCESS != rval) {
      free( groups );
      return rval;
    }
  }
  
DEBUGOUT("Reading sets.\n");
  
  rval = read_sets();
  if (rval != MB_SUCCESS) {
    free( groups );
    return rval;
  }
  
DEBUGOUT("Reading adjacencies.\n");
  
  rval = read_adjacencies( nodeSet );
  if (rval != MB_SUCCESS) {
    free( groups );
    return rval;
  }
  for (el_itor = elemList.begin(); el_itor != elemList.end(); ++el_itor) {
    rval = read_adjacencies( *el_itor );
    if (MB_SUCCESS != rval) {
      free( groups );
      return rval;
    }
  }

DEBUGOUT("Reading tags.\n");
  
  tag_names = mhdf_getTagNames( filePtr, &num_tags, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ));
    free( groups );
    return MB_FAILURE;
  }
  
  for (int t = 0; t < num_tags; ++t)
  {
    rval = read_tag( tag_names[t] );
    free( tag_names[t] );
    if (MB_SUCCESS != rval) {
      for (; t < num_tags; ++t)
        free( tag_names[t] );
      free( tag_names );
      free( groups );
      return MB_FAILURE;
    }
  }
  free( tag_names );

  free( groups );
  if (file_set) {
    DEBUGOUT("Creating entity set for file contents\n")
    rval = iFace->add_entities( file_set, nodeSet.range );
    if (MB_SUCCESS != rval)
      return rval;
    for (el_itor = elemList.begin(); el_itor != elemList.end(); ++el_itor) {
      rval = iFace->add_entities( file_set, el_itor->range);
      if (MB_SUCCESS != rval)
        return rval;
    }
    rval = iFace->add_entities( file_set, setSet.range );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
DEBUGOUT("Finishing read.\n");
  rval = read_qa( file_set );
  return rval;
}

MBErrorCode ReadHDF5::read_nodes()
{
  MBErrorCode rval;
  mhdf_Status status;
  long count, first_id;
  int dim;
  MBRange range;
  
  int cdim;
  rval = iFace->get_dimension( cdim );
  if (MB_SUCCESS != rval)
    return rval;
  
  hid_t data_id = mhdf_openNodeCoords( filePtr, &count, &dim, &first_id, &status );
  if (mhdf_isError( &status ))
  {
    //readUtil->report_error( mhdf_message( &status ));
    // Failed because no node data in file?
    nodeSet.range.clear();
    nodeSet.first_id = std::numeric_limits<long>::max();
    nodeSet.type = MBVERTEX;
    nodeSet.type2 = mhdf_node_type_handle();
    return MB_FILE_WRITE_ERROR;
  }
  
  if (cdim < dim)
  {
    rval = iFace->set_dimension( dim );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
  MBEntityHandle handle;
  std::vector<double*> arrays(dim);
  rval = readUtil->get_node_arrays( dim, (int)count, (int)first_id, 
                                    handle, arrays );
  if (MB_SUCCESS != rval)
  {
    mhdf_closeData( filePtr, data_id, &status );
    return rval;
  }
  
  nodeSet.range.clear();
  nodeSet.range.insert( handle, handle + count - 1 );
  nodeSet.first_id = first_id;
  nodeSet.type = MBVERTEX;
  nodeSet.type2 = mhdf_node_type_handle();
  for (int i = 0; i < dim; i++)
  {
    mhdf_readNodeCoordWithOpt( data_id, 0, count, i, arrays[i], ioProp, &status );
    if (mhdf_isError( &status ))
    {
      readUtil->report_error( mhdf_message(&status) );
      mhdf_closeData( filePtr, data_id, &status );
      return MB_FAILURE;
    }
  }
  for (int j = dim; j < cdim; j++)
    memset( arrays[j], 0, count * sizeof(double) );
  
  mhdf_closeData( filePtr, data_id, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message(&status) );
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

MBErrorCode ReadHDF5::read_elems( const char* elem_group )
{
  MBErrorCode rval;
  mhdf_Status status;
  char name[64];
  
    // Put elem set in list early so clean up code can 
    // get rid of them if we fail.
  ElemSet empty_set;
  empty_set.type2 = elem_group;
  elemList.push_back( empty_set );
  std::list<ElemSet>::iterator it = elemList.end();
  --it;
  ElemSet& elems = *it;
  
  mhdf_getElemTypeName( filePtr, elem_group, name, sizeof(name), &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  elems.type = MBCN::EntityTypeFromName( name );
  if (elems.type == MBMAXTYPE)
  {
    readUtil->report_error( "Unknown element type: \"%s\".\n", name );
    return MB_FAILURE;
  }
  
  int nodes_per_elem;
  long count, first_id;
  hid_t data_id = mhdf_openConnectivity( filePtr, 
                                         elem_group,
                                         &nodes_per_elem,
                                         &count,
                                         &first_id,
                                         &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }
  elems.first_id = first_id;
  
  MBEntityHandle handle;
  MBEntityHandle* array;
  rval = readUtil->get_element_array( (int)count,
                                       nodes_per_elem,
                                       elems.type,
                                      (int)first_id,
                                       handle, 
                                       array );
  if (MB_SUCCESS != rval)
  {
    mhdf_closeData( filePtr, data_id, &status );
    return rval;
  }
  
  elems.range.insert( handle, handle + count - 1 );
  mhdf_readConnectivityWithOpt( data_id, 0, count, handleType, array, ioProp, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    mhdf_closeData( filePtr, data_id, &status );
    return MB_FAILURE;
  }
  
  mhdf_closeData( filePtr, data_id, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  if (elems.type == MBPOLYHEDRON)
    rval = convert_id_to_handle( array, (size_t)(nodes_per_elem*count) );
  else
    rval = convert_id_to_handle( nodeSet, array, (size_t)(nodes_per_elem*count) );

  if (MB_SUCCESS != rval) return rval;
  
    // notify MOAB of the new elements
  MBErrorCode result = readUtil->update_adjacencies(handle, count,
                                                    nodes_per_elem, array);
  if (MB_SUCCESS != result) return result;

  return rval;
}

MBErrorCode ReadHDF5::read_poly( const char* elem_group )
{
  MBErrorCode rval;
  mhdf_Status status;
  char name[64];
   
  mhdf_getElemTypeName( filePtr, elem_group, name, sizeof(name), &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }
  MBEntityType type = MBCN::EntityTypeFromName( name );

  long count, first_id, data_len;
  hid_t handles[2];
  mhdf_openPolyConnectivity( filePtr, elem_group, &count, &data_len,
                             &first_id, handles, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
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
      return MB_FAILURE;
    }
    
    connectivity.resize( connend - prevend );
    mhdf_readPolyConnIDsWithOpt( handles[1], prevend+1, connectivity.size(), handleType,
                          &connectivity[0], ioProp, &status );
    if (mhdf_isError( &status ))
    {
      readUtil->report_error( mhdf_message( &status ) );
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      return MB_FAILURE;
    }
    
    rval= convert_id_to_handle( &connectivity[0], connectivity.size() );
    if (MB_SUCCESS != rval) 
    {
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      return rval;
    }
    
    rval = iFace->create_element( type, &connectivity[0], connectivity.size(), h );
    if (MB_SUCCESS != rval) 
    {
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      return rval;
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
    result = MB_FAILURE;
  }
  mhdf_closeData( filePtr, handles[1], &status );
  if (mhdf_isError( &status )) {
    readUtil->report_error( mhdf_message( &status ));
    result = MB_FAILURE;
  }
  return result;
}

template <typename T>
class auto_array {
private:
  T* data;
public:
  auto_array( size_t s ) : data(new T[s]) {}
  auto_array()           : data(0)        {}
  ~auto_array()                         { delete [] data; }
  T*       get       ()                 { return data;    }
  const T* get       ()           const { return data;    }
  T&       operator[]( size_t i )       { return data[i]; }
  const T& operator[]( size_t i ) const { return data[i]; }
  T&       operator* ()                 { return *data;   }
  const T* operator* ()           const { return *data;   }
  T*       operator->()                 { return data;    }
  const T* operator->()           const { return data;    }
// auto_ptr-style destrutive assigment
  auto_array( auto_array<T>& c ) {
    data = c.data;
    const_cast<auto_array<T>&>(c).data = 0;
  }
  auto_array<T>& operator=( auto_array<T>& c ) {
    data = c.data;
    const_cast<auto_array<T>&>(c).data = 0;
  }
};
  

MBErrorCode ReadHDF5::read_set_contents( hid_t meta_id, hid_t data_id,
                                         const unsigned long data_len )
{
  MBErrorCode rval;
  mhdf_Status status;
  
    // create extra buffers for storing end indices and flags
  const unsigned long offset_size = bufferSize / sizeof(long);
  auto_array<long> offsets(offset_size);
  auto_array<unsigned short> flags(offset_size);
  
  MBEntityHandle* buffer = (MBEntityHandle*)dataBuffer;
  size_t chunk_size = bufferSize / sizeof(MBEntityHandle);
  if (chunk_size % 2)
    --chunk_size; // makes reading range data easier.
  
  unsigned long set_offset = 0;  /* running offset into description table */
  unsigned long sets_remaining = setSet.range.size();
  unsigned long file_offset = 0; /* running offset into child table */
  MBRange::const_iterator set_iter = setSet.range.begin();
  while (sets_remaining) {
      // read end indices from meta data
    unsigned long set_count = sets_remaining < offset_size ? sets_remaining : offset_size;
    mhdf_readSetContentEndIndicesWithOpt( meta_id, set_offset, set_count, 
                                   H5T_NATIVE_LONG, offsets.get(), ioProp, &status );
    if (mhdf_isError( &status )) {
      readUtil->report_error( mhdf_message( &status ) );
      return MB_FAILURE;
    }
    mhdf_readSetFlagsWithOpt( meta_id, set_offset, set_count, H5T_NATIVE_USHORT, 
                       flags.get(), ioProp, &status );
    if (mhdf_isError( &status )) {
      readUtil->report_error( mhdf_message( &status ) );
      return MB_FAILURE;
    }
    
      // now iterate over sets, reading set content lists
    unsigned long i, r = 0; // number of sets for which data has been read.
    while (r < set_count) {
        // figure out how many we can read at once
      for (i = r; i < set_count; ++i)
        if ((unsigned long)(offsets[i] + 1 - file_offset) > chunk_size)
          break;
    
        // special case: set content list for single set greater than buffer
      if (i == r) {
           // Check if set contents are stored as ranges or a simple list
        bool ranged = (0 != (flags[r] & (unsigned short)mhdf_SET_RANGE_BIT));

        assert( set_iter != setSet.range.end() );
        MBEntityHandle h = *set_iter; ++set_iter;
        size_t remaining = offsets[r] + 1 - file_offset;
        if (remaining > data_len - file_offset) {
          readUtil->report_error( "Invalid set contents offset read from file." );
          return MB_FAILURE;
        }
        while (remaining)
        {
          size_t count = remaining > chunk_size ? chunk_size : remaining;
          remaining -= count;
          mhdf_readSetDataWithOpt( data_id, file_offset, count, handleType, buffer, ioProp, &status );
          if (mhdf_isError( &status )) {
            readUtil->report_error( mhdf_message( &status ) );
            return MB_FAILURE;
          }
          file_offset += count;

            // convert data from file ids to MBEntityHandles and add to set
          if (ranged)
          {
            if (count % 2 != 0) {
              readUtil->report_error( "Invalid ranged set contents spec." );
              return MB_FAILURE;
            }
            MBRange range;
            rval = convert_range_to_handle( buffer, count / 2, range );
            if (MB_SUCCESS != rval) {
              readUtil->report_error( "Invalid entities in set contents" );
              return rval;
            }
            rval = iFace->add_entities( h, range );
            if (MB_SUCCESS != rval)
              return rval;
          }
          else
          {
            rval = convert_id_to_handle( buffer, count );
            if (MB_SUCCESS != rval) {
              readUtil->report_error( "Invalid entities in set contents" );
              return rval;
            }
            rval = iFace->add_entities( h, buffer, count );
            if (MB_SUCCESS != rval)
              return rval;
          }
        } // while(remaining)
        ++r;
      } // end special case (really big set)
      
        // normal case - read data for several sets
      else {
        
          // read data for sets in [r,i)
          
        size_t count = offsets[i-1] + 1 - file_offset;
        if (count > data_len - file_offset) {
          readUtil->report_error( "Invalid set contents offset read from file." );
          return MB_FAILURE;
        }
          
        mhdf_readSetDataWithOpt( data_id, file_offset, count, handleType, buffer, ioProp, &status );
        if (mhdf_isError( &status )) {
          readUtil->report_error( mhdf_message( &status ) );
          return MB_FAILURE;
        }
          
          // add contents to each set
        size_t mem_offset = 0;
        for (; r < i; ++r) {
          bool ranged = (0 != (flags[r] & (long)mhdf_SET_RANGE_BIT));
          count = offsets[r] + 1 - file_offset;
          assert( set_iter != setSet.range.end() );
          MBEntityHandle h = *set_iter; ++set_iter;

          if (ranged)
          {
            if (count % 2 != 0) {
              readUtil->report_error( "Invalid ranged set contenst spec." );
              return MB_FAILURE;
            }
            MBRange range;
            rval = convert_range_to_handle( buffer+mem_offset, count / 2, range );
            if (MB_SUCCESS != rval) {
              readUtil->report_error( "Invalid entities in set contents" );
              return rval;
            }
            rval = iFace->add_entities( h, range );
            if (MB_SUCCESS != rval)
              return rval;
          }
          else
          {
            rval = convert_id_to_handle( buffer+mem_offset, count );
            if (MB_SUCCESS != rval) {
              readUtil->report_error( "Invalid entities in set contents" );
              return rval;
            }
            rval = iFace->add_entities( h, buffer+mem_offset, count );
            if (MB_SUCCESS != rval)
              return rval;
          }

          file_offset += count;
          mem_offset += count;
        }
      }
    } // while(r < num_sets)
    
    set_offset += set_count;
    sets_remaining -= set_count;
  } // while (sets_remaining)

  assert( set_iter == setSet.range.end() );
  return MB_SUCCESS;
}

MBErrorCode ReadHDF5::read_parents_children( bool parents,
                                             hid_t meta_id, 
                                             hid_t data_id,
                                             const unsigned long data_len )
{
  MBErrorCode rval;
  mhdf_Status status;
  
  // create an extra buffer for storing offsets
  const unsigned long offset_size = bufferSize / sizeof(long);
  auto_array<long> offsets( offset_size );
  // use the existing buffer for storing set child lists
  MBEntityHandle* buffer = (MBEntityHandle*)dataBuffer;
  size_t chunk_size = bufferSize / sizeof(MBEntityHandle);
  const size_t total_sets = setSet.range.size();
  
  unsigned long set_offset = 0;  /* running offset into description table */
  unsigned long sets_remaining = setSet.range.size();
  unsigned long file_offset = 0; /* running offset into child table */
  MBRange::const_iterator set_iter = setSet.range.begin();
  while (sets_remaining) {
      // read end indices from meta data
    unsigned long set_count = sets_remaining < offset_size ? sets_remaining : offset_size;
    if (parents)
      mhdf_readSetParentEndIndicesWithOpt( meta_id, set_offset, set_count, 
                                    H5T_NATIVE_LONG, offsets.get(), 
                                    ioProp, &status );
    else
      mhdf_readSetChildEndIndicesWithOpt(  meta_id, set_offset, set_count, 
                                    H5T_NATIVE_LONG, offsets.get(), 
                                    ioProp, &status );
    if (mhdf_isError( &status )) {
      readUtil->report_error( mhdf_message( &status ) );
      return MB_FAILURE;
    }
    
      // now iterate over sets, reading parent/child lists
    unsigned long i, r = 0; // number of sets for which data has been read.
    while (r < set_count) {
        // figure out how many we can read at once
      for (i = r; i < set_count; ++i)
        if ((unsigned long)(offsets[i] + 1 - file_offset) > chunk_size)
          break;
    
        // special case: children of one set greater than buffer
      if (i == r) {
        assert( set_iter != setSet.range.end() );
        MBEntityHandle h = *set_iter; ++set_iter;
        size_t remaining = offsets[r] + 1 - file_offset;
        if (remaining > data_len - file_offset) {
          readUtil->report_error( "Invalid set %s offset read from file.",
                                  parents ? "parent" : "child" );
          return MB_FAILURE;
        }
        while (remaining)
        {
          size_t count = remaining > chunk_size ? chunk_size : remaining;
          remaining -= count;
          mhdf_readSetParentsChildrenWithOpt( data_id, file_offset, count, handleType, 
                                              buffer, ioProp, &status );
          if (mhdf_isError( &status )) {
            readUtil->report_error( mhdf_message( &status ) );
            return MB_FAILURE;
          }
          file_offset += count;

            // convert from file_ids to set handles
          for (size_t j = 0; j < count; ++j) {
            buffer[j] -= setSet.first_id;
            if (buffer[j] >= total_sets) { 
              readUtil->report_error("Invalid set %s ID", parents ? "parent" : "child" );
              return MB_FAILURE;
            }
            buffer[j] += setSet.range.front();
          }

          if (parents)
            rval = iFace->add_parent_meshsets( h, buffer, count );
          else
            rval = iFace->add_child_meshsets( h, buffer, count );
          if (MB_SUCCESS != rval)
            return rval;
        } // while(remaining)
        ++r;
      } // end special case (really big set)
      
        // normal case - read data for several sets
      else {
        
          // read data for sets in [r,i)
        size_t count = offsets[i-1] + 1 - file_offset;
        if (count > data_len - file_offset) {
          readUtil->report_error( "Invalid set %s offset read from file.",
                                  parents ? "parent" : "child" );
          return MB_FAILURE;
        }
        mhdf_readSetParentsChildrenWithOpt( data_id, file_offset, count, handleType, 
                                            buffer, ioProp, &status );
        if (mhdf_isError( &status )) {
          readUtil->report_error( mhdf_message( &status ) );
          return MB_FAILURE;
        }
        
          // convert from file_ids to set handles
        for (size_t j = 0; j < count; ++j) {
          buffer[j] -= setSet.first_id;
          if (buffer[j] >= total_sets) { 
            readUtil->report_error("Invalid set %s ID", parents ? "parent" : "child" );
            return MB_FAILURE;
          }
          buffer[j] += setSet.range.front();
        }
          
          // add children to each set
        size_t mem_offset = 0;
        for (; r < i; ++r) {
          assert( set_iter != setSet.range.end() );
          MBEntityHandle h = *set_iter; ++set_iter;
          count = offsets[r] + 1 - file_offset;
          if (parents)
            rval = iFace->add_parent_meshsets( h, buffer+mem_offset, count );
          else
            rval = iFace->add_child_meshsets(  h, buffer+mem_offset, count );
          if (MB_SUCCESS != rval)
            return rval;

          file_offset += count;
          mem_offset += count;
        }
      }
    } // while(r < num_sets)
    
    set_offset += set_count;
    sets_remaining -= set_count;
  } // while (sets_remaining)

  assert( set_iter == setSet.range.end() );
  return MB_SUCCESS;
}


MBErrorCode ReadHDF5::read_sets()
{
  MBErrorCode rval;
  mhdf_Status status;
  
    // Check what data is in the file for sets
  int have_sets, have_data, have_children, have_parents;
  have_sets = mhdf_haveSets( filePtr, &have_data, &have_children, &have_parents, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }

  if (!have_sets)
    return MB_SUCCESS;
  
    // Open the list of sets
  long num_sets, first_id;
  hid_t meta_id = mhdf_openSetMeta( filePtr, &num_sets, &first_id, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }

    // Initialize internal bookkeeping
  setSet.first_id = first_id;
  setSet.type = MBENTITYSET;
  setSet.type2 = mhdf_set_type_handle();
  if (!num_sets) { // shouldn't happen if have_sets == true, but...
    mhdf_closeData( filePtr, meta_id, &status );
    return MB_SUCCESS;
  }
  
    // Get last used MeshSet handle
  MBRange junk;
  iFace->get_entities_by_type( 0, MBENTITYSET, junk );
  MBEntityHandle start_handle = junk.empty() ? 1 : junk.back() + 1;
  junk.clear();
  bool first = true;
  
  
    // Iterate over set metadata, creating all the sets.
    // Don't read any contents or parent/child links yet.
    // Create all the sets first, so all the contents
    // and parent/child links are valid.
  MBEntityHandle handle;
  size_t chunk_size = bufferSize / sizeof(unsigned);
  unsigned * buffer = reinterpret_cast<unsigned*>(dataBuffer);
  size_t remaining = num_sets, offset = 0;
  MBRange::iterator insert_iter = setSet.range.begin();
  while (remaining) {
      // Get a block of set flags
    size_t count = remaining > chunk_size ? chunk_size : remaining;
    mhdf_readSetFlagsWithOpt( meta_id, offset, count, H5T_NATIVE_UINT, 
                              buffer, ioProp, &status );
    if (mhdf_isError( &status )) {
      readUtil->report_error( mhdf_message( &status ) );
      mhdf_closeData( filePtr, meta_id, &status );
      return MB_FAILURE;
    }
    offset += count;
    remaining -= count;
    
      // clear ranged-storage bit.  Its internal data for the file
      // format, not one of MOAB's set flags.
    for (size_t i = 0;i < count; ++i)
      buffer[i] &= ~(unsigned)mhdf_SET_RANGE_BIT;
      
      // create block of sets
    rval = readUtil->create_entity_sets( count,
                                         buffer, 
                                         ID_FROM_HANDLE(start_handle),
                                         handle );
    if (MB_SUCCESS != rval) {
      mhdf_closeData( filePtr, meta_id, &status );
      return rval;
    }
    
      // We are careful to request start IDs such that the
      // resulting MBEntityHandles are always increasing.
      // We are relying on increasing handles, so make sure
      // that's what we get.
    if (!first && handle < start_handle) {
      readUtil->report_error( "Non-increasing handle space for mesh sets" );
      mhdf_closeData( filePtr, meta_id, &status );
      return MB_FAILURE;
    }
    first = false;
    start_handle = handle + count;
    
    insert_iter = setSet.range.insert( insert_iter, handle, handle + count - 1 );
  }
  assert( setSet.range.size() == (size_t)num_sets );
  
  MBErrorCode result = MB_SUCCESS;
  if (have_data) {
    long data_len; 
    hid_t data_id = mhdf_openSetData( filePtr, &data_len, &status );
    if (mhdf_isError(&status)) {
      readUtil->report_error( mhdf_message( &status ) );
      result = MB_FAILURE;
    }
    else {
      rval = read_set_contents( meta_id, data_id, data_len );
      mhdf_closeData( filePtr, data_id, &status );
      if (MB_SUCCESS != rval)
        result = rval;
    }
  }

  if (have_children) {
    long data_len; 
    hid_t data_id = mhdf_openSetChildren( filePtr, &data_len, &status );
    if (mhdf_isError(&status)) {
      readUtil->report_error( mhdf_message( &status ) );
      result = MB_FAILURE;
    }
    else {
      rval = read_parents_children( false, meta_id, data_id, data_len );
      mhdf_closeData( filePtr, data_id, &status );
      if (MB_SUCCESS != rval)
        result = rval;
    }
  }

  if (have_parents) {
    long data_len; 
    hid_t data_id = mhdf_openSetParents( filePtr, &data_len, &status );
    if (mhdf_isError(&status)) {
      readUtil->report_error( mhdf_message( &status ) );
      result = MB_FAILURE;
    }
    else {
      rval = read_parents_children( true, meta_id, data_id, data_len );
      mhdf_closeData( filePtr, data_id, &status );
      if (MB_SUCCESS != rval)
        result = rval;
    }
  }
  
  mhdf_closeData( filePtr, meta_id, &status );
  return result;
}

MBErrorCode ReadHDF5::read_adjacencies( ElemSet& elems )
{
  MBErrorCode rval;
  mhdf_Status status;
  hid_t table;
  long data_len;
  
  int adj = mhdf_haveAdjacency( filePtr, elems.type2, &status );
  if (mhdf_isError(&status))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  if (!adj)
    return MB_SUCCESS;
  
  table = mhdf_openAdjacency( filePtr, elems.type2, &data_len, &status );
  if (mhdf_isError(&status))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  MBEntityHandle* buffer = (MBEntityHandle*)dataBuffer;
  size_t chunk_size = bufferSize / sizeof(MBEntityHandle);
  size_t remaining = data_len;
  size_t leading = 0;
  size_t offset = 0;
  while (remaining)
  {
    size_t count = remaining > chunk_size ? chunk_size : remaining;
    count -= leading;
    remaining -= count;
    
    mhdf_readAdjacencyWithOpt( table, offset, count, handleType, buffer + leading,
                               ioProp, &status );
    if (mhdf_isError(&status))
    {
      readUtil->report_error( mhdf_message( &status ) );
      mhdf_closeData( filePtr, table, &status );
      return MB_FAILURE;
    }
    
    MBEntityHandle* iter = buffer;
    MBEntityHandle* end = buffer + count + leading;
    while (end - iter >= 3)
    {
      rval = convert_id_to_handle( elems, iter, 1 );
      MBEntityHandle entity = *iter;
      MBEntityHandle count = *++iter;
      if (MB_SUCCESS != rval || count < 1)
      {
        assert(0);
        mhdf_closeData( filePtr, table, &status );
        return rval == MB_SUCCESS ? MB_FAILURE : rval;
      }
      ++iter;
      
      if (end < count + iter)
      {
        iter -= 2;
        break;
      }
      
      rval = convert_id_to_handle( iter, count );
      if (MB_SUCCESS != rval)
      {
        assert(0);
        mhdf_closeData( filePtr, table, &status );
        return rval;
      }
      
      rval = iFace->add_adjacencies( entity, iter, count, false );
      if (MB_SUCCESS != rval)
      {
        assert(0);
        mhdf_closeData( filePtr, table, &status );
        return rval;
      }
      
      iter += count;
    }
    
    leading = end - iter;
    memmove( buffer, iter, leading );
  }
  
  assert(!leading);  // unexpected truncation of data
  
  mhdf_closeData( filePtr, table, &status );
  if (mhdf_isError(&status))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;  
}

MBErrorCode ReadHDF5::read_tag( const char* name )
{
  MBErrorCode rval;
  mhdf_Status status;
  std::string tag_type_name;
  
  mhdf_TagDataType mhdf_type;  // Enum for tag data type
  int tag_size;                // Size of tag
  int mhdf_storage;            // TSTT storage type (dense vs. sparse)
  int have_default;            // File contains default value for tag
  int have_global;             // File contains global value for tag
  int have_sparse;             // File contains sparse data table for tag
  hid_t hdf_type = 0;          // Type to use when reading tag data.
  int elem_size;               // Bytes required for one elem of array
  int array_size;              // If tag is not opaque, the number of data per entity
  MBTag handle;                // The handle for the tag
  MBDataType mb_type;          // The MOAB data type for the data
  MBTagType storage;
  
    // Get description of tag
  mhdf_getTagInfo( filePtr, name, &mhdf_type, &tag_size, &mhdf_storage,
                   &have_default, &have_global, &have_sparse, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }
  switch (mhdf_storage) {
    case mhdf_DENSE_TYPE : storage = MB_TAG_DENSE ; break;
    case mhdf_SPARSE_TYPE: storage = MB_TAG_SPARSE; break;
    case mhdf_BIT_TYPE   : storage = MB_TAG_BIT;    break;
    case mhdf_MESH_TYPE  : storage = MB_TAG_MESH;   break;
    default:
      readUtil->report_error( "Invalid storage type for tag '%s': %d\n", name, mhdf_storage );
      return MB_FAILURE;
  }

    // Type-specific stuff
  switch (mhdf_type)
  {
    case mhdf_BITFIELD:
    
      if (!tag_size || tag_size > 8)
      {
        readUtil->report_error( "Invalid bit tag:  class is MB_TAG_BIT, num bits = %d\n", tag_size );
        return MB_FAILURE;
      }
    
      elem_size = tag_size;
      tag_size = 1;
      array_size = 1;
      hdf_type = H5T_NATIVE_B8;
      mb_type = MB_TYPE_BIT;
      break;
      
    case mhdf_INTEGER:
    
      elem_size = sizeof(int);
      array_size = tag_size;
      hdf_type = H5T_NATIVE_INT;
      tag_size = elem_size * array_size;
      mb_type = MB_TYPE_INTEGER;
      break;
      
    case mhdf_FLOAT:
      
      elem_size = sizeof(double);
      array_size = tag_size;
      hdf_type = H5T_NATIVE_DOUBLE;
      tag_size = elem_size * array_size;
      mb_type = MB_TYPE_DOUBLE;
      break;
    
    case mhdf_BOOLEAN:
      
      elem_size = sizeof(int);
      array_size = tag_size;
      hdf_type = H5T_NATIVE_UINT;
      tag_size = elem_size * array_size;
      mb_type = MB_TYPE_INTEGER;
      break;
    
    case mhdf_ENTITY_ID:
      
      elem_size = sizeof(MBEntityHandle);
      array_size = tag_size;
      hdf_type = handleType;
      tag_size = elem_size * array_size;
      mb_type = MB_TYPE_HANDLE;
      break;
      
    case mhdf_OPAQUE:
    default:

      if (tag_size < 0) { // variable-length
        elem_size = 1;
        array_size = -1;
      }
      else {
        elem_size = tag_size;
        array_size = 1;
      }
      mb_type = MB_TYPE_OPAQUE;
      hdf_type = (hid_t)0;
      
        // Check for user-provided type
      MBTag type_handle;
      tag_type_name = "__hdf5_tag_type_";
      tag_type_name += name;
      rval = iFace->tag_get_handle( tag_type_name.c_str(), type_handle );
      if (MB_SUCCESS == rval)
      {
        rval = iFace->tag_get_data( type_handle, 0, 0, &hdf_type );
        if (MB_SUCCESS == rval)
          elem_size = H5Tget_size( hdf_type );
        else if (MB_TAG_NOT_FOUND != rval)
          return rval;
      }
      else if (MB_TAG_NOT_FOUND != rval)
        return rval;

      break;
  }
  if (tag_size < 0)  // variable length
    tag_size = MB_VARIABLE_LENGTH;

  
    // Create array type from base type if array
  if (array_size > 1)
  {
    hsize_t tmpsize = array_size;
#if defined(H5Tarray_create_vers) && H5Tarray_create_vers > 1  
    hdf_type = H5Tarray_create2( hdf_type, 1, &tmpsize );
#else
    hdf_type = H5Tarray_create( hdf_type, 1, &tmpsize, NULL );
#endif
    if (hdf_type < 0)
      return MB_FAILURE;
  }  

  
    // If default or global/mesh value in file, read it.
  void *default_ptr = 0, *global_ptr = 0;
  int default_size = 0, global_size = 0;
  if (have_default || have_global)
  {
    if (array_size == -1) { // variable-length tag
      default_size = have_default;
      global_size = have_global;
    }
    else {
      default_size = global_size = array_size;
    }
    
    assert( (default_size + global_size) * elem_size <= bufferSize );
    default_ptr = dataBuffer;
    global_ptr = dataBuffer + default_size*elem_size;

    mhdf_getTagValues( filePtr, name, hdf_type, default_ptr, global_ptr, &status );
    if (mhdf_isError( &status ))
    {
      readUtil->report_error( mhdf_message( &status ) );
      return MB_FAILURE;
    }
    
    if (MB_TYPE_HANDLE == mb_type) {
      if (have_default) {
        rval = convert_id_to_handle( (MBEntityHandle*)default_ptr, default_size );
        if (MB_SUCCESS != rval)
          have_default = 0;
      }
      if (have_global) {
        rval = convert_id_to_handle( (MBEntityHandle*)global_ptr, global_size );
        if (MB_SUCCESS != rval)
          have_global = 0;
      }
    }
  }
  global_size *= elem_size;
  default_size *= elem_size;
  
  
    // Check if tag already exists
  rval = iFace->tag_get_handle( name, handle );
  if (MB_SUCCESS == rval)
  {
    // If tag exists, make sure it is consistant with the type in the file
    int curr_size;
    MBDataType curr_type;
    MBTagType curr_store;
    
    rval = iFace->tag_get_size( handle, curr_size );
    if (MB_VARIABLE_DATA_LENGTH == rval)
      curr_size = -1;
    else if (MB_SUCCESS != rval)
      return rval;
      
    rval = iFace->tag_get_data_type( handle, curr_type );
    if (MB_SUCCESS != rval)
      return rval;
    
    rval = iFace->tag_get_type( handle, curr_store );
    if (MB_SUCCESS != rval)
      return rval;
    
    if ((curr_store != MB_TAG_BIT && curr_size != tag_size) || curr_type != mb_type ||
        ((curr_store == MB_TAG_BIT || storage == MB_TAG_BIT) && 
          curr_store != storage))
    {
      readUtil->report_error( "Tag type in file does not match type in "
                              "database for \"%s\"\n", name );
      return MB_FAILURE;
    }
  }
    // Create the tag if it doesn't exist
  else if (MB_TAG_NOT_FOUND == rval)
  {
    if (tag_size == MB_VARIABLE_LENGTH)
      rval = iFace->tag_create_variable_length( name, storage, mb_type,
                                                handle, default_ptr, default_size );
    else
      rval = iFace->tag_create( name, tag_size, storage, mb_type,
                                handle, default_ptr );
    if (MB_SUCCESS != rval)
      return rval;
  }
    // error
  else
    return rval;
    
  if (have_global) {
    rval = iFace->tag_set_data( handle, 0, 0, &global_ptr, &global_size );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
    // Read tag data
  MBErrorCode tmp = MB_SUCCESS;
  if (have_sparse) {
    if (tag_size == MB_VARIABLE_LENGTH)
      tmp = read_var_len_tag( handle, hdf_type, mb_type == MB_TYPE_HANDLE );
    else 
      tmp = read_sparse_tag( handle, hdf_type, tag_size, mb_type == MB_TYPE_HANDLE );
  }
  rval = read_dense_tag( handle, hdf_type, tag_size, mb_type == MB_TYPE_HANDLE );
  
  
    // If an array type, need to release the type object created above.
  if (array_size > 1)
    H5Tclose( hdf_type );
  
  return MB_SUCCESS == tmp ? rval : tmp;
}
  

MBErrorCode ReadHDF5::read_dense_tag( MBTag tag_handle,
                                      hid_t hdf_read_type,
                                      size_t read_size,
                                      bool is_handle_type )
{
  std::list<ElemSet>::iterator iter;
  const std::list<ElemSet>::iterator end = elemList.end();
  mhdf_Status status;
  std::string name;
  MBErrorCode rval;
  int have;
  
  rval = iFace->tag_get_name( tag_handle, name );
  if (MB_SUCCESS != rval)
    return rval;
  
  have = mhdf_haveDenseTag( filePtr, name.c_str(), nodeSet.type2, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  if (have)
  {
    rval = read_dense_tag( nodeSet, tag_handle, hdf_read_type, read_size, is_handle_type );
    if (!rval)
      return rval;
  }
  
   
  have = mhdf_haveDenseTag( filePtr, name.c_str(), setSet.type2, &status );
  if (mhdf_isError( &status ))
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  if (have)
  {
    rval = read_dense_tag( setSet, tag_handle, hdf_read_type, read_size, is_handle_type );
    if (!rval)
      return rval;
  }
  
  
  for (iter = elemList.begin(); iter != end; ++iter)
  {
    have = mhdf_haveDenseTag( filePtr, name.c_str(), iter->type2, &status );
    if (mhdf_isError( &status ))
    {
      readUtil->report_error( mhdf_message( &status ) );
      return MB_FAILURE;
    }
    
    if (have)
    {
      rval = read_dense_tag( *iter, tag_handle, hdf_read_type, read_size, is_handle_type );
      if (!rval)
        return rval;
    }
  }

  return MB_SUCCESS;
}
  
  
MBErrorCode ReadHDF5::read_dense_tag( ElemSet& set,
                                      MBTag tag_handle,
                                      hid_t hdf_read_type,
                                      size_t read_size, 
                                      bool is_handle_type )
{
  mhdf_Status status;
  std::string name;
  MBErrorCode rval;
  long num_values;
  
  rval = iFace->tag_get_name( tag_handle, name );
  if (MB_SUCCESS != rval)
    return rval;
  
  hid_t data = mhdf_openDenseTagData( filePtr, name.c_str(), 
                                      set.type2, &num_values, &status );
  if (mhdf_isError( &status ) )
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  if ((unsigned long)num_values != set.range.size())
  {
    assert( 0 );
    return MB_FAILURE;
  }
  
  MBRange::const_iterator iter = set.range.begin();
  MBRange subrange;
  
  assert ((hdf_read_type == 0) || (H5Tget_size(hdf_read_type) == read_size));
  size_t chunk_size = bufferSize / read_size;
  size_t remaining = set.range.size();
  size_t offset = 0;
  while (remaining)
  {
    size_t count = remaining > chunk_size ? chunk_size : remaining;
    remaining -= count;
    
    MBRange::const_iterator stop = iter;
    stop += count;
    subrange.clear();
    subrange.merge( iter, stop );
    iter = stop;
    
    mhdf_readDenseTagWithOpt( data, offset, count, hdf_read_type, dataBuffer, 
                              ioProp, &status );
    offset += count;
    if (mhdf_isError( &status ))
    {
      readUtil->report_error( mhdf_message( &status ) );
      mhdf_closeData( filePtr, data, &status );
      return MB_FAILURE;
    }
    
    if (is_handle_type)
    {
      rval = convert_id_to_handle( (MBEntityHandle*)dataBuffer, count * read_size / sizeof(MBEntityHandle) );
      if (MB_SUCCESS != rval)
      {
        mhdf_closeData( filePtr, data, &status );
        return rval;
      }
    }
    
    rval = iFace->tag_set_data( tag_handle, subrange, dataBuffer );
    if (MB_SUCCESS != rval)
    {
      mhdf_closeData( filePtr, data, &status );
      return MB_FAILURE;
    }
  }
  
  mhdf_closeData( filePtr, data, &status );
  if (mhdf_isError( &status ) )
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}


MBErrorCode ReadHDF5::read_sparse_tag( MBTag tag_handle,
                                       hid_t hdf_read_type,
                                       size_t read_size,
                                       bool is_handle_type )
{
  mhdf_Status status;
  std::string name;
  MBErrorCode rval;
  long num_values, data_size;
  hid_t data[3];
  MBTagType mbtype;
  assert ((hdf_read_type == 0) || (H5Tget_size(hdf_read_type) == read_size));
  
  rval = iFace->tag_get_name( tag_handle, name );
  if (MB_SUCCESS != rval)
    return rval;
  
  rval = iFace->tag_get_type( tag_handle, mbtype );
  if (MB_SUCCESS != rval)
    return rval;
  
  mhdf_openSparseTagData( filePtr, name.c_str(), &num_values, &data_size, data, &status );
  if (mhdf_isError( &status ) )
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }
    // fixed-length tag
  assert( num_values == data_size );
  
    // Split buffer into two portions: one for handles and one for data.
    
    // We want to read the same number of handles as data values for each
    // iteration.  Calculate the total number of entries to read in each
    // pass as the size of the buffer over the sum of the size of a handle
    // and value.  Subtract off the size of one value so we reserve space
    // for adjusting for data alignment.
  size_t chunk_size = (bufferSize - read_size) / (sizeof(MBEntityHandle) + read_size);
  
    // Use the first half of the buffer for the handles.
  MBEntityHandle* idbuf = (MBEntityHandle*)dataBuffer;
    // Use the latter portion of the buffer for data
  char* databuf = dataBuffer + (chunk_size * sizeof(MBEntityHandle));
    // To be safe, align tag data to the size of an entire tag value
  if ((size_t)databuf % read_size)
    databuf += read_size - ((size_t)databuf % read_size);
      // Make sure the above calculations are correct
  assert( databuf + chunk_size*read_size < dataBuffer + bufferSize );
  
  size_t remaining = (size_t)num_values;
  size_t offset = 0;
  while (remaining)
  {
    size_t count = remaining > chunk_size ? chunk_size : remaining;
    remaining -= count;
    
    mhdf_readSparseTagEntitiesWithOpt( data[0], offset, count, handleType, idbuf,
                                       ioProp, &status );
    if (mhdf_isError( &status ))
    {
      readUtil->report_error( mhdf_message( &status ) );
      mhdf_closeData( filePtr, data[0], &status );
      mhdf_closeData( filePtr, data[1], &status );
      return MB_FAILURE;
    }
    
    mhdf_readSparseTagValuesWithOpt( data[1], offset, count, hdf_read_type, 
                                     databuf, ioProp, &status );
    if (mhdf_isError( &status ))
    {
      readUtil->report_error( mhdf_message( &status ) );
      mhdf_closeData( filePtr, data[0], &status );
      mhdf_closeData( filePtr, data[1], &status );
      return MB_FAILURE;
    }
    
    offset += count;
    
    rval = convert_id_to_handle( idbuf, count );
    if (MB_SUCCESS != rval)
    {
      mhdf_closeData( filePtr, data[0], &status );
      mhdf_closeData( filePtr, data[1], &status );
      return rval;
    }
    
    if (is_handle_type)
    {
      rval = convert_id_to_handle( (MBEntityHandle*)databuf, count * read_size / sizeof(MBEntityHandle) );
      if (MB_SUCCESS != rval)
      {
        mhdf_closeData( filePtr, data[0], &status );
        mhdf_closeData( filePtr, data[1], &status );
        return rval;
      }
    }

/*** FIX ME - need to do one at a time for BIT tags!  This is stupid. ***/
    if (mbtype == MB_TAG_BIT)
    {
      rval = MB_SUCCESS;
      for (size_t i = 0; MB_SUCCESS == rval && i < count; ++i)
        rval = iFace->tag_set_data( tag_handle, idbuf + i, 1, databuf + i );
    }
    else
    {
      rval = iFace->tag_set_data( tag_handle, idbuf, count, databuf );
    }
    if (MB_SUCCESS != rval)
    {
      mhdf_closeData( filePtr, data[0], &status );
      mhdf_closeData( filePtr, data[1], &status );
      return rval;
    }
  }
  
  mhdf_closeData( filePtr, data[0], &status );
  if (mhdf_isError( &status ) )
    readUtil->report_error( mhdf_message( &status ) );
  mhdf_closeData( filePtr, data[1], &status );
  if (mhdf_isError( &status ) )
    readUtil->report_error( mhdf_message( &status ) );

  return MB_SUCCESS;
}

#define assert_range( ARRAY, BYTES ) \
  assert( (char*)(ARRAY) >= dataBuffer && ((char*)(ARRAY)) + (BYTES) <= dataBuffer + bufferSize )

MBErrorCode ReadHDF5::read_var_len_tag( MBTag tag_handle,
                                        hid_t hdf_read_type,
                                        bool is_handle_type )
{
  mhdf_Status status;
  std::string name;
  MBErrorCode rval;
  long num_values, num_data;
  hid_t data[3];
  MBTagType mbtype;
    // hdf_read_type is NULL (zero) for opaque tag data.
  long elem_size = hdf_read_type ? H5Tget_size( hdf_read_type ) : 1;
  if (elem_size < 1) // invalid type handle?
    return MB_FAILURE;
  
  rval = iFace->tag_get_name( tag_handle, name );
  if (MB_SUCCESS != rval)
    return rval;
  
  rval = iFace->tag_get_type( tag_handle, mbtype );
  if (MB_SUCCESS != rval)
    return rval;
  
  mhdf_openSparseTagData( filePtr, name.c_str(), &num_values, &num_data, data, &status );
  if (mhdf_isError( &status ) )
  {
    readUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }
    
    // Subdivide buffer into 5 chunks:
    // Be careful of order so alignment is valid
    // 1) pointer array (input for MOAB)
    // 2) offset array (read from file)
    // 3) entity handles
    // 4) tag sizes (calculated from file data)
    // 5) tag data
  const size_t avg_data_size = (num_data * elem_size) / num_values;
  const size_t per_ent_size = sizeof(void*) + sizeof(long) + sizeof(MBEntityHandle) + sizeof(int);
  long num_ent = bufferSize / (per_ent_size + 2*avg_data_size);
  if (num_ent == 0) 
    num_ent = bufferSize / 4 / sizeof(void*);
  const size_t data_buffer_size = (bufferSize - num_ent * per_ent_size) / elem_size;
  const void** const pointer_buffer = reinterpret_cast<const void**>(dataBuffer);
  long* const end_idx_buffer = reinterpret_cast<long*>(pointer_buffer + num_ent);
  char* handle_buffer_start = reinterpret_cast<char*>(end_idx_buffer + num_ent);
  if (((size_t)handle_buffer_start) % sizeof(MBEntityHandle))
    handle_buffer_start += sizeof(MBEntityHandle) - ((size_t)handle_buffer_start);
  MBEntityHandle* const handle_buffer = reinterpret_cast<MBEntityHandle*>(handle_buffer_start);
  int* const size_buffer = reinterpret_cast<int*>(handle_buffer + num_ent);
  char* const data_buffer = reinterpret_cast<char*>(size_buffer + num_ent);
  
  
    // do num_ent blocks of entities
  long remaining = num_values;
  long offset = 0;
  long data_offset = 0;
  while (remaining) {
    const long count = remaining < num_ent ? remaining : num_ent;
    remaining -= count;

      // read entity IDs
    assert_range( handle_buffer, count * sizeof(MBEntityHandle) );
    mhdf_readSparseTagEntitiesWithOpt( data[0], offset, count, handleType, 
                                       handle_buffer, ioProp, &status );
    if (mhdf_isError( &status ))
    {
      readUtil->report_error( mhdf_message( &status ) );
      mhdf_closeData( filePtr, data[0], &status );
      mhdf_closeData( filePtr, data[1], &status );
      mhdf_closeData( filePtr, data[2], &status );
      return MB_FAILURE;
    }
      // convert entity ID to MBEntityHandle
    rval = convert_id_to_handle( handle_buffer, count );
    if (MB_SUCCESS != rval)
    {
      mhdf_closeData( filePtr, data[0], &status );
      mhdf_closeData( filePtr, data[1], &status );
      mhdf_closeData( filePtr, data[2], &status );
      return rval;
    }
      // read end index of tag value
    assert_range( end_idx_buffer, count * sizeof(long) );
    mhdf_readSparseTagIndicesWithOpt( data[2], offset, count, H5T_NATIVE_LONG, 
                                      end_idx_buffer, ioProp, &status );
    if (mhdf_isError( &status ))
    {
      readUtil->report_error( mhdf_message( &status ) );
      mhdf_closeData( filePtr, data[0], &status );
      mhdf_closeData( filePtr, data[1], &status );
      mhdf_closeData( filePtr, data[2], &status );
      return MB_FAILURE;
    }
      // update entity offset
    offset += count;
    
    char* pointer = data_buffer;
    long finished = 0;
    while (finished < count) {
      
        // iterate over entities until we have enough to fill the data buffer
      long i = finished;
      long prev_end_idx = data_offset - 1;
      for (; i < count && end_idx_buffer[i] - data_offset < (long)data_buffer_size; ++i) {
          // calculate tag size for entity
        const int size = (end_idx_buffer[i] - prev_end_idx);
        if (size < 0 || (long)size != (end_idx_buffer[i] - prev_end_idx)) {
          readUtil->report_error( "Invalid end index in variable length tag data for tag: \"%s\"", name.c_str() );
          mhdf_closeData( filePtr, data[0], &status );
          mhdf_closeData( filePtr, data[1], &status );
          mhdf_closeData( filePtr, data[2], &status );
          return MB_FAILURE;
        }
        assert_range( size_buffer + i - finished, sizeof(int) );
        size_buffer[i-finished] = size * elem_size;
        assert_range( pointer_buffer + i - finished, sizeof(void*) );
        pointer_buffer[i-finished] = pointer;
        assert_range( pointer_buffer[i-finished], size_buffer[i - finished] );
        pointer += size_buffer[i-finished];
        prev_end_idx = end_idx_buffer[i];
      }
      
      const size_t num_val = prev_end_idx - data_offset + 1;
      if (num_val) {
          // read data
        assert( num_val <= data_buffer_size );
        assert_range( data_buffer, num_val * elem_size );
        mhdf_readSparseTagValuesWithOpt( data[1], data_offset, num_val, 
                                         hdf_read_type, data_buffer,
                                         ioProp, &status );
        if (mhdf_isError( &status )) {
          readUtil->report_error( mhdf_message( &status ) );
          mhdf_closeData( filePtr, data[0], &status );
          mhdf_closeData( filePtr, data[1], &status );
          mhdf_closeData( filePtr, data[2], &status );
          return MB_FAILURE;
        }
        data_offset += num_val;

        if (is_handle_type) {
          rval = convert_id_to_handle( reinterpret_cast<MBEntityHandle*>(data_buffer), num_val );
          if (MB_SUCCESS != rval) {
            mhdf_closeData( filePtr, data[0], &status );
            mhdf_closeData( filePtr, data[1], &status );
            mhdf_closeData( filePtr, data[2], &status );
            return rval;
          }
        }
          // put data in tags
        rval = iFace->tag_set_data( tag_handle, handle_buffer + finished, i - finished, pointer_buffer, size_buffer );
        if (MB_SUCCESS != rval) {
          mhdf_closeData( filePtr, data[0], &status );
          mhdf_closeData( filePtr, data[1], &status );
          mhdf_closeData( filePtr, data[2], &status );
          return rval;
        }

        finished = i;
      }
      else if (finished < i) { // zero-length tag values???
        finished = i;
      }
        // if tag value to big for buffer
      else {
        const int size = (end_idx_buffer[i] - prev_end_idx);
        std::vector<char*> tmp_buffer( size * elem_size );
        mhdf_readSparseTagValuesWithOpt( data[1], data_offset, size, 
                                         hdf_read_type, &tmp_buffer[0],
                                         ioProp, &status );
        if (mhdf_isError( &status )) {
          readUtil->report_error( mhdf_message( &status ) );
          mhdf_closeData( filePtr, data[0], &status );
          mhdf_closeData( filePtr, data[1], &status );
          mhdf_closeData( filePtr, data[2], &status );
          return MB_FAILURE;
        }
        data_offset += size;
        
        if (is_handle_type) {
          rval = convert_id_to_handle( reinterpret_cast<MBEntityHandle*>(&tmp_buffer[0]), size );
          if (MB_SUCCESS != rval) {
            mhdf_closeData( filePtr, data[0], &status );
            mhdf_closeData( filePtr, data[1], &status );
            mhdf_closeData( filePtr, data[2], &status );
            return rval;
          }
        }
        
        pointer_buffer[0] = &tmp_buffer[0];
        size_buffer[0] = size * elem_size;
        rval = iFace->tag_set_data( tag_handle, handle_buffer + i, 1, pointer_buffer, size_buffer );
        if (MB_SUCCESS != rval) {
          mhdf_closeData( filePtr, data[0], &status );
          mhdf_closeData( filePtr, data[1], &status );
          mhdf_closeData( filePtr, data[2], &status );
          return rval;
        }
        
        prev_end_idx = end_idx_buffer[i];
        ++finished;
      }
    }
  }
  
  mhdf_closeData( filePtr, data[0], &status );
  if (mhdf_isError( &status ) )
    readUtil->report_error( mhdf_message( &status ) );
  mhdf_closeData( filePtr, data[1], &status );
  if (mhdf_isError( &status ) )
    readUtil->report_error( mhdf_message( &status ) );
  mhdf_closeData( filePtr, data[2], &status );
  if (mhdf_isError( &status ) )
    readUtil->report_error( mhdf_message( &status ) );

  return MB_SUCCESS;
}

MBErrorCode ReadHDF5::convert_id_to_handle( const ElemSet& elems,
                                            MBEntityHandle* array,
                                            size_t size )
{
  MBEntityHandle offset = elems.first_id;
  MBEntityHandle last = offset + elems.range.size();
  for (MBEntityHandle *const end = array + size; array != end; ++array)
  {
    if (*array >= last || *array < (MBEntityHandle)offset)
      return MB_FAILURE;
    MBRange:: const_iterator iter = elems.range.begin();
    iter += *array - offset;
    *array = *iter;
  }
  
  return MB_SUCCESS;
}

MBErrorCode ReadHDF5::convert_id_to_handle( MBEntityHandle* array, 
                                            size_t size )
{
  MBEntityHandle offset = 1;
  MBEntityHandle last = 0;
  ElemSet* set = 0;
  std::list<ElemSet>::iterator iter;
  const std::list<ElemSet>::iterator i_end = elemList.end();

  for (MBEntityHandle *const end = array + size; array != end; ++array)
  {
      // special case for ZERO
    if (!*array)
      continue;
    
    if (nodeSet.first_id && (*array < offset || *array >= last))
    {
      offset = nodeSet.first_id;
      last = offset + nodeSet.range.size();
      set = &nodeSet;
    }
    if (setSet.first_id && (*array < offset || *array >= last))
    {
      offset = setSet.first_id;
      last = offset + setSet.range.size();
      set = &setSet;
    }
    iter = elemList.begin();
    while (*array < offset || *array >= last)
    {
      if (iter == i_end)
      {
        return MB_FAILURE;
      }
      
      set = &*iter;
      offset = set->first_id;
      last = offset + set->range.size();
      ++iter;
    }
  
    MBRange:: const_iterator riter = set->range.begin();
    riter += *array - offset;
    *array = *riter;
  }
  
  return MB_SUCCESS;
}

MBErrorCode ReadHDF5::convert_range_to_handle( const MBEntityHandle* array,
                                               size_t num_ranges,
                                               MBRange& range )
{
  const MBEntityHandle *const end = array + 2*num_ranges;
  MBEntityHandle offset = 1;
  MBEntityHandle last = 0;
  ElemSet* set = 0;
  std::list<ElemSet>::iterator iter;
  const std::list<ElemSet>::iterator i_end = elemList.end();
  MBEntityHandle start = *(array++);
  MBEntityHandle count = *(array++);
  
  while (true)
  {
    if (nodeSet.first_id && (start < offset || start >= last))
    {
      offset = nodeSet.first_id;
      last = offset + nodeSet.range.size();
      set = &nodeSet;
    }
    if (setSet.first_id && (start < offset || start >= last))
    {
      offset = setSet.first_id;
      last = offset + setSet.range.size();
      set = &setSet;
    }
    iter = elemList.begin();
    while (start < offset || start >= last)
    {
      if (iter == i_end)
      {
        return MB_FAILURE;
      }
      
      set = &*iter;
      offset = set->first_id;
      last = offset + set->range.size();
      ++iter;
    }
  
    MBEntityHandle s_rem = set->range.size() - (start - offset);
    MBEntityHandle num = count > s_rem ? s_rem : count;
    MBRange::const_iterator riter = set->range.begin();
    riter += (start - offset);
    MBRange::const_iterator rend = riter;
    rend += num;
    assert( riter != rend );
    MBEntityHandle h_start = *riter++;
    MBEntityHandle h_prev = h_start;
    
    while (riter != rend)
    {
      if (h_prev + 1 != *riter)
      {
        range.insert( h_start, h_prev );
        h_start = *riter;
      }
      h_prev = *riter;
      ++riter;
    }
    range.insert( h_start, h_prev );
    
    count -= num;
    start += num;
    if (count == 0)
    {
      if (array == end)
        break;
      
      start = *(array++);
      count = *(array++);
    }
  }
  
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
    return MB_FAILURE;
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

  
    
