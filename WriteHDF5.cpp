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

#include <assert.h>
#include <sys/time.h>
#include <time.h>
#include <H5Tpublic.h>
#include "MBInterface.hpp"
#include "MBInternals.hpp"
#include "WriteHDF5.hpp"
#include "RangeTree.hpp"
#include "mhdf.h"

#undef DEBUG

#ifdef DEBUG
#  define DEBUGOUT(A) fputs( A, stderr )
#  include <stdio.h>
#else
#  define DEBUGOUT(A)
#endif

#ifdef DEBUG
# include <H5Epublic.h>
  extern "C" herr_t hdf_error_handler( void*  )
  {
    H5Eprint( stderr );
    assert( 0 );
  }
#endif

#define WRITE_HDF5_BUFFER_SIZE (40*1024*1024)

const hid_t WriteHDF5::id_type = H5T_NATIVE_INT;


MBWriterIface* WriteHDF5::factory( MBInterface* iface )
  { return new WriteHDF5( iface ); }

WriteHDF5::WriteHDF5( MBInterface* iface )
  : bufferSize( WRITE_HDF5_BUFFER_SIZE ),
    dataBuffer( 0 ),
    iFace( iface ), 
    writeUtil( 0 ), 
    filePtr( 0 ), 
    createdIdTag(false)
{
}

MBErrorCode WriteHDF5::init()
{
  MBErrorCode rval;
  id_t zero_int = 0;

  if (writeUtil) 
    return MB_SUCCESS;
 
#ifdef DEBUG
  H5Eset_auto( &hdf_error_handler, writeUtil );
#endif
 
  register_known_tag_types( iFace );
 
  rval = iFace->query_interface( "MBWriteUtilIface", (void**)&writeUtil );
  if (MB_SUCCESS != rval)
    return rval;
  
  rval = iFace->tag_get_handle( GLOBAL_ID_TAG_NAME, idTag );
  if (MB_TAG_NOT_FOUND == rval)
  {
    rval = iFace->tag_create( GLOBAL_ID_TAG_NAME, sizeof(id_t), 
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
                                   const MBEntityHandle* set_array,
                                   const int num_sets,
                                   std::vector<std::string>& qa_records,
                                   int user_dimension )
{
  MBErrorCode result;
  mhdf_Status rval;
  std::string tagname;
  std::vector<MBTag> tag_list;
  std::vector<MBTag>::iterator t_itor;

  std::list<ExportSet>::iterator ex_itor;
  
  if (MB_SUCCESS != init())
    return MB_FAILURE;

DEBUGOUT("Gathering Mesh\n");
  
    // Gather mesh to export
  exportList.clear();
  if (0 == num_sets)
  {
    result = gather_all_mesh( );
    if (MB_SUCCESS != result)
      return result;
  }
  else
  {
    std::vector<MBEntityHandle> passed_export_list(num_sets);
    memcpy( &passed_export_list[0], set_array, sizeof(MBEntityHandle)*num_sets );
    result = gather_mesh_info( passed_export_list );
    if (MB_SUCCESS != result)
      return result;
    
      // Mark all entities invalid.  Later the ones we are
      // exporting will be marked valid.  This way we can
      // distinguish adjacent entities and such that aren't
      // being exported.
    result = clear_all_id_tags();
    if (MB_SUCCESS != result)
      return result;
  }
  
  if (nodeSet.range.size() == 0)
    return MB_FAILURE;
  
DEBUGOUT("Checking ID space\n");

    // Make sure ID space is sufficient
  MBEntityHandle elem_count = nodeSet.range.size() + setSet.range.size();
  for (ex_itor = exportList.begin(); ex_itor != exportList.end(); ++ex_itor)
    elem_count += ex_itor->range.size();
  MBEntityHandle max_id = (MBEntityHandle)1 << (8*sizeof(id_t)-1);
  if (elem_count > max_id)
  {
    writeUtil->report_error("ID space insufficient for mesh size.\n");
    return MB_FAILURE;
  }

DEBUGOUT( "Creating File\n" );  
  
  const char* type_names[MBMAXTYPE];
  bzero( type_names, MBMAXTYPE * sizeof(char*) );
  for (MBEntityType i = MBEDGE; i < MBENTITYSET; ++i)
    type_names[i] = MBCN::EntityTypeName( i );
 
    // Create the file
  filePtr = mhdf_createFile( filename, overwrite, type_names, MBMAXTYPE, &rval );
  if (!filePtr)
  {
    fprintf(stderr, mhdf_message( &rval ));
    return MB_FAILURE;
  }
  
  result = write_qa( qa_records );
  if (MB_SUCCESS != result)
    return result;
  
  dataBuffer = (char*)malloc( bufferSize );
  if (!dataBuffer)
    goto write_fail;

DEBUGOUT("Writing Nodes.\n");
  
    // Write nodes
  if (write_nodes( user_dimension ) != MB_SUCCESS)
    goto write_fail;

DEBUGOUT("Writing connectivity.\n");
  
    // Write connectivity
  for (ex_itor = exportList.begin(); ex_itor != exportList.end(); ++ex_itor)
  {
    if (ex_itor->type == MBPOLYGON || ex_itor->type == MBPOLYHEDRON)
      result = write_poly( *ex_itor );
    else
      result = write_elems( *ex_itor );
    if (MB_SUCCESS != result)
      goto write_fail;
  }

DEBUGOUT("Writing sets.\n");
  
    // Write meshsets
  if (write_sets() != MB_SUCCESS)
    goto write_fail;

DEBUGOUT("Writing adjacencies.\n");
  
    // Write adjacencies
  // Tim says don't save node adjacencies!
  //if (write_adjacencies( nodeSet ) != MB_SUCCESS)
  //  goto write_fail;
  for (ex_itor = exportList.begin(); ex_itor != exportList.end(); ++ex_itor)
    if (write_adjacencies( *ex_itor ) != MB_SUCCESS)
      goto write_fail;

DEBUGOUT("Writing tags.\n");
  
    // Get list of Tags to write
  result = iFace->tag_get_tags( tag_list );
  if (MB_SUCCESS != result)
    goto write_fail;

    // Write tags
  for (t_itor = tag_list.begin(); t_itor != tag_list.end(); ++t_itor)
  {
      // Don't write global ID tag
    if (*t_itor == idTag)
      continue;
    
      // Don't write tags that have name beginning with "__"
    if (iFace->tag_get_name( *t_itor, tagname ) != MB_SUCCESS)
      goto write_fail;
    if (tagname[0] == '_' && tagname[1] == '_')
      continue;
    
    if (write_tag( *t_itor ) != MB_SUCCESS)
      goto write_fail;
  }

    // Clean up and exit.
  free( dataBuffer );
  dataBuffer = 0;
  mhdf_closeFile( filePtr, &rval );
  filePtr = 0;
  if (mhdf_isError( &rval ))
  {
    writeUtil->report_error("%s", mhdf_message( &rval ));
    return MB_FAILURE;
  }
    
  return MB_SUCCESS;
  
write_fail:
  
  if (dataBuffer)
  {
    free( dataBuffer );
    dataBuffer = 0;
  }
  mhdf_closeFile( filePtr, &rval );
  filePtr = 0;
  return MB_FAILURE;
}

MBErrorCode WriteHDF5::clear_all_id_tags()
{
  id_t cleared_value = -1;
  MBRange range;
  for (MBEntityType type = MBVERTEX; type < MBMAXTYPE; ++type)
  {
    range.clear();
    MBErrorCode rval = iFace->get_entities_by_type( 0, type, range, false );
    if (MB_SUCCESS != rval)
      return rval;
    
    for (MBRange::iterator itor = range.begin(); itor != range.end(); ++itor)
    {
      MBEntityHandle handle = *itor;
      rval = iFace->tag_set_data( idTag, &handle, 1, &cleared_value );
      if (MB_SUCCESS != rval)
        return rval;
    }
  }
  return MB_SUCCESS;
}
   

MBErrorCode WriteHDF5::midnode_combinations( MBEntityType type,
                                    std::vector<int>& combinations )
{
  combinations.clear();

  int dimension = MBCN::Dimension( type );
  int num_faces = 0;
  int num_edges = 0;
  switch (dimension) 
  {
    case 3:
      num_faces = MBCN::NumSubEntities( type, 2 );
      num_edges = MBCN::NumSubEntities( type, 1 );
      combinations.push_back(num_faces + num_edges + 1);
      combinations.push_back(num_faces + num_edges);
      combinations.push_back(num_faces + 1);
      combinations.push_back(num_faces);
    case 2:
      num_edges = MBCN::NumSubEntities( type, 1 );
      combinations.push_back(num_edges + 1);
      combinations.push_back(num_edges);
    case 1:
      combinations.push_back(1);
      combinations.push_back(0);
      break;
    default:
      assert(0);
      return MB_FAILURE;
  }

  return MB_SUCCESS;
}

MBErrorCode WriteHDF5::subrange_by_type_and_conn( const MBRange& input,
                                                  MBEntityType type,
                                                  int nodecount,
                                                  MBRange& output )
{
  MBRange::const_iterator itor;
  const MBEntityHandle* junk;
  MBErrorCode rval;
  int num_nodes;
  output.clear();
  for (itor = input.begin(); itor != input.end(); ++itor)
  {
    rval = iFace->get_connectivity( *itor, junk, num_nodes, false );
    if (MB_SUCCESS != rval)
      return MB_FAILURE;
    if (num_nodes == nodecount && 
        iFace->type_from_handle(*itor) == type)
      output.insert( *itor );
  }
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5::subrange_by_type( const MBRange& input,
                                         MBEntityType type,
                                         MBRange& output )
{
  int err;
  MBEntityHandle first, last;
  first = CREATE_HANDLE(  type, 0, err ); assert(!err);
  last = CREATE_HANDLE( type+1, 0, err ); assert(!err);
  MBRange::const_iterator start, end;
  start = input.lower_bound( input.begin(), input.end(), first );
  end = input.lower_bound( start, input.end(), last );
  
  output.clear();
  output.merge( start, end );
  return MB_SUCCESS;
}

                                         

MBErrorCode WriteHDF5::gather_mesh_info( 
                           const std::vector<MBEntityHandle>& export_sets )
{
  MBErrorCode rval;
  std::list<MBEntityHandle>::iterator export_itor;
  
  nodeSet.range.clear();
  setSet.range.clear();
  exportList.clear();
  
  int dim;
  MBRange range;
  MBRange ranges[4];
  
    // Gather list of all related sets
  std::vector<MBEntityHandle> stack(export_sets);
  std::copy( export_sets.begin(), export_sets.end(), stack.begin() );
  std::vector<MBEntityHandle> set_children;
  while( !stack.empty() )
  {
    MBEntityHandle meshset = stack.back(); stack.pop_back();
    setSet.range.insert( meshset );
  
      // Get contained sets
    range.clear();
    rval = iFace->get_entities_by_type( meshset, MBENTITYSET, range, true );
    if (MB_SUCCESS != rval)
      return MB_FAILURE;
    for (MBRange::iterator ritor = range.begin(); ritor != range.end(); ++ritor)
      if (setSet.range.find( *ritor ) == setSet.range.end())
        stack.push_back( *ritor );
    
      // Get child sets
    set_children.clear();
    rval = iFace->get_child_meshsets( meshset, set_children, 1 );
    if (MB_SUCCESS != rval)
      return MB_FAILURE;
    for (std::vector<MBEntityHandle>::iterator vitor = set_children.begin();
         vitor != set_children.end(); ++vitor )
      if (setSet.range.find( *vitor ) == setSet.range.end())
        stack.push_back( *vitor );
  }
  
    // Gather list of all mesh entities from list of sets,
    // grouped by dimension.
  for (MBRange::iterator setitor = setSet.range.begin();
       setitor != setSet.range.end(); ++setitor)
  {
    for (dim = 0; dim < 4; ++dim)
    {
      range.clear();
      rval = iFace->get_entities_by_dimension( *setitor, dim, range, false );
      if (MB_SUCCESS != rval)
        return MB_FAILURE;

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
      if (MB_SUCCESS != rval)
        return MB_FAILURE;
      ranges[cdim].merge( range );
    }  
    range.clear();
    rval = writeUtil->gather_nodes_from_elements( ranges[dim], 0, range );
    if (MB_SUCCESS != rval)
      return MB_FAILURE;
    ranges[0].merge( range );      
  }
  
    // Split lists by element type and number of nodes
  nodeSet.range = ranges[0];
  nodeSet.type = MBVERTEX;
  nodeSet.num_nodes = 1;
  setSet.type = MBENTITYSET;
  setSet.num_nodes = 0;
  
  std::vector<int> node_counts;
  for (MBEntityType type = MBEDGE; type < MBENTITYSET; ++type)
  {
    ExportSet set;
    dim = MBCN::Dimension(type);

    if (ranges[dim].empty())
      continue;
    
    if (type == MBPOLYGON || type == MBPOLYHEDRON)
    {
      rval = subrange_by_type( ranges[dim], type, set.range );
      if (MB_SUCCESS != rval)
        return MB_FAILURE;
      
      if (!set.range.empty())
      {
        set.type = type;
        set.num_nodes = 0;
        exportList.push_back( set );
      }
    }
    else
    {      
      node_counts.clear();
      rval = midnode_combinations( type, node_counts );
      if (MB_SUCCESS != rval)
        return MB_FAILURE;

      const int num_corners = MBCN::VerticesPerEntity( type );
      for (std::vector<int>::iterator itor = node_counts.begin();
           itor != node_counts.end(); ++itor)
      {
        set.range.clear();
        int num_nodes = *itor + num_corners;
        rval = subrange_by_type_and_conn( ranges[dim], type, num_nodes, set.range );
        if (MB_SUCCESS != rval)
          return MB_FAILURE;

        if (!set.range.empty())
        {
          set.type = type;
          set.num_nodes = num_nodes;
          exportList.push_back( set );
        }
      }
    }
  }
    
  return MB_SUCCESS;  
}


MBErrorCode WriteHDF5::gather_all_mesh( )
{
  MBErrorCode rval;
  
    // Get all nodes
  nodeSet.range.clear();
  rval = iFace->get_entities_by_type( 0, MBVERTEX, nodeSet.range, true );
  if (MB_SUCCESS != rval)
    return MB_FAILURE;
  nodeSet.type = MBVERTEX;
  nodeSet.num_nodes = 1;
  
    // Get all sets
  setSet.range.clear();
  rval = iFace->get_entities_by_type( 0, MBENTITYSET, setSet.range, true );
  if (MB_SUCCESS != rval)
    return MB_FAILURE;
  setSet.type = MBENTITYSET;
  setSet.num_nodes = 0;
  
    // MOAB always returns the default set / global mesh
    // Don't want it.
  int size;
  MBRange::iterator gs = setSet.range.begin();
  assert( gs != setSet.range.end() &&
          MB_SUCCESS == iFace->get_number_entities_by_handle( *gs, size ) && 
          0 == size &&
          MB_SUCCESS == iFace->num_child_meshsets( *gs, &size ) && 
          0 == size &&
          MB_SUCCESS == iFace->num_parent_meshsets( *gs, &size ) && 
          0 == size );
  setSet.range.erase( gs );

    // Get all elements, grouped by type and number of higher-order nodes
  exportList.clear();
  std::vector<int> node_counts;
  for (MBEntityType type = MBEDGE; type < MBENTITYSET; ++type)
  {
    ExportSet set;
    
    MBRange range;
    rval = iFace->get_entities_by_type( 0, type, range, true );
    if (MB_SUCCESS != rval)
      return MB_FAILURE;
      
    if (range.empty())
      continue;
    
    if (type == MBPOLYGON || type == MBPOLYHEDRON)
    {
      set.range = range;
      set.type = type;
      set.num_nodes = 0;
      exportList.push_back( set );
      continue;
    }
    
    rval = midnode_combinations( type, node_counts );
    if (MB_SUCCESS != rval)
      return MB_FAILURE;
      
    const int num_corners = MBCN::VerticesPerEntity( type );
    for (std::vector<int>::iterator itor = node_counts.begin();
         itor != node_counts.end(); ++itor)
    {
      set.range.clear();
      int num_nodes = *itor + num_corners;
      rval = subrange_by_type_and_conn( range, type, num_nodes, set.range );
      if (MB_SUCCESS != rval)
        return MB_FAILURE;
      
      if (!set.range.empty())
      {
        set.type = type;
        set.num_nodes = num_nodes;
        exportList.push_back( set );
      }
    }
  }
  
  return MB_SUCCESS;
}
  
MBErrorCode WriteHDF5::write_nodes( int req_dim )
{
  mhdf_Status status;
  int dim, mesh_dim;
  MBErrorCode rval;
  hid_t node_table;
  long first_id;
  nodeSet.type2 = mhdf_node_type_handle();
  
  rval = iFace->get_dimension( mesh_dim );
  if (MB_SUCCESS != rval)
    return rval;
  
  if (req_dim < 1) 
    req_dim = mesh_dim;
  dim = req_dim > mesh_dim ? mesh_dim : req_dim;
  
  node_table = mhdf_createNodeCoords( filePtr, req_dim, nodeSet.range.size(), &first_id, &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
  writeUtil->assign_ids( nodeSet.range, idTag, (id_t)first_id );
  
  double* buffer = (double*)dataBuffer;
  int chunk_size = bufferSize / sizeof(double);
  
  long remaining = nodeSet.range.size();
  long offset = 0;
  MBRange::const_iterator iter = nodeSet.range.begin();
  while (remaining)
  {
    long count = chunk_size < remaining ? chunk_size : remaining;
    remaining -= count;
    MBRange::const_iterator end = iter;
    end += count;
    
    for (int d = 0; d < req_dim; d++)
    {
      if (d < dim)
      {
        rval = writeUtil->get_node_array( d, iter, end, count, buffer );
        if (MB_SUCCESS != rval)
        {
          mhdf_closeData( filePtr, node_table, &status );
          return rval;
        }
      }
      else
      {
        bzero( buffer, count * sizeof(double) );
      }
    
      mhdf_writeNodeCoord( node_table, offset, count, d, buffer, &status );
      if (mhdf_isError( &status ))
      {
        writeUtil->report_error( "%s\n", mhdf_message( &status ) );
        mhdf_closeData( filePtr, node_table, &status );
        return MB_FAILURE;
      }
    }
    
    iter = end;
    offset += count;
  }
  
  mhdf_closeData( filePtr, node_table, &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
 
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5::write_elems( ExportSet& elems )
{
  mhdf_Status status;
  MBErrorCode rval;
  char name[64];
  long first_id;
  
  sprintf( name, "%s%d", MBCN::EntityTypeName(elems.type), elems.num_nodes );
  elems.type2 = mhdf_addElement( filePtr, name, elems.type, &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }

  hid_t elem_table = mhdf_createConnectivity( filePtr,
                                              elems.type2,
                                              elems.num_nodes,
                                              elems.range.size(),
                                              &first_id,
                                              &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
  writeUtil->assign_ids( elems.range, idTag, (id_t)first_id );
  
  
  id_t* buffer = (id_t*)dataBuffer;
  int chunk_size = bufferSize / (elems.num_nodes * sizeof(id_t));
  long offset = 0;
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
    iter = next;
    if (MB_SUCCESS != rval)
    {
      mhdf_closeData( filePtr, elem_table, &status );
      return rval;
    }    
    
    mhdf_writeConnectivity( elem_table, offset, count, 
                            id_type, buffer, &status );
    if (mhdf_isError( &status ))
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      mhdf_closeData( filePtr, elem_table, &status );
      return MB_FAILURE;
    }
    
    offset += count;
  }

  mhdf_closeData( filePtr, elem_table, &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
 
  return MB_SUCCESS;
}


MBErrorCode WriteHDF5::write_poly( ExportSet& elems )
{
  mhdf_Status status;
  MBErrorCode rval;
  char name[64];
  long first_id;
  
  assert( elems.type == MBPOLYGON || elems.type == MBPOLYHEDRON );
  
  // Create the element group in the file
  sprintf( name, "%s", MBCN::EntityTypeName(elems.type) );
  elems.type2 = mhdf_addElement( filePtr, name, elems.type, &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }

  // Determine how much data must be written
  int mb_count;
  rval = writeUtil->get_poly_array_size( elems.range.begin(), 
                                         elems.range.end(),
                                         mb_count );
  if (MB_SUCCESS != rval)
    return rval;
  
  long count = mb_count;
  const MBRange::const_iterator end = elems.range.end();
  MBRange::const_iterator iter;

    // Create the tables in the file and assign IDs to polys
  hid_t handles[2];
  mhdf_createPolyConnectivity( filePtr, elems.type2, elems.range.size(),
                               count, &first_id, handles, &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
  writeUtil->assign_ids( elems.range, idTag, (id_t)first_id );
  
    // Split the data buffer into two chunks, one for the indices
    // and one for the IDs.  Assume average of 4 IDs per poly.
  
  size_t chunk_size = bufferSize / (5*sizeof(id_t));
  id_t* idx_buffer = (id_t*)dataBuffer;
  id_t* conn_buffer = idx_buffer + chunk_size;
  
  long offset[2] = {0,0};
  id_t index_offset = 0;
  iter = elems.range.begin();
  
  while (iter != end)
  {
    size_t num_idx = chunk_size;
    size_t num_conn = 4*chunk_size;
    rval = writeUtil->get_poly_arrays( iter, end, idTag,
                                       num_conn, conn_buffer,
                                       num_idx,  idx_buffer,
                                       index_offset );
    if (MB_SUCCESS != rval)
    {
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      return rval;
    }
    
    mhdf_writePolyConnIndices( handles[0], offset[0], num_idx, id_type, idx_buffer, &status );
    if (mhdf_isError( &status ))
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      return MB_FAILURE;
    }
    offset[0] += num_idx;
    
    mhdf_writePolyConnIDs( handles[1], offset[1], num_conn, id_type, conn_buffer, &status );
    if (mhdf_isError( &status ))
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      mhdf_closeData( filePtr, handles[0], &status );
      mhdf_closeData( filePtr, handles[1], &status );
      return MB_FAILURE;
    }
    offset[1] += num_conn;
  }
  
  assert((unsigned)offset[0] == elems.range.size());
  assert(offset[1] == mb_count);

  mhdf_closeData( filePtr, handles[0], &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }

  mhdf_closeData( filePtr, handles[1], &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
 
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5::get_set_info( MBEntityHandle set,
                                     long& num_entities,
                                     long& num_children,
                                     unsigned long& flags )
{
  MBErrorCode rval;
  int i;
  unsigned int u;
  
  rval = iFace->get_number_entities_by_handle( set, i, false );
  if (MB_SUCCESS != rval)
    return rval;
  num_entities = i;

  rval = iFace->num_child_meshsets( set, &i );
  if (MB_SUCCESS != rval)
    return rval;
  num_children = i;

  rval = iFace->get_meshset_options( set, u );
  if (MB_SUCCESS != rval)
    return rval;
  flags = u;
  
  return MB_SUCCESS;
}


MBErrorCode WriteHDF5::write_sets( )
{
  mhdf_Status status;
  MBRange& sets = setSet.range;
  MBErrorCode rval;
  setSet.type2 = mhdf_set_type_handle();
  long first_id;
  
  /* If no sets, just return success */
  if (sets.empty())
    return MB_SUCCESS;
  
  /* Write set description table */
  
  /* Create the table */
  hid_t table = mhdf_createSetMeta( filePtr, sets.size(), &first_id, &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
  writeUtil->assign_ids( sets, idTag, (id_t)first_id );

  /* Set up I/O buffer */
  int chunk_size = bufferSize / (3*sizeof(long));
  long* meta_buf = (long*)dataBuffer;
  
  /* Write set descriptions */
  MBRange::const_iterator itor = sets.begin();
  std::vector<MBEntityHandle> handle_list;
  MBRange set_contents;
  long offset = 0, count = 0, child_count = 0, data_count = 0, remaining = sets.size();
  while (remaining)
  {
    count = chunk_size < remaining ? chunk_size : remaining;
    remaining -= count;
    long* buf_iter = meta_buf;
    unsigned long flags;
    for (int i = 0; i < count; i++, itor++, buf_iter += 3)
    {
      rval = get_set_info( *itor, buf_iter[0], buf_iter[1], flags );
      if (MB_SUCCESS != rval)
      {
        mhdf_closeData( filePtr, table, &status );
        return rval;
      }
      child_count += buf_iter[1];

      /* Check if set can be written as ranges of ids */
      buf_iter[2] = flags;
      if ((flags&MESHSET_SET) && !(flags&MESHSET_ORDERED))
      {
        set_contents.clear();
        rval = iFace->get_entities_by_handle( *itor, set_contents, false );
        if (MB_SUCCESS != rval)
        {
          mhdf_closeData( filePtr, table, &status );
          return rval;
        }
        
        int length;
        rval = range_to_id_list( set_contents, NULL, length );
        if (MB_SUCCESS != rval)
        {
          mhdf_closeData( filePtr, table, &status );
          return rval;
        }
        
        if (length < buf_iter[0])
        {
          buf_iter[0] = length;
          buf_iter[2] |= mhdf_SET_RANGE_BIT;
        }
        data_count += length;
      }
      else
      {
        data_count += buf_iter[0];
      }
    }
    
    mhdf_writeSetMeta( table, offset, count, H5T_NATIVE_LONG, meta_buf, &status );
    if (mhdf_isError( &status ))
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      mhdf_closeData( filePtr, table, &status );
      return MB_FAILURE;
    }
    offset += count;
  }
  
  mhdf_closeData( filePtr, table, &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  /* Write set contents */
  
  if (data_count > 0)
  {

    /* Create the table */  
    table = mhdf_createSetData( filePtr, data_count, &status );
    if (mhdf_isError( &status ))
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      return MB_FAILURE;
    }

    /* Write contents of each set */
    offset = 0;
    std::vector<id_t> list;
    MBRange range;
    for (itor = sets.begin(); itor != sets.end(); ++itor)
    {
      long junk1, junk2;
      id_t junk3;
      unsigned long flags;
      rval = get_set_info( *itor, junk1, junk2, flags );
      if (MB_SUCCESS != rval)
      {
        mhdf_closeData( filePtr, table, &status );
        return rval;
      }

      /* Ranges of IDs or just list of IDs? */
      if ((flags&MESHSET_SET) && !(flags&MESHSET_ORDERED))
      {
        list.clear();
        range.clear();

        rval = iFace->get_entities_by_handle( *itor, range, false );
        if (MB_SUCCESS != rval)
        {
          mhdf_closeData( filePtr, table, &status );
          return rval;
        }

        if (range.size() == 0)
          continue;

        rval = range_to_id_list( range, &list, junk3 );
        if (MB_SUCCESS != rval)
        {
          mhdf_closeData( filePtr, table, &status );
          return rval;
        }
      }
      else
      {
        list.clear();
        rval = iFace->get_entities_by_handle( *itor, handle_list, false );
        if (MB_SUCCESS != rval)
        {
          mhdf_closeData( filePtr, table, &status );
          return rval;
        }

        if (handle_list.size() == 0)
          continue;

        rval = vector_to_id_list( handle_list, list );
        if (MB_SUCCESS != rval)
        {
          mhdf_closeData( filePtr, table, &status );
          return rval;
        }
      }

      mhdf_writeSetData( table, offset, list.size(), id_type, 
                         &list[0], &status );
      if (mhdf_isError( &status ))
      {
        writeUtil->report_error( "%s\n", mhdf_message( &status ) );
        mhdf_closeData( filePtr, table, &status );
        return MB_FAILURE;
      }

      offset += list.size();
    }  

    mhdf_closeData( filePtr, table, &status );
    if (mhdf_isError( &status ))
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      return MB_FAILURE;
    }
  }
  
   
  /* Write set children */
  offset = 0;
  if (child_count > 0)
  {  
    std::vector<id_t> list;
    table = mhdf_createSetChildren( filePtr, child_count, &status );
    if (mhdf_isError( &status ))
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      return MB_FAILURE;
    }

    count = 0;
    for (itor = sets.begin(); itor != sets.end(); ++itor)
    {
      handle_list.clear();
      list.clear();
      rval = iFace->get_child_meshsets( *itor, handle_list, 1 );
      if (MB_SUCCESS != rval)
      {
        mhdf_closeData( filePtr, table, &status );
        return rval;
      }

      if (handle_list.size() == 0)
        continue;

      rval = vector_to_id_list( handle_list, list );
      if (MB_SUCCESS != rval)
      {
        mhdf_closeData( filePtr, table, &status );
        return rval;
      }


      mhdf_writeSetChildren( table, offset, list.size(), id_type, 
                             &list[0], &status );
      if (mhdf_isError( &status ))
      {
        writeUtil->report_error( "%s\n", mhdf_message( &status ) );
        mhdf_closeData( filePtr, table, &status );
        return MB_FAILURE;
      }
      offset += list.size();
    }

    mhdf_closeData( filePtr, table, &status );
    if (mhdf_isError( &status ))
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      return MB_FAILURE;
    }
  }
   
  return MB_SUCCESS;
}


MBErrorCode WriteHDF5::range_to_id_list( const MBRange& input_range,
                                         std::vector<id_t>* output_id_list,
                                         id_t& output_length )
{
  RangeTree<id_t> idtree;
  std::vector<id_t>::iterator v_iter;
  MBRange::const_iterator r_iter;
  MBRange::const_iterator const r_end = input_range.end();
  id_t linear_size = input_range.size();
  id_t ranged_size;
  id_t id;
  MBErrorCode rval;
  MBEntityHandle handle;
  
  for (r_iter = input_range.begin(); r_iter != r_end; ++r_iter)
  {
    handle = *r_iter;
    rval = iFace->tag_get_data( idTag, &handle, 1, &id );
    if (MB_SUCCESS != rval)
      return rval;
    idtree.insert( id );
  }
  ranged_size = idtree.num_ranges() * 2;
  
  if (ranged_size < linear_size)
  {
    output_length = ranged_size;
    if (!output_id_list)
      return MB_SUCCESS;
      
    RangeTree<id_t>::span_iterator s_iter;
    RangeTree<id_t>::span_iterator const s_end = idtree.span_end();
    output_id_list->resize( ranged_size );
    v_iter = output_id_list->begin();
    for (s_iter = idtree.span_begin(); s_iter != s_end; ++s_iter)
    {
      *v_iter = (*s_iter).first;
      ++v_iter;
      *v_iter = (*s_iter).second - (*s_iter).first + 1;
      ++v_iter;
    }
  }
  else
  {
    output_length = linear_size;
    if (!output_id_list)
      return MB_SUCCESS;
      
    output_id_list->resize( linear_size );
    v_iter = output_id_list->begin();
    r_iter = input_range.begin();
    while (r_iter != r_end)
    {
      handle = *r_iter;
      rval = iFace->tag_get_data( idTag, &handle, 1, &id );
      if (MB_SUCCESS != rval)
        return rval;
      *v_iter = id;
      ++v_iter;
      ++r_iter;
    }
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
    if (MB_SUCCESS != rval)
      return rval;
    *o_iter = id;
    ++o_iter;
    ++i_iter;
  }

  return MB_SUCCESS;
}


template <class T> static void 
erase_from_vector( std::vector<T> vector, T value )
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
  erase_from_vector( adj, (id_t)0 );
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
  for (iter = elements.range.begin(); iter != end; ++iter)
  {
    adj_list.clear();
    rval = get_adjacencies( *iter, adj_list);
    if (MB_SUCCESS != rval)
      return rval;

    if (adj_list.size() > 0)
      count += adj_list.size() + 2;
  }
  
  if (count == 0)
    return MB_SUCCESS;
  
  /* Create data list */
  hid_t table = mhdf_createAdjacency( filePtr, elements.type2, count, &status );
  if (mhdf_isError(&status))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  /* Write data */
  long offset = 0;
  id_t* buffer = (id_t*)dataBuffer;
  long chunk_size = bufferSize / sizeof(id_t); 
  count = 0;
  for (iter = elements.range.begin(); iter != end; ++iter)
  {
    adj_list.clear();
    rval = get_adjacencies( *iter, adj_list );
    if(MB_SUCCESS != rval)
    {
      mhdf_closeData( filePtr, table, &status );
      return rval;
    }
    if (adj_list.size() == 0)
      continue;
    
    if (count + adj_list.size() + 2 > (unsigned long)chunk_size)
    {
      mhdf_writeAdjacency( table, offset, count, id_type, buffer, &status );
      if (mhdf_isError(&status))
      {
        writeUtil->report_error( "%s\n", mhdf_message( &status ) );
        mhdf_closeData( filePtr, table, &status );
        return MB_FAILURE;
      }
      
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
    if (mhdf_isError(&status))
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      mhdf_closeData( filePtr, table, &status );
      return MB_FAILURE;
    }

    offset += count;
    count = 0;
  }
  
  mhdf_closeData( filePtr, table, &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}


MBErrorCode WriteHDF5::write_tag( MBTag tag_handle )
{
  MBErrorCode rval;
  MBTagType tag_type;
  MBTag type_handle;
  int tag_size, mem_size;
  mhdf_Status status;
  hid_t hdf_tag_type;
  
  rval = iFace->tag_get_type( tag_handle, tag_type );
  if (MB_SUCCESS != rval) 
    return rval;
 
  rval = iFace->tag_get_size( tag_handle, tag_size );
  if (MB_SUCCESS != rval)
    return rval;
  
  bool sparse = true;
  bool have_type = false;
  std::string tag_type_name = "__hdf5_tag_type_";
  std::string tag_name;
  rval = iFace->tag_get_name( tag_handle, tag_name );
  if (MB_SUCCESS != rval)
    return rval;
  
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

  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  if (sparse)
    rval = write_sparse_tag( tag_handle, hdf_tag_type );
  else
    rval = write_dense_tag( tag_handle, hdf_tag_type );
    
  return rval;
}


MBErrorCode WriteHDF5::write_dense_tag( MBTag handle,
                                        hid_t type )
{
  MBErrorCode rval = MB_SUCCESS;
  
  if (!nodeSet.range.empty())
    rval = write_dense_tag( nodeSet, handle, type );
  if (MB_SUCCESS != rval)
    return rval;
  
  std::list<ExportSet>::iterator iter, end = exportList.end();
  for (iter = exportList.begin(); iter != end; ++iter)
  {
    MBErrorCode rval = write_dense_tag( *iter, handle, type );
    if (MB_SUCCESS != rval)
      return rval;
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
  if (mhdf_isError(&status))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
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
    if (MB_SUCCESS != rval)
    {
      mhdf_closeData( filePtr, data_handle, &status );
      return MB_FAILURE;
    }
    
    mhdf_writeDenseTag( data_handle, offset, count, type, dataBuffer, &status );
    if (mhdf_isError( &status ))
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      mhdf_closeData( filePtr, data_handle, &status );
      return MB_FAILURE;
    }
    
    offset += count;
  }
  
  mhdf_closeData( filePtr, data_handle, &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5::write_sparse_tag( MBTag handle,
                                         hid_t type )
{
  MBErrorCode rval;
  mhdf_Status status;
  hid_t tables[2];
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
  
    // Get one big range of all entities that have tag.
  MBRange range, e_range;
  
  e_range = nodeSet.range;
  rval = iFace->get_entities_by_type_and_tag( 0, nodeSet.type,
    &handle, NULL, 1, e_range, MBInterface::INTERSECT, false );
  if (MB_SUCCESS != rval)
    return rval;
  range.merge( e_range );
  
  std::list<ExportSet>::iterator e_iter, e_end = exportList.end();
  for (e_iter = exportList.begin(); e_iter != e_end; ++e_iter)
  {
    e_range = e_iter->range;
    rval = iFace->get_entities_by_type_and_tag( 0, e_iter->type, 
      &handle, NULL, 1, e_range, MBInterface::INTERSECT, false );
    if (MB_SUCCESS != rval)
      return rval;
    range.merge( e_range );
  }
  
  e_range = setSet.range;
  rval = iFace->get_entities_by_type_and_tag( 0, setSet.type,
    &handle, NULL, 1, e_range, MBInterface::INTERSECT, false );
  if (MB_SUCCESS != rval)
    return rval;
  range.merge( e_range );
  
  if (range.empty())
    return MB_SUCCESS;
  
  
  
    // Create data table to write to
  mhdf_createSparseTagData( filePtr, name.c_str(), range.size(), tables, &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }

    // Set up data buffer for writing IDs
  size_t chunk_size = bufferSize / sizeof(id_t);
  id_t* id_buffer = (id_t*)dataBuffer;
  
    // Write IDs of tagged entities.
  long remaining = range.size();
  long offset = 0;
  MBRange::iterator iter = range.begin();
  while (remaining)
  {
      // write "chunk_size" blocks of data
    long count = (unsigned long)remaining > chunk_size ? chunk_size : remaining;
    remaining -= count;
    MBRange::iterator stop = iter;
    stop += count;
    e_range.clear();
    e_range.merge( iter, stop );
    iter = stop;
    assert(e_range.size() == (unsigned)count);
    
    rval = iFace->tag_get_data( idTag, e_range, id_buffer );
    if (MB_SUCCESS != rval)
    {
      mhdf_closeData( filePtr, tables[0], &status );
      mhdf_closeData( filePtr, tables[1], &status );
      return rval;
    }
    
      // write the data
    mhdf_writeSparseTagEntities( tables[0], offset, count, id_type, 
                                 id_buffer, &status );
    if (mhdf_isError( &status ))
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      mhdf_closeData( filePtr, tables[0], &status );
      mhdf_closeData( filePtr, tables[1], &status );
      return MB_FAILURE;
    }
   
    offset += count;
  } // while (remaining)
  mhdf_closeData( filePtr, tables[0], &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
    // Set up data buffer for writing tag values
  chunk_size = bufferSize / tag_size;
  assert( chunk_size > 0 );
  char* tag_buffer = (char*)dataBuffer;
  
    // Write the tag values
  remaining = range.size();
  offset = 0;
  iter = range.begin();
  while (remaining)
  {
      // write "chunk_size" blocks of data
    long count = (unsigned long)remaining > chunk_size ? chunk_size : remaining;
    remaining -= count;
    bzero( tag_buffer, count * tag_size );
    MBRange::iterator stop = iter;
    stop += count;
    e_range.clear();
    e_range.merge( iter, stop );
    iter = stop;
    assert(e_range.size() == (unsigned)count);
 
 /** Fix me - stupid API requires we get these one at a time for BIT tags */
    if (mb_type == MB_TAG_BIT)
    {
      rval = MB_SUCCESS;
      char* buf_iter = tag_buffer;
      for (MBRange::iterator it = e_range.begin(); 
           MB_SUCCESS == rval && it != e_range.end(); 
           ++it, buf_iter += tag_size)
        rval = iFace->tag_get_data( handle, &*it, 1, buf_iter );
    }
    else
    {
      rval = iFace->tag_get_data( handle, e_range, tag_buffer );
    }
    if (MB_SUCCESS != rval)
    {
      mhdf_closeData( filePtr, tables[1], &status );
      return rval;
    }
    
      // write the data
    mhdf_writeSparseTagValues( tables[1], offset, count,
                               type, tag_buffer, &status );
    if (mhdf_isError( &status ))
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      mhdf_closeData( filePtr, tables[1], &status );
      return MB_FAILURE;
    }
   
    offset += count;
  } // while (remaining)
  
  mhdf_closeData( filePtr, tables[1], &status );
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
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
  if (mhdf_isError( &status ))
  {
    writeUtil->report_error( mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}


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
