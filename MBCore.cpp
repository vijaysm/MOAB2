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

#ifdef WIN32
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "MBCore.hpp"
#include "TagServer.hpp"
#include "MBMeshSet.hpp"
#include "assert.h"
#include "AEntityFactory.hpp"
#include "MBReadUtil.hpp"
#include "MBWriteUtil.hpp"
#include "MBCN.hpp"
#include "HigherOrderFactory.hpp"
#include "EntitySequenceManager.hpp"
#include "MBError.hpp"
#include "MBReaderWriterSet.hpp"
#include "MBReaderIface.hpp"
#include "MBWriterIface.hpp"
#include "WriteNCDF.hpp"
#include "MBTagConventions.hpp"
#ifdef LINUX
# include <dlfcn.h>
# include <dirent.h>
#endif


#ifdef XPCOM_MB
#include "nsMemory.h"
#endif

using namespace std;

//! Constructor
MBCore::MBCore() 
{
#ifdef XPCOM_MB
  NS_INIT_ISUPPORTS();
#endif

  if (initialize() != MB_SUCCESS)
  {
    printf("Error initializing MB\n");
    exit(1);
  }

}

//! destructor
MBCore::~MBCore()
{
  if(mMBWriteUtil)
    delete mMBWriteUtil;
  if(mMBReadUtil) 
    delete mMBReadUtil;
  
  mMBWriteUtil = NULL;
  mMBReadUtil = NULL;

  deinitialize();
}


MBErrorCode MBCore::initialize()
{
  geometricDimension = 3;
  materialTag      = 0;
  neumannBCTag     = 0;
  dirichletBCTag   = 0;
  geomDimensionTag = 0;
  globalIdTag      = 0;

    //Initialize MBset Members
  cachedMsPtr = NULL;
  cachedEntityHandle = 0;
  maxMeshSetid = 1;

  tagServer = new TagServer();
  if (!tagServer)
    return MB_MEMORY_ALLOCATION_FAILED;
  
  sequenceManager = new EntitySequenceManager();
  if (!sequenceManager)
    return MB_MEMORY_ALLOCATION_FAILED;

  aEntityFactory = new AEntityFactory(this);
  if (!aEntityFactory)
    return MB_MEMORY_ALLOCATION_FAILED;

  mError = new MBError;

  mMBWriteUtil = NULL;
  mMBReadUtil = NULL;
    
    // Readers and writers try to get pointers to above utils.
    // Do this after pointers are initialized. (Pointers should
    // really be initialized in constructor to avoid this kind
    // of thing -- j.kraftcheck.)
  readerWriterSet = new MBReaderWriterSet( this, mError );
  if (!readerWriterSet)
    return MB_MEMORY_ALLOCATION_FAILED;

  MBErrorCode result = create_meshset(0, myMeshSet);
  if (MB_SUCCESS != result) return result;
  
  material_tag();
  neumannBC_tag();
  dirichletBC_tag();
  geom_dimension_tag();
  globalId_tag();

  return MB_SUCCESS;
}

void MBCore::deinitialize()
{
  if (aEntityFactory)
    delete aEntityFactory;

  aEntityFactory = 0;

  if (tagServer)
    delete tagServer;

  tagServer = 0;
  
  if (sequenceManager)
    delete sequenceManager;

  sequenceManager = 0;

  std::map<MBEntityHandle, MBMeshSet*>::iterator iter;
  for(iter = global_mesh_set_list.begin();
      iter != global_mesh_set_list.end();
      ++iter)
  {
    iter->second->set_adj_factory(NULL);
    delete iter->second;
  }
  
  global_mesh_set_list.clear();

  if(mError)
    delete mError;
  mError = 0;
}

MBErrorCode MBCore::query_interface(const std::string& iface_name, void** iface)
{
  if(iface_name == "MBReadUtilIface")
  {
    if(mMBReadUtil)
      *iface = (MBReadUtilIface*)mMBReadUtil;
    else
      *iface = (MBReadUtilIface*)(mMBReadUtil = new MBReadUtil(this, mError));
    return MB_SUCCESS;
  }
  else if(iface_name == "MBWriteUtilIface")
  {
    if(mMBWriteUtil)
      *iface = (MBWriteUtilIface*)mMBWriteUtil;
    else
      *iface = (MBWriteUtilIface*)(mMBWriteUtil = new MBWriteUtil(this, mError));
    return MB_SUCCESS;
  }
  else if(iface_name == "MBReaderWriterSet")
  {
    *iface = reader_writer_set();
    return MB_SUCCESS;
  }
  return MB_FAILURE;
}


MBErrorCode MBCore::release_interface(const std::string& iface_name, void* iface)
{
  if(iface == NULL)
    return MB_FAILURE;

  if(iface_name == "MBReadUtilIface")
  {
      // Is it possible to crash here?  We should fail gracefully instead.
    return MB_SUCCESS;
  }
  else if(iface_name == "MBWriteUtilIface")
  {
    return MB_SUCCESS;
  }
  else if(iface_name == "MBReaderWriterSet")
  {
    return MB_SUCCESS;
  }
  
  return MB_FAILURE;
}


#ifdef XPCOM_MB
// provides basic implementation of nsISupports methods
NS_IMPL_ISUPPORTS1_CI(MBCore, MBInterface);
#endif

int MBCore::QueryInterface(const MBuuid& uuid, MBUnknownInterface** iface)
{
  *iface = 0;
  if(uuid == IDD_MBUnknown)
    *iface = this;
  if(uuid == IDD_MBCore)
    *iface = this;
  else
    return 0;
  return 1;
}

//! get the type from a handle, returns type
MBEntityType MBCore::type_from_handle(const MBEntityHandle handle) const
{
  return TYPE_FROM_HANDLE(handle);
}
  
//! get the id from a handle, returns id
unsigned int MBCore::id_from_handle(const MBEntityHandle handle) const
{
  return ID_FROM_HANDLE(handle);
}

//! get a handle from an id and type
MBErrorCode MBCore::handle_from_id(const MBEntityType type, 
                                     const unsigned int id, 
                                     MBEntityHandle& handle) const
{
  static int err;
  handle = CREATE_HANDLE(type, id, err);

    //check to see if handle exists 
  MBEntitySequence *dummy_seq = 0;
  MBErrorCode error_code = sequence_manager()->find(handle, dummy_seq);
  if(error_code == MB_SUCCESS && type != MBVERTEX)
    error_code = dummy_seq->is_valid_entity(handle) ? MB_SUCCESS : MB_ENTITY_NOT_FOUND;

  return error_code; 
}

int MBCore::dimension_from_handle(const MBEntityHandle handle) const
{
  return MBCN::Dimension(TYPE_FROM_HANDLE(handle));
}

//! load mesh from data in file
//! NOTE: if there is mesh already present, the new mesh will be added
MBErrorCode  MBCore::load_mesh( const char *file_name,
                                const int* block_id_list,
                                const int num_blocks )
{
  MBErrorCode rval;
  const MBReaderWriterSet* set = reader_writer_set();
  
    // Try using the file extension to select a reader
  MBReaderIface* reader = set->get_file_extension_reader( file_name );
  if (reader)
  { 
    rval = reader->load_file( file_name, block_id_list, num_blocks );
    delete reader;
    return rval;
  }
  
    // Try all the readers
  MBReaderWriterSet::iter_type iter;
  for (iter = set->begin(); iter != set->end(); ++iter)
  {
    MBReaderIface* reader = iter->make_reader( this );
    if (NULL != reader)
    {
      rval = reader->load_file( file_name, block_id_list, num_blocks );
      delete reader;
      if (MB_SUCCESS == rval)
        return MB_SUCCESS;
    }
  }

  return MB_FAILURE; 
}

MBErrorCode  MBCore::write_mesh(const char *file_name,
                                  const MBEntityHandle *output_list,
                                  const int num_sets)
{
  MBErrorCode rval;
  const MBReaderWriterSet* set = reader_writer_set();
  std::vector<std::string> qa_records;
  const bool overwrite = true;

  MBWriterIface* writer = set->get_file_extension_writer( file_name );
  if (writer == NULL)
  {
    WriteNCDF exowriter(this);
    rval = exowriter.write_file(file_name, overwrite, output_list, num_sets, qa_records, 0);
  }
  else
  {
    rval = writer->write_file(file_name, overwrite, output_list, num_sets, qa_records );
    delete writer;
  }
  
  return rval; 
}




//! deletes all mesh entities from this datastore
MBErrorCode MBCore::delete_mesh()
{

  MBErrorCode result = MB_SUCCESS;

    // perform all deinitialization procedures to clean up
  if (aEntityFactory)
    delete aEntityFactory;
  aEntityFactory = new AEntityFactory(this);

  tagServer->reset_all_data();
  
  if (sequenceManager)
    delete sequenceManager;
  sequenceManager = new EntitySequenceManager();

  std::map<MBEntityHandle, MBMeshSet*>::iterator iter;
  for(iter = global_mesh_set_list.begin();
      iter != global_mesh_set_list.end();
      ++iter)
  {
    iter->second->set_adj_factory(NULL);
    delete iter->second;
  }

  global_mesh_set_list.clear();

  cachedMsPtr = NULL;
  cachedEntityHandle = 0;
  maxMeshSetid = 1;
  
  result = create_meshset(0, myMeshSet);

  return result;
}

  //! get overall geometric dimension
MBErrorCode MBCore::get_dimension(int &dim) const
{
  dim = geometricDimension;
  return MB_SUCCESS;
}

  //! set overall geometric dimension
  /** Returns error if setting to 3 dimensions, mesh has been created, and 
   *  there are only 2 dimensions on that mesh
   */
MBErrorCode MBCore::set_dimension(const int dim) 
{
    // check to see if current dimension is smaller
  if (geometricDimension < dim) 
  {
      // need to check the number of entities 
    int num;
    /*MBErrorCode result = */ get_number_entities_by_dimension(0, geometricDimension, num);
    
      // test written to be more readable but possibly less efficient
      //if (MB_SUCCESS != result) return MB_FAILURE;
      //else if (0 != num && dim == 2 && ycoordTag == 0) return MB_FAILURE;
      //else if (0 != num && dim == 3 && (ycoordTag == 0 || zcoordTag == 0)) return MB_FAILURE;
      //TODO -- replace this with not using xcoordTag, etc...
  }
    
    // if we got here, it's ok to set dimension
  geometricDimension = dim;
  return MB_SUCCESS;
}

  //! get blocked vertex coordinates for all vertices
  /** Blocked = all x, then all y, etc. 
   */
MBErrorCode MBCore::get_vertex_coordinates(std::vector<double> &coords) const
{
    // INEFFICIENT implementation for now, until we get blocked tag access
  MBRange vertices;
  MBErrorCode result = get_entities_by_type(0, MBVERTEX, vertices);
  if (MB_SUCCESS != result) return result;
  
    // the least we can do is resize the vector and only go through the 
    // vertex list once
  int num_verts = vertices.size();
  int vec_pos = 0;
  double xyz[3];
  coords.resize(geometricDimension*num_verts);
  for (MBRange::iterator it = vertices.begin(); it != vertices.end(); it++) 
  {
    result = get_coords(&(*it), 1, xyz);
    if (MB_SUCCESS != result) return result;

    coords[vec_pos] = xyz[0];
    coords[num_verts+vec_pos] = xyz[1];
    coords[2*num_verts+vec_pos] = xyz[2];

    vec_pos++;
  }
  
  return result;
}

MBErrorCode  MBCore::get_coords(const MBRange& entities, double *coords) const
{

  int node_index = 0;

  // lets get ready to iterate the entity sequences and copy data
  // we'll iterate the map in the sequence manager and
  // the entities range at the same time

  MBRange::const_iterator range_iter = entities.begin();
  MBRange::const_iterator range_iter_end = entities.end();

  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter
    = sequence_manager()->entity_map(MBVERTEX)->begin();
  
  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter_end
    = sequence_manager()->entity_map(MBVERTEX)->end();

  // lets find the entity sequence which holds the first entity
  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator seq_iter_lookahead = seq_iter;
  seq_iter_lookahead++;
  for( ; seq_iter_lookahead != seq_iter_end && 
      seq_iter_lookahead->second->get_start_handle() < *range_iter; )
  {
    ++seq_iter;
    ++seq_iter_lookahead;
  }

  // a look ahead iterator
  MBRange::const_iterator range_iter_lookahead = range_iter;

  // our main loop
  for(; range_iter != range_iter_end && seq_iter != seq_iter_end; /* ++ is handled in loop*/ )
  {
    // find a range that fits in the current entity sequence
    for(; range_iter_lookahead != range_iter_end && 
        *range_iter_lookahead <= seq_iter->second->get_end_handle(); 
        ++range_iter_lookahead)
    {}

    double* coord_array[3];
    static_cast<VertexEntitySequence*>(seq_iter->second)->get_coordinate_arrays(
        coord_array[0], coord_array[1], coord_array[2]);
    MBEntityHandle start_ent = seq_iter->second->get_start_handle();

    // for each of the entities in this entity sequence, copy data
    for(MBRange::const_iterator tmp_iter = range_iter; 
        tmp_iter != range_iter_lookahead;
        ++tmp_iter)
    {
      coords[node_index] = coord_array[0][*tmp_iter - start_ent];
      node_index++;
      coords[node_index] = coord_array[1][*tmp_iter - start_ent];
      node_index++;
      coords[node_index] = coord_array[2][*tmp_iter - start_ent];
      node_index++;
    }

    // go to the next entity sequence
    ++seq_iter;
    // start with the next entities
    range_iter = range_iter_lookahead;
  }
  return MB_SUCCESS;
}

MBErrorCode  MBCore::get_coords(const MBEntityHandle* entities, 
                                  const int num_entities, 
                                  double *coords) const
{
  MBErrorCode status = MB_SUCCESS;

  const MBEntityHandle* end = entities + num_entities;
  MBEntitySequence* seq = 0;

  for(const MBEntityHandle* iter = entities; iter != end; ++iter)
  {
    if(TYPE_FROM_HANDLE(*iter) != MBVERTEX)
      return MB_TYPE_OUT_OF_RANGE;

    status = sequence_manager()->find(*iter, seq);
    if(seq == NULL || status != MB_SUCCESS || !seq->is_valid_entity(*iter) )
      return MB_ENTITY_NOT_FOUND;
    
    status = static_cast<VertexEntitySequence*>(seq)->get_coordinates(*iter, coords);

    coords += 3;

    if(status != MB_SUCCESS)
      return status;
  }

  return MB_SUCCESS; 
}


MBErrorCode  MBCore::get_coords(const MBEntityHandle entity_handle, 
                                  const double *& x, const double *& y, const double *& z) const
{

  MBErrorCode status = MB_TYPE_OUT_OF_RANGE;

  if ( TYPE_FROM_HANDLE(entity_handle) == MBVERTEX )
  {
    MBEntitySequence* seq = 0;
    status = sequence_manager()->find(entity_handle, seq);

    if (seq == 0 || status != MB_SUCCESS || !seq->is_valid_entity(entity_handle) )
      return MB_ENTITY_NOT_FOUND;

    status = static_cast<VertexEntitySequence*>(seq)->get_coordinates_ref(entity_handle, 
                                                                          x, y, z);

  }

  return status; 

}

//! set the coordinate information for this handle if it is of type Vertex
//! otherwise, return an error
MBErrorCode  MBCore::set_coords(MBEntityHandle *entity_handles, const int num_entities, 
                                  const double *coords)
{

  MBErrorCode status = MB_SUCCESS;

  int i, j = 0;

  for (i = 0; i < num_entities; i++) {
    if ( TYPE_FROM_HANDLE(entity_handles[i]) == MBVERTEX )
    {
      MBEntitySequence* seq = 0;
      status = sequence_manager()->find(entity_handles[i], seq);

      if (seq != 0 && status == MB_SUCCESS) {
        status = static_cast<VertexEntitySequence*>(seq)->set_coordinates(entity_handles[i], coords[j], coords[j+1], coords[j+2]);
        j += 3;
      }
    }
    else if (status == MB_SUCCESS)
      status = MB_TYPE_OUT_OF_RANGE;
  }

  return status; 

}

  //! get global connectivity array for specified entity type
  /**  Assumes just vertices, no higher order nodes
   */
MBErrorCode MBCore::get_connectivity_by_type(const MBEntityType type, 
                                               std::vector<MBEntityHandle> &connect) const
{
    // inefficient implementation until we get blocked tag access
  
    // get the range of entities of this type
  MBRange this_range;
  MBErrorCode result = get_entities_by_type(0, type, this_range);
  
  int num_ents = this_range.size();
  connect.reserve(num_ents*MBCN::VerticesPerEntity(type));
  
    // now loop over these entities, getting connectivity for each
  for (MBRange::iterator this_it = this_range.begin(); 
       this_it != this_range.end();
       this_it++)
  {
    const MBEntityHandle *connect_vec;
    result = get_connectivity(*this_it, connect_vec, num_ents, true);
    if (MB_SUCCESS != result) 
      return result;
    connect.insert(connect.end(), &connect_vec[0], &connect_vec[num_ents]); 
  }
  
  return MB_SUCCESS;
}
  

//! get the connectivity for element /handles.  For non-element handles, return an error
MBErrorCode  MBCore::get_connectivity(const MBEntityHandle *entity_handles, 
                                        const int num_handles,
                                        std::vector<MBEntityHandle> &connectivity,
                                        bool topological_connectivity) const
{
  if (num_handles == 0) return MB_FAILURE;
  
    // Make sure the entity should have a connectivity.
  MBEntityType type = TYPE_FROM_HANDLE(entity_handles[0]);

  int i;
  connectivity.clear();

  // handle the special case where someone asks for the connectivity
  // of a vertex.  This actually kind of makes sense for a sphere element.
  if (type == MBVERTEX)
  {
    connectivity.reserve(num_handles);
    std::copy(entity_handles, entity_handles+num_handles, 
              std::back_inserter(connectivity));
    return MB_SUCCESS;
  }

    // WARNING: This is very dependent on the ordering of the MBEntityType enum
  if(type <= MBVERTEX || type >= MBENTITYSET)
    return MB_TYPE_OUT_OF_RANGE;

  MBErrorCode result = MB_SUCCESS, temp_result;
  for (i = 0; i < num_handles; i++) {
    MBEntitySequence* seq = 0;

      // We know that connectivity is stored in an EntitySequence so jump straight
      // to the entity sequence
    temp_result = sequence_manager()->find(entity_handles[i], seq);
    if (seq == NULL || !seq->is_valid_entity(entity_handles[i])) {
      result = MB_ENTITY_NOT_FOUND;
      continue;
    }
    else if (temp_result != MB_SUCCESS) {
      result = temp_result;
      continue;
    }

    ElementEntitySequence *elem_seq = static_cast<ElementEntitySequence*>(seq);
      // let's be smart about this...
    temp_result = elem_seq->get_connectivity(entity_handles[i], connectivity,
                                             topological_connectivity);
    if (MB_SUCCESS != temp_result) {
      result = temp_result;
      continue;
    }
  }
  
  return result;
}

//! get the connectivity for element handles.  For non-element handles, return an error
MBErrorCode MBCore::get_connectivity(const MBEntityHandle entity_handle, 
                                     const MBEntityHandle*& connectivity,
                                     int& number_nodes,
                                     bool topological_connectivity) const
{

  MBErrorCode status;

    // Make sure the entity should have a connectivity.
  MBEntityType type = TYPE_FROM_HANDLE(entity_handle);
  
    // WARNING: This is very dependent on the ordering of the MBEntityType enum
  if(type <= MBVERTEX || type >= MBENTITYSET)
    return MB_TYPE_OUT_OF_RANGE;
  
  MBEntitySequence* seq = 0;

    // We know that connectivity is stored in an EntitySequence so jump straight
    // to the entity sequence
  status = sequence_manager()->find(entity_handle, seq);
  if (seq == 0 || status != MB_SUCCESS || !seq->is_valid_entity(entity_handle)) 
    return MB_ENTITY_NOT_FOUND;

  return static_cast<ElementEntitySequence*>(seq)->get_connectivity(entity_handle, connectivity,
                                                                    number_nodes,
                                                                    topological_connectivity);
}

//! set the connectivity for element handles.  For non-element handles, return an error
MBErrorCode  MBCore::set_connectivity(const MBEntityHandle entity_handle, 
                                        std::vector<MBEntityHandle> &connectivity)
{
  MBErrorCode status = MB_FAILURE;

    // Make sure the entity should have a connectivity.
    // WARNING: This is very dependent on the ordering of the MBEntityType enum
  MBEntityType type = TYPE_FROM_HANDLE(entity_handle);
  
  static std::vector<MBEntityHandle> tmp(31);
  
  MBEntitySequence* seq = 0;

  if (type < MBVERTEX || type > MBENTITYSET)
    return MB_TYPE_OUT_OF_RANGE;
  
  status = sequence_manager()->find(entity_handle, seq);
  if (seq == 0 || status != MB_SUCCESS || !seq->is_valid_entity(entity_handle))
    return (status != MB_SUCCESS ? status : MB_ENTITY_NOT_FOUND);

  status = static_cast<ElementEntitySequence*>(seq)->get_connectivity(entity_handle, tmp);
  if (status != MB_SUCCESS) return status;
  
  status = static_cast<ElementEntitySequence*>(seq)->set_connectivity(entity_handle, 
                                                                      &connectivity[0], 
                                                                      connectivity.size());
  if (status != MB_SUCCESS) return status;

  aEntityFactory->notify_change_connectivity(
    entity_handle, &tmp[0], &connectivity[0], connectivity.size());

  return status;
}

MBErrorCode MBCore::get_adjacencies(const MBEntityHandle *from_entities,
                                      const int num_entities,
                                      const int to_dimension,
                                      const bool create_if_missing,
                                      std::vector<MBEntityHandle> &adj_entities,
                                      const int operation_type)
{
  MBErrorCode result;
  if (num_entities == 1) {
    if(to_dimension == 0 && TYPE_FROM_HANDLE(from_entities[0]) != MBPOLYHEDRON)
      result = get_connectivity(&from_entities[0], 1, adj_entities);
    else
      result = aEntityFactory->get_adjacencies(from_entities[0], to_dimension, 
                                               create_if_missing, adj_entities);
    std::vector<MBEntityHandle>::iterator iter = 
      std::remove(adj_entities.begin(), adj_entities.end(), MBEntityHandle(0));
    if(iter != adj_entities.end())
      adj_entities.erase(iter, adj_entities.end());
    return result;
  }
  
  MBRange temp_range, temp_range2;
  std::copy(from_entities, from_entities+num_entities,
            mb_range_inserter(temp_range));
  
  result = get_adjacencies(temp_range, to_dimension, create_if_missing, temp_range2,
                           operation_type);
  if (MB_SUCCESS != result) 
    return result;
  
  adj_entities.clear();
  adj_entities.reserve(temp_range2.size());
  for (MBRange::const_iterator it = temp_range2.begin();
       it != temp_range2.end();
       ++it)
    adj_entities.push_back(*it);
  
  return MB_SUCCESS;
}

MBErrorCode MBCore::get_adjacencies(const MBRange &from_entities,
                                      const int to_dimension,
                                      const bool create_if_missing,
                                      MBRange &adj_entities,
                                      const int operation_type)
{
  if (operation_type != MBInterface::INTERSECT &&
      operation_type != MBInterface::UNION) return MB_FAILURE;

  if(from_entities.size() == 0)
    return MB_SUCCESS;

  MBRange temp_range;
  std::vector<MBEntityHandle> temp_vec;
  MBErrorCode result = MB_SUCCESS, tmp_result;

  for (MBRange::const_iterator from_it = from_entities.begin(); 
       from_it != from_entities.end(); from_it++) 
  {
      // running results kept in adj_entities; clear temp_vec and temp_range, which are working space
    temp_vec.clear();
    temp_range.clear();

      // get the next set of adjacencies
    if(to_dimension == 0 && TYPE_FROM_HANDLE(*from_it) != MBPOLYHEDRON)
      tmp_result = get_connectivity(&(*from_it), 1, temp_vec);
    else
      tmp_result = aEntityFactory->get_adjacencies(*from_it, to_dimension, 
                                                   create_if_missing, temp_vec);
    
    if (MB_SUCCESS != tmp_result) result = tmp_result;
    
      // if we're on the first iteration and we didn't come in with entities,
      // just get the first results and move on
    if (adj_entities.empty() && from_it == from_entities.begin()) {
      std::copy(temp_vec.begin(), temp_vec.end(), mb_range_inserter(adj_entities));
      continue;
    }

      // operate on the vectors
    if (operation_type == MBInterface::INTERSECT) {
        // only have to sort if we're doing intersection
      std::sort(temp_vec.begin(), temp_vec.end());
      std::set_intersection(adj_entities.begin(), adj_entities.end(), 
                            temp_vec.begin(), temp_vec.end(),
                            mb_range_inserter(temp_range));
      adj_entities = temp_range;
    }
    else if (operation_type == MBInterface::UNION) {
      std::copy(temp_vec.begin(), temp_vec.end(), mb_range_inserter(adj_entities));
    }
  }

  return result;
}

MBErrorCode MBCore::add_adjacencies(const MBEntityHandle entity_handle, 
                                      const MBEntityHandle *adjacencies,
                                      const int num_handles,
                                      bool both_ways)
{
  MBErrorCode result = MB_SUCCESS, temp_result;
  
  for (const MBEntityHandle *it = adjacencies; 
       it != adjacencies+num_handles; it++) {
    temp_result = aEntityFactory->add_adjacency(entity_handle, *it, both_ways);
    if (MB_SUCCESS != temp_result) result = temp_result;
  }

  return result;
}

MBErrorCode MBCore::remove_adjacencies(const MBEntityHandle entity_handle,
                                         const MBEntityHandle *adjacencies,
                                         const int num_handles)
{
  MBErrorCode result = MB_SUCCESS, temp_result;
  
  for (const MBEntityHandle *it = adjacencies; 
       it != adjacencies+num_handles; it++) {
    temp_result = aEntityFactory->remove_adjacency(entity_handle, *it);
    if (MB_SUCCESS != temp_result) result = temp_result;
    temp_result = aEntityFactory->remove_adjacency(*it, entity_handle);
    if (MB_SUCCESS != temp_result) result = temp_result;
  }

  return result;
}

MBErrorCode MBCore::get_entities_by_dimension(const MBEntityHandle meshset,
                                                const int dimension, 
                                                MBRange &entities,
                                                const bool recursive) const
{
  MBErrorCode result = MB_SUCCESS;
  if (0 != meshset) {
    MBMeshSet *ms_ptr = update_cache(meshset);
    if(NULL == ms_ptr) return MB_ENTITY_NOT_FOUND;
    MBRange dum_range;
    result = ms_ptr->get_entities(dum_range, recursive);
    if (MB_SUCCESS != result) return result;
    for (MBRange::reverse_iterator it = dum_range.rbegin(); it != dum_range.rend(); it++)
      if (MBCN::Dimension(TYPE_FROM_HANDLE(*it)) == dimension) entities.insert(*it);
  }
  else {
    for (MBEntityType this_type = MBCN::TypeDimensionMap[dimension].first;
         this_type <= MBCN::TypeDimensionMap[dimension].second;
         this_type++) {
      if (this_type == MBENTITYSET) {
        for (std::map< MBEntityHandle, MBMeshSet*>::const_iterator it = 
               global_mesh_set_list.begin(); it != global_mesh_set_list.end(); it++)
          entities.insert((*it).first);
        result = MB_SUCCESS;
      }
      else
        result = sequence_manager()->get_entities( this_type, entities );
    }
  }

  return result;
}

MBErrorCode MBCore::get_entities_by_type(const MBEntityHandle meshset,
                                           const MBEntityType type, 
                                           MBRange &entities,
                                           const bool recursive) const
{
  MBErrorCode result;
  if (0 != meshset) {
    MBMeshSet *ms_ptr = update_cache(meshset);
    if(NULL == ms_ptr) return MB_ENTITY_NOT_FOUND;
    MBRange dum_range;
    result = ms_ptr->get_entities(dum_range, recursive);
    if (MB_SUCCESS != result) return result;
    for (MBRange::reverse_iterator it = dum_range.rbegin(); it != dum_range.rend(); it++)
      if (TYPE_FROM_HANDLE(*it) == type) entities.insert(*it);
  }
  else if (type == MBENTITYSET) {
    for (std::map< MBEntityHandle, MBMeshSet*>::const_iterator it = 
           global_mesh_set_list.begin(); it != global_mesh_set_list.end(); it++)
      entities.insert((*it).first);
    result = MB_SUCCESS;
  }
  else
    result = sequence_manager()->get_entities( type, entities );

  return result;
}

MBErrorCode MBCore::get_entities_by_type_and_tag(const MBEntityHandle meshset,
                                                   const MBEntityType type,
                                                   const MBTag *tags,
                                                   const void** values,
                                                   const int num_tags,
                                                   MBRange &entities,
                                                   const int condition,
                                                   const bool recursive) const
{
  MBErrorCode result;
  if (0 != meshset) {
    MBMeshSet *ms_ptr = update_cache(meshset);
    if(NULL == ms_ptr) return MB_ENTITY_NOT_FOUND;
    MBRange dum_range;
    result = ms_ptr->get_entities(dum_range, recursive);
    if (MB_SUCCESS != result) return result;
    result = tagServer->get_entities_with_tag_values(dum_range, type, 
                                                     tags, values, num_tags, 
                                                     entities, condition);  
  }
  
  else 
    result = tagServer->get_entities_with_tag_values( type, tags, values, num_tags, 
                                                      entities, condition);
  
  return result;
}


MBErrorCode MBCore::get_entities_by_handle(const MBEntityHandle meshset,
                                             MBRange &entities,
                                             const bool recursive) const
{
  MBErrorCode result;
  if (0 != meshset) {
    MBMeshSet *ms_ptr = update_cache(meshset);
    if(NULL == ms_ptr) return MB_ENTITY_NOT_FOUND;
    MBRange dum_range;
    result = ms_ptr->get_entities(entities, recursive);
    if (MB_SUCCESS != result) return result;
  }
  else {
    result = MB_SUCCESS;
    for (MBEntityType tp = MBVERTEX; tp < MBMAXTYPE; tp++) {
      MBErrorCode tmp_result = sequence_manager()->get_entities( tp, entities );
      if (tmp_result != MB_SUCCESS) result = tmp_result;
    }
  }

  return result;
}


MBErrorCode MBCore::get_entities_by_handle(const MBEntityHandle meshset,
                                   std::vector<MBEntityHandle> &entities,
                                   const bool recursive) const
{
  MBErrorCode result;
  if (0 != meshset) {
    MBMeshSet *ms_ptr = update_cache(meshset);
    if(NULL == ms_ptr) return MB_ENTITY_NOT_FOUND;
    MBRange dum_range;
    result = ms_ptr->get_entities(entities, recursive);
    if (MB_SUCCESS != result) return result;
  }
  else {
    MBRange dum_range;
    result = get_entities_by_handle(meshset, dum_range, recursive);
    if (MB_SUCCESS != result) 
      return result;
    entities.reserve(entities.size() + dum_range.size());
    for (MBRange::const_iterator it = dum_range.begin();
         it != dum_range.end();
         ++it)
      entities.push_back(*it);
  }
  
  return result;
}

  //! get # entities of a given dimension
MBErrorCode MBCore::get_number_entities_by_dimension(const MBEntityHandle meshset,
                                                       const int dim, 
                                                       int &number,
                                                       const bool recursive) const
{
  MBErrorCode result;
  
  number = 0;

  if (0 != meshset) {
    MBMeshSet *ms_ptr = update_cache(meshset);
    if(NULL == ms_ptr) return MB_ENTITY_NOT_FOUND;
    MBRange dum_range;
    result = ms_ptr->get_entities(dum_range, recursive);
    for (MBRange::iterator it = dum_range.begin(); it != dum_range.end(); it++)
      if (MBCN::Dimension(TYPE_FROM_HANDLE(*it)) == dim) number++;
    return result;
  }

  for (MBEntityType this_type = MBCN::TypeDimensionMap[dim].first;
       this_type <= MBCN::TypeDimensionMap[dim].second;
       this_type++) {
    int dummy = 0;
    result = get_number_entities_by_type(0, this_type, dummy);
    if (result != MB_SUCCESS) {
      number = 0;
      return result;
    }
    number += dummy;
  }
  
  return MB_SUCCESS;
}

//! returns the number of entities with a given type and tag
MBErrorCode MBCore::get_number_entities_by_type(const MBEntityHandle meshset,
                                                  const MBEntityType type, 
                                                  int& num_ent,
                                                  const bool recursive) const
{
  MBRange dum_ents;
  MBErrorCode result = get_entities_by_type(meshset, type, dum_ents, recursive);
  num_ent = dum_ents.size();
  return result;
}

MBErrorCode MBCore::get_number_entities_by_type_and_tag(const MBEntityHandle meshset,
                                                          const MBEntityType type,
                                                          const MBTag *tag_handles,
                                                          const void** values,
                                                          const int num_tags,
                                                          int &num_entities,
                                                          const bool recursive) const
{
  MBRange dum_ents;
  MBErrorCode result = get_entities_by_type_and_tag(meshset, type, tag_handles, values, num_tags, 
                                                     dum_ents, recursive);
  num_entities = dum_ents.size();
  return result;
}

MBErrorCode MBCore::get_number_entities_by_handle(const MBEntityHandle meshset,
                                          int& num_ent,
                                          const bool recursive) const
{
  MBErrorCode result;
  if (0 != meshset) {
    MBMeshSet *ms_ptr = update_cache(meshset);
    if(NULL == ms_ptr) return MB_ENTITY_NOT_FOUND;
    MBRange dum_range;
    result = ms_ptr->get_entities(dum_range, recursive);
    num_ent = dum_range.size();
    return result;
  }

  num_ent = 0;
  for (MBEntityType this_type = MBVERTEX;
       this_type < MBMAXTYPE;
       this_type++) {
    int dummy = 0;
    result = get_number_entities_by_type(0, this_type, dummy);
    if (result != MB_SUCCESS) {
      num_ent = 0;
      return result;
    }
    num_ent += dummy;
  }

  return MB_SUCCESS;
}

//! return the tag data for a given EntityHandle and MBTag
MBErrorCode  MBCore::tag_get_data(const MBTag tag_handle, 
                                    const MBEntityHandle* entity_handles, 
                                    const int num_entities,
                                    void *tag_data) const
{
  MBTagType tag_type;
  
  if (NULL == entity_handles && 0 == num_entities && 
      MB_SUCCESS == tag_get_type(tag_handle, tag_type))
    return tagServer->get_data(tag_handle, &myMeshSet, 1, tag_data);

  else return tagServer->get_data(tag_handle, entity_handles, num_entities, tag_data);
}

//! return the tag data for a given EntityHandle and MBTag
MBErrorCode  MBCore::tag_get_data(const MBTag tag_handle, 
                                    const MBRange& entity_handles,
                                    void *tag_data) const
{
  return tagServer->get_data(tag_handle, entity_handles, tag_data);
}

//! set the data  for given EntityHandles and MBTag
MBErrorCode  MBCore::tag_set_data(const MBTag tag_handle, 
                                    const MBEntityHandle* entity_handles, 
                                    const int num_entities,
                                    const void *tag_data)
{
  MBTagType tag_type;
  
  if (NULL == entity_handles && 0 == num_entities && 
      MB_SUCCESS == tag_get_type(tag_handle, tag_type))
    return tagServer->set_data(tag_handle, &myMeshSet, 1, tag_data);

  //verify handles
  MBEntitySequence* seq;
  const MBEntityHandle* iter;
  const MBEntityHandle* end = entity_handles + num_entities;
  for(iter = entity_handles; iter != end; ++iter)
  {
    if (TYPE_FROM_HANDLE(*iter) == MBENTITYSET) continue;
    
    else if(sequenceManager->find(*iter, seq) != MB_SUCCESS)
      return MB_ENTITY_NOT_FOUND;
  }

  return tagServer->set_data(tag_handle, entity_handles, num_entities, tag_data);
}

//! set the data  for given EntityHandles and MBTag
MBErrorCode  MBCore::tag_set_data(const MBTag tag_handle, 
                                    const MBRange& entity_handles, 
                                    const void *tag_data)
{
  //verify handles
  MBRange::const_iterator iter;
  MBErrorCode result;
  for(iter = entity_handles.begin(); 
      iter != entity_handles.end() && TYPE_FROM_HANDLE(*iter) != MBENTITYSET;
      ++iter)
  {
    MBEntitySequence* seq = NULL;
    result = sequenceManager->find(*iter, seq);
    if(result != MB_SUCCESS)
      return result;
  }

  return tagServer->set_data(tag_handle, entity_handles, tag_data);
}

//! adds a sparse tag for this specific EntityHandle/tag_name combination
MBErrorCode MBCore::tag_create(const char *tag_name,
                                 const int tag_size, 
                                 const MBTagType tag_type,
                                 MBTag &tag_handle, 
                                 const void *default_value)
{
    // Don't know what to do with the default value yet.
  default_value = default_value;

  return tagServer->add_tag(tag_name, tag_size, tag_type, tag_handle, default_value);
}

//! removes the tag from the entity
MBErrorCode  MBCore::tag_delete_data(const MBTag tag_handle, 
                                       const MBEntityHandle *entity_handles,
                                       const int num_handles)
{
  if (reinterpret_cast<long>(tag_handle) & TagInfo::TagBitProperties[MB_TAG_DENSE])
    return MB_FAILURE;

  MBErrorCode status = MB_SUCCESS, temp_status;
  for (int i = 0; i < num_handles; i++) {
    temp_status = tagServer->remove_data(tag_handle, entity_handles[i]);
    if (temp_status != MB_SUCCESS) status = temp_status;
  }

  return status;
}

//! removes the tag from the entity
MBErrorCode  MBCore::tag_delete_data(const MBTag tag_handle, 
                                       const MBRange &entity_handles)
{
  if (reinterpret_cast<long>(tag_handle) & TagInfo::TagBitProperties[MB_TAG_DENSE])
    return MB_FAILURE;

  MBErrorCode status = MB_SUCCESS, temp_status;
  for (MBRange::const_iterator it = entity_handles.begin(); it != entity_handles.end(); it++) {
    temp_status = tagServer->remove_data(tag_handle, *it);
    if (temp_status != MB_SUCCESS) status = temp_status;
  }

  return status;
}

//! removes the tag from MB
MBErrorCode  MBCore::tag_delete(MBTag tag_handle)
{
  return tag_server()->remove_tag(tag_handle);
}

//! gets the tag name string for the tag_handle
MBErrorCode  MBCore::tag_get_name(const MBTag tag_handle, 
                                    std::string& tag_name) const
{
  const TagInfo* tag_info = tagServer->get_tag_info( tag_handle );
  if(!tag_info)
    return MB_TAG_NOT_FOUND;
  
  tag_name = tag_info->get_name();
  return MB_SUCCESS;

}

//! gets tag handle from its name.
//! the type must be specified because the same name is valid for multiple types
MBErrorCode  MBCore::tag_get_handle(const char *tag_name, 
                                      MBTag &tag_handle) const
{
  MBErrorCode status = MB_TAG_NOT_FOUND;

  tag_handle = tagServer->get_handle( tag_name );
  
  if (tag_handle != 0)
  {
    status = MB_SUCCESS;
  }

  return status;
}

  //! get size of tag in bytes
MBErrorCode MBCore::tag_get_size(const MBTag tag_handle, int &tag_size) const
{
  const TagInfo* tag_info = tagServer->get_tag_info( tag_handle );
  if(!tag_info)
    return MB_TAG_NOT_FOUND;
  
  tag_size = tag_info->get_size();
  return MB_SUCCESS;
}

  //! get default value of the tag
MBErrorCode MBCore::tag_get_default_value(const MBTag tag_handle, void *def_value) const
{
  const TagInfo* tag_info = tagServer->get_tag_info( tag_handle );
  if(!tag_info)
    return MB_TAG_NOT_FOUND;

  if (NULL == def_value) return MB_FAILURE;
  
  if (tag_info->default_value() == NULL)
    return MB_ENTITY_NOT_FOUND;
  
  memcpy(def_value, tag_info->default_value(), tag_info->get_size());

  return MB_SUCCESS;
}

  //! get type of tag (sparse, dense, etc.; 0 = dense, 1 = sparse, 2 = bit, 3 = static)
MBErrorCode MBCore::tag_get_type(const MBTag tag_handle, MBTagType &tag_type) const
{
  if( reinterpret_cast<long>(tag_handle) & TagInfo::TagBitProperties[MB_TAG_DENSE])
    tag_type = MB_TAG_DENSE;
  else if( reinterpret_cast<long>(tag_handle) & TagInfo::TagBitProperties[MB_TAG_SPARSE])
    tag_type = MB_TAG_SPARSE;
  else if( reinterpret_cast<long>(tag_handle) & TagInfo::TagBitProperties[MB_TAG_BIT])
    tag_type = MB_TAG_BIT;
  else if( reinterpret_cast<long>(tag_handle) & TagInfo::TagBitProperties[MB_TAG_MESH])
    tag_type = MB_TAG_MESH;
  else {
    tag_type = MB_TAG_LAST;
    return MB_TAG_NOT_FOUND;
  }
  
  return MB_SUCCESS;
}

  //! get handles for all tags defined
MBErrorCode MBCore::tag_get_tags(std::vector<MBTag> &tag_handles) const
{
  return tagServer->get_tags(tag_handles);
}

  //! Get handles for all tags defined on this entity
MBErrorCode MBCore::tag_get_tags_on_entity(const MBEntityHandle entity,
                                            std::vector<MBTag> &tag_handles) const 
{
  if (0 == entity)
    return tagServer->get_tags(myMeshSet, tag_handles);
  else return tagServer->get_tags(entity, tag_handles);
}

MBTag MBCore::material_tag()
{
  if (0 == materialTag)
    tagServer->add_tag(MATERIAL_SET_TAG_NAME, sizeof(int), 
                       MB_TAG_SPARSE, materialTag);
  return materialTag;
}

MBTag MBCore::neumannBC_tag()
{
  if (0 == neumannBCTag)
    tagServer->add_tag(NEUMANN_SET_TAG_NAME, sizeof(int), 
                       MB_TAG_SPARSE, neumannBCTag);
  return neumannBCTag;
}

MBTag MBCore::dirichletBC_tag()
{
  if (0 == dirichletBCTag)
    tagServer->add_tag(DIRICHLET_SET_TAG_NAME, sizeof(int), 
                       MB_TAG_SPARSE, dirichletBCTag);
  return dirichletBCTag;
}

MBTag MBCore::globalId_tag()
{
  if (0 == globalIdTag)
    tagServer->add_tag(GLOBAL_ID_TAG_NAME, sizeof(int), 
                       MB_TAG_DENSE, globalIdTag);
  return globalIdTag;
}

MBTag MBCore::geom_dimension_tag()
{
  if (0 == geomDimensionTag)
    tagServer->add_tag(GEOM_DIMENSION_TAG_NAME, sizeof(int), 
                       MB_TAG_SPARSE, geomDimensionTag);
  return geomDimensionTag;
}

//! creates an element based on the type and connectivity.  returns a handle and error code
MBErrorCode MBCore::create_element(const MBEntityType type, 
                                   const MBEntityHandle *connectivity,
                                   const int num_nodes, 
                                   MBEntityHandle &handle)
{

    // make sure we have enough vertices for this entity type
  if(num_nodes < MBCN::VerticesPerEntity(type))
    return MB_FAILURE;
  
  MBErrorCode status = sequence_manager()->create_element(type, connectivity, num_nodes, handle);
  if (MB_SUCCESS == status)
    status = aEntityFactory->notify_create_entity( handle, connectivity, num_nodes); 

  return status;

}

//! creates a vertex based on coordinates, returns a handle and error code
MBErrorCode MBCore::create_vertex(const double coords[3], MBEntityHandle &handle )
{
    // get an available vertex handle
  return sequence_manager()->create_vertex(coords, handle);
}

// merge entities without doing a down check because it's already been done
MBErrorCode MBCore::merge_entities_up(MBEntityHandle entity_to_keep,
                                        MBEntityHandle entity_to_remove,
                                        bool auto_merge,
                                        bool delete_removed_entity)
{
    // The two entities to merge must be of the same type
  MBEntityType type_to_keep = TYPE_FROM_HANDLE(entity_to_keep);

  if (type_to_keep != TYPE_FROM_HANDLE(entity_to_remove))
    return MB_TYPE_OUT_OF_RANGE;

  int ent_dim = MBCN::Dimension(type_to_keep);
  
    // Make sure both entities exist before trying to merge.
  if(ent_dim > 0)
  {
    MBEntitySequence* seq = 0;
    MBErrorCode status;
    status = sequence_manager()->find(entity_to_keep, seq);
    if(seq == 0 || status != MB_SUCCESS ||
       !seq->is_valid_entity(entity_to_keep))
      return MB_ENTITY_NOT_FOUND;
    status = sequence_manager()->find(entity_to_remove, seq);
    if(seq == 0 || status != MB_SUCCESS ||
       !seq->is_valid_entity(entity_to_remove))
      return MB_ENTITY_NOT_FOUND;
  }

  std::vector<MBEntityHandle> adjacencies, tmp, conn;
  std::vector<MBEntityHandle>::iterator iter;
  MBErrorCode result;

    // gather all the adjacencies
  result = aEntityFactory->get_adjacencies(entity_to_remove, 1, false, tmp);
  if(result != MB_SUCCESS)
    return result;

  adjacencies.insert(adjacencies.end(), tmp.begin(), tmp.end());

  result = aEntityFactory->get_adjacencies(entity_to_remove, 2, false, tmp);
  if(result != MB_SUCCESS)
    return result;

  adjacencies.insert(adjacencies.end(), tmp.begin(), tmp.end());

  result = aEntityFactory->get_adjacencies(entity_to_remove, 3, false, tmp);
  if(result != MB_SUCCESS)
    return result;

  adjacencies.insert(adjacencies.end(), tmp.begin(), tmp.end());

  result = aEntityFactory->get_adjacencies(entity_to_remove, 4, false, tmp);
  if (result != MB_SUCCESS)
        return result;
    
  std::copy(tmp.begin(), tmp.end(),
            std::back_inserter<std::vector<MBEntityHandle> >(adjacencies));

  for(iter = adjacencies.begin(); iter != adjacencies.end(); iter++)
  {
      // attempt to update the connectivity
    if(ent_dim == 0)
    {
      result = get_connectivity(&(*iter), 1, conn);

      if(result == MB_SUCCESS)
      {
        std::replace(conn.begin(), conn.end(), entity_to_remove, entity_to_keep);
        set_connectivity(*iter, conn);
      }
    }
    else
    {
        // update the adjacency with the new entity
      result = aEntityFactory->remove_adjacency(entity_to_remove, *iter);
      if(result != MB_SUCCESS)
        return result;
      result = aEntityFactory->add_adjacency(entity_to_keep, *iter, false);
      if(result != MB_SUCCESS)
        return result;
    }

      // TODO: figure out how to update the meshsets
  }

    // delete the merged out entity
  if (delete_removed_entity) delete_entities(&entity_to_remove, 1);

    // Does this merge cause other entities to be coincident?
  if(auto_merge && ent_dim < 3)
  {
    std::vector<MBEntityHandle> conn2;
    std::vector<MBEntityHandle>::iterator jter;
    int i = ent_dim + 1;
    for(; i<4; i++)
    {
      result = aEntityFactory->get_adjacencies(entity_to_keep, i, false, tmp);
      if(result != MB_SUCCESS)
        return result;

      for(iter=tmp.begin(); iter<tmp.end(); iter++)
      {
        result = get_connectivity(&(*iter), 1, conn);
        if(result != MB_SUCCESS)
          continue;

        for(jter=iter+1; jter<tmp.end(); jter++)
        {
          result = get_connectivity(&(*jter), 1, conn2);
          if(result != MB_SUCCESS)
            continue;

          if(conn == conn2)
          {
            result = merge_entities_up(*iter, *jter, false, delete_removed_entity);
            if(result != MB_SUCCESS && result != MB_ENTITY_NOT_FOUND)
              return result;
          }
        }
      }
    }
  }

  return MB_SUCCESS;
}

//! merges two  entities
MBErrorCode MBCore::merge_entities( MBEntityHandle entity_to_keep, 
                                      MBEntityHandle entity_to_remove,
                                      bool auto_merge,
                                      bool delete_removed_entity)
{
    // The two entities to merge must be of the same type
  MBEntityType type_to_keep = TYPE_FROM_HANDLE(entity_to_keep);

  if (type_to_keep != TYPE_FROM_HANDLE(entity_to_remove))
    return MB_TYPE_OUT_OF_RANGE;

    // Make sure both entities exist before trying to merge.
  MBEntitySequence* seq = 0;
  MBErrorCode result, status;
  status = sequence_manager()->find(entity_to_keep, seq);
  if(seq == 0 || status != MB_SUCCESS ||
     !seq->is_valid_entity(entity_to_keep))
    return MB_ENTITY_NOT_FOUND;
  status = sequence_manager()->find(entity_to_remove, seq);
  if(seq == 0 || status != MB_SUCCESS ||
     !seq->is_valid_entity(entity_to_remove))
    return MB_ENTITY_NOT_FOUND;
  
    // If auto_merge is not set, all sub-entities should
    // be merged if the entities are to be merged.
  int ent_dim = MBCN::Dimension(type_to_keep);
  if(ent_dim > 0)
  {
    std::vector<MBEntityHandle> adjacencies, tmp, conn;
    std::vector<MBEntityHandle> adjacencies2, tmp2, conn2;
    std::vector<MBEntityHandle>::iterator iter, jter;

    result = get_connectivity(&entity_to_keep, 1, conn);
    if(result != MB_SUCCESS)
      return result;
    result = get_connectivity(&entity_to_remove, 1, conn2);
    if(result != MB_SUCCESS)
      return result;

      // Check to see if we can merge before pulling adjacencies.
    if(!auto_merge && conn != conn2)
      return MB_FAILURE;

    int i = 1;
    for(; i<ent_dim; i++)
    {
      result = aEntityFactory->get_adjacencies(entity_to_keep, i, false, tmp);
      if(result != MB_SUCCESS)
        return result;
      result = aEntityFactory->get_adjacencies(entity_to_remove, i, false, tmp2);
      if(result != MB_SUCCESS)
        return result;
      if(tmp.size() != tmp2.size())
        return MB_FAILURE;

        // Sort the entities for comparison
      std::sort(tmp.begin(), tmp.end());
      std::sort(tmp2.begin(), tmp2.end());

        // Put these sets in with the other ones.
      adjacencies.insert(adjacencies.end(), tmp.begin(), tmp.end());
      adjacencies.insert(adjacencies.end(), tmp2.begin(), tmp2.end());
    }

      // Make sure we can merge before doind so.
    if(auto_merge)
    {
        // Merge all the vertecies for the two entities.
      for(iter=conn.begin(), jter=conn2.begin();
          iter != conn.end() || jter != conn2.end(); iter++, jter++)
      {
        if(*iter != *jter)
        {
          result = merge_entities_up(*iter, *jter, true, delete_removed_entity);
          if(result != MB_SUCCESS)
            return result;
        }
      }

        // Merge the remainder of the sub-entities.
      for(iter=adjacencies.begin(), jter=adjacencies2.begin();
          iter != adjacencies.end() || jter != adjacencies2.end(); iter++, jter++)
      {
        if(*iter != *jter)
        {
            // Check to see if entities still exist. They may have
            // been merged when merging a lower entity.
          status = sequence_manager()->find(*iter, seq);
          if(seq != 0 && status == MB_SUCCESS && seq->is_valid_entity(*iter))
          {
            status = sequence_manager()->find(*jter, seq);
            if(seq != 0 && status == MB_SUCCESS && seq->is_valid_entity(*jter))
            {
              result = merge_entities_up(*iter, *jter, true, delete_removed_entity);
              if(result != MB_SUCCESS && result != MB_ENTITY_NOT_FOUND)
                return result;
            }
          }
        }
      }
    }
    else
    {
        // Don't merge if the adjacencies don't match.
      std::sort(adjacencies.begin(), adjacencies.end());
      std::sort(adjacencies2.begin(), adjacencies2.end());
      if(adjacencies != adjacencies2)
        return MB_FAILURE;
    }
  }

  return merge_entities_up(entity_to_keep, entity_to_remove, auto_merge, 
                           delete_removed_entity);
}

//! deletes an entity vector
MBErrorCode MBCore::delete_entities(const MBEntityHandle *entities,
                                      const int num_entities)
{
  MBRange entity_range;
  std::copy(entities, entities+num_entities,
            mb_range_inserter(entity_range));
  return delete_entities(entity_range);
}


//! deletes an entity range
MBErrorCode MBCore::delete_entities(const MBRange &range)
{
  std::vector<MBEntityHandle> adj_list;
  MBErrorCode result = MB_SUCCESS, temp_result;
  
  for (MBRange::const_iterator it = range.begin(); it != range.end(); it++) {
    MBEntityType type = TYPE_FROM_HANDLE(*it);

      // delete sets a different way
    if (MBENTITYSET == type) {
      std::map<MBEntityHandle, MBMeshSet*>::iterator iter;
      iter = global_mesh_set_list.find(*it);
      if(iter != global_mesh_set_list.end())
      {
        delete iter->second;
        tagServer->reset_data(*it);
        global_mesh_set_list.erase(iter);
      }
      continue;
    }
    
      // tell AEntityFactory that this element is going away
    temp_result = aEntityFactory->notify_delete_entity(*it);
    if (MB_SUCCESS != temp_result) {
      result = temp_result;
      continue;
    }

      // reset and/or clean out data associated with this entity handle
    temp_result = tagServer->reset_data(*it);
    if (MB_SUCCESS != temp_result) {
      result = temp_result;
      continue;
    }

      // now delete the entity
    temp_result = sequence_manager()->delete_entity(*it);
    if (MB_SUCCESS != temp_result) {
      result = temp_result;
      continue;
    }
  }

  return result;
}

MBErrorCode MBCore::list_entities(const MBEntityHandle *entities,
                                    const int num_entities) const
{
  MBRange temp_range;
  MBErrorCode result;
  if (NULL == entities && num_entities <= 0) {
      // just list the numbers of each entity type
    int num_ents;
    std::cout << std::endl;
    std::cout << "Number of entities per type: " << std::endl;
    for (MBEntityType this_type = MBVERTEX; this_type < MBMAXTYPE; this_type++) {
      result = get_number_entities_by_type(0, this_type, num_ents);
      std::cout << MBCN::EntityTypeName(this_type) << ": " << num_ents << std::endl;
    }
    std::cout << std::endl;

      // if negative num_entities, list the set hierarchy too
    if (0 > num_entities) {
      MBRange sets;
      result = this->get_entities_by_type(0, MBENTITYSET, sets);
      if (MB_SUCCESS != result) return result;
      for (MBRange::iterator rit = sets.begin(); rit != sets.end(); rit++) {
        this->print(*rit, "", false);
        result = this->get_number_entities_by_handle(*rit, num_ents);
        std::cout << "(" << num_ents << " total entities)" << std::endl;
      }
    }
    
    return MB_SUCCESS;
  }
      
  else if (NULL == entities) {

      // list all entities of all types
    std::cout << std::endl;
    for (MBEntityType this_type = MBVERTEX; this_type < MBMAXTYPE; this_type++) {
      temp_range.clear();
      std::cout << "Entities of type: " << MBCN::EntityTypeName(this_type) << ": " << std::endl;
      result = get_entities_by_type(0, this_type, temp_range);
      result = list_entities(temp_range);
    }
    std::cout << std::endl;
    
    return MB_SUCCESS;
  }

  else {
    std::copy(entities, entities+num_entities, mb_range_inserter(temp_range));
    return list_entities(temp_range);
  }
}
  
MBErrorCode MBCore::list_entities(const MBRange &entities) const
{
  MBErrorCode result;
  MBHandleVec adj_vec;

    // list entities
  for (MBRange::const_iterator iter = entities.begin(); iter != entities.end(); iter++) {
    MBEntityType this_type = TYPE_FROM_HANDLE(*iter);
    std::cout << MBCN::EntityTypeName(this_type) << " " << ID_FROM_HANDLE(*iter) << ":" << endl;

    if (this_type == MBVERTEX) {
      double coords[3];
      result = get_coords(&(*iter), 1, coords);
      if (MB_SUCCESS != result) return result;
      std::cout << "Coordinates: (" << coords[0] << ", " << coords[1] << ", " << coords[2] 
           << ")" << std::endl;
    }
    else if (this_type == MBENTITYSET)
      this->print(*iter, "");
    
    std::cout << "Adjacencies:" << std::endl;
    bool some = false;
    for (int dim = 0; dim <= 3; dim++) {
      if (dim == MBCN::Dimension(this_type)) continue;
      adj_vec.clear();
        // use const_cast here 'cuz we're in a const function and we're passing 'false' for
        // create_if_missing, so we know we won't change anything
      result = (const_cast<MBCore*>(this))->get_adjacencies(&(*iter), 1, dim, false, adj_vec);
      if (MB_FAILURE == result) continue;
      for (MBHandleVec::iterator adj_it = adj_vec.begin(); adj_it != adj_vec.end(); adj_it++) {
        if (adj_it != adj_vec.begin()) std::cout << ", ";
        else std::cout << std::endl;
        std::cout << MBCN::EntityTypeName(TYPE_FROM_HANDLE(*adj_it)) << " " << ID_FROM_HANDLE(*adj_it);
      }
      if (!adj_vec.empty()) {
        std::cout << std::endl;
        some = true;
      }
    }
    if (!some) std::cout << "(none)" << std::endl;

    std::cout << std::endl;
  }

  return MB_SUCCESS;
}

MBErrorCode MBCore::convert_entities( const MBEntityHandle meshset, 
                                        const bool mid_side, const bool mid_face, const bool mid_volume,
                                        MBInterface::HONodeAddedRemoved* function_object )
{
  HigherOrderFactory fact(this, function_object);
  return fact.convert(meshset, mid_side, mid_face, mid_volume);
}

  //! function to get the side number given two elements; returns
  //! MB_FAILURE if child not related to parent; does *not* create adjacencies
  //! between parent and child
MBErrorCode MBCore::side_number(const MBEntityHandle parent,
                                  const MBEntityHandle child,
                                  int &side_number,
                                  int &sense,
                                  int &offset) const
{
    // get the connectivity of parent and child
  const MBEntityHandle *parent_conn, *child_conn;
  int num_parent_vertices, num_child_vertices;
  MBErrorCode result = get_connectivity(parent, parent_conn, num_parent_vertices, true);
  if (MB_SUCCESS != result) return result;
  result = get_connectivity(child, child_conn, num_child_vertices, true);
  if (MB_SUCCESS != result) return result;

    // call handle vector-based function
  int temp_result = MBCN::SideNumber(parent_conn, TYPE_FROM_HANDLE(parent),
                                       child_conn, num_child_vertices, 
                                       MBCN::Dimension(TYPE_FROM_HANDLE(child)), 
                                       side_number, sense, offset);
  return (0 == temp_result ? MB_SUCCESS : MB_FAILURE);
}

  //! given an entity and the connectivity and type of one of its subfacets, find the
  //! high order node on that subfacet, if any
MBErrorCode MBCore::high_order_node(const MBEntityHandle parent_handle,
                                      const MBEntityHandle *subfacet_conn,
                                      const MBEntityType subfacet_type,
                                      MBEntityHandle &high_order_node) const
{
  high_order_node = 0;

  MBEntityType parent_type = TYPE_FROM_HANDLE(parent_handle);

    // get the parent's connectivity
  const MBEntityHandle *parent_conn;
  int num_parent_vertices;
  MBErrorCode result = get_connectivity(parent_handle, parent_conn, 
                                         num_parent_vertices, false);
  if (result != MB_SUCCESS) return result;

    // find whether this entity has ho nodes
  bool mid_nodes[3];
  MBCN::HasMidNodes(parent_type, num_parent_vertices, mid_nodes);

    // check whether this entity has mid nodes on this dimension subfacet; 
    // use dimension-1 because vertices don't have mid nodes
  if (!mid_nodes[MBCN::Dimension(subfacet_type)-1]) return MB_SUCCESS;

    // ok, we have mid nodes; now must compute expected index in connectivity array; 
    // ho nodes stored for edges, faces then entity

    // offset starts with # corner vertices
  int offset = MBCN::VerticesPerEntity(parent_type);
  int i;

  for (i = 0; i < MBCN::Dimension(subfacet_type)-1; i++)
      // for each dimension lower than that of the subfacet we're looking for, 
      // if this entity has midnodes in that dimension, increment offset by # 
      // of subfacets of that dimension; use dimension-1 in loop because 
      // canon numbering table only has 2 positions, for edges and faces;
    if (mid_nodes[i]) offset += MBCN::mConnectivityMap[parent_type][i].num_sub_elements;

    // now add the index of this subfacet; only need to if it's not the highest dimension
  if (subfacet_type != parent_type) {
    int dum, side_no, temp_offset;
    int temp_result = 
      MBCN::SideNumber(parent_conn, parent_type, subfacet_conn, 
                         MBCN::VerticesPerEntity(subfacet_type), subfacet_type,
                         side_no, dum, temp_offset);
    if(temp_result != 0) return MB_FAILURE;

    offset += side_no;
  }

    // offset shouldn't be off the end of the connectivity vector
  if (offset >= num_parent_vertices) return MB_INDEX_OUT_OF_RANGE;

  high_order_node = parent_conn[offset];

  return MB_SUCCESS;
}

  //! given an entity and a target dimension & side number, get that entity
MBErrorCode MBCore::side_element(const MBEntityHandle source_entity,
                                   const int dim, 
                                   const int side_number,
                                   MBEntityHandle &target_entity) const
{
    // get a handle on the connectivity
  const MBEntityHandle *verts;
  int num_verts;
  MBErrorCode result = get_connectivity(source_entity, verts, num_verts);
  if (MB_SUCCESS != result) return result;
  
    // get the vertices comprising the target entity
  MBRange side_verts, target_ents;
  const MBEntityType source_type = TYPE_FROM_HANDLE(source_entity);
    // first get the indices
  std::vector<int> vertex_indices;
  int temp_result = 
    MBCN::AdjacentSubEntities(source_type, &side_number, 1, dim, 0, vertex_indices);
  if (0 != temp_result) return MB_FAILURE;
    // now get the actual vertices
  for (unsigned int i = 0; i < vertex_indices.size(); i++)
    side_verts.insert(verts[vertex_indices[i]]);
  
    // now look for an entity of the correct type
    // use const_cast here 'cuz we're in a const function and we're passing 'false' for
    // create_if_missing, so we know we won't change anything
  result = (const_cast<MBCore*>(this))->get_adjacencies(side_verts, dim, false, target_ents);
  if (MB_SUCCESS != result && MB_MULTIPLE_ENTITIES_FOUND != result) return result;
  
  if (!target_ents.empty() &&
         TYPE_FROM_HANDLE(*(target_ents.begin())) != 
         MBCN::mConnectivityMap[source_type][dim-1].target_type[side_number])
    return MB_ENTITY_NOT_FOUND;

  if (!target_ents.empty()) target_entity = *(target_ents.begin());
  
  return result;
}

//-------------------------MBSet Functions---------------------//


MBErrorCode MBCore::create_meshset(const unsigned int options, 
                                     MBEntityHandle &ms_handle)
{
  int error;
  ms_handle = CREATE_HANDLE(MBENTITYSET, maxMeshSetid, error );

  //create the mesh set
  MBMeshSet *mesh_set = NULL;
  if(options & MESHSET_SET)
    mesh_set = new MBMeshSet_MBRange( ms_handle, this, a_entity_factory(), 
                                         options & MESHSET_TRACK_OWNER  );
  else if(options & MESHSET_ORDERED)
    mesh_set = new MBMeshSet_Vector(ms_handle, this, a_entity_factory(), 
                                      options & MESHSET_TRACK_OWNER );
  else {
      // unspecified - create as a set
    mesh_set = new MBMeshSet_MBRange( ms_handle, this, a_entity_factory(), 
                                         options & MESHSET_TRACK_OWNER  );
  }

  if( mesh_set == NULL )
    return MB_FAILURE;
  
  //add it to global map & update cached iter
  global_mesh_set_list.insert( std::pair<MBEntityHandle,MBMeshSet*>(ms_handle,mesh_set)); 
      
  maxMeshSetid++;

  //update cached values
  cachedMsPtr = mesh_set;
  cachedEntityHandle = ms_handle;

  return MB_SUCCESS;
}

MBErrorCode MBCore::get_meshset_options( const MBEntityHandle ms_handle, 
                                          unsigned int& options) const
{
  MBMeshSet *ms_ptr = update_cache( ms_handle );
  if (!ms_ptr) 
    return MB_ENTITY_NOT_FOUND;
  
  if (dynamic_cast<MBMeshSet_MBRange*>(ms_ptr))
    options = MESHSET_SET;
  else if (dynamic_cast<MBMeshSet_Vector*>(ms_ptr))
    options = MESHSET_ORDERED;
  else
    { assert(0); return MB_FAILURE; }
  
  if (ms_ptr->tracking())
    options |= MESHSET_TRACK_OWNER;
  
  return MB_SUCCESS;
}

MBErrorCode MBCore::clear_meshset(MBEntityHandle *ms_handles,
                                    const int num_meshsets)
{
  MBRange range;
  std::copy(ms_handles, ms_handles+num_meshsets, mb_range_inserter(range));
  return clear_meshset(range);
}

MBErrorCode MBCore::clear_meshset(MBRange &ms_handles)
{
  MBErrorCode result = MB_SUCCESS, temp_result;
  
  for (MBRange::iterator it = ms_handles.begin(); it != ms_handles.end(); it++) {
    
    MBMeshSet *ms_ptr = update_cache(*it);
    if(ms_ptr)
      temp_result = ms_ptr->clear();
    else
      temp_result = MB_ENTITY_NOT_FOUND;
    if (MB_SUCCESS != temp_result) result = temp_result;
  }

  return result;
}

MBErrorCode MBCore::subtract_meshset(MBEntityHandle meshset1, const MBEntityHandle meshset2)
{ MBMeshSet *ms_ptr1 = update_cache( meshset1 );
  if( !ms_ptr1 )
    return MB_ENTITY_NOT_FOUND;

  std::map<MBEntityHandle, MBMeshSet*>::iterator ms_iter;
  ms_iter = global_mesh_set_list.find( meshset2 );
  if( ms_iter == global_mesh_set_list.end() )
    return MB_ENTITY_NOT_FOUND;

  return ms_ptr1->subtract(ms_iter->second); 
}


MBErrorCode MBCore::intersect_meshset(MBEntityHandle meshset1, const MBEntityHandle meshset2)
{

  MBMeshSet *ms_ptr1 = update_cache( meshset1 );
  if( !ms_ptr1 )
    return MB_ENTITY_NOT_FOUND;

  std::map<MBEntityHandle, MBMeshSet*>::iterator ms_iter;
  ms_iter = global_mesh_set_list.find( meshset2 );
  if( ms_iter == global_mesh_set_list.end() )
    return MB_ENTITY_NOT_FOUND;

  return ms_ptr1->intersect(ms_iter->second); 
}

MBErrorCode MBCore::unite_meshset(MBEntityHandle meshset1, const MBEntityHandle meshset2)
{
  MBMeshSet *ms_ptr1 = update_cache( meshset1 );
  if( !ms_ptr1 )
    return MB_ENTITY_NOT_FOUND;

  std::map<MBEntityHandle, MBMeshSet*>::iterator ms_iter;
  ms_iter = global_mesh_set_list.find( meshset2 );
  if( ms_iter == global_mesh_set_list.end() )
    return MB_ENTITY_NOT_FOUND;

  return ms_ptr1->unite(ms_iter->second ); 
}

MBErrorCode MBCore::add_entities(MBEntityHandle meshset, 
                                   const MBRange &entities)
{
  if(entities.empty())
    return MB_SUCCESS;

  MBMeshSet *ms_ptr = update_cache( meshset );
  if( !ms_ptr )
    return MB_ENTITY_NOT_FOUND;

  return ms_ptr->add_entities( entities);
}

MBErrorCode MBCore::add_entities(MBEntityHandle meshset, 
                                   const MBEntityHandle *entities,
                                   const int num_entities)
{
  if(num_entities == 0)
    return MB_SUCCESS;

  MBMeshSet *ms_ptr = update_cache( meshset );
  if( !ms_ptr )
    return MB_ENTITY_NOT_FOUND;

  return ms_ptr->add_entities( entities, num_entities);
}


//! remove a range of entities from a meshset
MBErrorCode MBCore::remove_entities(MBEntityHandle meshset, 
                                      const MBRange &entities)
{
  if(entities.empty())
    return MB_SUCCESS;

  MBMeshSet *ms_ptr = update_cache( meshset );
  if( !ms_ptr )
    return MB_ENTITY_NOT_FOUND;

  return ms_ptr->remove_entities( entities);
}

//! remove a vector of entities from a meshset
MBErrorCode MBCore::remove_entities( MBEntityHandle meshset, 
                                       const MBEntityHandle *entities,
                                       const int num_entities)
{
  if(num_entities == 0)
    return MB_SUCCESS;

  MBMeshSet *ms_ptr = update_cache( meshset );
  if( !ms_ptr )
    return MB_ENTITY_NOT_FOUND;

  return ms_ptr->remove_entities( entities, num_entities );
}

MBErrorCode MBCore::get_parent_meshsets(const MBEntityHandle meshset,
                                          std::vector<MBEntityHandle> &parents,
                                          const int num_hops) const
{
  MBMeshSet *ms_ptr = update_cache( meshset );
  if( !ms_ptr )
    return MB_ENTITY_NOT_FOUND;
  
  return ms_ptr->get_parents( num_hops, parents); 
}

MBErrorCode MBCore::get_child_meshsets(const MBEntityHandle meshset,
                                         std::vector<MBEntityHandle> &children,
                                         const int num_hops) const
{
  MBMeshSet *ms_ptr = update_cache( meshset );
  if( !ms_ptr )
    return MB_ENTITY_NOT_FOUND;
  
  return ms_ptr->get_children( num_hops, children); 
}


MBErrorCode MBCore::num_parent_meshsets(const MBEntityHandle meshset, int* number) const
{

  MBMeshSet *ms_ptr = update_cache( meshset );
  if( !ms_ptr )
    return MB_ENTITY_NOT_FOUND;
 
  int error; 
  *number = ms_ptr->num_parents( &error );

  return MB_SUCCESS;
}

MBErrorCode MBCore::num_child_meshsets(const MBEntityHandle meshset, int* number ) const
{

  MBMeshSet *ms_ptr = update_cache( meshset );
  if( !ms_ptr )
    return MB_ENTITY_NOT_FOUND;
 
  int error; 
  *number = ms_ptr->num_children( &error );

  return MB_SUCCESS;
}


MBErrorCode MBCore::add_parent_meshset(MBEntityHandle meshset, 
                                         const MBEntityHandle parent_meshset)
{

  MBMeshSet *parent_ms_ptr = update_cache( parent_meshset );
  if( !parent_ms_ptr )
    return MB_ENTITY_NOT_FOUND;

  MBMeshSet *ms_ptr = update_cache( meshset );
  if( !ms_ptr )
    return MB_ENTITY_NOT_FOUND;

  ms_ptr->add_parent( parent_ms_ptr );

  return MB_SUCCESS;
}

MBErrorCode MBCore::add_child_meshset(MBEntityHandle meshset, 
                                        const MBEntityHandle child_meshset)
{

  MBMeshSet *child_ms_ptr = update_cache( child_meshset );
  if( !child_ms_ptr )
    return MB_ENTITY_NOT_FOUND;

  MBMeshSet *ms_ptr = update_cache( meshset );
  if( !ms_ptr )
    return MB_ENTITY_NOT_FOUND;

  ms_ptr->add_child( child_ms_ptr );

  return MB_SUCCESS;
}


MBErrorCode MBCore::add_parent_child(MBEntityHandle parent, 
                                       MBEntityHandle child)
{
  MBMeshSet *parent_ms_ptr = update_cache( parent );
  if( !parent_ms_ptr )
    return MB_ENTITY_NOT_FOUND;

  MBMeshSet *child_ms_ptr = update_cache( child );
  if( !child_ms_ptr )
    return MB_ENTITY_NOT_FOUND;
  
  return MBMeshSet::add_parent_child( parent_ms_ptr, child_ms_ptr);

}

MBErrorCode MBCore::remove_parent_child(MBEntityHandle parent, 
                                          MBEntityHandle child)
{
  MBMeshSet *parent_ms_ptr = update_cache( parent );
  if( !parent_ms_ptr )
    return MB_ENTITY_NOT_FOUND;

  MBMeshSet *child_ms_ptr = update_cache( child );
  if( !child_ms_ptr )
    return MB_ENTITY_NOT_FOUND;
  
  return MBMeshSet::remove_parent_child( parent_ms_ptr, child_ms_ptr);
}


MBErrorCode MBCore::remove_parent_meshset(MBEntityHandle meshset, 
                                            const MBEntityHandle parent_meshset)
{

  MBMeshSet *parent_ms_ptr = update_cache( parent_meshset );
  if( !parent_ms_ptr )
    return MB_ENTITY_NOT_FOUND;

  MBMeshSet *ms_ptr = update_cache( meshset );
  if( !ms_ptr )
    return MB_ENTITY_NOT_FOUND;

  ms_ptr->remove_parent( parent_ms_ptr );

  return MB_SUCCESS;
}

MBErrorCode MBCore::remove_child_meshset(MBEntityHandle meshset, 
                                           const MBEntityHandle child_meshset)
{

  MBMeshSet *child_ms_ptr = update_cache( child_meshset );
  if( !child_ms_ptr )
    return MB_ENTITY_NOT_FOUND;

  MBMeshSet *ms_ptr = update_cache( meshset );
  if( !ms_ptr )
    return MB_ENTITY_NOT_FOUND;

  ms_ptr->remove_child( child_ms_ptr );

  return MB_SUCCESS;
}


MBErrorCode MBCore::get_last_error(std::string& info) const
{
  return mError->get_last_error(info);
}

void MBCore::print(const MBEntityHandle ms_handle, const char *prefix,
                   bool first_call) const
{
    // get the entities
  MBRange entities;
  if (0 != ms_handle) {
    MBMeshSet *ms_ptr = update_cache( ms_handle );
    if( !ms_ptr )
      return;

    ms_ptr->get_entities(entities, false);
    if (!first_call)
      std::cout << prefix << "MBENTITYSET " << ID_FROM_HANDLE(ms_handle) 
                << std::endl;
  }
  else {
    get_entities_by_dimension(0, 3, entities);
    if (entities.empty()) get_entities_by_dimension(0, 2, entities);
    if (entities.empty()) get_entities_by_dimension(0, 1, entities);
    get_entities_by_dimension(0, 0, entities);
    get_entities_by_type(0, MBENTITYSET, entities);
    std::cout << prefix << "--MB: " << std::endl;
  }
    
  std::string indent_prefix = prefix;
  indent_prefix += "  ";
  for (MBRange::iterator it = entities.begin(); 
       it != entities.end(); it++) 
  {
    if (TYPE_FROM_HANDLE(*it) == MBENTITYSET) {
      print(*it, indent_prefix.c_str(), false);
    }
    else if (first_call) {
      std::cout << prefix << MBCN::EntityTypeName(TYPE_FROM_HANDLE(*it)) << " " 
                << ID_FROM_HANDLE(*it) << std::endl;
    }
  }
}

MBMeshSet* MBCore::update_cache( const MBEntityHandle ms_handle ) const
{
  
  //update cached values
  if( ms_handle != cachedEntityHandle )
  {
    std::map<MBEntityHandle, MBMeshSet*>::const_iterator ms_iter
      = global_mesh_set_list.find( ms_handle ); 

    if( ms_iter == global_mesh_set_list.end() )
    {
      MBMeshSet *ms_ptr = NULL;
      return ms_ptr;
    }

    cachedEntityHandle = ms_iter->first;
    cachedMsPtr = ms_iter->second;
    return cachedMsPtr;
  }

  return  cachedMsPtr;

}

//! define non-inline versions of these functions for debugging
MBEntityHandle ifh(MBEntityHandle handle) 
{
  return ID_FROM_HANDLE(handle);
}

MBEntityType tfh(MBEntityHandle handle)
{
  return TYPE_FROM_HANDLE(handle);
}


