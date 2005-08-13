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
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif


#include "EntitySequence.hpp"
#include "EntitySequenceManager.hpp"
#include "MBCN.hpp"
#include <assert.h>
#include "AEntityFactory.hpp" 

using namespace std;

//Constructor
MBEntitySequence::MBEntitySequence(EntitySequenceManager* seq_manager,
  MBEntityHandle start_handle, int num_entities )
  : mSequenceManager(seq_manager)
{
  assert(MB_START_ID <= ID_FROM_HANDLE(start_handle));
  assert(TYPE_FROM_HANDLE(start_handle) < MBMAXTYPE);

  // set first entity handle in sequence
  mStartEntityHandle = start_handle;

  // set number of entity handles in sequence
  mNumAllocated = num_entities;

  // no entities created yet
  mNumEntities = 0;

  mLastDeletedIndex = -1;
}


MBEntitySequence::~MBEntitySequence()
{ }



VertexEntitySequence::VertexEntitySequence(EntitySequenceManager* seq_manager,
                                           MBEntityHandle start_handle, int num_entities,
                                           bool all_handles_used)
: MBEntitySequence(seq_manager, start_handle, num_entities)
{
  seq_manager->entity_sequence_created(this);

  // allocate the arrays
  mCoords = new double* [3];
  mCoords[0] = new double[num_entities];
  mCoords[1] = new double[num_entities];
  mCoords[2] = new double[num_entities];
  
  if(all_handles_used)
  {
    mNumEntities = num_entities;
    mFirstFreeIndex = -1;
    std::vector<bool>(mNumAllocated, false).swap(mFreeEntities);
  }
  else
  {
    seq_manager->notify_not_full(this);
    std::vector<bool>(mNumAllocated, true).swap(mFreeEntities);
    mNumEntities = 0;
    mFirstFreeIndex = 0;
    for(int i=0; i<num_entities; i++)
    {
      reinterpret_cast<int&>(mCoords[0][i]) = i+1;
    }
    reinterpret_cast<int&>(mCoords[0][num_entities-1]) = -1;
  }
}

VertexEntitySequence::~VertexEntitySequence()
{
  mSequenceManager->entity_sequence_deleted(this);
  delete [] mCoords[0];
  delete [] mCoords[1];
  delete [] mCoords[2];
  delete [] mCoords;
}


MBEntityHandle VertexEntitySequence::get_unused_handle()
{
  if(mFirstFreeIndex == -1)
    return 0;

  MBEntityHandle new_handle = mStartEntityHandle + mFirstFreeIndex;
  mFreeEntities[mFirstFreeIndex] = false;

  mFirstFreeIndex = reinterpret_cast<int&>(mCoords[0][mFirstFreeIndex]);

  mNumEntities++;

  if(mNumEntities == mNumAllocated)
    mSequenceManager->notify_full(this);

  if( mLastDeletedIndex == (int)( new_handle - mStartEntityHandle ))
    mLastDeletedIndex = -1;

  return new_handle;
}

void VertexEntitySequence::free_handle(MBEntityHandle handle)
{

  if(!is_valid_entity(handle))
    return;

  if(mNumEntities == mNumAllocated)
    mSequenceManager->notify_not_full(this);

  mFreeEntities[handle - mStartEntityHandle] = true;

  int prev_free_handle = -1;
  const int handle_index = handle - mStartEntityHandle;
  
  int index = mFirstFreeIndex;
  if( mLastDeletedIndex != -1 && handle_index > mLastDeletedIndex)
    index = mLastDeletedIndex;

  mLastDeletedIndex = handle_index;

  for(;
      (index != -1) && (index < handle_index); 
      index = reinterpret_cast<int&>(mCoords[0][index]) )
  {
    prev_free_handle = index;
  }

  // was full and now has one free handle
  if(prev_free_handle == -1 && mFirstFreeIndex == -1)
  {
    mFirstFreeIndex = handle_index;
    reinterpret_cast<int&>(mCoords[0][handle_index]) = -1;
  }
  // insert before all free handles
  else if(prev_free_handle == -1)
  {
    reinterpret_cast<int&>(mCoords[0][handle_index]) = mFirstFreeIndex;
    mFirstFreeIndex = handle_index;
  }
  // insert in middle or end
  else
  {
    reinterpret_cast<int&>(mCoords[0][handle_index]) = 
      reinterpret_cast<int&>(mCoords[0][prev_free_handle]);
    reinterpret_cast<int&>(mCoords[0][prev_free_handle]) = handle_index;
  }

  mNumEntities--;

}

bool VertexEntitySequence::is_valid_entity(MBEntityHandle entity) const
{
  return !mFreeEntities[entity - mStartEntityHandle];
}


void VertexEntitySequence::get_entities(MBRange& entities) const
{
  MBRange::iterator iter = entities.insert(mStartEntityHandle, mStartEntityHandle+mNumAllocated-1);
  
  for(int index = mFirstFreeIndex; (index != -1); 
      index = reinterpret_cast<int&>(mCoords[0][index]))
  {
    for(; *iter != (mStartEntityHandle+index); ++iter);
    entities.erase(iter);
  }

}



ElementEntitySequence::ElementEntitySequence(EntitySequenceManager* seq_manager,
                        MBEntityHandle start_handle, int num_entities,
                        int nodes_per_element, bool all_handles_used,
                                             bool allocate_connect)
: MBEntitySequence(seq_manager, start_handle, num_entities)
{
  assert(nodes_per_element < 28);
  if (allocate_connect) {
    mElements = new MBEntityHandle[num_entities*nodes_per_element];
    memset(mElements, 0, sizeof(MBEntityHandle)*num_entities*nodes_per_element);
  }
  else {
    mElements = NULL;
  }
  
  mNodesPerElement = nodes_per_element;

  seq_manager->entity_sequence_created(this);

  if(all_handles_used)
  {
    mNumEntities = num_entities;
    std::vector<bool>(mNumAllocated, false).swap(mFreeEntities);
    mFirstFreeIndex = -1;
  }
  else
  {
    seq_manager->notify_not_full(this);
    std::vector<bool>(mNumAllocated, true).swap(mFreeEntities);
    mNumEntities = 0;
    mFirstFreeIndex = 0;
    if (nodes_per_element)
    {
      int max = (num_entities-1)*nodes_per_element;
      for(int i=0; i<max ; i+=nodes_per_element)
      {
        mElements[i] = (i/nodes_per_element)+1;
      }
      reinterpret_cast<int&>(mElements[max]) = -1;
    }
  }
}


ElementEntitySequence::~ElementEntitySequence()
{
  mSequenceManager->entity_sequence_deleted(this);
  if (NULL != mElements)
    delete [] mElements;
}

MBEntityHandle ElementEntitySequence::get_unused_handle()
{
  if(mFirstFreeIndex == -1)
    return 0;

  mFreeEntities[mFirstFreeIndex] = false;
  MBEntityHandle new_handle = mStartEntityHandle + mFirstFreeIndex;

  mFirstFreeIndex = reinterpret_cast<int&>(mElements[mFirstFreeIndex*mNodesPerElement]);

  mNumEntities++;

  if(mNumEntities == mNumAllocated)
    mSequenceManager->notify_full(this);
  
  if( mLastDeletedIndex == (int)( new_handle - mStartEntityHandle ))
    mLastDeletedIndex = -1;

  return new_handle;
}

void ElementEntitySequence::free_handle(MBEntityHandle handle)
{
  if(!is_valid_entity(handle))
    return;

  if(mNumEntities == mNumAllocated)
    mSequenceManager->notify_not_full(this);

  mFreeEntities[handle-mStartEntityHandle] = true;
  
  int prev_free_index = -1;
  const int handle_index = handle - mStartEntityHandle;

  int index = mFirstFreeIndex;
  if( mLastDeletedIndex != -1 && handle_index > mLastDeletedIndex)
    index = mLastDeletedIndex;

  mLastDeletedIndex = handle_index;
  
  for(;
      (index != -1) && (index < handle_index); 
      index = reinterpret_cast<int&>(mElements[index*mNodesPerElement]) )
  {
    prev_free_index = index;
  }
  
  memset(&mElements[handle_index*mNodesPerElement], 0, sizeof(MBEntityHandle)*mNodesPerElement);

  // was full and now has one free handle
  if(prev_free_index == -1 && mFirstFreeIndex == -1)
  {
    mFirstFreeIndex = handle_index;
    reinterpret_cast<int&>(mElements[handle_index*mNodesPerElement]) = -1;
  }
  // insert before all free handles
  else if(prev_free_index == -1)
  {
    reinterpret_cast<int&>(mElements[handle_index*mNodesPerElement]) = mFirstFreeIndex;
    mFirstFreeIndex = handle_index;
  }
  // insert in middle or end
  else
  {
    mElements[handle_index*mNodesPerElement] = mElements[prev_free_index*mNodesPerElement];
    reinterpret_cast<int&>(mElements[prev_free_index*mNodesPerElement]) = handle_index;
  }

  mNumEntities--;

}



void ElementEntitySequence::get_entities(MBRange& entities) const
{
  MBRange::iterator iter = entities.insert(mStartEntityHandle, mStartEntityHandle+mNumAllocated-1);
  
  for(int index = mFirstFreeIndex; index != -1;
      index = reinterpret_cast<int&>(mElements[index*mNodesPerElement]) )
  {
    for(; *iter != (mStartEntityHandle+index); ++iter);
    entities.erase(iter);
  }
}



MBErrorCode ElementEntitySequence::split(MBEntityHandle split_location, 
    MBEntitySequence*& new_sequence)
{
  // only split (begin, end)
  assert(mStartEntityHandle < split_location);
  assert(get_end_handle() >= split_location);

  // make a new sequence
  ElementEntitySequence* seq = new ElementEntitySequence( mSequenceManager, split_location, 
      get_end_handle() - split_location + 1 , mNodesPerElement, false);
  new_sequence = seq;

  // copy data into new sequence
  memcpy(seq->mElements, &mElements[mNodesPerElement*(split_location - mStartEntityHandle)],
         seq->mNumAllocated*mNodesPerElement*sizeof(MBEntityHandle));

  //copy free handles over too
  std::copy(mFreeEntities.begin()+(split_location-mStartEntityHandle), mFreeEntities.end(), 
      seq->mFreeEntities.begin());
 
  // make a new shorter array for this sequence and copy data over 
  mNumAllocated = split_location - mStartEntityHandle;
  MBEntityHandle* tmp = new MBEntityHandle[mNumAllocated*mNodesPerElement];
  memcpy(tmp, mElements, mNumAllocated*mNodesPerElement*sizeof(MBEntityHandle));
  delete [] mElements;
  mElements = tmp;

  // shrink capacity to what we need
  mFreeEntities.resize(mNumAllocated);
  std::vector<bool>(mFreeEntities).swap(mFreeEntities);


  // need to recompute : mNumEntities, mFirstFreeIndex for both sequences
  
  mNumEntities = 0;
  mFirstFreeIndex = -1;
  
  std::vector<bool>::reverse_iterator iter = mFreeEntities.rbegin();
  std::vector<bool>::reverse_iterator end_iter = mFreeEntities.rend();
  int index = mFreeEntities.size() - 1;
  int last_index = -1;

  for(; iter != end_iter; )
  {
    if(*iter == true)
    {
      reinterpret_cast<int&>(mElements[index*mNodesPerElement]) = last_index;
      last_index = index;
    }
    else
    {
      ++mNumEntities;
    }

    ++iter;
    --index;
  }
  mFirstFreeIndex = last_index;

  if(mNumEntities == mNumAllocated)
    mSequenceManager->notify_full(this);
  else
    mSequenceManager->notify_not_full(this);
  
  
  seq->mNumEntities = 0;
  seq->mFirstFreeIndex = -1;
  
  iter = seq->mFreeEntities.rbegin();
  end_iter = seq->mFreeEntities.rend();
  index = seq->mFreeEntities.size() - 1;
  last_index = -1;

  for(; iter != end_iter; )
  {
    if(*iter == true)
    {
      reinterpret_cast<int&>(seq->mElements[index*mNodesPerElement]) = last_index;
      last_index = index;
    }
    else
    {
      ++seq->mNumEntities;
    }

    ++iter;
    --index;
  }
  seq->mFirstFreeIndex = last_index;
  
  if(seq->mNumEntities == seq->mNumAllocated)
    mSequenceManager->notify_full(seq);
  else
    mSequenceManager->notify_not_full(seq);

  return MB_SUCCESS;
}



MBErrorCode ElementEntitySequence::convert_realloc(bool& mid_edge_nodes,
    bool& mid_face_nodes, bool& mid_volume_nodes, MBCore* MB, 
    MBTag delete_mark_bit )
{

  MBEntityType this_type = get_type();

  // figure out how many nodes per element we'll end up with
  int num_corner_nodes = MBCN::VerticesPerEntity(this_type);
  int new_nodes_per_element = num_corner_nodes;
  if(mid_edge_nodes)
    new_nodes_per_element += MBCN::mConnectivityMap[this_type][0].num_sub_elements;
  if(mid_face_nodes)
    new_nodes_per_element += MBCN::mConnectivityMap[this_type][1].num_sub_elements;
  if(mid_volume_nodes)
    new_nodes_per_element++;

  if(new_nodes_per_element == mNodesPerElement)
    return MB_SUCCESS;

  // make the new array
  MBEntityHandle* new_array = new MBEntityHandle[mNumAllocated*new_nodes_per_element];
  memset(new_array, 0, sizeof(MBEntityHandle)*mNumAllocated*new_nodes_per_element);

  // copy the corner nodes
  MBEntityHandle* old_end_handle = mElements + mNodesPerElement*mNumAllocated;
  MBEntityHandle* old_iter = mElements;
  MBEntityHandle* new_iter = new_array;
  int amount_to_copy = num_corner_nodes*sizeof(MBEntityHandle);
  for(; old_iter < old_end_handle; )
  {
    memcpy(new_iter, old_iter, amount_to_copy);
    old_iter += mNodesPerElement;
    new_iter += new_nodes_per_element;
  }

  // copy mid edge nodes if they existed and we want to keep them
  if(mid_edge_nodes && has_mid_edge_nodes())
  {
    amount_to_copy = MBCN::mConnectivityMap[this_type][0].num_sub_elements * sizeof(MBEntityHandle);
    old_end_handle = mElements + mNodesPerElement*mNumAllocated;
    old_iter = mElements + MBCN::VerticesPerEntity(this_type);
    new_iter = new_array + MBCN::VerticesPerEntity(this_type);
    for(; old_iter < old_end_handle; )
    {
      memcpy(new_iter, old_iter, amount_to_copy);
      old_iter += mNodesPerElement;
      new_iter += new_nodes_per_element;
    }
  }
  // need to delete them
  else if(has_mid_edge_nodes())
  {
    int number_to_delete = MBCN::mConnectivityMap[this_type][0].num_sub_elements;
    old_end_handle = mElements + mNodesPerElement*mNumAllocated;
    old_iter = mElements + MBCN::VerticesPerEntity(this_type);

    std::vector<MBEntityHandle> nodes_processed; 
    std::vector<MBEntityHandle>::iterator tmp_iter; 
    
    for(; old_iter < old_end_handle; )
    {
      //tag each node-to-delete with parent element entity handle and number 
      //of elements it is found on
      for(int i=0; i<number_to_delete; i++)
      {
        //see if node has been processed yet
        tmp_iter = std::find( nodes_processed.begin(), nodes_processed.end(), old_iter[i] );
        if( tmp_iter == nodes_processed.end() && old_iter[i] != 0 )
        {
          nodes_processed.push_back( old_iter[i] );

          //determines if node should be deleted or not
          //(makes sure it's not on other entities besides those in this sequence) 
          if( tag_for_deletion( old_iter[i], MB ) )
          {
            //tag node as deletable
            unsigned char bit = 0x1;
            MB->tag_set_data(delete_mark_bit, &(old_iter[i]), 1, &bit);
          }
        }
      }
      old_iter += mNodesPerElement;
    }
  }

  // copy mid face nodes if they existed and we want to keep them
  if(mid_face_nodes && has_mid_face_nodes())
  {
    amount_to_copy = MBCN::mConnectivityMap[this_type][1].num_sub_elements * sizeof(MBEntityHandle);
    old_end_handle = mElements + mNodesPerElement*mNumAllocated;
    old_iter = mElements + MBCN::VerticesPerEntity(this_type);
    if(has_mid_edge_nodes())
      old_iter += MBCN::mConnectivityMap[this_type][0].num_sub_elements;
    new_iter = new_array + MBCN::VerticesPerEntity(this_type);
    if(mid_edge_nodes)
      new_iter += MBCN::mConnectivityMap[this_type][0].num_sub_elements;

    for(; old_iter < old_end_handle; )
    {
      memcpy(new_iter, old_iter, amount_to_copy);
      old_iter += mNodesPerElement;
      new_iter += new_nodes_per_element;
    }
  }
  // need to delete them
  else if(has_mid_face_nodes())
  {
    int number_to_delete = MBCN::mConnectivityMap[this_type][1].num_sub_elements;
    old_end_handle = mElements + mNodesPerElement*mNumAllocated;
    old_iter = mElements + MBCN::VerticesPerEntity(this_type);
    if(has_mid_edge_nodes())
      old_iter+=MBCN::mConnectivityMap[this_type][0].num_sub_elements;

    std::vector<MBEntityHandle> nodes_processed; 
    std::vector<MBEntityHandle>::iterator tmp_iter; 

    for(; old_iter < old_end_handle; )
    {
      //tag each node-to-delete with parent element entity handle and number 
      //of elements it is found on
      for(int i=0; i<number_to_delete; i++)
      {
        //see if node has been processed yet
        tmp_iter = std::find( nodes_processed.begin(), nodes_processed.end(), old_iter[i] );
        if( tmp_iter == nodes_processed.end() && old_iter[i] != 0 )
        {
          nodes_processed.push_back( old_iter[i] );

          //determines if node should be deleted or not
          //(makes sure it's not on other entities besides those in this sequence) 
          if( tag_for_deletion( old_iter[i], MB ) )
          {
            //tag node as deletable
            unsigned char bit = 0x1;
            MB->tag_set_data(delete_mark_bit, &(old_iter[i]), 1, &bit);
          }
        }
      }
      old_iter += mNodesPerElement;
    }

  }

  // copy mid volume nodes if they existed and we want to keep them
  if(mid_volume_nodes && has_mid_volume_nodes())
  {
    old_end_handle = mElements + mNodesPerElement*mNumAllocated;
    old_iter = mElements + MBCN::VerticesPerEntity(this_type);
    if(has_mid_edge_nodes())
      old_iter += MBCN::mConnectivityMap[this_type][0].num_sub_elements;
    if(has_mid_face_nodes())
      old_iter += MBCN::mConnectivityMap[this_type][1].num_sub_elements;
    new_iter = new_array + (new_nodes_per_element - 1);
   
    for(; old_iter < old_end_handle; )
    {
      *new_iter = *old_iter;
      old_iter += mNodesPerElement;
      new_iter += new_nodes_per_element;
    }
  }
  // need to delete them
  else if(has_mid_volume_nodes())
  {
    old_end_handle = mElements + mNodesPerElement*mNumAllocated;
    old_iter = mElements + MBCN::VerticesPerEntity(this_type);
    if(has_mid_edge_nodes())
      old_iter += MBCN::mConnectivityMap[this_type][0].num_sub_elements;
    if(has_mid_face_nodes())
      old_iter += MBCN::mConnectivityMap[this_type][1].num_sub_elements;
   
    for(; old_iter < old_end_handle; )
    {
      if( *old_iter != 0 && tag_for_deletion( *old_iter, MB ) )
      {
        //tag node as deletable
        unsigned char bit = 0x1;
        MB->tag_set_data(delete_mark_bit, &(*old_iter), 1, &bit);
      }
      old_iter += mNodesPerElement;
    }

  }

  // did we actually allocate new space for these?
  if(mid_edge_nodes)
    mid_edge_nodes = !has_mid_edge_nodes();
  if(mid_face_nodes)
    mid_face_nodes = !has_mid_face_nodes();
  if(mid_volume_nodes)
    mid_volume_nodes = !has_mid_volume_nodes();
  
  // swap arrays
  delete [] mElements;
  mElements = new_array;

  mNodesPerElement = new_nodes_per_element;
  
  return MB_SUCCESS;
}

bool ElementEntitySequence::has_mid_edge_nodes() const
{
  return MBCN::HasMidEdgeNodes(get_type(), mNodesPerElement);
}

bool ElementEntitySequence::has_mid_face_nodes() const
{
  return MBCN::HasMidFaceNodes(get_type(), mNodesPerElement);
}

bool ElementEntitySequence::has_mid_volume_nodes() const
{
  return MBCN::HasMidRegionNodes(get_type(), mNodesPerElement);
}

bool ElementEntitySequence::tag_for_deletion( MBEntityHandle &node, MBCore *MB ) 
{
  //get type of this sequence
  MBEntityType this_type = get_type();

  //find 'parent' element (in this sequence)
  MBEntityHandle parent_handle = 0;
  parent_handle = ( ( &node - mElements ) / mNodesPerElement ) + mStartEntityHandle;

  //get dimension of 'parent' element
  int this_dimension = MB->dimension_from_handle( parent_handle );

  //get connectivity of 'parent' element
  std::vector<MBEntityHandle> connectivity;
  MB->get_connectivity( &parent_handle, 1, connectivity);

  //tells us if higher order node is on 
  int dimension, side_number; 
  MBCN::HONodeParent( &(connectivity[0]), this_type, connectivity.size(),
                        &node, dimension, side_number );  

  //it MUST be a higher-order node
  bool delete_node = false;

  assert( dimension != -1 );
  assert( side_number != -1 );

  //could be a mid-volume/face/edge node on a hex/face/edge respectively
  //if so...delete it bc/ no one else owns it too
  if( dimension == this_dimension && side_number == 0 )
    delete_node = true;
  else //the node could also be on a lower order entity of 'tmp_entity' 
  {
    //get 'side' of 'parent_handle' that node is on 
    MBEntityHandle target_entity = 0;
    MB->side_element( parent_handle, dimension, side_number, target_entity );

    if( target_entity )
    {
      AEntityFactory *a_fact = MB->a_entity_factory();
      MBEntityHandle low_meshset;
      int dum;
      low_meshset = CREATE_HANDLE(MBENTITYSET, 0, dum);

      //just get corner nodes of target_entity
      connectivity.clear();
      MB->get_connectivity(&( target_entity), 1, connectivity, true  );

      //for each node, get all common adjacencies of nodes in 'parent_handle' 
      std::vector<MBEntityHandle> adj_list_1, adj_list_2, adj_entities;
      a_fact->get_adjacencies(connectivity[0], adj_list_1);

      // remove meshsets
      adj_list_1.erase(std::remove_if(adj_list_1.begin(), adj_list_1.end(), 
           std::bind2nd(std::greater<MBEntityHandle>(),low_meshset)), adj_list_1.end());

      unsigned int i; 
      for( i=1; i<connectivity.size(); i++)
      {
        adj_list_2.clear();
        a_fact->get_adjacencies(connectivity[i], adj_list_2);

        // remove meshsets
        adj_list_2.erase(std::remove_if(adj_list_2.begin(), adj_list_2.end(), 
             std::bind2nd(std::greater<MBEntityHandle>(),low_meshset)), adj_list_2.end());
       
        //intersect the 2 lists 
        adj_entities.clear();
        std::set_intersection(adj_list_1.begin(), adj_list_1.end(), 
                              adj_list_2.begin(), adj_list_2.end(), 
                              std::back_inserter< std::vector<MBEntityHandle> >(adj_entities));
        adj_list_1.clear();
        adj_list_1 = adj_entities;
      } 

      assert( adj_entities.size() );  //has to have at least one adjacency 

      //see if node is in other elements, not in this sequence...if so, delete it 
      for( i=0; i<adj_entities.size(); i++)
      {
        if( adj_entities[i] >= get_start_handle() &&
            adj_entities[i] <= get_end_handle() )
        {
          delete_node = false;
          break;
        }
        else 
          delete_node = true;
      }             
    }
    else //there is no lower order entity that also contains node 
      delete_node = true;
  }

  return delete_node;

}
