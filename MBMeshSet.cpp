//
//-------------------------------------------------------------------------
// Filename      : MBMeshSet.cpp 
// Creator       : Tim Tautges
//
// Date          : 02/01/02
//
// Owner         : Tim Tautges
//
// Description   : MBMeshSet is the MB implementation of MBSet
//-------------------------------------------------------------------------

#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif



#include <algorithm>
#include <utility>
#include <set>
#include "assert.h"
#include "MBInternals.hpp"
#include "MBMeshSet.hpp"
#include "AEntityFactory.hpp"

using namespace std;

MBMeshSet::MBMeshSet( MBEntityHandle handle, MBInterface* mdb, AEntityFactory* a_fact, bool track_ownership ) 
{
  mEntityHandle = handle; 
  mAdjFact = a_fact;
  mTracking = track_ownership;
  mMB = mdb;
}

MBMeshSet::~MBMeshSet() 
{
}


MBErrorCode MBMeshSet::get_children(const int num_hops, 
    std::vector<MBEntityHandle> &children) const
{
    // compute new value of num_hops, either decremented or equal to -1
  int this_hops = (-1 == num_hops ? -1 : num_hops-1);

  std::vector<MBMeshSet*>::const_iterator i = childMeshSets.begin();
    // go through child sets; if we find a unique child, and we haven't exceeded
    // the hops, recurse
  while (i != childMeshSets.end()) 
  {
    std::vector<MBEntityHandle>::iterator j = std::find(children.begin(), children.end(), (*i)->mEntityHandle);
    if (j == children.end()) 
    {
      children.push_back((*i)->mEntityHandle);
      if (0 != num_hops)
        (*i)->get_children(this_hops, children);
    }
    i++;
  }

  return MB_SUCCESS;
}


MBErrorCode MBMeshSet::get_parents(const int num_hops, 
    std::vector<MBEntityHandle> &parents) const
{
    // compute new value of num_hops, either decremented or equal to -1
  int this_hops = (-1 == num_hops ? -1 : num_hops-1);

  std::vector<MBMeshSet*>::const_iterator i = parentMeshSets.begin();
    // go through parent sets; if we find a unique parent, and we haven't exceeded
    // the hops, recurse
  while (i != parentMeshSets.end()) 
  {
    std::vector<MBEntityHandle>::iterator j = std::find(parents.begin(), parents.end(), (*i)->mEntityHandle );
    if (j == parents.end()) 
    {
      parents.push_back( (*i)->mEntityHandle );
      if (0 != num_hops)
        (*i)->get_parents(this_hops, parents);
    }
    i++;
  }

  return MB_SUCCESS;
}

  //! add a parent/child link between the meshsets; returns error if entities are already
  //! related or if child is already a parent of parent
MBErrorCode MBMeshSet::add_parent_child(MBMeshSet* parent_meshset, 
                                 MBMeshSet* child_meshset) 
{ 
  parent_meshset->add_child(child_meshset);
  child_meshset->add_parent(parent_meshset);
  return MB_SUCCESS;
}

  //! remove a parent/child link between the meshsets; returns error if entities
  //! are not related
MBErrorCode MBMeshSet::remove_parent_child(MBMeshSet* parent_meshset, 
                                    MBMeshSet* child_meshset) 
{ 
  parent_meshset->remove_child(child_meshset);
  child_meshset->remove_parent(parent_meshset);
  return MB_SUCCESS;
}


int MBMeshSet::add_parent(MBMeshSet *parent) 
{
  //add parent to this's parent-mesh-set-list
  if (std::find(parentMeshSets.begin(),
                parentMeshSets.end(), parent) == parentMeshSets.end()) 
  {
    parentMeshSets.push_back(parent);
    return 1; 
  }

  else 
    return 0;

}

int MBMeshSet::add_child(MBMeshSet *child) 
{
  //add child to this's child-mesh-set-list
  if (std::find(childMeshSets.begin(),
                childMeshSets.end(), child) == childMeshSets.end()) 
  {
    childMeshSets.push_back(child);
    return 1; 
  }

  else 
    return 0;
}

int MBMeshSet::remove_parent(MBMeshSet *parent) 
{
  //erace position from this's parent-mesh-set-list
  std::vector<MBMeshSet*>::iterator position = 
    std::find(parentMeshSets.begin(), parentMeshSets.end(), parent);

  if ( position != parentMeshSets.end()) 
  {
    parentMeshSets.erase(position);
    return 1; 
  }

  else 
    return 0;
}

int MBMeshSet::remove_child(MBMeshSet *child) 
{
  //erace position from this's child-mesh-set-list
  std::vector<MBMeshSet*>::iterator position = 
    std::find(childMeshSets.begin(), childMeshSets.end(), child);
  
  if ( position != childMeshSets.end()) 
  {
    childMeshSets.erase(position);
    return 1; 
  }

  else 
    return 0;
}

  //! return the number of child/parent relations for this meshset
int MBMeshSet::num_children(int *) const
{ 
  return childMeshSets.size();
}

int MBMeshSet::num_parents(int *) const
{ 
  return parentMeshSets.size();
}


MBMeshSet_MBRange::~MBMeshSet_MBRange()
{
  // clean up usage tracking  
  if(mTracking && mAdjFact)
  {
    for(MBRange::iterator iter = mRange.begin();
        iter != mRange.end(); ++iter)
    {
      mAdjFact->remove_adjacency(*iter, mEntityHandle);
    }
  }
}

MBErrorCode MBMeshSet_MBRange::clear()
{
  if(mTracking && mAdjFact)
  {
    for(MBRange::iterator iter = mRange.begin();
        iter != mRange.end(); ++iter)
    {
      mAdjFact->remove_adjacency(*iter, mEntityHandle);
    }
  }
  mRange.clear();
  return MB_SUCCESS;
}

MBErrorCode MBMeshSet_MBRange::get_entities(std::vector<MBEntityHandle>& entity_list,
                                                const bool recursive) const
{
  for(MBRange::const_iterator iter = mRange.begin();
      iter != mRange.end(); ++iter)
  {
    if (recursive && TYPE_FROM_HANDLE(*iter) == MBENTITYSET) {
      if (NULL != mMB) mMB->get_entities_by_handle(*iter, entity_list, recursive);
    }
    else
      entity_list.push_back(*iter);
  }
  return MB_SUCCESS;
}

MBErrorCode MBMeshSet_MBRange::get_entities(MBRange& entity_list,
                                                const bool recursive) const
{
    
  for(MBRange::const_iterator iter = mRange.begin();
      iter != mRange.end(); ++iter)
  {
    if (recursive && TYPE_FROM_HANDLE(*iter) == MBENTITYSET) {
      if (NULL != mMB) mMB->get_entities_by_handle(*iter, entity_list, recursive);
    }
    else
      entity_list.insert(*iter);
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBMeshSet_MBRange::get_entities_by_type(MBEntityType type,
    std::vector<MBEntityHandle>& entity_list) const
{
  MBRange::pair_iterator iter = mRange.begin(), end_iter = mRange.end();
  
  for(; iter != end_iter; ++iter)
    if(TYPE_FROM_HANDLE(iter->first) == type)
      break;

  for(; iter != end_iter; ++iter)
  {
    if (TYPE_FROM_HANDLE(iter->first) != type) break;
    for(MBEntityHandle curr_handle = iter->first;
        curr_handle <= iter->second; ++curr_handle)
    {
      assert(TYPE_FROM_HANDLE(curr_handle) == type);
      entity_list.push_back(curr_handle);
    }
  }

  return MB_SUCCESS;

}

MBErrorCode MBMeshSet_MBRange::get_entities_by_type(MBEntityType type,
    MBRange& entity_list) const
{

  MBRange::pair_iterator iter = mRange.begin(), end_iter = mRange.end();
  
  for(; iter != end_iter; ++iter)
    if(TYPE_FROM_HANDLE(iter->first) == type)
      break;

  for(; iter != end_iter; ++iter)
    if(TYPE_FROM_HANDLE(iter->first) == type)
      entity_list.insert(iter->first, iter->second);
    else
      break;
  
  return MB_SUCCESS;

}

MBErrorCode MBMeshSet_MBRange::add_entities(const MBEntityHandle *entities,
                                                const int num_entities)
{
  int i;
  for(i = 0; i < num_entities; i++)
    mRange.insert(entities[i]);

  if(mTracking && mAdjFact)
  {
    for(i = 0; i < num_entities; i++)
      mAdjFact->add_adjacency(entities[i], mEntityHandle);
  }
  return MB_SUCCESS;
}


MBErrorCode MBMeshSet_MBRange::add_entities(
    const MBRange& entities)
{

  mRange.merge(entities);
  
  if(mTracking && mAdjFact)
  {
    for(MBRange::const_iterator iter = entities.begin();
        iter != entities.end(); ++iter)
      mAdjFact->add_adjacency(*iter, mEntityHandle);
  }

  return MB_SUCCESS;
}

MBErrorCode MBMeshSet_MBRange::remove_entities(
    const MBRange& entities)
{
  MBRange::const_iterator iter = entities.begin();
  for(; iter != entities.end(); ++iter)
  {
    MBRange::iterator found = mRange.find(*iter);
    if(found != mRange.end())
    {
      mRange.erase(found);
      if(mTracking && mAdjFact)
        mAdjFact->remove_adjacency(*iter, mEntityHandle);
    }
  }

  return MB_SUCCESS;

}

MBErrorCode MBMeshSet_MBRange::remove_entities(const MBEntityHandle *entities,
                                                   const int num_entities)
{
  for(int i = 0; i < num_entities; i++)
  {
    MBRange::iterator found = mRange.find(entities[i]);
    if(found != mRange.end())
    {
      mRange.erase(found);
      if(mTracking && mAdjFact)
        mAdjFact->remove_adjacency(entities[i], mEntityHandle);
    }
  }

  return MB_SUCCESS;

}


unsigned int MBMeshSet_MBRange::num_entities(int*) const
{
  return mRange.size();
}

unsigned int MBMeshSet_MBRange::num_entities_by_type(MBEntityType type) const
{
  type = type;
  return 0;
}


MBErrorCode MBMeshSet_MBRange::subtract(const MBMeshSet *meshset_2)
{
  MBRange other_range;
  meshset_2->get_entities(other_range, false);

  return remove_entities(other_range);
}

MBErrorCode MBMeshSet_MBRange::intersect(const MBMeshSet *meshset_2)
{
  MBRange other_range;
  meshset_2->get_entities(other_range, false);

  std::set<MBEntityHandle> tmp;

  std::set_intersection(mRange.begin(), mRange.end(), other_range.begin(), other_range.end(), 
      std::inserter< std::set<MBEntityHandle> >(tmp, tmp.end()));

  mRange.clear();

  for(std::set<MBEntityHandle>::reverse_iterator iter = tmp.rbegin();
      iter != tmp.rend(); ++iter)
  {
    mRange.insert(*iter);
  }

  //track owner
  if(mTracking && mAdjFact)
  {
    tmp.clear();
    std::set_difference(mRange.begin(), mRange.end(), other_range.begin(), other_range.end(),
        std::inserter< std::set<MBEntityHandle> >(tmp, tmp.end()));

    for(std::set<MBEntityHandle>::iterator iter = tmp.begin();
        iter != tmp.end(); ++iter)
      mAdjFact->remove_adjacency(*iter, mEntityHandle);
  }

  return MB_SUCCESS;
}

MBErrorCode MBMeshSet_MBRange::unite(const MBMeshSet *meshset_2)
{
  MBRange other_range;
  meshset_2->get_entities(other_range, false);
  return add_entities(other_range);
}




MBMeshSet_Vector::~MBMeshSet_Vector()
{
  // clean up usage tracking  
  if(mTracking && mAdjFact)
  {
    for(std::vector<MBEntityHandle>::iterator iter = mVector.begin();
        iter != mVector.end(); ++iter)
    {
      mAdjFact->remove_adjacency(*iter, mEntityHandle);
    }
  }
}

MBErrorCode MBMeshSet_Vector::clear()
{
  if(mTracking && mAdjFact)
  {
    for(std::vector<MBEntityHandle>::iterator iter = mVector.begin();
        iter != mVector.end(); ++iter)
    {
      mAdjFact->remove_adjacency(*iter, mEntityHandle);
    }
  }
  mVector.clear();
  return MB_SUCCESS;
}

MBErrorCode MBMeshSet_Vector::get_entities(std::vector<MBEntityHandle>& entity_list,
                                              const bool recursive) const
{
  for(std::vector<MBEntityHandle>::const_iterator iter = mVector.begin();
      iter != mVector.end(); ++iter)
  {
    if (recursive && TYPE_FROM_HANDLE(*iter) == MBENTITYSET) {
      if (NULL != mMB) mMB->get_entities_by_handle(*iter, entity_list, recursive);
    }
    else
      entity_list.push_back(*iter);
  }
  return MB_SUCCESS;
}

MBErrorCode MBMeshSet_Vector::get_entities(MBRange& entity_list,
                                              const bool recursive) const
{
  for(std::vector<MBEntityHandle>::const_iterator iter = mVector.begin();
      iter != mVector.end(); ++iter)
  {
    if (recursive && TYPE_FROM_HANDLE(*iter) == MBENTITYSET) {
      if (NULL != mMB) mMB->get_entities_by_handle(*iter, entity_list, recursive);
    }
    else
      entity_list.insert(*iter);
  }
  return MB_SUCCESS;
}

MBErrorCode MBMeshSet_Vector::get_entities_by_type(MBEntityType type,
    std::vector<MBEntityHandle>& entity_list) const
{
  for (std::vector<MBEntityHandle>::const_iterator iter = mVector.begin();
    iter != mVector.end(); ++iter)
  {
    if(TYPE_FROM_HANDLE(*iter) == type)
      entity_list.push_back(*iter);
  }  

  return MB_SUCCESS;
}

MBErrorCode MBMeshSet_Vector::get_entities_by_type(MBEntityType type,
    MBRange& entity_list) const
{

  for(std::vector< MBEntityHandle>::const_iterator iter = mVector.begin();
    iter != mVector.end(); ++iter)
  {
    if(TYPE_FROM_HANDLE(*iter) == type)
      entity_list.insert(*iter);
  }  

  return MB_SUCCESS;
}

MBErrorCode MBMeshSet_Vector::add_entities(const MBEntityHandle *entities,
                                              const int num_entities)
{
  unsigned int prev_size = mVector.size();
  mVector.resize(prev_size + num_entities);
  std::copy(entities, entities+num_entities, &mVector[prev_size]);

  if(mTracking && mAdjFact)
  {
    for(int i = 0; i < num_entities; i++)
      mAdjFact->add_adjacency(entities[i], mEntityHandle);
  }

  return MB_SUCCESS;
}


MBErrorCode MBMeshSet_Vector::add_entities(
    const MBRange& entities)
{

  unsigned int prev_size = mVector.size();
  mVector.resize(prev_size + entities.size());
  std::copy(entities.begin(), entities.end(), &mVector[prev_size]);
  
  if(mTracking && mAdjFact)
  {
    for(MBRange::const_iterator iter = entities.begin();
        iter != entities.end(); ++iter)
      mAdjFact->add_adjacency(*iter, mEntityHandle);
  }

  return MB_SUCCESS;
}


MBErrorCode MBMeshSet_Vector::remove_entities(
    const MBRange& entities)
{
  std::vector<MBEntityHandle>::iterator iter = mVector.begin();
  for(; iter < mVector.end(); )
  {
    if(entities.find(*iter) != entities.end())
    {
      if(mTracking && mAdjFact)
        mAdjFact->remove_adjacency(*iter, mEntityHandle);
      iter = mVector.erase(iter);
    }
    else
      ++iter;
  }

  return MB_SUCCESS;

}

MBErrorCode MBMeshSet_Vector::remove_entities(const MBEntityHandle *entities, 
                                                 const int num_entities)
{
  std::vector<MBEntityHandle>::iterator temp_iter;
  for(int i = 0; i != num_entities; i++ )
  {
    temp_iter = std::find( mVector.begin(), mVector.end(), entities[i]); 
    if( temp_iter != mVector.end() )
    {
      if(mTracking && mAdjFact)
        mAdjFact->remove_adjacency(entities[i], mEntityHandle);
      mVector.erase(temp_iter);
    }
  }
  return MB_SUCCESS;
}


unsigned int MBMeshSet_Vector::num_entities(int*) const
{
  return mVector.size();
}

unsigned int MBMeshSet_Vector::num_entities_by_type(MBEntityType type) const
{
  type = type;
  return 0;
}


MBErrorCode MBMeshSet_Vector::subtract(const MBMeshSet *meshset_2)
{
  std::vector<MBEntityHandle> other_vector;
  meshset_2->get_entities(other_vector, false);

  return remove_entities(&other_vector[0], other_vector.size());
}

MBErrorCode MBMeshSet_Vector::intersect(const MBMeshSet *meshset_2)
{
  MBRange other_range;
  meshset_2->get_entities(other_range, false);

  std::sort(mVector.begin(), mVector.end());

  std::set<MBEntityHandle> tmp;

  std::set_intersection(mVector.begin(), mVector.end(), other_range.begin(), other_range.end(), 
      std::inserter< std::set<MBEntityHandle> >(tmp, tmp.end()));

  mVector.resize(tmp.size());

  std::copy(tmp.begin(), tmp.end(), mVector.begin());

  //track owner
  tmp.clear();
  std::set_difference(mVector.begin(), mVector.end(), other_range.begin(), other_range.end(),
      std::inserter< std::set<MBEntityHandle> >(tmp, tmp.end()));

  if(mTracking && mAdjFact)
  {
    for(std::set<MBEntityHandle>::iterator iter = tmp.begin();
        iter != tmp.end(); ++iter)
      mAdjFact->remove_adjacency(*iter, mEntityHandle);
  }

  return MB_SUCCESS;
}

MBErrorCode MBMeshSet_Vector::unite(const MBMeshSet *meshset_2)
{
  std::vector<MBEntityHandle> other_vector;
  meshset_2->get_entities(other_vector, false);
  return add_entities(&other_vector[0], other_vector.size());
}

