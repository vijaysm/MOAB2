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




/*  Bit tools for MB
 *
 *  File   :      MBBits.hpp
 *  Creator:      Clinton Stimpson
 *  Date   :      10-10-2002
 */

#include "MBBits.hpp"
#include "MBRange.hpp"


/*! bit masks for how many bits to mask */
const unsigned int MBBitManipulator::masks[9] = 
{
  0x0,
  0x01,
  0x03,
  0x07,
  0x0F,
  0x1F,
  0x3F,
  0x7F,
  0xFF
};

/*! page size of 512 bytes for storing bits */
const int MBBitPage::mPageSize = 512;

/*! returns an available tag id to use when getting and setting bits */
MBErrorCode MBBitServer::reserve_tag_id(int num_bits, MBTagId& tag_id)
{
  tag_id = 0;

  // make sure we get a good number of bits
  if(num_bits <= 0 || num_bits >8)
    return MB_FAILURE;

  // see if we have any bit page groups that aren't being used
  for(std::vector<MBBitPageGroup*>::iterator iter = (*mBitPageGroups).begin();
      iter != (*mBitPageGroups).end(); ++iter)
  {
    if(*iter == NULL)
    {
      tag_id = iter - (*mBitPageGroups).begin() + 1;
      break;
    }
  }

  // if they are all being used, make space for a new one
  if(!tag_id)
  {
    for(int i=0; i<(int)MBMAXTYPE; i++)
      mBitPageGroups[i].push_back( NULL );
    tag_id = (*mBitPageGroups).size();
  }

  for(int i=0; i<(int)MBMAXTYPE; i++)
    mBitPageGroups[i][tag_id-1] = new MBBitPageGroup(num_bits);

  mBitPageGroupsSize = (*mBitPageGroups).size();

  // can we really fail?
  return MB_SUCCESS;
  
}


/*! give back a tag id that was used to set and get bits */
MBErrorCode MBBitServer::release_tag_id(MBTagId tag_id)
{
  // tag ids begin with 1
  --tag_id;
  
  // make sure tag_id is good
  if(tag_id >= (*mBitPageGroups).size())
    return MB_TAG_NOT_FOUND;

  // mark it unused
  for(int i=0; i<(int)MBMAXTYPE; i++)
  {
    delete mBitPageGroups[i][tag_id];
    mBitPageGroups[i][tag_id] = NULL;
  }

  // clean up a bit if this is the last one
  if((tag_id+1) == (unsigned short)(*mBitPageGroups).size())
  {
    for(int i=0; i<(int)MBMAXTYPE; i++)
    {
      mBitPageGroups[i].resize(mBitPageGroups[i].size()-1);
    }
  }
  
  mBitPageGroupsSize = (*mBitPageGroups).size();
  
  return MB_SUCCESS;
}


void MBBitServer::reset_data()
{
  for(int i = 0; i<(int)MBMAXTYPE; i++)
  {
    for (std::vector<MBBitPageGroup*>::iterator iter = mBitPageGroups[i].begin();
        iter != mBitPageGroups[i].end();
        ++iter)
    {
      if(*iter == NULL)
        continue;
      int tag_size = (*iter)->tag_size();
      delete *iter;
      *iter = new MBBitPageGroup(tag_size);
    }
  }
}

MBErrorCode MBBitServer::get_memory_use( MBTagId tag_id,
                                         unsigned long& total,
                                         unsigned long& per_entity ) const
{
  per_entity = 1; // cannot return fraction of bytes
  
  unsigned max_num_tags = 0;
  for (unsigned i = 0; i < MBMAXTYPE; ++i) {
    if (mBitPageGroups[i].size() > max_num_tags)
      max_num_tags = mBitPageGroups[i].size();

    if (tag_id >= mBitPageGroups[i].size())
      continue;
      
    total += mBitPageGroups[i][tag_id]->get_memory_use();
    total += sizeof(MBBitPageGroup*) * mBitPageGroups[i].capacity() / mBitPageGroups[i].size();
  }
  total += sizeof( std::vector<MBBitPageGroup*> );
  total += sizeof(*this) / max_num_tags;

  return MB_SUCCESS;
}

unsigned long MBBitPageGroup::get_memory_use() const
{
  unsigned long result = sizeof(*this);
  result += sizeof(MBBitPage*) * mBitPages.capacity();
  for (unsigned long i = 0; i < mBitPages.size(); ++i)
    if (mBitPages[i] && mBitPages[i]->has_data())
      result += MBBitPage::mPageSize;
  return result;
}

//! get the entities
MBErrorCode MBBitPageGroup::get_entities(MBEntityType type, MBRange& entities)
{
  std::vector<MBBitPage*>::iterator iter;
  int dum =0;
  MBEntityHandle handle = CREATE_HANDLE(type, MB_START_ID, dum);
  bool first_time = true;
  for(iter = mBitPages.begin(); iter < mBitPages.end(); ++iter)
  {
    if(*iter)
    {
      if((*iter)->has_data())
      {
        entities.insert(handle, handle + mOffsetFactor);
      }        
    }
    if(first_time)
    {
      first_time = false;
      handle = handle + mOffsetFactor - MB_START_ID;
    }
    else
      handle += mOffsetFactor;
  }
  return MB_SUCCESS;
}

MBErrorCode MBBitServer::get_entities(const MBRange &range, MBTagId tag_id, MBEntityType type, 
                                        MBRange& entities)
{
  --tag_id;
  if(tag_id >= mBitPageGroupsSize || (*mBitPageGroups)[tag_id] == NULL)
    return MB_FAILURE;

  MBRange dum_range;
  MBErrorCode result = mBitPageGroups[type][tag_id]->get_entities(type, dum_range);
  if (MB_FAILURE == result) return result;

  std::set_intersection(dum_range.begin(), dum_range.end(),
                        range.begin(), range.end(),
                        mb_range_inserter(entities));
  
  return result;
}

MBErrorCode MBBitServer::get_entities_with_tag_value(MBTagId tag_id, 
                                                       MBEntityType type,
                                                       MBRange& entities, 
                                                       const unsigned char bits)
{
  MBRange possibles;
  MBErrorCode result = get_entities(tag_id, type, possibles);
  if (MB_SUCCESS != result || possibles.empty()) return result;
  MBErrorCode tmp_result;
  unsigned char dum;
  for (MBRange::iterator it = possibles.begin(); it != possibles.end(); it++) {
    tmp_result = get_bits(tag_id, *it, dum);
    if (dum == bits) entities.insert(*it);
    if (tmp_result != MB_SUCCESS) result = tmp_result;
  }
  
  return result;
}

MBErrorCode MBBitServer::get_entities_with_tag_value(const MBRange &range,
                                                       MBTagId tag_id, MBEntityType type, 
                                                       MBRange& entities,
                                                       const unsigned char bits)
{
  MBRange temp1, temp2;
  MBErrorCode result = get_entities(tag_id, type, temp1);
  if (MB_SUCCESS != result) return result;
  std::set_intersection(range.begin(), range.end(),
                        temp1.begin(), temp1.end(), mb_range_inserter(temp2));
  if (temp2.empty()) return result;
  
  unsigned char dum;
  MBErrorCode tmp_result;
  for (MBRange::iterator it = temp2.begin(); it != temp2.end(); it++) {
    tmp_result = get_bits(tag_id, *it, dum);
    if (dum == bits) entities.insert(*it);
    if (tmp_result != MB_SUCCESS) result = tmp_result;
  }
  
  return result;
}

MBErrorCode MBBitServer::get_number_entities( const MBTagId tag_id,
                                                const MBEntityType type,
                                                int& num_entities) 
{
  MBRange dum_range;
  MBErrorCode result = get_entities(tag_id, type, dum_range);
  num_entities = dum_range.size();
  return result;
}

  
MBErrorCode MBBitServer::get_number_entities( const MBRange &range,
                                                const MBTagId tag_id, 
                                                const MBEntityType type,
                                                int& num_entities)
{
  MBRange dum_range;
  MBErrorCode result = get_entities(range, tag_id, type, dum_range);
  num_entities = dum_range.size();
  return result;
}

#ifdef BIT_MANIP_TEST

int main()
{

  unsigned char arr[10] = { 0 };

  unsigned char bits=0;

  MBBitManipulator::set_bits(7, 2, 0x3, arr);

  bits = MBBitManipulator::get_bits(6, 3, arr);
  

  return 0;
}


#endif


// unit test
#ifdef TEST

#include <stdio.h>
#include <assert.h>

int main()
{

  MBBitServer bit_server;
  MBTagId tag_id;
  bit_server.reserve_tag_id(5, tag_id);
  unsigned char bits;
  bit_server.get_bits(tag_id, 400, bits);
  assert(bits == 0x0);

  bit_server.set_bits(tag_id, 600, 0x9);
  bit_server.get_bits(tag_id, 600, bits);
  assert(bits == 0x9);

  //bit_server.release_tag_id(tag_id);
    
  srand(0xb);

  for(int i=100; i--; )
  {
    int num_bits = rand() % 9;
    if(bit_server.reserve_tag_id(num_bits, tag_id) == MB_SUCCESS)
    {
      for(int j=1000; j--;)
      {
        unsigned short handle = rand();
        unsigned char some_bits = rand();
        some_bits <<= 8-num_bits;
        some_bits >>= 8-num_bits;
        bit_server.set_bits(tag_id, handle, some_bits);
        unsigned char some_bits_again = rand();
        bit_server.get_bits(tag_id, handle, some_bits_again);
        if(some_bits_again != some_bits)
        {
          bit_server.get_bits(tag_id, handle, some_bits_again);
          printf("ERROR\n  num bits = %i\n", num_bits);
          printf( "   tag id    %i\n", tag_id);
          printf( "   handle    %u\n", handle);
          printf( "   set value %i\n", (int)some_bits);
          printf( "   get value %i\n", (int)some_bits_again);
          assert(some_bits_again == some_bits);
        }
      }
      
      bit_server.release_tag_id(tag_id);
    }
  }
  return 0;
}

#endif


