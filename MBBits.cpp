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
#include <algorithm>


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
MBErrorCode MBBitServer::reserve_tag_id(int num_bits, MBTagId tag_id)
{
  for(int i=0; i<(int)MBMAXTYPE; i++)
    if (tag_id >= mBitPageGroups[i].size())
      mBitPageGroups[i].resize( tag_id + 1, 0 );
  
  for(int i=0; i<(int)MBMAXTYPE; i++)
    if (mBitPageGroups[i][tag_id-1])
      return MB_FAILURE;

  for(int i=0; i<(int)MBMAXTYPE; i++)
    mBitPageGroups[i][tag_id-1] = new MBBitPageGroup(num_bits);

  mBitPageGroupsSize = (*mBitPageGroups).size();

  return MB_SUCCESS;
}

  
MBErrorCode MBBitPage::get_entities_with_value( unsigned char value, 
                                                int offset, 
                                                int count, 
                                                int num_bits_per_flag,
                                                MBEntityHandle first,
                                                MBRange& results ) const
{
  if (mBitArray) {
    MBRange::iterator hint = results.begin();
    for (int i = 0; i < count; ++i) 
      if (value == MBBitManipulator::get_bits( (offset+i)*num_bits_per_flag, 
                                               num_bits_per_flag, mBitArray))
        hint = results.insert( hint, first + i );
  }
  return MB_SUCCESS;
}

void MBBitPage::alloc_array( int num_bits_per_flag, 
                             const unsigned char* default_value )
{
  assert(!mBitArray);
  
  mBitArray = new unsigned char[mPageSize];
    // Modifed by J.Kraftcheck : 31 Jan, 2008:
    //  Need to initialize to default value to ensure that we return
    //  the default value for unset entities.  

    // Zero memory if no default value.  Also, if default value is
    // zero, we can zero all the memory w/out worring about the
    // number of bits per entity.
  if (!default_value || !*default_value)
    memset(mBitArray, 0, mPageSize);
    // Otherwise initialize memory using default value
  else {
      // Mask unused bits of default value so that we can set
      // individual bits using bitwise-OR w/out having to worry
      // about masking unwanted stuff.
    unsigned char defval = (*default_value) & ((1u << num_bits_per_flag) - 1);

    switch (num_bits_per_flag) {
      // If number of bits is a power of two (a byte contains a whole
      // number of tag bit values) then use memset to initialize the memory.
      // Note fall-through for switch cases:  for 1-bit tags we first
      // copy the lsb into the adjacent bit, then fall through to 2-bit
      // case, copying last two bits into next two, and so on.
      case 1: defval |= (defval << 1);
      case 2: defval |= (defval << 2);
      case 4: defval |= (defval << 4);
      case 8: memset( mBitArray, defval, mPageSize );
        break;
      // If num_bits_per_flag is not a power of two, then values do
      // not align with byte boundaries.  Need to initialize values
      // individually.
      default:
        memset(mBitArray, 0, mPageSize);
        // Subtract 1 from mPageSize because last byte is unused, allowing
        // questionable MBBitManipulator code to read/write 1 past end 
        // of array.
        for (int i = 0; i < 8 * (mPageSize-1); i += num_bits_per_flag)
          MBBitManipulator::set_bits( i, num_bits_per_flag, defval, mBitArray );
    }
  }
}


MBErrorCode MBBitPageGroup::get_bits( MBEntityHandle start, 
                                      MBEntityID count,
                                      unsigned char* bits,
                                      const unsigned char* def_val)
{
  MBErrorCode result = MB_SUCCESS, tmp_result;
  
  
  MBEntityID id = ID_FROM_HANDLE(start);
  MBEntityHandle page = id/mOffsetFactor;
  MBEntityID offset = id%mOffsetFactor;
  MBEntityID pcount = mOffsetFactor - offset;
  while (count) {
    if (pcount > count)
      pcount = count;
    
    if (page >= mBitPagesSize) {
      memset( bits, def_val ? *def_val : 0, count );
      break;
    }
    
    tmp_result = mBitPages[page]->get_bits( offset, pcount, mBitsPerFlag, bits, def_val);
    if (MB_SUCCESS != tmp_result) {
      memset( bits, 0, pcount );
      result = tmp_result;
    }
    
    count -= pcount;
    bits += pcount;
    ++page;
    offset = 0;
    pcount = mOffsetFactor;
  }
  
  return result;
}

MBErrorCode MBBitPageGroup::set_bits( MBEntityHandle start, 
                                      MBEntityID count,
                                      const unsigned char* bits,
                                      const unsigned char* def_val)
{
  MBErrorCode result = MB_SUCCESS, tmp_result;
  MBEntityID id = ID_FROM_HANDLE(start);
  MBEntityHandle page = id/mOffsetFactor;
  MBEntityID offset = id%mOffsetFactor;
  MBEntityID pcount = mOffsetFactor - offset;
  while (count) {
    if (pcount > count)
      pcount = count;
    
    if (page >= mBitPagesSize) {
      for(int j = page - mBitPagesSize +1; j--;)
        mBitPages.push_back(new MBBitPage());
      mBitPagesSize = mBitPages.size();
    }

    tmp_result = mBitPages[page]->set_bits( offset, pcount, mBitsPerFlag, bits, def_val);
    if (MB_SUCCESS != tmp_result) {
      result = tmp_result;
    }
    
    count -= pcount;
    bits += pcount;
    ++page;
    offset = 0;
    pcount = mOffsetFactor;
  }
  
  return result;
}

  
MBErrorCode MBBitPageGroup::get_entities_with_value( unsigned char value,
                                                     MBEntityHandle first,
                                                     MBEntityHandle last,
                                                     MBRange& results )
{
  MBErrorCode rval;
  assert(last >= first);
  MBEntityID count = last - first + 1;
  
  MBEntityID id = ID_FROM_HANDLE(first);
  MBEntityHandle page = id/mOffsetFactor;
  MBEntityID offset = id%mOffsetFactor;
  MBEntityID pcount = mOffsetFactor - offset;
  while (count && page < mBitPagesSize) {
    if (pcount > count)
      pcount = count;
    
    rval = mBitPages[page]->get_entities_with_value( value, offset, pcount, mBitsPerFlag, first, results );
    if (MB_SUCCESS != rval)
      return rval;
    
    first += pcount;
    count -= pcount;
    ++page;
    offset = 0;
    pcount = mOffsetFactor;
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBBitServer::set_bits( MBTagId tag_id, 
                                   const MBRange& handles,
                                   const unsigned char* data,
                                   const unsigned char* default_val)
{
  --tag_id; // First ID is 1.
  if(tag_id >= mBitPageGroups[0].size() || mBitPageGroups[0][tag_id] == NULL)
    return MB_TAG_NOT_FOUND;

  MBErrorCode rval;
  MBRange::const_pair_iterator i;
  for (i = handles.const_pair_begin(); i != handles.const_pair_end(); ++i) {
    MBEntityType type = TYPE_FROM_HANDLE(i->first);
    assert(TYPE_FROM_HANDLE(i->second) == type); // should be true because id of zero is never used
    MBEntityID count = i->second - i->first + 1;
    rval = mBitPageGroups[type][tag_id]->set_bits(i->first, count, data, default_val);
    if (MB_SUCCESS != rval)
      return rval;
    data += count;
  }
  return MB_SUCCESS;
}

MBErrorCode MBBitServer::get_bits( MBTagId tag_id, 
                                   const MBRange& handles,
                                   unsigned char* data,
                                   const unsigned char* default_val)
{
  --tag_id; // First ID is 1.
  if(tag_id >= mBitPageGroups[0].size() || mBitPageGroups[0][tag_id] == NULL)
    return MB_TAG_NOT_FOUND;

  MBErrorCode rval;
  MBRange::const_pair_iterator i;
  for (i = handles.const_pair_begin(); i != handles.const_pair_end(); ++i) {
    MBEntityType type = TYPE_FROM_HANDLE(i->first);
    assert(TYPE_FROM_HANDLE(i->second) == type); // should be true because id of zero is never used
    MBEntityID count = i->second - i->first + 1;
    rval = mBitPageGroups[type][tag_id]->get_bits(i->first, count, data, default_val);
    if (MB_SUCCESS != rval)
      return rval;
    data += count;
  }
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

MBErrorCode MBBitServer::get_entities_with_tag_value( MBTagId tag_id, 
                                                      MBEntityType type,
                                                      MBRange& entities, 
                                                      const unsigned char bits)
{
  --tag_id; // First ID is 1.
  if(tag_id >= mBitPageGroups[type].size() || mBitPageGroups[type][tag_id] == NULL)
    return MB_TAG_NOT_FOUND;

  return mBitPageGroups[type][tag_id]->
    get_entities_with_value( bits, FIRST_HANDLE(type), LAST_HANDLE(type), entities );
}

MBErrorCode MBBitServer::get_entities_with_tag_value( const MBRange &range,
                                                      MBTagId tag_id, 
                                                      MBEntityType type, 
                                                      MBRange& entities,
                                                      const unsigned char bits)
{
  --tag_id; // First ID is 1.
  if(tag_id >= mBitPageGroups[0].size() || mBitPageGroups[0][tag_id] == NULL)
    return MB_TAG_NOT_FOUND;

  MBErrorCode rval;
  MBRange::const_pair_iterator i;
  for (i = range.const_pair_begin(); i != range.const_pair_end(); ++i) {
    MBEntityType this_type = TYPE_FROM_HANDLE(i->first);
    assert(TYPE_FROM_HANDLE(i->second) == this_type); // should be true because id of zero is never used
    if (type < this_type)
      continue;
    if (type > this_type)
      break;
      
    rval = mBitPageGroups[type][tag_id]->
      get_entities_with_value( bits, i->first, i->second, entities );
    if (MB_SUCCESS != rval)
      return rval;
  }
  return MB_SUCCESS;
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


