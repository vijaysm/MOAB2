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

#ifndef MB_BITS_HPP
#define MB_BITS_HPP

#ifndef IS_BUILDING_MB
#error "MBBits.hpp isn't supposed to be included into an application"
#endif

#include "MBInterface.hpp"
#include "MBInternals.hpp"

#include <assert.h>

//! bit manipulator class
class MBBitManipulator
{
public:

  //! set bits from value into bit array given an offset and bits per offset
  static MBErrorCode set_bits(int bit_offset, int num_bits, 
                       unsigned char value, unsigned char* bit_array );

  //! gets bits from bit_array given an offset and bits per offset
  static unsigned char get_bits(int bit_offset, int num_bits, 
                      unsigned char* bit_array );

  // bit_array example
  //
  //                 010101010101001100100111010010101010010111001010
  // char            \1 byte/\1 byte/\1 byte/\1 byte/\1 byte/\1 byte/
  // 3 bits per flag \1/\2/\3/\4/\5/\6/\7/\8/\9/\./\./\./\./\./\./
  //
  // lets offset 5 and get those bits
  // int x = get_bits( 5, 3, char_array )
  // bits returned are  001
  //
  // Assumptions made by this class:
  //   num_bits is 8 or less, if it is greater, there is no guarantee 
  //   you'll get what you want.

private:

  //! quick mod 8 operation for finding out how far some bits
  //! are from a byte boundary
  static int mod_8(int val)
  {
    return val & 0x7;
  }  

  // stores the masks indexed by how many bits to mask
  static const unsigned int masks[9];
};

/*! set the bits in an array
    takes bit_offset which is the ith bit in the array
    takes num_bits which is number of bits from the ith bit to set
    takes value which is the value to set
    takes_byte array which is the array of bits
*/
inline MBErrorCode MBBitManipulator::set_bits(
    int bit_offset,
    int num_bits,
    unsigned char value,
    unsigned char* byte_array 
    )
{

  // check the value to make sure it is good
  if((value | masks[num_bits]) != masks[num_bits])
    return MB_FAILURE;

  // offset our pointer to where we want it
  byte_array += bit_offset/8;

  // copy data from two bytes into our unsigned short
  unsigned short int old_data = *(byte_array +1) << 8;
  old_data += *byte_array;

  // shift the value so it is lined up with word boundary of bit array
  unsigned short int aligned_value = value << mod_8(bit_offset);

  // shift the mask where we need it
  unsigned short int aligned_mask = masks[num_bits] << mod_8(bit_offset);

  // AND aligned_ask and byte_bound to preserve bits that are on
  unsigned short int x = (~aligned_mask) & old_data;

  // OR to set value
  old_data = x | aligned_value;

  // set the data
  *byte_array = (old_data & 0xFF); 
  *(byte_array + 1) = (old_data >> 8);

  return MB_SUCCESS;
}

/*! get the bits in an array
    takes bit_offset which is the ith bit in the array
    takes num_bits which is number of bits from the ith bit to get
    takes_byte array which is the array of bits
    return the bits
*/
inline unsigned char MBBitManipulator::get_bits(
    int bit_offset,
    int num_bits,
    unsigned char* byte_array
    )
{
  // offset our pointer to where we want it
  byte_array += bit_offset/8;
  
  // copy data from two bytes into our unsigned short
  unsigned short int data = *(byte_array +1) << 8;
  data += *byte_array;

  // shift data where we need it
  unsigned int aligned_data = data >> mod_8(bit_offset);

  // AND caputures value of the flag
  return aligned_data & masks[num_bits];
}


//! bit page class
/*! This class stores a page of memory for storing bits
*/
class MBBitPage
{

public:

  //! the page size
  static const int mPageSize;
  
  // page sizes are 512 bytes   
  //   for 1 bit per entity, that is 4096 bit flags
  //   for 2 bits per entity, that is 2048 bit flags
  //   for 3 bits per entity, that is 1365.333 bit flags
  //   for 4 bits per entity, that is 1024 bit flags
  //   for 5 bits per entity, that is 819.2 bit flags
  //   for 6 bits per entity, that is 682.666 bit flags
  //   for 7 bits per entity, that is 585.142857
  //   for 8 bits per entity, that is 512 bit flags
  //
  //   for some, we'll have left over bits on the end.

  // default constructor
  MBBitPage() : mBitArray(NULL) {}
  
  // default destructor
  ~MBBitPage()
  {  
    if(mBitArray)
      delete [] mBitArray;
  }
  
  //! get the bits from a bit page
  MBErrorCode get_bits(int offset, int num_bits_per_flag, unsigned char& bits);
  
  //! set the bits in a bit page
  MBErrorCode set_bits(int offset, int num_bits_per_flag, unsigned char bits);
  
  //! set the bits in a bit page only if space has been allocated
  MBErrorCode weak_set_bits(int offset, int num_bits_per_flag, unsigned char bits);

  bool has_data() const { return mBitArray != NULL; }

private:
  //! bit array  uses lazy allocation
  unsigned char* mBitArray;
  
  //! don't allow copying of these bit pages
  MBBitPage(const MBBitPage&) 
    : mBitArray(NULL) 
  { 
    // not even by self
    assert(0);
  }
  
  //! don't allow copying of these bit pages
  MBBitPage& operator=(const MBBitPage&) 
  { 
    // not even by self
    assert(0);
    return *this;
  }
  
};

/*! get the bits from a bit page
    takes bit offset into the page
    takes how many bits to get
    return bits
*/
inline MBErrorCode MBBitPage::get_bits(int offset, int num_bits_per_flag, unsigned char& bits)
{
  // because the default bits are 0x0, 
  // we'll return that if no memory has been allocated
  if(!mBitArray)
  {
    bits = 0;
    return MB_SUCCESS;
  }

  // get bits using bit manipulator
  bits = MBBitManipulator::get_bits(
      offset*num_bits_per_flag, num_bits_per_flag, mBitArray);

  return MB_SUCCESS;
}

/*! set the bits in a bit page
    takes bit offset into the page
    takes how many bits to set
    takes the bits to set
*/
inline MBErrorCode MBBitPage::set_bits(int offset, 
    int num_bits_per_flag, unsigned char bits)
{
  // if memory hasn't been allocated, allocate it and zero the memory
  if(!mBitArray)
  {
    mBitArray = new unsigned char[mPageSize];
    memset(mBitArray, 0, mPageSize);
  }

  // set the bits using bit manipulator
  return MBBitManipulator::set_bits(
      offset*num_bits_per_flag, num_bits_per_flag, bits, mBitArray);
}

/*! weak set bits only sets bits if memory has been allocated
    takes bit offset
    takes number of bits to set
    takes the bits to set
*/
inline MBErrorCode MBBitPage::weak_set_bits(int offset, 
    int num_bits_per_flag, unsigned char bits)
{
  return mBitArray ? MBBitManipulator::set_bits(
      offset*num_bits_per_flag, num_bits_per_flag, bits, mBitArray) : MB_SUCCESS;
}

//! class which is a collection of bit pages
class MBBitPageGroup
{
public:
  // default constructor
  MBBitPageGroup(int bits_per_flag)
    : mBitsPerFlag(bits_per_flag)
  {
    //compute the offset factor based on the number of bits for each entity
    //this offset will most likely leave some unused bits at the end of a page
    mOffsetFactor = compute_offset_factor(bits_per_flag);
    mBitPagesSize = 0;
  }

  // default destructor
  ~MBBitPageGroup() 
  { 
    // delete each bit page
    for (std::vector<MBBitPage*>::iterator iter = mBitPages.begin(); iter != mBitPages.end(); ++iter)
      delete *iter;

    // clean out the vector of pointers
    mBitPages.clear(); 
  }

  //! get bits from bit pages
  MBErrorCode get_bits(MBEntityHandle handle, unsigned char& bits);

  //! set bits in bit pages
  MBErrorCode set_bits(MBEntityHandle handle, unsigned char bits);

  //! set bits in bit pages only if the bit page allocated memory
  MBErrorCode weak_set_bits(MBEntityHandle handle, unsigned char bits);

  MBErrorCode get_entities(MBEntityType type, MBRange& entities);

  //! if this page group contains this entity, return true, otherwise false
  bool contains(const MBEntityHandle handle) const;

  int tag_size() const { return mBitsPerFlag; }

private:

  // compute offset factor to use when computing which page to index into
  int compute_offset_factor(int num_bits)
  {
    // subtract one from page size to prevent reading past
    // the page one byte (BitManipulator does this).
    // this results in one byte at the end of the page
    // that isn't used, but that's no big deal.
    return (MBBitPage::mPageSize - 1) * 8 / num_bits;
  }
  
  //!  number of bits for each entity  
  unsigned short mBitsPerFlag;
  
  //! offset factor used when computing which page to jump to
  unsigned short mOffsetFactor;

  unsigned int  mBitPagesSize;

  //! vector of bit pages
  std::vector<MBBitPage*> mBitPages;
  
  //! don't allow copy of this class
  MBBitPageGroup& operator=(const MBBitPageGroup&)
  {
    assert(0);
    return *this;
  }

  //! don't allow copy of this class
  MBBitPageGroup(const MBBitPageGroup&)
  {
    assert(0);
  }

};

/*! get the bits from bit pages
    takes entity handle
    return the bits
*/
inline MBErrorCode MBBitPageGroup::get_bits(MBEntityHandle handle, unsigned char& bits)
{
  // strip off the entity type
  handle = ID_FROM_HANDLE(handle);
  // figure out which page to jump to
  unsigned int which_page = handle / mOffsetFactor;
  
  // if the page isn't there, just return 0x0
  if(which_page >= mBitPagesSize)
  {
    bits = 0;
    return MB_SUCCESS;
  }

  // return bits from bit page
  return mBitPages[which_page]->get_bits( (handle - ( which_page * mOffsetFactor )), mBitsPerFlag, bits);
}

/*! set the bits in bit pages
    takes entity handle
    takes the bits to set
*/
inline MBErrorCode MBBitPageGroup::set_bits(MBEntityHandle handle, unsigned char bits)
{
  // strip off the entity type
  handle = ID_FROM_HANDLE(handle);
  // figure out which page to jump to
  unsigned int which_page = handle / mOffsetFactor;

  // if the page doesn't exist, make one
  if(which_page >= mBitPagesSize)
  {
    for(int j= which_page - mBitPagesSize +1; j--;)
      mBitPages.push_back(new MBBitPage());
    mBitPagesSize = mBitPages.size();
  }

  // return set of bits in page
  return mBitPages[which_page]->set_bits( (handle - ( which_page * mOffsetFactor )), mBitsPerFlag, bits);
}


/*! set the bits in bit pages
    if the page didn't allocate memory yet, the bits aren't set
    takes entity handle
    takes the bits to set
*/
inline MBErrorCode MBBitPageGroup::weak_set_bits(MBEntityHandle handle, unsigned char bits)
{
  // strip off entity type
  handle = ID_FROM_HANDLE(handle);
  // find out which page to jump to
  unsigned int which_page = handle / mOffsetFactor;
  
  // if the page doesn't exist, return
  if(which_page >= mBitPagesSize)
    return MB_SUCCESS;
 
  // try to set bits 
  return mBitPages[which_page]->weak_set_bits( (handle - ( which_page * mOffsetFactor )), mBitsPerFlag, bits);
}

inline bool MBBitPageGroup::contains(const MBEntityHandle handle) const
{
  // strip off the entity type
  unsigned int entity_id = ID_FROM_HANDLE(handle);
  // figure out which page to jump to
  unsigned int which_page = entity_id / mOffsetFactor;
  
  // if the page isn't there, just return 0x0
  return (which_page >= mBitPagesSize ? false : true);
}

//! MBBitServer class provides constant time 
//! lookup for bit flags tagged on entities
class MBBitServer
{
public:
  //! default constructor
  MBBitServer()
  {
    mBitPageGroupsSize = 0;
  }
  //! default destructor
  ~MBBitServer()
  {
    // clean things out
    for(int i = 0; i<(int)MBMAXTYPE; i++)
    {
      for (std::vector<MBBitPageGroup*>::iterator iter = mBitPageGroups[i].begin();
          iter != mBitPageGroups[i].end();
          ++iter)
      {
        if(*iter == NULL)
          continue;
        delete *iter;
      }

      mBitPageGroups[i].clear();
    }
  }

  void reset_data();

  //! return an available tag id for use
  MBErrorCode reserve_tag_id(int num_bits, MBTagId& tag_id);
  //! release a tag id for reuse
  MBErrorCode release_tag_id(MBTagId tag_id);

  //! get the bits associated with an entity handle
  MBErrorCode get_bits(MBTagId tag_id, MBEntityHandle handle, unsigned char& bits);
  //! set the bits associated with an entity handle
  MBErrorCode set_bits(MBTagId tag_id, MBEntityHandle handle, unsigned char bits);
  //! set the bits associated with an entity handle, only if memory has been allocated
  MBErrorCode weak_set_bits(MBTagId tag_id, MBEntityHandle handle, unsigned char bits);

  MBErrorCode get_entities(MBTagId tag_id, MBEntityType type, MBRange& entities);

  MBErrorCode get_entities_with_tag_value(MBTagId tag_id, 
                                           MBEntityType type, 
                                           MBRange& entities,
                                           const unsigned char bits);

  MBErrorCode get_entities(const MBRange &range,
                            MBTagId tag_id, 
                            MBEntityType type, 
                            MBRange& entities);

  MBErrorCode get_entities_with_tag_value(const MBRange &range,
                                           MBTagId tag_id, MBEntityType type, 
                                           MBRange& entities,
                                           const unsigned char bits);

    //! get all tags defined on an entity
  MBErrorCode get_tags(const MBEntityHandle entity,
                        std::vector<MBTag> &tags);

  MBErrorCode get_number_entities(const MBTagId tag_id, 
                                   const MBEntityType type,
                                   int& num_entities);
  
  MBErrorCode get_number_entities( const MBRange &range,
                                    const MBTagId tag_id, 
                                    const MBEntityType type,
                                    int& num_entities);
  
private:
  //! bit pages are indexed by tag id and entity type
  std::vector< MBBitPageGroup* > mBitPageGroups[MBMAXTYPE];
  unsigned long mBitPageGroupsSize;

};

/*! get some bits based on a tag id and handle */
inline MBErrorCode MBBitServer::get_bits(MBTagId tag_id, 
    MBEntityHandle handle, unsigned char& bits)
{
  if(tag_id >= (*mBitPageGroups).size() || (*mBitPageGroups)[tag_id] == NULL)
    return MB_FAILURE;

  return mBitPageGroups[TYPE_FROM_HANDLE(handle)][tag_id]->get_bits(handle, bits);

}

/*! set some bits based on a tag id and handle */
inline MBErrorCode MBBitServer::set_bits(MBTagId tag_id, 
    MBEntityHandle handle, unsigned char bits)
{
  if(tag_id >= mBitPageGroupsSize || (*mBitPageGroups)[tag_id] == NULL)
    return MB_FAILURE;

  return mBitPageGroups[TYPE_FROM_HANDLE(handle)][tag_id]->set_bits(handle, bits);
}

/*! set some bits based on a tag id and handle only if memory has been allocated*/
inline MBErrorCode MBBitServer::weak_set_bits(MBTagId tag_id, 
    MBEntityHandle handle, unsigned char bits)
{
  if(tag_id >= mBitPageGroupsSize || (*mBitPageGroups)[tag_id] == NULL)
    return MB_SUCCESS;

  return mBitPageGroups[TYPE_FROM_HANDLE(handle)][tag_id]->weak_set_bits(handle, bits);
}

inline MBErrorCode MBBitServer::get_entities(MBTagId tag_id, MBEntityType type, MBRange& entities)
{
  if(tag_id >= mBitPageGroupsSize || (*mBitPageGroups)[tag_id] == NULL)
    return MB_FAILURE;

  return mBitPageGroups[type][tag_id]->get_entities(type, entities);
}

    //! get all tags defined on an entity
inline MBErrorCode MBBitServer::get_tags(const MBEntityHandle entity,
                                           std::vector<MBTag> &tags) 
{
    // get the tags defined for this type
  MBEntityType this_type = TYPE_FROM_HANDLE(entity);

  for (long i = 0; i < (long) mBitPageGroups[this_type].size(); i++) {
    if (mBitPageGroups[this_type][i] != NULL &&
        mBitPageGroups[this_type][i]->contains(entity))
      tags.push_back(TAG_HANDLE_FROM_ID(i, MB_TAG_BIT));
  }
  return MB_SUCCESS;
}

#endif


