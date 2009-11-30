#ifndef BIT_TAG_SERVER_HPP
#define BIT_TAG_SERVER_HPP

#include "MBTypes.h"
#include "MBInternals.hpp"
#include <algorithm>
#include <vector>
#include <assert.h>

class MBRange;
class BitTag;
class BitPage;

/**\brief Manage all bit tag data 
 *
 * This class manages data for all bit tags.
 */
class BitTagServer
{
public:

  /**\brief Clear all tag data
   *
   * Clear all tag values, but keep the actual tags 
   */
  void reset_data();

    /** Allocate space for a tag with the specified ID.
     * Returns MB_TAG_NOT_FOUND if ID is already in use.
     *\param num_bits Number of bits in each per-entity tag value
     */
  MBErrorCode reserve_tag_id( int num_bits, MBTagId tag_id );

    /**\brief Deallocate tag
     *
     * Mask specified tag id as unused and release any memory
     * associated with the tag.
     */
  MBErrorCode release_tag_id( MBTagId tag_id );

    /**\brief Get tag values for an array of entity handles */
  MBErrorCode get_bits( MBTagId tag_id, 
                        const MBEntityHandle* handles, 
                        int num_handles, 
                        unsigned char* data,
                        const unsigned char* default_value ) const;
    /**\brief Get tag values for a range of entity handles */
  MBErrorCode get_bits( MBTagId tag_id, 
                        const MBRange& handles, 
                        unsigned char* data,
                        const unsigned char* default_value ) const;

    /**\brief Set tag values for an array of entity handles */
  MBErrorCode set_bits( MBTagId tag_id, 
                        const MBEntityHandle* handles, 
                        int num_handles,
                        const unsigned char* data, 
                        const unsigned char* default_value );
    /**\brief Set tag values for a range of entity handles */
  MBErrorCode set_bits( MBTagId tag_id, 
                        const MBRange& handles, 
                        const unsigned char* data, 
                        const unsigned char* default_value );

    /**\brief Clear tag values for an array of entity handles */
  MBErrorCode clear_bits( MBTagId tag_id, 
                          const MBEntityHandle* handles, 
                          int num_handles,
                          const unsigned char* default_value );
    /**\brief Clear tag values for a range of entity handles */
  MBErrorCode clear_bits( MBTagId tag_id, 
                          const MBRange& handles, 
                          const unsigned char* default_value );
    
    /**\brief Get entities for which an explicit tag value is stored */
  MBErrorCode get_entities( MBTagId tag_id, 
                            MBRange& entities ) const;

    /**\brief Get entities for which an explicit tag value is stored */
  MBErrorCode get_entities( MBTagId tag_id, 
                            MBEntityType type, 
                            MBRange& entities ) const;

    /**\brief Get entities for which an explicit tag of the specified value is stored */
  MBErrorCode get_entities_with_tag_value( MBTagId tag_id, 
                                           MBEntityType type, 
                                           MBRange& entities,
                                           unsigned char bits ) const;

    /**\brief Get entities for which an explicit tag value is stored */
  MBErrorCode get_entities( const MBRange &range,
                            MBTagId tag_id, 
                            MBEntityType type, 
                            MBRange& entities ) const;

    /**\brief Get entities for which an explicit tag of the specified value is stored */
  MBErrorCode get_entities_with_tag_value( const MBRange &range,
                                           MBTagId tag_id, 
                                           MBEntityType type, 
                                           MBRange& entities,
                                           unsigned char bits ) const;

    //! get all tags defined on an entity
  MBErrorCode get_tags( MBEntityHandle entity,
                        std::vector<MBTag> &tags ) const;

    /**\brief Get number of entities for which an explicit tag value is stored */
  MBErrorCode get_number_entities( MBTagId tag_id, 
                                   MBEntityType type,
                                   int& num_entities ) const;
  
    /**\brief Get number of entities for which an explicit tag value is stored */
  MBErrorCode get_number_entities( const MBRange &range,
                                   MBTagId tag_id, 
                                   MBEntityType type,
                                   int& num_entities ) const;
  
  MBErrorCode get_memory_use( MBTagId tag,
                              unsigned long& total,
                              unsigned long& per_entity ) const;

private:
  
    /**\brief Get BitTag instance for specified tag ID, or NULL if none */
  inline BitTag* get_tag( MBTagId id );
  
    /**\brief Get BitTag instance for specified tag ID, or NULL if none */
  inline const BitTag* get_tag( MBTagId id ) const;
  
    /**\brief Array of BitTag instances, indexed by (tag id - 1) */
  std::vector<BitTag> tagList;
  
};

/**\brief All data for a single bit tag */
class BitTag 
{
  public:
  
  BitTag() : inUse(false) {}
  ~BitTag() { release(); }
  
    // Do destructive copying (assiging moves data from
    // source to destination, such that source is invalid)
    // so that resizing the vector doesn't copy data
  BitTag( const BitTag& other ) 
    { const_cast<BitTag&>(other).swap(*this); }
  BitTag& operator=( const BitTag& other ) 
    { const_cast<BitTag&>(other).swap(*this); return *this; }
    
    /**\brief Get tag values for an array of entity handles */
  MBErrorCode get_bits( const MBEntityHandle* handles, 
                        int num_handles, 
                        unsigned char* data,
                        const unsigned char* default_value ) const;
    /**\brief Get tag values for a range of entity handles */
  MBErrorCode get_bits( const MBRange& handles, 
                        unsigned char* data,
                        const unsigned char* default_value ) const;

    /**\brief Set tag values for an array of entity handles */
  MBErrorCode set_bits( const MBEntityHandle* handles, 
                        int num_handles,
                        const unsigned char* data, 
                        const unsigned char* default_value );
    /**\brief Set tag values for a range of entity handles */
  MBErrorCode set_bits( const MBRange& handles, 
                        const unsigned char* data, 
                        const unsigned char* default_value );

    /**\brief Clear tag values for an array of entity handles */
  MBErrorCode clear_bits( const MBEntityHandle* handles, 
                          int num_handles,
                          const unsigned char* default_value );
    /**\brief Clear tag values for a range of entity handles */
  MBErrorCode clear_bits( const MBRange& handles, 
                          const unsigned char* default_value );
    

    /**\brief Get entities for which an explicit tag value is stored */
  MBErrorCode get_entities( MBRange& entities ) const;

    /**\brief Get entities for which an explicit tag value is stored */
  MBErrorCode get_entities( MBEntityType type, 
                            MBRange& entities ) const;

    /**\brief Get entities for which an explicit tag value is stored */
  MBErrorCode get_entities( const MBRange &range,
                            MBEntityType type, 
                            MBRange& entities ) const;

    /**\brief Get entities for which an explicit tag of the specified value is stored */
  MBErrorCode get_entities_with_bits( MBEntityType type, 
                                      MBRange& entities,
                                      unsigned char bits ) const;

    /**\brief Get entities for which an explicit tag of the specified value is stored */
  MBErrorCode get_entities_with_bits( const MBRange &range,
                                      MBEntityType type, 
                                      MBRange& entities,
                                      unsigned char bits ) const;

    /**\brief Get number of entities for which an explicit tag value is stored */
  MBErrorCode get_number_entities( MBEntityType type,
                                   int& num_entities ) const;
  
    /**\brief Get number of entities for which an explicit tag value is stored */
  MBErrorCode get_number_entities( const MBRange &range,
                                   MBEntityType type,
                                   int& num_entities ) const;
  
  MBErrorCode get_memory_use( unsigned long& total,
                              unsigned long& per_entity ) const;
                              
  enum { Ln2PageSize = 12,              //!< Constant: log2(PageSize)
         PageSize = (1u << Ln2PageSize) //!< Constant: Bytes per BitPage (power of 2)
      };
  
    /**\return false if instance is not in use, true otherwise */
  bool in_use() const { return inUse; }
  
    /**\brief Release all stored tag data
     *
     * Releases all tag data such that tag is still defined, but stores
     * no value for any entity.
     */
  void reset_data();
  
    /**\brief Delete tag.
     *
     * Releases all tag data such that tag is no longer defined.  
     * Equivalent to reset_data() and setting inUse to false.
     */
  void release();
  
    /**\brief Mark this instance is in use.
     *
     * Initialize this BitTag instance such that inUse 
     * returns true.
     *\param bits_per_entity Desired bits per entity.  Note:
     *         implementation is free to allocate more than the
     *         requested number of bits for each entity
     */
  MBErrorCode reserve( unsigned int bits_per_entity );
  
    /**\brief Helper function for destructive assignment */
  void swap( BitTag& other ) {
    std::swap( inUse, other.inUse );
    std::swap( storedBitsPerEntity, other.storedBitsPerEntity );
    std::swap( requestedBitsPerEntity, other.requestedBitsPerEntity );
    for (MBEntityType t = (MBEntityType)0; t != MBMAXTYPE; ++t)
      pageList[t].swap( other.pageList[t] );
  }

  private:
  
  bool inUse; //!< This instance is in use (represents data for some tag)
  std::vector<BitPage*> pageList[MBMAXTYPE]; //!< Array of BitPage instances storing actual data.
  unsigned int requestedBitsPerEntity; //!< user-requested bits per entity
  unsigned int storedBitsPerEntity;    //!< allocated bits per entity (power of 2)
  unsigned int pageShift;              //!< Ln2PageSize + log2(storedBitsPerEntity)
  
  /**\brief Get indices from handle 
   *
   *\param type   Output: entity type
   *\param page   Output: index into pageList[type]
   *\param offset Output: index into pageList[type][page]
   */
  void unpack( MBEntityHandle h, MBEntityType& type, size_t& page, int& offset ) const
    { 
      type = TYPE_FROM_HANDLE(h);
      h = ID_FROM_HANDLE(h);
      page = ((size_t)h) >> pageShift;   // h / (offset*storedBitsPerEntity)
      offset = h & ((1u<<pageShift)-1u); // h % (offset*storedBitsPerEntity)
    }
    
  /**\brief Get the number of tag values that are stored in each BitPage */
  int ents_per_page() const { return 8*PageSize/storedBitsPerEntity; }
};

/**\brief bit tag data
 *
 * This class represents a fixed-size block of memory in which bit tag
 * values are stored.  
 */
class BitPage
{
public:
  /**\brief Initialize memory
   *
   *\param bits_per_ent  Number of bits in each tag value.
   *                     MUST BE A POWER OF TWO.
   *\param init_val      The lower bits_per_ent bits of this byte are
   *                     used to initialize each tag value.
   */
  BitPage( int bits_per_ent, unsigned char init_val );
  
  /**\brief Get tag values
   *
   * Get 'count' tag values, beginning with the one at 'offset'.
   *\param offset Offset into list of values, where a value of zero indicates
   *              the first tag value, a value of one indicates the second
   *              tag value, etc.  NOTE:  This is the value offset, not the
   *              bit offset.
   *\param count  Number of consecutive tag values to get.
   *\param bits_per_ent  Number of bits composing each tag value.  
   *                     NOTE: Must be a power of two.
   *\param data   Memory into which to copy tag values.  Each value is copied
   *              into a separate byte, such that the lower bits of the bit
   *              contain the tag value and any unused higher bits are zero.
   */
  void get_bits( int offset, int count, int bits_per_ent, unsigned char* data ) const;
  
  /**\brief Set tag values
   *
   * Set 'count' tag values, beginning with the one at 'offset'.
   *\param offset Offset into list of values, where a value of zero indicates
   *              the first tag value, a value of one indicates the second
   *              tag value, etc.  NOTE:  This is the value offset, not the
   *              bit offset.
   *\param count  Number of consecutive tag values to set.
   *\param bits_per_ent  Number of bits composing each tag value.  
   *                     NOTE: Must be a power of two.
   *\param data   Memory from which to copy tag values.  Each value is copied
   *              from a separate byte.  The lower 'bits_per_ent' of each
   *              byte are used as the tag value.  Any additional higher bits
   *              are ignored.
   */
  void set_bits( int offset, int count, int bits_per_ent, const unsigned char* data );
  
  /**\brief Set several tag values to the same value.
   *
   * Set 'count' tag values to specified value.
   *\param offset Offset into list of values, where a value of zero indicates
   *              the first tag value, a value of one indicates the second
   *              tag value, etc.  NOTE:  This is the value offset, not the
   *              bit offset.
   *\param count  Number of consecutive tag values to set.
   *\param bits_per_ent  Number of bits composing each tag value.  
   *                     NOTE: Must be a power of two.
   *\param value  The lower 'bits_per_ent' of this
   *              byte are used as the tag value.  Any additional higher bits
   *              are ignored.
   */
  void set_bits( int offset, int count, int bits_per_ent, unsigned char value );

  /**\brief Get tag value
   *
   * Get one tag value.
   *\param offset Offset into list of values, where a value of zero indicates
   *              the first tag value, a value of one indicates the second
   *              tag value, etc.  NOTE:  This is the value offset, not the
   *              bit offset.
   *\param bits_per_ent  Number of bits composing each tag value.  
   *                     NOTE: Must be a power of two.
   *\return       A byte containing the tag value in the lower bits with
   *              any unused higher bits zeroed.
   */
  unsigned char get_bits( int offset, int bits_per_ent ) const;
  
  /**\brief Set tag value
   *
   * Set tag value.
   *\param offset Offset into list of values, where a value of zero indicates
   *              the first tag value, a value of one indicates the second
   *              tag value, etc.  NOTE:  This is the value offset, not the
   *              bit offset.
   *\param bits_per_ent  Number of bits composing each tag value.  
   *                     NOTE: Must be a power of two.
   *\param value  The lower 'bits_per_ent' of this
   *              byte are used as the tag value.  Any additional higher bits
   *              are ignored.
   */
  void set_bits( int offset, int bits_per_ent, unsigned char data );
  
  /**\brief Search stored values for specified value.
   *
   * Find the offsets n in the data at which the specified value occurs,
   * and for each one insert 'start + n' into the passed MBRange.
   *\param value   The value to look for
   *\param offset  The offset at which to begin searching
   *\param count   The number of values to search
   *\param bits_per_ent Number of bits composing each tag value.
   *\param results Result list.
   *\param start   The handle of the entity corresponding to the 
   *               tag value stored at 'offset'
   */
  void search( unsigned char value, int offset, int count, 
               int bits_per_ent, MBRange& results, MBEntityHandle start ) const;

private:

  /**\brief The actual array of bytes */
  char byteArray[BitTag::PageSize];
};

inline unsigned char BitPage::get_bits( int offset, int per_ent ) const
{
    // Assume per_ent is a power of two, which should be guaranteed
    // by higher-level code.
  unsigned char mask = (1<<per_ent)-1; // 2^per_ent - 1
  int byte = (offset * per_ent) >> 3; // shifting 3 is dividing by eight
  int bit =  (offset * per_ent) & 7;  // masking with 7 is modulo eight
  assert(byte < BitTag::PageSize);
  return (byteArray[byte] >> bit) & mask;
} 

inline void BitPage::set_bits( int offset, int per_ent, unsigned char bits )
{
  int byte = (offset * per_ent) >> 3; // shifting 3 is dividing by eight
  int bit =  (offset * per_ent) & 7;  // masking with 7 is modulo eight
  assert(byte < BitTag::PageSize);
    // Assume per_ent is a power of two, which should be guaranteed
    // by higher-level code.
  unsigned char mask = ((1<<per_ent)-1) << bit;
  byteArray[byte] = (byteArray[byte] & ~mask) | ((bits << bit) & mask);
} 

inline void BitPage::get_bits( int offset, int count, int per_ent, unsigned char* data ) const
{
  unsigned char* end = data+count;
  while (data != end)
    *(data++) = get_bits( offset++, per_ent );
}

inline void BitPage::set_bits( int offset, int count, int per_ent, const unsigned char* data )
{
  const unsigned char* end = data+count;
  while (data != end)
    set_bits( offset++, per_ent, *(data++) );
}

inline void BitPage::set_bits( int offset, int count, int per_ent, unsigned char value )
{
  int end = offset + count;
  while (offset < end)
    set_bits( offset++, per_ent, value );
}

  
inline BitTag* BitTagServer::get_tag( MBTagId id )
  { return id-1 >= tagList.size() || !tagList[id-1].in_use() ? 0 : &tagList[id-1]; }
  
inline const BitTag* BitTagServer::get_tag( MBTagId id ) const
  { return id-1 >= tagList.size() || !tagList[id-1].in_use() ? 0 : &tagList[id-1]; }

inline MBErrorCode 
BitTagServer::get_bits( MBTagId tag_id, 
                        const MBEntityHandle* handles, 
                        int num_handles, 
                        unsigned char* data,
                        const unsigned char* default_value ) const
{
  if (const BitTag* ptr = get_tag(tag_id))
    return ptr->get_bits( handles, num_handles, data, default_value );
  else
    return MB_TAG_NOT_FOUND;
}

inline MBErrorCode 
BitTagServer::get_bits( MBTagId tag_id, 
                        const MBRange& handles, 
                        unsigned char* data,
                        const unsigned char* default_value ) const
{
  if (const BitTag* ptr = get_tag(tag_id))
    return ptr->get_bits( handles, data, default_value );
  else
    return MB_TAG_NOT_FOUND;
}

inline MBErrorCode 
BitTagServer::set_bits( MBTagId tag_id, 
                        const MBEntityHandle* handles, 
                        int num_handles,
                        const unsigned char* data, 
                        const unsigned char* default_value )
{
  if (BitTag* ptr = get_tag(tag_id))
    return ptr->set_bits( handles, num_handles, data, default_value );
  else
    return MB_TAG_NOT_FOUND;
}
inline MBErrorCode 
BitTagServer::set_bits( MBTagId tag_id, 
                        const MBRange& handles, 
                        const unsigned char* data, 
                        const unsigned char* default_value )
{
  if (BitTag* ptr = get_tag(tag_id))
    return ptr->set_bits( handles, data, default_value );
  else
    return MB_TAG_NOT_FOUND;
}

inline MBErrorCode 
BitTagServer::clear_bits( MBTagId tag_id, 
                          const MBEntityHandle* handles, 
                          int num_handles,
                          const unsigned char* default_value )
{
  if (BitTag* ptr = get_tag(tag_id))
    return ptr->clear_bits( handles, num_handles, default_value );
  else
    return MB_TAG_NOT_FOUND;
}
inline MBErrorCode 
BitTagServer::clear_bits( MBTagId tag_id, 
                          const MBRange& handles, 
                          const unsigned char* default_value )
{
  if (BitTag* ptr = get_tag(tag_id))
    return ptr->clear_bits( handles, default_value );
  else
    return MB_TAG_NOT_FOUND;
}
    

inline MBErrorCode 
BitTagServer::get_entities( MBTagId tag_id, MBRange& entities ) const
{
  if (const BitTag* ptr = get_tag(tag_id))
    return ptr->get_entities( entities );
  else
    return MB_TAG_NOT_FOUND;
}

inline MBErrorCode 
BitTagServer::get_entities( MBTagId tag_id, 
                            MBEntityType type, 
                            MBRange& entities ) const
{
  if (const BitTag* ptr = get_tag(tag_id))
    return ptr->get_entities( type, entities );
  else
    return MB_TAG_NOT_FOUND;
}

inline MBErrorCode 
BitTagServer::get_entities_with_tag_value( MBTagId tag_id, 
                                           MBEntityType type, 
                                           MBRange& entities,
                                           unsigned char bits ) const
{
  if (const BitTag* ptr = get_tag(tag_id))
    return ptr->get_entities_with_bits( type, entities, bits );
  else
    return MB_TAG_NOT_FOUND;
}

inline MBErrorCode 
BitTagServer::get_entities( const MBRange &range,
                            MBTagId tag_id, 
                            MBEntityType type, 
                            MBRange& entities ) const
{
  if (const BitTag* ptr = get_tag(tag_id))
    return ptr->get_entities( range, type, entities );
  else
    return MB_TAG_NOT_FOUND;
}

inline MBErrorCode 
BitTagServer::get_entities_with_tag_value( const MBRange &range,
                                           MBTagId tag_id, 
                                           MBEntityType type, 
                                           MBRange& entities,
                                           unsigned char bits ) const
{
  if (const BitTag* ptr = get_tag(tag_id))
    return ptr->get_entities_with_bits( range, type, entities, bits );
  else
    return MB_TAG_NOT_FOUND;
}

inline MBErrorCode 
BitTagServer::get_number_entities( const MBTagId tag_id, 
                                   const MBEntityType type,
                                   int& num_entities ) const
{
  if (const BitTag* ptr = get_tag(tag_id))
    return ptr->get_number_entities( type, num_entities );
  else
    return MB_TAG_NOT_FOUND;
}
  
inline MBErrorCode 
BitTagServer::get_number_entities( const MBRange &range,
                                   const MBTagId tag_id, 
                                   const MBEntityType type,
                                   int& num_entities ) const
{
  if (const BitTag* ptr = get_tag(tag_id))
    return ptr->get_number_entities( range, type, num_entities );
  else
    return MB_TAG_NOT_FOUND;
}
  
inline MBErrorCode 
BitTagServer::get_memory_use( MBTagId tag,
                              unsigned long& total,
                              unsigned long& per_entity ) const
{
  if (const BitTag* ptr = get_tag(tag)) {
    MBErrorCode result = ptr->get_memory_use( total, per_entity );
    total += sizeof(*this);
    return result;
  }
  else
    return MB_TAG_NOT_FOUND;
}

#endif
