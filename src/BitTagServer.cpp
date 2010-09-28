#include "BitTagServer.hpp"
#include "moab/Range.hpp"
#include "Internals.hpp"
#include <stdlib.h>
#include <string.h>

namespace moab {

void BitPage::search( unsigned char value, int offset, int count, 
                      int per_ent, Range& results, EntityHandle start ) const
{
  const int end = offset + count;
  Range::iterator hint = results.begin();
  while (offset != end) {
    if (get_bits( offset, per_ent ) == value)
      hint = results.insert( hint, start );
    ++offset;
    ++start;
  }
}

BitPage::BitPage( int per_ent, unsigned char init_val )
{
  unsigned char mask = (1<<per_ent)-1; // 2^per_ent - 1
  init_val &= mask;
  switch (per_ent) {
    default: assert(false); abort(); break; // must be power of two
      // Note: no breaks. fall through such that all bits in init_val are set
    case 1: init_val |= (init_val << 1);
    case 2: init_val |= (init_val << 2);
    case 4: init_val |= (init_val << 4);
    case 8: ;
  }
  memset( byteArray, init_val, BitTag::PageSize );
}

ErrorCode BitTagServer::reserve_tag_id( int num_bits, TagId tag_id )
{
  if (tag_id > tagList.size())
    tagList.resize( tag_id );
  else if (tagList[tag_id-1].in_use())
    return MB_ALREADY_ALLOCATED;
  
  return tagList[tag_id-1].reserve( num_bits );
}

ErrorCode BitTagServer::release_tag_id( TagId tag_id )
{
  if (tag_id > tagList.size() && !tagList[tag_id-1].in_use())
    return MB_TAG_NOT_FOUND;
  
  tagList[tag_id-1].release();
  return MB_SUCCESS;
}

void BitTagServer::reset_data()
{
  for (std::vector<BitTag>::iterator i = tagList.begin(); i != tagList.end(); ++i) 
    i->reset_data();
}

ErrorCode BitTagServer::get_tags( EntityHandle entity,  
                                  std::vector<Tag> &tags ) const
{
  for (size_t i = 0; i < tagList.size(); ++i) 
    if (tagList[i].in_use() && tagList[i].is_tagged(entity))
      tags.push_back( (Tag)(i+1) );
  return MB_SUCCESS;
}

ErrorCode BitTag::reserve( unsigned bits )
{
  if (in_use() || bits > 8)
    return MB_FAILURE;
  
  inUse = true;
  requestedBitsPerEntity = bits;
    // store smallest power of two greater than or 
    // equal to the number of bits
  storedBitsPerEntity = 1;
  unsigned ln2storedbits = 0;
  while (storedBitsPerEntity < bits) {
    storedBitsPerEntity *= 2;
    ++ln2storedbits;
  }
  
  // pageShift = log2( ents_per_page() )
  //           = log2( 8 * pageSize / storedBitsPerEntity )
  //           = log2(8) + log2(pageSize) - log2(storedBitsPerEntity)
  //           = 3 + Ln2PageSize - ln2storedbits;
  pageShift = 3 + Ln2PageSize - ln2storedbits;

  return MB_SUCCESS;
}

void BitTag::reset_data()
{
  for (EntityType t = (EntityType)0; t != MBMAXTYPE; ++t) {
    for (size_t i = 0; i < pageList[t].size(); ++i)
      delete pageList[t][i];
    pageList[t].clear();
  }
}

void BitTag::release()
{
  reset_data();
  inUse = false;
}

ErrorCode BitTag::get_bits( const EntityHandle* handles, 
                              int num_handles, 
                              unsigned char* data,
                              const unsigned char* default_value ) const
{
  EntityType type;
  size_t page;
  int offset;
  unsigned char def = default_value ? *default_value : 0;
  for (int i = 0; i < num_handles; ++i) {
    unpack( handles[i], type, page, offset );
    if (pageList[type].size() <= page || !pageList[type][page]) 
      data[i] = def;
    else
      data[i] = pageList[type][page]->get_bits( offset, storedBitsPerEntity );
  }
  return MB_SUCCESS;
}

ErrorCode BitTag::set_bits( const EntityHandle* handles, 
                              int num_handles, 
                              const unsigned char* data,
                              const unsigned char* default_value )
{
  EntityType type;
  size_t page;
  int offset;
  for (int i = 0; i < num_handles; ++i) {
    unpack( handles[i], type, page, offset );
    if (pageList[type].size() <= page)
      pageList[type].resize(page+1, 0);
    if (!pageList[type][page])
      pageList[type][page] = new BitPage( storedBitsPerEntity, 
                              default_value ? *default_value : 0 );
    pageList[type][page]->set_bits( offset, storedBitsPerEntity, data[i] );
  }
  return MB_SUCCESS;
}

ErrorCode BitTag::set_bits( const EntityHandle* handles, 
                            int num_handles, 
                            unsigned char value,
                            const unsigned char* default_value )
{
  EntityType type;
  size_t page;
  int offset;
  for (int i = 0; i < num_handles; ++i) {
    unpack( handles[i], type, page, offset );
    if (pageList[type].size() <= page)
      pageList[type].resize(page+1, 0);
    if (!pageList[type][page])
      pageList[type][page] = new BitPage( storedBitsPerEntity, 
                              default_value ? *default_value : 0 );
    pageList[type][page]->set_bits( offset, storedBitsPerEntity, value );
  }
  return MB_SUCCESS;
}

ErrorCode BitTag::clear_bits( const EntityHandle* handles, 
                                int num_handles, 
                                const unsigned char* default_value )
{
  EntityType type;
  size_t page;
  int offset;
  const unsigned char val = default_value ? *default_value : 0;
  for (int i = 0; i < num_handles; ++i) {
    unpack( handles[i], type, page, offset );
    if (pageList[type].size() > page && pageList[type][page])
      pageList[type][page]->set_bits( offset, storedBitsPerEntity, val );
  }
  return MB_SUCCESS;
}

ErrorCode BitTag::get_bits( const Range& handles, 
                              unsigned char* data,
                              const unsigned char* default_value ) const
{
  EntityType type;
  EntityID count;
  size_t page;
  int offset, per_page = ents_per_page();
  unsigned char def = default_value ? *default_value : 0;
  Range::const_pair_iterator i;
  for (i = handles.const_pair_begin(); i != handles.const_pair_end(); ++i) {
    unpack( i->first, type, page, offset );
    assert(TYPE_FROM_HANDLE(i->second) == type); // should be true because id of zero is never used
    count = i->second - i->first + 1;
    if (page >= pageList[type].size()) {
      memset( data, def, count );
      data += count;
      continue;
    }
    
    while (count) {
      size_t pcount = std::min( (EntityID)(per_page - offset), count );
      if (pageList[type][page])
        pageList[type][page]->get_bits( offset, pcount, storedBitsPerEntity, data );
      else
        memset( data, def, pcount );
      data += pcount;
      count -= pcount; 
      offset = 0;
      ++page;
    }
  }
  return MB_SUCCESS;
}

ErrorCode BitTag::set_bits( const Range& handles, 
                              const unsigned char* data, 
                              const unsigned char* default_value )
{
  EntityType type;
  EntityID count;
  size_t page;
  int offset, per_page = ents_per_page();
  unsigned char def = default_value ? *default_value : 0;
  Range::const_pair_iterator i;
  for (i = handles.const_pair_begin(); i != handles.const_pair_end(); ++i) {
    unpack( i->first, type, page, offset );
    assert(TYPE_FROM_HANDLE(i->second) == type); // should be true because id of zero is never used
    count = i->second - i->first + 1;
    
    while (count) {
      if (page >= pageList[type].size())
        pageList[type].resize( page+1, 0 );
      if (!pageList[type][page])
        pageList[type][page] = new BitPage( storedBitsPerEntity, def );

      size_t pcount = std::min( (EntityID)(per_page - offset), count );
      pageList[type][page]->set_bits( offset, pcount, storedBitsPerEntity, data );
      data += pcount;
      count -= pcount; 
      offset = 0;
      ++page;
    }
  }
  return MB_SUCCESS;
}

ErrorCode BitTag::set_bits( const Range& handles, 
                            unsigned char value, 
                            const unsigned char* default_value )
{
  EntityType type;
  EntityID count;
  size_t page;
  int offset, per_page = ents_per_page();
  unsigned char def = default_value ? *default_value : 0;
  Range::const_pair_iterator i;
  for (i = handles.const_pair_begin(); i != handles.const_pair_end(); ++i) {
    unpack( i->first, type, page, offset );
    assert(TYPE_FROM_HANDLE(i->second) == type); // should be true because id of zero is never used
    count = i->second - i->first + 1;
    
    while (count) {
      if (page >= pageList[type].size())
        pageList[type].resize( page+1, 0 );
      if (!pageList[type][page])
        pageList[type][page] = new BitPage( storedBitsPerEntity, def );

      size_t pcount = std::min( (EntityID)(per_page - offset), count );
      pageList[type][page]->set_bits( offset, pcount, storedBitsPerEntity, value );
      count -= pcount; 
      offset = 0;
      ++page;
    }
  }
  return MB_SUCCESS;
}

ErrorCode BitTag::clear_bits( const Range& handles, 
                                const unsigned char* default_value )
{
  EntityType type;
  EntityID count;
  size_t page;
  int offset, per_page = ents_per_page();
  unsigned char val = default_value ? *default_value : 0;
  Range::const_pair_iterator i;
  for (i = handles.const_pair_begin(); i != handles.const_pair_end(); ++i) {
    unpack( i->first, type, page, offset );
    assert(TYPE_FROM_HANDLE(i->second) == type); // should be true because id of zero is never used
    count = i->second - i->first + 1;
    
    while (count) {
      size_t pcount = std::min( (EntityID)(per_page - offset), count );
      if (page < pageList[type].size() && pageList[type][page])
        pageList[type][page]->set_bits( offset, pcount, storedBitsPerEntity, val );
      count -= pcount; 
      offset = 0;
      ++page;
    }
  }
  return MB_SUCCESS;
}

ErrorCode BitTag::get_entities( Range& entities ) const
{
  ErrorCode rval = MB_SUCCESS;
  EntityType type = MBMAXTYPE;
  while (type--) {
    rval = get_entities( type, entities );
    if (MB_SUCCESS != rval)
      break;
  }
  return rval;
}
    

ErrorCode BitTag::get_entities( EntityType type, Range& entities ) const
{
  if (MBMAXTYPE == type)
    return get_entities(entities);

  const int per_page = ents_per_page();
  Range::iterator hint = entities.begin();
  for (size_t i = 0; i < pageList[type].size(); ++i) {
    if (pageList[type][i]) {
      EntityID id = i * per_page;
      EntityHandle h = CREATE_HANDLE( type, id );
      EntityHandle last = h + per_page - 1;
        // never zero ID
      if (!id) ++h;
      hint = entities.insert( hint, h, last );
    }
  }
  return MB_SUCCESS;
}

ErrorCode BitTag::get_entities( const Range& range, 
                                  EntityType in_type,
                                  Range& entities ) const
{
  EntityType type;
  EntityID count;
  size_t page;
  int offset, per_page = ents_per_page();
  Range::const_iterator j, i = in_type == MBMAXTYPE ? range.begin() : range.lower_bound( in_type );
  EntityHandle h;
  while (i != range.end()) {
    h = *i;
    unpack( h, type, page, offset );
    if (in_type != MBMAXTYPE && type != in_type)
      break;
    
    i = i.end_of_block();
    count = *i - h + 1;
    ++i;
    while (count > 0) {
      EntityID pcount = std::min( count, (EntityID)(per_page - offset) );
      if (page < pageList[type].size() && pageList[type][page]) 
        entities.insert( h, h+pcount-1 );
    
      count -= pcount;
      h += pcount;
      assert(TYPE_FROM_HANDLE(h) == type);
      offset = 0;
      ++page;
    }
  }
  return MB_SUCCESS;
}

ErrorCode BitTag::get_entities_with_bits( EntityType type, 
                                            Range& entities,
                                            unsigned char bits ) const
{
  if (MBMAXTYPE == type) {
    ErrorCode rval;
    for (--type; type >= MBVERTEX; --type) {
      rval = get_entities_with_bits( type, entities, bits );
      if (MB_SUCCESS != rval)
        return rval;
    }
    return MB_SUCCESS;
  }

  const int per_page = ents_per_page();
  for (size_t i = 0; i < pageList[type].size(); ++i) {
    if (pageList[type][i]) {
      EntityID id = i * per_page;
      EntityHandle h = CREATE_HANDLE( type, id );
      int off = !i; // never zero ID
      pageList[type][i]->search( bits, off, per_page-off, storedBitsPerEntity, 
                                 entities, h+off );
    }
  }
  return MB_SUCCESS;
}

ErrorCode BitTag::get_entities_with_bits( const Range &range,
                                            EntityType in_type, 
                                            Range& entities,
                                            unsigned char bits ) const
{
  if (MBMAXTYPE == in_type) {
    ErrorCode rval;
    for (--in_type; in_type >= MBVERTEX; --in_type) {
      rval = get_entities_with_bits( range, in_type, entities, bits );
      if (MB_SUCCESS != rval)
        return rval;
    }
    return MB_SUCCESS;
  }

  EntityType type;
  EntityID count;
  size_t page;
  int offset, per_page = ents_per_page();
  Range::const_iterator j, i = range.lower_bound( in_type );
  EntityHandle h;
  while (i != range.end()) {
    h = *i;
    unpack( h, type, page, offset );
    if (type != in_type)
      break;
    
    i = i.end_of_block();
    count = *i - h + 1;
    ++i;
    while (count > 0) {
      EntityID pcount = std::min( count, (EntityID)(per_page - offset) );
      if (page < pageList[type].size() && pageList[type][page]) 
        pageList[type][page]->search( bits, offset, pcount, 
                                      storedBitsPerEntity,
                                      entities, h );
    
      count -= pcount;
      h += pcount;
      assert(TYPE_FROM_HANDLE(h) == type);
      offset = 0;
      ++page;
    }
  }
  return MB_SUCCESS;
}

ErrorCode BitTag::get_number_entities( EntityType type,
                                         int& num_entities ) const
{
  num_entities = 0;
  if (MBMAXTYPE == type) {
    int tmp;
    ErrorCode rval;
    for (--type; type >= MBVERTEX; --type) {
      rval = get_number_entities( type, tmp );
      if (MB_SUCCESS != rval)
        return rval;
      num_entities += tmp;
    }
    return MB_SUCCESS;
  }

  const int per_page = ents_per_page();
  if (pageList[type].empty())
    return MB_SUCCESS;
  if (pageList[type][0])
    num_entities = per_page - 1; // never zero ID
  for (size_t i = 1; i < pageList[type].size(); ++i) 
    if (pageList[type][i])
      num_entities += per_page;
  return MB_SUCCESS;
}

ErrorCode BitTag::get_number_entities( const Range &range,
                                         EntityType type,
                                         int& num_entities ) const
{
  Range tmp;
  ErrorCode result = get_entities( range, type, tmp );
  num_entities = tmp.size();
  return result;
}

ErrorCode BitTag::get_memory_use( unsigned long& total,
                                  unsigned long& per_entity ) const
{
  per_entity = (storedBitsPerEntity > 4); // cannot return fraction of bytes, so round
  total = 0;
  for (EntityType t = (EntityType)0; t < MBMAXTYPE; ++t) {
    total += pageList[t].capacity() * sizeof(BitPage*);
    for (size_t i = 0; i < pageList[t].size(); ++i)
      if (pageList[t][i])
        total += sizeof(BitPage);
  }
  return MB_SUCCESS;
}
  
} // namespace moab
