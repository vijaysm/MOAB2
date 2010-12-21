/** \file   MeshTag.cpp
 *  \author Jason Kraftcheck 
 *  \date   2010-12-14
 */

#include "moab/Interface.hpp"
#include "MeshTag.hpp"
#include "SysUtil.hpp"

namespace moab {

static inline bool all_root_set( const EntityHandle* array, size_t len )
{
  for (size_t i = 0; i < len; ++i)
    if (array[i])
      return false;
  return true;
}

MeshTag::MeshTag( const char * name, 
                  int size, 
                  DataType type, 
                  const void * default_value,
                  int default_value_size)
 : TagInfo( name, size, type, default_value, default_value_size )
 {}
  
MeshTag::~MeshTag() {}

TagType MeshTag::get_storage_type() const 
  { return MB_TAG_MESH; }

ErrorCode MeshTag::release_all_data( SequenceManager*, bool )
  { return MB_SUCCESS; }

ErrorCode MeshTag::get_data( const SequenceManager*,
                             const EntityHandle* entities,
                             size_t num_entities,
                             void* data ) const
{
  if (!all_root_set( entities, num_entities ))
    return MB_TAG_NOT_FOUND;

  const void* ptr;
  int len;

  if (!mValue.empty()) {
    ptr = &mValue[0];
    len = mValue.size();
  }
  else if (get_default_value()) {
    ptr = get_default_value();
    len = get_default_value_size();
  }
  else {
    return MB_TAG_NOT_FOUND;
  }

  SysUtil::setmem( data, ptr, len, num_entities );
  return MB_SUCCESS;
}
  

ErrorCode MeshTag::get_data( const SequenceManager*,
                             const Range&,
                             void* ) const
{
  if (variable_length())
    return MB_VARIABLE_DATA_LENGTH;
  else
    return MB_TAG_NOT_FOUND;
}
                      
ErrorCode MeshTag::get_data( const SequenceManager*,
                             const EntityHandle* entities,
                             size_t num_entities,
                             const void** data_ptrs,
                             int* data_lengths ) const 
{
  const void* ptr;
  int len;

  if (!mValue.empty()) {
    ptr = &mValue[0];
    len = mValue.size();
  }
  else if (get_default_value()) {
    ptr = get_default_value();
    len = get_default_value_size();
  }
  else {
    return MB_TAG_NOT_FOUND;
  }
    
  for (size_t i = 0; i < num_entities; ++i) {
    if (entities[i]) return MB_TAG_NOT_FOUND; // not root set
    data_ptrs[i] = ptr;
    if (data_lengths)
      data_lengths[i] = len;
  }
  return MB_SUCCESS;
}
                      
                      
ErrorCode MeshTag::get_data( const SequenceManager*,
                             const Range&,
                             const void**,
                             int* ) const
{
  return MB_TAG_NOT_FOUND;
}
  
ErrorCode MeshTag::set_data( SequenceManager*,
                             const EntityHandle* entities,
                             size_t num_entities,
                             const void* data )
{
  if (variable_length())
    return MB_VARIABLE_DATA_LENGTH;
  if (!all_root_set( entities, num_entities ))
    return MB_TAG_NOT_FOUND;
  
  if (num_entities > 0) {
    mValue.resize( get_size() );
    const unsigned char* bytes = reinterpret_cast<const unsigned char*>(data);
    memcpy( &mValue[0], bytes + get_size() * (num_entities - 1), get_size() );
  }
  return MB_SUCCESS;
}
 
ErrorCode MeshTag::set_data( SequenceManager*,
                             const Range&,
                             const void* )
{
  if (variable_length())
    return MB_VARIABLE_DATA_LENGTH;
  else
    return MB_TYPE_OUT_OF_RANGE;
}

ErrorCode MeshTag::set_data( SequenceManager*,
                             const EntityHandle* entities,
                             size_t num_entities,
                             void const* const* data_ptrs,
                             const int* data_lengths )
{
  if (!all_root_set( entities, num_entities ))
    return MB_TAG_NOT_FOUND;
  
  ErrorCode valid = validate_lengths( data_lengths, num_entities );
  if (MB_SUCCESS != valid)
    return valid;
  
  if (num_entities > 0) {
    mValue.resize( data_lengths[num_entities-1] );
    memcpy( &mValue[0], data_ptrs[num_entities-1], mValue.size() );
  }
  return MB_SUCCESS;
}
                      
                      
ErrorCode MeshTag::set_data( SequenceManager*,
                             const Range&,
                             void const* const*,
                             const int* )
{
  return MB_TYPE_OUT_OF_RANGE;
}

ErrorCode MeshTag::clear_data( SequenceManager*,
                               const EntityHandle* entities,
                               size_t num_entities,
                               const void* value_ptr,
                               int value_len )
{
  if (!all_root_set( entities, num_entities ))
    return MB_TAG_NOT_FOUND;
  
  ErrorCode valid = validate_lengths( value_len ? &value_len : 0, 1 );
  if (MB_SUCCESS != valid)
    return valid;
  
  if (num_entities > 0) {
    mValue.resize( value_len );
    memcpy( &mValue[0], value_ptr, value_len );
  }

  return MB_SUCCESS;
}

ErrorCode MeshTag::clear_data( SequenceManager*,
                               const Range&,
                               const void*,
                               int )
{
  return MB_TYPE_OUT_OF_RANGE;
}

ErrorCode MeshTag::remove_data( SequenceManager*,
                                const EntityHandle* entities,
                                size_t num_entities )
{
  if (!all_root_set( entities, num_entities ))
    return MB_TAG_NOT_FOUND;
  
  if (num_entities)
    mValue.clear();;
  return MB_SUCCESS;
}

ErrorCode MeshTag::remove_data( SequenceManager*,
                                const Range& )
{
  return MB_TAG_NOT_FOUND;
}

ErrorCode MeshTag::remove_all_entity_data( SequenceManager* )
{
  return MB_SUCCESS;
}

ErrorCode MeshTag::tag_iterate( SequenceManager*,
                                Range::iterator&,
                                const Range::iterator&,
                                void*& )
{
  return MB_TAG_NOT_FOUND;
}

ErrorCode MeshTag::get_tagged_entities( const SequenceManager*,
                                        Range&,
                                        EntityType,
                                        const Range* ) const
{
  return MB_SUCCESS;
}

ErrorCode MeshTag::num_tagged_entities( const SequenceManager*,
                                        size_t&,
                                        EntityType,
                                        const Range* ) const
{
  return MB_SUCCESS;
}

ErrorCode MeshTag::find_entities_with_value( const SequenceManager*,
                                             Range&,
                                             const void*,
                                             int,
                                             EntityType,
                                             const Range* ) const
{
  return MB_SUCCESS;
}

bool MeshTag::is_tagged( const SequenceManager*, EntityHandle h ) const
  { return !h && !mValue.empty(); }

ErrorCode MeshTag::get_memory_use( const SequenceManager*,
                                   unsigned long& total,
                                   unsigned long& per_entity ) const
{
  total = TagInfo::get_memory_use() + sizeof(*this) + mValue.size();
  per_entity = 0;
  return MB_SUCCESS;
}


} // namespace moab
