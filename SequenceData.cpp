#include "SequenceData.hpp"
#include <assert.h>

SequenceData::~SequenceData()
{
  for (int i = -numSequenceData; i <= numTagData; ++i)
    free( arraySet[i] );
  free( arraySet - numSequenceData );
}

void* SequenceData::create_data( int index, int bytes_per_ent, void* initial_value )
{  
  char* array = (char*)malloc( bytes_per_ent * size() );
  if (initial_value) {
    memcpy( array, initial_value, bytes_per_ent );
    const size_t last = size() / 2;
    size_t count;
    for (count = 1; count < last; count *= 2) 
      memcpy( array + count * bytes_per_ent, array, count * bytes_per_ent );
    memcpy( array + count * bytes_per_ent, array, (size() - count) * bytes_per_ent ); 
  }
  else {
    memset( array, 0, bytes_per_ent * size() );
  }
  
  arraySet[index] = array;
  return array;
}

void* SequenceData::create_sequence_data( int array_num,
                                          int bytes_per_ent,
                                          void* initial_value )
{
  const int index = -1 - array_num;
  assert( array_num < numSequenceData );
  assert( !arraySet[index] );
  return create_data( index, bytes_per_ent, initial_value );
}


void* SequenceData::create_custom_data( int array_num, size_t total_bytes )
{
  const int index = -1 - array_num;
  assert( array_num < numSequenceData );
  assert( !arraySet[index] );

  void* array = malloc( total_bytes );
  memset( array, 0, total_bytes );
  arraySet[index] = array;
  return array;
}

SequenceData::AdjacencyDataType* SequenceData::allocate_adjacency_data()
{
  assert( !arraySet[0] );
  const size_t s = sizeof(AdjacencyDataType*) * size();
  arraySet[0] = malloc( s );
  memset( arraySet[0], 0, s );
  return reinterpret_cast<AdjacencyDataType*>(arraySet[0]);
}

void* SequenceData::create_tag_data( int tag_num,
                                     int bytes_per_ent,
                                     void* initial_val )
{
  const int index = tag_num + 1;
  if (tag_num >= numTagData) {
    void** list = arraySet - numSequenceData;
    const size_t size = sizeof(void*) * (numSequenceData + tag_num + 1);
    list = (void**)realloc( list, size );
    arraySet = list + numSequenceData;
    memset( arraySet + numTagData, 0, sizeof(void*) * (tag_num + 1 - numTagData) );
  }
  
  assert( !arraySet[index] );
  return create_data( index, bytes_per_ent, initial_val );
}

SequenceData* SequenceData::subset( MBEntityHandle start,
                                    MBEntityHandle end,
                                    const int* sequence_data_sizes,
                                    const int* tag_data_sizes ) const
{
  return new SequenceData( this, start, end, sequence_data_sizes, tag_data_sizes );
}

SequenceData::SequenceData( const SequenceData* from,
                            MBEntityHandle start, 
                            MBEntityHandle end,
                            const int* sequence_data_sizes,
                            const int* tag_data_sizes )
  : numSequenceData( from->numSequenceData ),
    numTagData( from->numTagData ),
    startHandle( start ),
    endHandle( end )
{
  assert( start <= end );
  assert( from != 0 );
  assert( from->start_handle() <= start );
  assert( from->end_handle() >= end );

  void** array = (void**)malloc( sizeof(void*) * (numSequenceData + numTagData + 1) );
  arraySet = array + numSequenceData;
  const size_t offset = start - from->start_handle();
  const size_t count = end - start + 1;
  
  for (int i = 0; i < numSequenceData; ++i)
    copy_data_subset( -1 - i, sequence_data_sizes[i], from->get_sequence_data(i), offset, count );
  copy_data_subset( 0, sizeof(AdjacencyDataType*), from->get_adjacency_data(), offset, count );
  for (int i = 0; i< numTagData; ++i)
    copy_data_subset( 1 + i, tag_data_sizes[i], from->get_tag_data(i), offset, count );
}

void SequenceData::copy_data_subset( int index, 
                                     int size_per_ent, 
                                     const void* source, 
                                     size_t offset, 
                                     size_t count )
{
  if (!source)
    arraySet[index] = 0;
  else {
    arraySet[index] = malloc( count * size_per_ent );
    memcpy( arraySet[index], 
            (const char*)source + offset * size_per_ent, 
            count * size_per_ent );
  }
}

