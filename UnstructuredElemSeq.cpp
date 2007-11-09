#include "UnstructuredElemSeq.hpp"
#include "SequenceData.hpp"
#include "MBCN.hpp"

UnstructuredElemSeq::UnstructuredElemSeq( MBEntityHandle start_handle, 
                                          MBEntityID entity_count, 
                                          unsigned nodes_per_entity,
                                          SequenceData* data )
  : ElementSequence( start_handle, entity_count, nodes_per_entity, data )
  {}
  

UnstructuredElemSeq::UnstructuredElemSeq( MBEntityHandle start_handle, 
                                          MBEntityID entity_count, 
                                          unsigned nodes_per_entity,
                                          MBEntityID data_size )
  : ElementSequence( start_handle, entity_count, nodes_per_entity,
                      new SequenceData( 1, start_handle, start_handle + data_size - 1))
{
  data()->create_sequence_data( 0, nodes_per_entity * sizeof(MBEntityHandle) );
}


UnstructuredElemSeq::~UnstructuredElemSeq()
{}


int UnstructuredElemSeq::values_per_entity() const
  { return nodes_per_element(); }


EntitySequence*
UnstructuredElemSeq::split( MBEntityHandle here )
{
  if (here <= start_handle() || here > end_handle())
    return 0;
  
  return new UnstructuredElemSeq( *this, here );
}


SequenceData*
UnstructuredElemSeq::create_data_subset( MBEntityHandle start, MBEntityHandle end ) const
{
  int esize = nodes_per_element() * sizeof(MBEntityHandle);
  return data()->subset(start, end, &esize, 0 );
}


void
UnstructuredElemSeq::get_const_memory_use( unsigned long& bytes_per_entity,
                                           unsigned long& size_of_sequence ) const
{
  bytes_per_entity = nodes_per_element() * sizeof(MBEntityHandle);
  size_of_sequence = sizeof(this);
}

MBErrorCode
UnstructuredElemSeq::get_connectivity( MBEntityHandle handle,
                                       std::vector<MBEntityHandle>& connect,
                                       bool topological ) const
{
  MBEntityHandle const* conn = get_array() + nodes_per_element() * (handle - start_handle());
  int len = topological ? MBCN::VerticesPerEntity(type()) : nodes_per_element();
  connect.reserve( connect.size() + len );
  std::copy( conn, conn+len, std::back_inserter( connect ) );
  return MB_SUCCESS;
}


MBErrorCode
UnstructuredElemSeq::get_connectivity( MBEntityHandle handle,
                                       MBEntityHandle const*& conn_ptr,
                                       int& len,
                                       bool topological,
                                       std::vector<MBEntityHandle>* ) const
{
  conn_ptr = get_array() + nodes_per_element() * (handle - start_handle());
  len = topological ? MBCN::VerticesPerEntity(type()) : nodes_per_element();
  return MB_SUCCESS;
}


MBErrorCode
UnstructuredElemSeq::set_connectivity( MBEntityHandle handle,
                                       MBEntityHandle const* connect,
                                       int connect_length )
{
  if ((unsigned)connect_length != nodes_per_element())
    return MB_INDEX_OUT_OF_RANGE;
  MBEntityHandle* conn_ptr = get_array() + nodes_per_element() * (handle - start_handle());
  std::copy( connect, connect+connect_length, conn_ptr );
  return MB_SUCCESS;
}


MBEntityHandle* UnstructuredElemSeq::get_connectivity_array()
  { return get_array(); }

MBErrorCode UnstructuredElemSeq::push_back( MBEntityID count )
  { return EntitySequence::append_entities(count); }

MBErrorCode UnstructuredElemSeq::push_front( MBEntityID count )
  { return EntitySequence::prepend_entities(count); }

