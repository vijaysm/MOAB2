#include "PolyElementSeq.hpp"

PolyElementSeq::~PolyElementSeq() {}
  
EntitySequence* PolyElementSeq::split( MBEntityHandle here )
  { return new PolyElementSeq( *this, here ); }
                       

MBErrorCode
PolyElementSeq::get_connectivity( MBEntityHandle handle,
                                  std::vector<MBEntityHandle>& connect,
                                  bool ) const
{
  MBEntityHandle const* conn = get_array() + nodes_per_element() * (handle - start_handle());
  int len = nodes_per_element();
  connect.reserve( connect.size() + len );
  std::copy( conn, conn+len, std::back_inserter( connect ) );
  return MB_SUCCESS;
}


MBErrorCode
PolyElementSeq::get_connectivity( MBEntityHandle handle,
                                  MBEntityHandle const*& conn_ptr,
                                  int& len,
                                  bool,
                                  std::vector<MBEntityHandle>* ) const
{
  conn_ptr = get_array() + nodes_per_element() * (handle - start_handle());
  len = nodes_per_element();
  return MB_SUCCESS;
}
