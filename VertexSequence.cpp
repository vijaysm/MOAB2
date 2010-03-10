#include "VertexSequence.hpp"

VertexSequence::~VertexSequence() {}

EntitySequence* VertexSequence::split( MBEntityHandle here )
  { return new VertexSequence( *this, here ); }

SequenceData* VertexSequence::create_data_subset( MBEntityHandle start,
                                                  MBEntityHandle end ) const
{
  const int sizes[] = { sizeof(double), sizeof(double), sizeof(double) };
  return data()->subset(start, end, sizes );
}
  
MBErrorCode VertexSequence::push_back( MBEntityID count )
  { return EntitySequence::append_entities(count); }

MBErrorCode VertexSequence::push_front( MBEntityID count )
  { return EntitySequence::prepend_entities(count); }

void VertexSequence::get_const_memory_use( unsigned long& per_ent, unsigned long& seq ) const
{
  per_ent = 3 * sizeof(double);
  seq = sizeof(*this);
}
