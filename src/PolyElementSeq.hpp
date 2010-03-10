#ifndef POLY_ELEMENT_SEQ_HPP
#define POLY_ELEMENT_SEQ_HPP

#include "UnstructuredElemSeq.hpp"


class PolyElementSeq : public UnstructuredElemSeq
{
public:
  PolyElementSeq( MBEntityHandle start_handle, 
                  MBEntityID entity_count, 
                  unsigned nodes_per_entity,
                  SequenceData* data )
    : UnstructuredElemSeq( start_handle, entity_count, nodes_per_entity, data )
    {}

  PolyElementSeq( MBEntityHandle start_handle, 
                  MBEntityID entity_count, 
                  unsigned nodes_per_entity,
                  MBEntityID sequence_data_size)
    : UnstructuredElemSeq( start_handle, entity_count, nodes_per_entity, sequence_data_size )
    {}

  virtual ~PolyElementSeq();
  
  virtual EntitySequence* split( MBEntityHandle here );
                       
  virtual MBErrorCode get_connectivity( MBEntityHandle handle,
                                        std::vector<MBEntityHandle>& connect,
                                        bool topological = false ) const;
  
  virtual MBErrorCode get_connectivity( MBEntityHandle handle,
                                        MBEntityHandle const*& connect,
                                        int &connect_length,
                                        bool topological = false,
                                        std::vector<MBEntityHandle>* storage = 0
                                       ) const;

protected:

  PolyElementSeq( PolyElementSeq& split_from, MBEntityHandle here )
    : UnstructuredElemSeq( split_from, here )
   {}
};

#endif
