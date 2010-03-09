#ifndef UNSTRUCTURED_ELEM_SEQ_HPP
#define UNSTRUCTURED_ELEM_SEQ_HPP

#include "ElementSequence.hpp"
#include "SequenceData.hpp"


class UnstructuredElemSeq : public ElementSequence
{
public:

  UnstructuredElemSeq( MBEntityHandle start_handle, 
                       MBEntityID entity_count, 
                       unsigned nodes_per_entity,
                       SequenceData* data );

  UnstructuredElemSeq( MBEntityHandle start_handle, 
                       MBEntityID entity_count, 
                       unsigned nodes_per_entity,
                       MBEntityID sequence_data_size);

  virtual ~UnstructuredElemSeq();

  int values_per_entity() const;
  
  virtual EntitySequence* split( MBEntityHandle here );
  
  SequenceData* create_data_subset( MBEntityHandle start, MBEntityHandle end ) const;
                       
  virtual MBErrorCode get_connectivity( MBEntityHandle handle,
                                        std::vector<MBEntityHandle>& connect,
                                        bool topological = false ) const;
  
  virtual MBErrorCode get_connectivity( MBEntityHandle handle,
                                        MBEntityHandle const*& connect,
                                        int &connect_length,
                                        bool topological = false,
                                        std::vector<MBEntityHandle>* storage = 0
                                       ) const;

  MBErrorCode set_connectivity( MBEntityHandle handle,
                                MBEntityHandle const* connect,
                                int connect_length );
  
  MBEntityHandle* get_connectivity_array();
  
  MBErrorCode push_front( MBEntityID count );
  MBErrorCode push_back ( MBEntityID count );
  
  
  void get_const_memory_use( unsigned long& bytes_per_entity,
                             unsigned long& size_of_sequence ) const;
protected:

  inline MBEntityHandle const* get_array() const
  {
    return reinterpret_cast<MBEntityHandle const*>(data()->get_sequence_data(0))
      + nodes_per_element() * (start_handle() - data()->start_handle());
  }

  inline MBEntityHandle* get_array()
  {
    return reinterpret_cast<MBEntityHandle*>(data()->get_sequence_data(0))
      + nodes_per_element() * (start_handle() - data()->start_handle());
  }

  UnstructuredElemSeq( UnstructuredElemSeq& split_from, MBEntityHandle here )
    : ElementSequence( split_from, here )
   {}
};

#endif

  
