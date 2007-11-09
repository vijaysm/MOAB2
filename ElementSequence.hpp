#ifndef ELEMENT_SEQUENCE_HPP
#define ELEMENT_SEQUENCE_HPP

#include "EntitySequence.hpp"
#include "SequenceData.hpp"
#include "MBCN.hpp"

class ElementSequence : public EntitySequence
{
public:
  
  ElementSequence( MBEntityHandle start,
                   MBEntityID count,
                   unsigned int nodes_per_element,
                   SequenceData* data )
    : EntitySequence( start, count, data ), 
      nodesPerElement(nodes_per_element)
    {}
                   
  virtual ~ElementSequence() {}
  
  inline unsigned int nodes_per_element() const { return nodesPerElement; }
  
  virtual MBErrorCode get_connectivity( MBEntityHandle handle,
                                        std::vector<MBEntityHandle>& connect,
                                        bool topological = false ) const = 0;
  
  virtual MBErrorCode get_connectivity( MBEntityHandle handle,
                                        MBEntityHandle const*& connect,
                                        int &connect_length,
                                        bool topological = false,
                                        std::vector<MBEntityHandle>* storage = 0
                                       ) const = 0;

  virtual MBErrorCode set_connectivity( MBEntityHandle handle,
                                        MBEntityHandle const* connect,
                                        int connect_length ) = 0;

  inline MBEntityHandle const* get_connectivity_array() const;
  
  virtual MBEntityHandle* get_connectivity_array() = 0;
  
  inline bool has_mid_edge_nodes() const;
  inline bool has_mid_face_nodes() const;
  inline bool has_mid_volume_nodes() const;

protected:

  ElementSequence( ElementSequence& split_from, MBEntityHandle here )
    : EntitySequence( split_from, here ),
      nodesPerElement( split_from.nodesPerElement )
    {}

private:
  
  unsigned nodesPerElement;
};

inline MBEntityHandle const*
ElementSequence::get_connectivity_array() const
  { return const_cast<ElementSequence*>(this)->get_connectivity_array(); }

inline bool
ElementSequence::has_mid_edge_nodes() const
  { return MBCN::HasMidEdgeNodes( type(), nodes_per_element() ); }

inline bool
ElementSequence::has_mid_face_nodes() const
  { return MBCN::HasMidFaceNodes( type(), nodes_per_element() ); }

inline bool
ElementSequence::has_mid_volume_nodes() const
  { return MBCN::HasMidRegionNodes( type(), nodes_per_element() ); }
  
#endif
