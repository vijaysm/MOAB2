/*!
 *  \class   PolyEntitySequence
 *  \authors Tim Tautges
 *  \date    2/04
 *  \brief   PolyEntitySequence is a sequence of polygons or polyhedra
 *        These entities represent connectivity as a list of
 *        handles, which represent vertices for polygons, or faces for polyhedra.  The
 *        sequence stores the number of vertices for each entity by storing the ending
 *        index for its connectivity array.
 *
 */ 

#ifndef POLY_ENTITY_SEQUENCE_HPP
#define POLY_ENTITY_SEQUENCE_HPP

#ifndef IS_BUILDING_MB
#error "PolyEntitySequence.hpp isn't supposed to be included into an application"
#endif

#include "EntitySequence.hpp"

class PolyEntitySequence : public ElementEntitySequence
{
public:
  
  PolyEntitySequence(EntitySequenceManager* seq_manager, MBEntityHandle start_handle, int num_entities,
                     int nodes_per_element, bool all_handles_used,
                     bool allocate_connect = true);
  virtual ~PolyEntitySequence();

  MBEntityHandle get_unused_handle();
  void free_handle(MBEntityHandle) {}

  virtual MBErrorCode get_connectivity(MBEntityHandle entity, 
                                       std::vector<MBEntityHandle>& connectivity,
                                       const bool topological_connectivity = false) const;
  virtual MBErrorCode get_connectivity(MBEntityHandle entity, 
                                       const MBEntityHandle*& connectivity,
                                       int &num_vertices,
                                       const bool topological_connectivity = false) const;

  MBErrorCode set_connectivity(MBEntityHandle entity, const MBEntityHandle *conn,
                               const int num_vertices);

  virtual MBErrorCode get_connectivity_array(MBEntityHandle*& conn_array);
  
  MBErrorCode PolyEntitySequence::get_index_array(int*& index_array);
  
  virtual MBErrorCode split(MBEntityHandle , MBEntitySequence*& )
    {return MB_FAILURE;}

  // reallocated the sequence to hold extra/less nodes, pass in what you want, and will return whether it needed
  // reallocate space for those nodes
  MBErrorCode convert_realloc(bool&, bool&, bool&, 
                              MBCore* , MBTag ) {return MB_FAILURE;}
  
  bool has_mid_edge_nodes() const {return false;}
  bool has_mid_face_nodes() const {return false;}
  bool has_mid_volume_nodes() const {return false;}

  virtual bool is_valid_entity(MBEntityHandle entity) const;

  MBErrorCode add_entity(const MBEntityHandle *conn,
                         const int num_conn, MBEntityHandle &handle);
  

  //! get entities in this range, will not add unused entities
  virtual void get_entities(MBRange& entities) const;

private:

  std::vector<MBEntityHandle> polyConn;

  std::vector<int> mLastIndex;
};

inline MBEntityHandle PolyEntitySequence::get_unused_handle() 
{
  return mNumEntities + mStartEntityHandle;
}

#endif
