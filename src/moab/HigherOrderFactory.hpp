/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/*!
 *  \class   HigherOrderFactory
 *  \authors Clinton Stimpson
 *  \date    11/25/02
 *  \brief   
 *          
 */ 

#ifndef HIGHER_ORDER_FACTORY_HPP
#define HIGHER_ORDER_FACTORY_HPP

#ifndef IS_BUILDING_MB
#error "HigherOrderFactory.hpp isn't supposed to be included into an application"
#endif

#include "MBInterface.hpp"
class ElementSequence;
class MBCore;

class HigherOrderFactory
{
public:

  HigherOrderFactory(MBCore* MB, MBInterface::HONodeAddedRemoved* function_object);
  ~HigherOrderFactory();

  MBErrorCode convert(const MBEntityHandle meshset, const bool mid_edge_nodes, 
                       const bool mid_face_nodes, const bool mid_volume_nodes);

  MBErrorCode convert( const MBRange& entities, const bool mid_edge_nodes, 
                       const bool mid_face_nodes, const bool mid_volume_nodes);

  unsigned char mNodeMap[MBMAXTYPE][8][8];

private:

  //static bool mMapInitialized;
  void initialize_map();

  MBCore* mMB;
  MBInterface::HONodeAddedRemoved* mHONodeAddedRemoved;

  MBErrorCode convert_sequence(ElementSequence* sequence, 
                               MBEntityHandle sequence_subset_start,
                               MBEntityHandle sequence_subset_end,
                               bool mid_edge_nodes, 
                               bool mid_face_nodes,
                               bool mid_volume_nodes);
  MBErrorCode add_mid_edge_nodes(ElementSequence*);
  MBErrorCode add_mid_face_nodes(ElementSequence*);
  MBErrorCode add_mid_volume_nodes(ElementSequence*);

  //! returns the handle of the first center node found between the two corner nodes.
  //! returns zero if none found
  //! entities that share those two corner nodes and have space allocated for mid-edge nodes are returned in a vector
  MBEntityHandle center_node_exist(MBEntityHandle corner1, MBEntityHandle corner2,
         std::vector<MBEntityHandle>& adj_entities);
  
  //! returns the handle of the first center node found between the 3-4 corner nodes.
  //! set the last node to zero if you want only 3 nodes
  //! returns zero if none found
  //! entities that share those corner nodes and have space allocated for mid face nodes are returned in a vector
  MBEntityHandle center_node_exist(MBEntityHandle corners[4], std::vector<MBEntityHandle>& adj_entities);

  //! adds a center node to element between corner nodes, returns success
  bool add_center_node(MBEntityType type, MBEntityHandle* element_conn, int conn_size, 
      MBEntityHandle corner_node1, MBEntityHandle corner_node2, MBEntityHandle center_node);


  MBErrorCode copy_corner_nodes( ElementSequence* src, ElementSequence* dst );
  MBErrorCode copy_mid_edge_nodes( ElementSequence* src, ElementSequence* dst ); 
  MBErrorCode copy_mid_face_nodes( ElementSequence* src, ElementSequence* dst ); 
  MBErrorCode copy_mid_volume_nodes( ElementSequence* src, ElementSequence* dst ); 
  MBErrorCode copy_nodes( ElementSequence* src, 
                          ElementSequence* dst,
                          unsigned nodes_per_elem_to_copy,
                          unsigned src_conn_offset,
                          unsigned dst_conn_offset ); 

  MBErrorCode remove_mid_edge_nodes( ElementSequence* seq, 
                                     MBEntityHandle start,
                                     MBEntityHandle stop,
                                     MBTag deletable_ndoes );
  MBErrorCode remove_mid_face_nodes( ElementSequence* seq, 
                                     MBEntityHandle start,
                                     MBEntityHandle stop,
                                     MBTag deletable_ndoes );
  MBErrorCode remove_mid_volume_nodes( ElementSequence* seq, 
                                       MBEntityHandle start,
                                       MBEntityHandle stop,
                                       MBTag deletable_ndoes );
  MBErrorCode remove_ho_nodes( ElementSequence* sequence,
                               MBEntityHandle subset_start_handle,
                               MBEntityHandle subset_end_handle,
                               int nodes_per_elem_to_remove,
                               int elem_conn_offset_to_remove,
                               MBTag deletable_nodes );
  bool tag_for_deletion( MBEntityHandle element_with_node,
                         int node_index_in_elem_connectivity,
                         ElementSequence* sequence );
};


#endif

