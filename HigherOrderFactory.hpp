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
#include "MBRange.hpp"
class ElementEntitySequence;
class MBCore;

class HigherOrderFactory
{
public:

  HigherOrderFactory(MBCore* MB, MBInterface::HONodeAddedRemoved* function_object);
  ~HigherOrderFactory();

  MBErrorCode convert(const MBEntityHandle meshset, const bool mid_edge_nodes, 
                       const bool mid_face_nodes, const bool mid_volume_nodes);

  unsigned char mNodeMap[MBMAXTYPE][8][8];

private:

  //static bool mMapInitialized;
  void initialize_map();

  MBCore* mMB;
  MBInterface::HONodeAddedRemoved* mHONodeAddedRemoved;

  MBErrorCode convert_sequence(ElementEntitySequence*, const bool mid_edge_nodes, 
                                const bool mid_face_nodes, const bool mid_volume_nodes);
  MBErrorCode add_mid_edge_nodes(ElementEntitySequence*);
  MBErrorCode add_mid_face_nodes(ElementEntitySequence*);
  MBErrorCode add_mid_volume_nodes(ElementEntitySequence*);

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

};


#endif

