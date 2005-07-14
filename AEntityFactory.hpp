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

#ifndef AENTITY_FACTORY_HPP
#define AENTITY_FACTORY_HPP

#ifndef IS_BUILDING_MB
#error "AEntityFactory.hpp isn't supposed to be included into an application"
#endif

#include "MBInterface.hpp"
#include "DenseTagCollections.hpp"
#include <vector>

typedef std::vector<MBEntityHandle> MBAdjacencyVector;

//! class AEntityFactory
class AEntityFactory 
{
public:

  //! require an MBInterface object in order to access tags on that interface
  AEntityFactory(MBInterface *mdb);

  //! destructor
  ~AEntityFactory();
  
//! add an adjacency from from_ent to to_ent; if both_ways is true, add one
//! in reverse too
//! NOTE: this function is defined even though we may only be implementing
//! vertex-based up-adjacencies
  MBErrorCode add_adjacency(MBEntityHandle from_ent,
                             MBEntityHandle to_ent,
                             const bool both_ways = false);

//! remove an adjacency from from the base_entity.
  MBErrorCode remove_adjacency(MBEntityHandle base_entity,
                                MBEntityHandle adjacency_to_remove);

//! remove all adjacencies from from the base_entity.
  MBErrorCode remove_all_adjacencies(MBEntityHandle base_entity,
                                     const bool delete_adj_list = false);
  
//! get the elements contained by source_entity, of
//! type target_type, passing back in target_entities; if create_if_missing
//! is true and no entity is found, one is created; if create_adjacency_option
//! is >= 0, adjacencies from entities of that dimension to each target_entity
//! are created (this function uses AEntityFactory::get_element for each element)
  MBErrorCode get_elements(MBEntityHandle source_entity,
                            const unsigned int target_dimension,
                            std::vector<MBEntityHandle> &target_entities,
                            const bool create_if_missing,
                            const int create_adjacency_option = -1);

    //! get the vertices for a polyhedron (special implementation because for polyhedra
    //! connectivity array stores faces)
  MBErrorCode get_polyhedron_vertices(const MBEntityHandle source_entity, 
                                      std::vector<MBEntityHandle> &target_entities);
  
//! get the meshsets that are in source_entitiy's adjacency vector
  MBErrorCode get_associated_meshsets( MBEntityHandle source_entity, 
                                        std::vector<MBEntityHandle> &target_entities );

//! get the element defined by the vertices in vertex_list, of the
//! type target_type, passing back in target_entity; if create_if_missing
//! is true and no entity is found, one is created; if create_adjacency_option
//! is >= 0, adjacencies from entities of that dimension to target_entity
//! are created (only create_adjacency_option=0 is supported right now,
//! so that never creates other ancillary entities); explicitly require
//! the vertex_list_size for consistency, even though we could probably get
//! it from target_type
  MBErrorCode get_element(const MBEntityHandle *vertex_list,
                           const int vertex_list_size,
                           const MBEntityType target_type,
                           MBEntityHandle &target_entity,
                           const bool create_if_missing,
                           const MBEntityHandle source_entity = 0,
                           const int create_adjacency_option = -1);

  MBErrorCode get_adjacencies(const MBEntityHandle entity,
                               const unsigned int to_dimension,
                               bool create_if_missing,
                               std::vector<MBEntityHandle> &adjacent_entities);

    //! return const array * for adjacencies
  MBErrorCode get_adjacencies(MBEntityHandle entity,
                               const MBEntityHandle *&adjacent_entities,
                               int &num_entities);
  
  //! returns the entities in sorted order
  MBErrorCode get_adjacencies(MBEntityHandle entity,
                               std::vector<MBEntityHandle>& adjacent_entities);
  

  //! creates vertex to element adjacency information
  MBErrorCode create_vert_elem_adjacencies();

  //! returns whether vertex to element adjacencies are being stored
  bool vert_elem_adjacencies() const { return mVertElemAdj; }

  //! calling code notifying this that an entity is getting deleted
  MBErrorCode notify_delete_entity(MBEntityHandle entity);

  //! calling code notifying this that to update connectivity of 'entity' 
  MBErrorCode notify_create_entity(const MBEntityHandle entity, 
                                    const MBEntityHandle *node_array,
                                    const int number_nodes);

  //! calling code notifying that an entity changed its connectivity
  MBErrorCode notify_change_connectivity(MBEntityHandle entity, 
                                          const MBEntityHandle* old_array,
                                          const MBEntityHandle* new_array, int number_nodes);

    //! return true if 2 entities are explicitly adjacent
  bool explicitly_adjacent(const MBEntityHandle ent1,
                           const MBEntityHandle ent2);

private:

  //! private constructor to prevent the construction of a default one
  AEntityFactory();

  //! interface associated with this tool
  MBInterface *thisMB;

  //! adjacencies collection
  std::vector<DensePageGroup*> mDensePageGroups;

  //! whether vertex to element adjacencies are begin done
  bool mVertElemAdj;
  
  //! compare vertex_list to the vertices in this_entity, 
  //!  and return true if they contain the same vertices
  bool entities_equivalent(const MBEntityHandle this_entity, 
                           const MBEntityHandle *vertex_list, 
                           const int vertex_list_size,
                           const MBEntityType target_type);

  MBErrorCode get_zero_to_n_elements(MBEntityHandle source_entity,
                            const unsigned int target_dimension,
                            std::vector<MBEntityHandle> &target_entities,
                            const bool create_if_missing,
                            const int create_adjacency_option = -1);

  MBErrorCode get_down_adjacency_elements(MBEntityHandle source_entity,
                            const unsigned int target_dimension,
                            std::vector<MBEntityHandle> &target_entities,
                            const bool create_if_missing,
                            const int create_adjacency_option = -1);

  MBErrorCode get_down_adjacency_elements_poly(MBEntityHandle source_entity,
                                               const unsigned int target_dimension,
                                               std::vector<MBEntityHandle> &target_entities,
                                               const bool create_if_missing,
                                               const int create_adjacency_option = -1);

  MBErrorCode get_up_adjacency_elements(MBEntityHandle source_entity,
                            const unsigned int target_dimension,
                            std::vector<MBEntityHandle> &target_entities,
                            const bool create_if_missing,
                            const int create_adjacency_option = -1);

};

/*! 
   calling code is notifying this that an entity is going to be deleted
   from the database
*/
inline MBErrorCode AEntityFactory::notify_delete_entity(MBEntityHandle entity)
{
  // remove any references to this entity from other entities
  return remove_all_adjacencies(entity, true);
}


#endif
