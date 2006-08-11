/*
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

/**\file MeshSetManager.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2006-08-10
 */

#ifndef MB_MESH_SET_MANAGER_HPP
#define MB_MESH_SET_MANAGER_HPP

#include "MBInterface.hpp"

#include <vector>
#include <deque>

class MBMeshSet;
class AEntityFactory;

/**\brief Manage global list of MBMeshSets 
 *
 * Maintain the global list of MBMeshSets.  Provide
 * handle->MBMeshSet* mapping.
 */
class MeshSetManager {
public:

  MeshSetManager( AEntityFactory* a_ent_fact );
  ~MeshSetManager();
  
  MBErrorCode create_mesh_set( unsigned options, 
                               MBEntityHandle start_id,
                               unsigned int proc_id,
                               MBEntityHandle& output_handle );
  
  MBErrorCode get_options( MBEntityHandle meshset, unsigned& options ) const;
  
  /**\brief Destroy mesh set
   *
   * Destroy the MBMeshSet instance corresponding to the passed 
   * handle.  
   *\return MB_SUCCESS, MB_ENTITY_NOT_FOUND or MB_TYPE_OUT_OF_RANGE
   */
  MBErrorCode destroy_mesh_set( MBEntityHandle handle );
  
  /**\brief Remove contents of mesh set */
  MBErrorCode clear_mesh_set( MBEntityHandle meshset );
  
  void get_all_mesh_sets( MBRange& output_list ) const;
  
  MBErrorCode get_entities( MBEntityHandle meshset, MBRange& entities, bool recursive ) const;
  MBErrorCode get_entities( MBEntityHandle meshset, std::vector<MBEntityHandle>& entities ) const;
  MBErrorCode get_dimension( MBEntityHandle meshset, int dim, MBRange& entities, bool recursive ) const;
  MBErrorCode get_type( MBEntityHandle meshset, MBEntityType type, MBRange& entities, bool recursive ) const;
  
  MBErrorCode num_entities( MBEntityHandle, int& number, bool recursive ) const;
  MBErrorCode num_dimension( MBEntityHandle, int dimension, int& number, bool recursive ) const;
  MBErrorCode num_type( MBEntityHandle, MBEntityType type, int& number, bool recursive ) const;

  MBErrorCode subtract ( MBEntityHandle from, MBEntityHandle subtr );
  MBErrorCode intersect( MBEntityHandle into, MBEntityHandle other );
  MBErrorCode unite    ( MBEntityHandle into, MBEntityHandle other );
  
  MBErrorCode add_entities( MBEntityHandle to, const MBRange& entities );
  MBErrorCode add_entities( MBEntityHandle to, const MBEntityHandle* array, int array_len );
  MBErrorCode remove_entities( MBEntityHandle from, const MBRange& entities );
  MBErrorCode remove_entities( MBEntityHandle form, const MBEntityHandle* array, int array_len );
  
  MBErrorCode get_parents ( MBEntityHandle of, std::vector<MBEntityHandle>& parents, int num_hops ) const;
  MBErrorCode get_children( MBEntityHandle of, std::vector<MBEntityHandle>& children, int num_hops ) const;
  MBErrorCode num_parents ( MBEntityHandle of, int& number, int num_hops ) const;
  MBErrorCode num_children( MBEntityHandle of, int& number, int num_hops ) const;
  
  MBErrorCode add_parent( MBEntityHandle to, MBEntityHandle parent );
  MBErrorCode add_child ( MBEntityHandle to, MBEntityHandle child );
  MBErrorCode add_parent_child( MBEntityHandle parent, MBEntityHandle child );
  
  MBErrorCode remove_parent( MBEntityHandle from, MBEntityHandle parent );
  MBErrorCode remove_child ( MBEntityHandle from, MBEntityHandle child  );
  MBErrorCode remove_parent_child( MBEntityHandle parent, MBEntityHandle child );

private:

  //! Recursively get all contained entity sets.
  MBErrorCode recursive_get_sets( MBEntityHandle start_set, 
                                std::vector<MBMeshSet*>& sets ) const;
  
  //! Do breadth-first search.  Get parent or child entity sets.
  //!\param meshset The set to start the query from
  //!\param results As output, a sorted list of entity set handles.
  //!               Any entity sets in this list when it is passed in
  //!               will be skipped during the traversal (except for 
  //!               the input set, which is always traversed.)
  //!\param depth   Get all sets within this many hops from the input
  //!               set.
  //!\param parents If true, traverse parent links.  If false, traverse
  //!               child links.
  MBErrorCode get_parent_child_meshsets( MBEntityHandle meshset,
                                    std::vector<MBEntityHandle>& results,
                                    int num_hops, bool parents ) const;

  /**\brief Get MBMeshSet* given handle 
   *
   * Get the MBMeshSet instance corresponding to the passed 
   * handle. 
   *\return Mesh set instance, or NULL if not found.
   */
  MBMeshSet* get_mesh_set( MBEntityHandle handle ) const;
  
  struct ProcData {
    std::vector<MBMeshSet*> list;
    std::deque<MBEntityHandle> free;
  };

    // Keep separate lists for each processor, otherwise
    // will end up with big holes (unused blocks) in the
    // std::vector.
  std::vector<ProcData> procSet; 
  
  AEntityFactory* aEntityFactory;
};


#endif
