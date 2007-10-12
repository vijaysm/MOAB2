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



#ifndef MB_SKINNER_HPP
#define MB_SKINNER_HPP

#include "MBForward.hpp"
#include <vector>

class MBSkinner 
{

  enum direction{FORWARD=1, REVERSE=-1};
protected:
  //! the MB instance that this works with
  MBInterface* thisMB;

  MBTag mDeletableMBTag;
  MBTag mAdjTag;
  int mTargetDim;

public:
  //! constructor, takes mdb instance
  MBSkinner(MBInterface* mdb) 
    : thisMB(mdb), mDeletableMBTag(0), mAdjTag(0){}

  //! destructor
  ~MBSkinner();

  MBErrorCode find_geometric_skin(MBRange &forward_target_entities);
  
  // will accept entities all of one dimension
  // and return entities of n-1 dimension
  MBErrorCode find_skin( const MBRange &entities,
                         MBRange &forward_lower_entities,
                         MBRange &reverse_lower_entities,
                         bool create_vert_elem_adjs = false);

    // get skin entities of prescribed dimension
  MBErrorCode find_skin(const MBRange &entities,
                        int dim,
                        MBRange &skin_entities,
                        bool create_vert_elem_adjs = false);

  MBErrorCode classify_2d_boundary( const MBRange &boundary,
                                     const MBRange &bar_elements,
                                     MBEntityHandle boundary_edges,
                                     MBEntityHandle inferred_edges,
                                     MBEntityHandle non_manifold_edges,
                                     MBEntityHandle other_edges,
                                     int &number_boundary_nodes);
  
  //!given a skin of dimension 2, will classify and return edges
  //! as boundary, inferred, and non-manifold, and the rest (other)
  MBErrorCode classify_2d_boundary( const MBRange  &boundary,
                                     const MBRange  &mesh_1d_elements,
                                     MBRange  &boundary_edges,
                                     MBRange  &inferred_edges,
                                     MBRange  &non_manifold_edges,
                                     MBRange  &other_edges,
                                     int &number_boundary_nodes);

protected:
  
  void initialize();
  
  void deinitialize();

  void add_adjacency(MBEntityHandle entity);
  
  void add_adjacency(MBEntityHandle entity, const MBEntityHandle *conn,
                     const int num_nodes);

  MBErrorCode remove_adjacency(MBEntityHandle entity);

  bool entity_deletable(MBEntityHandle entity);

  void find_match( MBEntityType type, 
                   const MBEntityHandle *conn, 
                   const int num_nodes,
                   MBEntityHandle& match,
                   MBSkinner::direction &direct);

  bool connectivity_match(const MBEntityHandle *conn1,
                          const MBEntityHandle *conn2,
                          const int num_verts,
                          MBSkinner::direction &direct);

  void find_inferred_edges(MBRange &skin_boundary,
                           MBRange &candidate_edges,
                           MBRange &inferred_edges,
                           double reference_angle_degrees);

  bool has_larger_angle(MBEntityHandle &entity1,
                       MBEntityHandle &entity2,
                       double reference_angle_cosine);

};


#endif

