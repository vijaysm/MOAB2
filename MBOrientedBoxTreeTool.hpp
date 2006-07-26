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

/**\file MBOrientedBoxTreeTool.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2006-07-18
 */

#ifndef MB_ORIENTED_BOX_TREE_TOOL_HPP
#define MB_ORIENTED_BOX_TREE_TOOL_HPP

#include "MBInterface.hpp"

#include <iosfwd>

/** If true, deviate from the VTk algorithm for chosing the
 *  split direction by also considering the number of entities entirely
 *  to each side of the plane (triangles that do not intersect the plane)
 */
#define MB_OOB_SPLIT_BY_NON_INTERSECTING 0

class MBRange;
class MBOrientedBox;


class MBOrientedBoxTreeTool
{
  public:
  
    /**\brief Misc. knobs controlling tree subdivision
     *
     * Available settings for controlling when and how nodes in the tree
     * are split.  The constructor will initialize to the default
     * settings.  All settings except best_split_ratio control when
     * a node is subdivied.  best_split_ratio influences the choice
     * of how the node is subdivied.
     *
     * A calculated ratio is used in the determination of when and how
     * to split a node.  The ratio is calculated as:
     * - \f$max(\frac{|n_L - n_R|}{n_L+n_R}, f*\frac{n_I}{n_L+n_R})\f$
     * - \f$n_L\f$ : num entities to be placed in left child
     * - \f$n_R\f$ : num entities to be placed in right child
     * - \f$f\f$ : Settings::intersect_ratio_factor
     * - \f$n_I\f$: num entities intersecting split plane
     *
     * ALL of the following conditions must be met for a node to be further
     * subdivied:
     *  - Depth must be less than max_depth
     *  - Node must contain more than max_leaf_entities entities.
     *  - The 'ratio' must be less than worst_split_ratio
     *
     * The node will be subdivied using a plane normal to one of the
     * box axis and containg the box center.  The planes are tested 
     * beginning with the one orthogonal to the longest box axis and
     * finishing with the one orthogonal to the shortest box axis.  The
     * search will stop at the first plane for which the 'ratio' is
     * at least Settings::best_split_ratio .  Giving Settings::best_split_ratio
     * a non-zero value gives preference to a split orthogonal to larger
     * box dimensions.
     */
    struct Settings {
      public:
        Settings();              //!< set defaults
        int max_leaf_entities;   //!< Average number of entities per leaf
        int max_depth;           //!< Maximum tree depth
        //! Must be in [best_split_ratio,1.0]
        //! A tree node will not be split if the ratio of children
        //! in the child nodes is greater than this value.
        double worst_split_ratio;
        //! Must be in [0.0,worst_split_ratio]
        //! The search for an optimal split plane for splitting a node
        //! will stop if at least this ratio is achieved for the number of
        //! entities on each side of the split plane.
        double best_split_ratio;
#if MB_OOB_SPLIT_BY_NON_INTERSECTING
        double intersect_ratio_factor;
#endif
        //! Check if settings are valid.
        bool valid() const;
    };
  
    MBOrientedBoxTreeTool( MBInterface* i, 
                           const char* tag_name = 0 ) ;
  
    /**\brief Build oriented bounding box tree
     *
     * Build an oriented bounding box tree.  
     *\param entities A list of either vertices or 2-D elements (not both)
     *                for which to build a tree.
     *\param set_handle_out A handle for the entity set representing the
     *                root of the tree.
     */
    MBErrorCode build( const MBRange& entities, 
                       MBEntityHandle& set_handle_out,
                       const Settings* settings = 0 );
    
    /**\brief Get oriented box at node in tree
     *
     * Get the oriented box for a node in an oriented bounding box tree.
     */
    MBErrorCode box( MBEntityHandle node_set,
                     MBOrientedBox& box );
    
    /**\brief Get oriented box at node in tree
     *
     * Get the oriented box for a node in an oriented bounding box tree.
     */
    MBErrorCode box( MBEntityHandle node_set,
                     double center[3],
                     double axis1[3],
                     double axis2[3],
                     double axis3[3] );
                         
    /**\brief Test for leaf node / get children
     *
     * Test if a node in an oriented bouding box tree is a leaf
     * node and if it is not, get the two child nodes.
     */
    MBErrorCode children( MBEntityHandle node_set,
                          bool& is_leaf,
                          MBEntityHandle* children = 0 );

    MBErrorCode delete_tree( MBEntityHandle root_set );

    /**\brief Print out tree
     *
     * Print the tree to an output stream in a human-readable form.
     *\param tree_root_set  Entity set representing tree root.
     *\param list_contents  If true, list entities in each tree node,
     *                      If false, just list number of entities.
     *\param id_tag_name    If specified, must be the name of an existing
     *                      integer tag containing an ID for the entities.
     *                      Not used if list_contents is false.
     */
    void print( MBEntityHandle tree_root_set, 
                std::ostream& stream,
                bool list_contents = false,
                const char* id_tag_name = 0 );
                
    /**\brief Print tree statistics
     *
     * Print misc. stats. describing tree
     */
    MBErrorCode stats( MBEntityHandle tree_root_set, std::ostream& stream );
  
    class Op {
      public:
        virtual MBErrorCode operator()( MBEntityHandle node,
                                        int depth,
                                        bool& descend ) = 0;
        virtual ~Op(); // probably isn't necessary in this case, and
                       // does nothing, but lots of compilers warn if
                       // virtual function but no virtual destructor.
    };
    
    /**\brief Visistor pattern - do operation for each tree node
     *
     * Do a preorder traversal of the tree, calling the method
     * in the passed operation instance for each node in the tree.
     * Parent node is visited before either child (pre-order traversal).
     * If operator method passes back the 'descend' argument as false,
     * traversal will not descend to the children of the current node.
     */
    MBErrorCode preorder_traverse( MBEntityHandle root_set,
                                   Op& operation );
  
    MBInterface* get_moab_instance() const { return instance; }
  
  private:
  
    MBErrorCode build_tree( const MBRange& entities, 
                            MBEntityHandle& set, 
                            int depth,
                            const Settings& settings );
  
    MBInterface* instance;
    MBTag tagHandle;
};

#endif
