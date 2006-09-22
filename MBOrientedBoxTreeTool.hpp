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
#include <list>

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
        int max_depth;           //!< Maximum tree depth - 0->no limit
        //! Must be in [best_split_ratio,1.0]
        //! A tree node will not be split if the ratio of children
        //! in the child nodes is greater than this value.
        double worst_split_ratio;
        //! Must be in [0.0,worst_split_ratio]
        //! The search for an optimal split plane for splitting a node
        //! will stop if at least this ratio is achieved for the number of
        //! entities on each side of the split plane.
        double best_split_ratio;
        //! Flags used to create entity sets representing tree nodes
        unsigned int set_options;
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
     
    /**\brief Build a tree of sets, where each set contains triangles.
     *
     * Build a tree of sets.  Each set must contain at least one triangle
     * to define it's geometry.  Each passed set will become a leaf of
     * the OBB tree.  Settings controlling tree depth are ignored by
     * this method.  The tree will be as deep as it needs to be for each
     * input set to be a leaf.
     *
     * To build a tree representing the surfaces of a geometric volume,
     * 1) Build and OBB tree for each surface using the 'build' method
     * 2) Add each surface to the contents of the resulting OBB tree root set
     * 3) Build a tree from all the surface OBB tree root sets using this
     *    method to get a combined tree for the volume.
     */
    MBErrorCode join_trees( const MBRange& tree_roots,
                            MBEntityHandle& root_set_out,
                            const Settings* settings = 0 );

    /**\brief Intersect a ray with the triangles contained within the tree
     *
     * Intersect a ray with the triangles contained in the tree and return
     * the distance at which the intersection occured.
     *\param distances_out The output list of intersection points on the ray.
     *\param root_set      The MBENTITYSET representing the root of the tree.
     *\param tolerance     The tolerance to use in intersection checks.
     *\param ray_point     The base point of the ray.
     *\param unit_ray_dir  The ray direction vector (must be unit length)
     *\param ray_length    Optional ray length (intersect segment instead of ray.)
     */
    MBErrorCode ray_intersect_triangles( std::vector<double>& distances_out,
                                         MBEntityHandle root_set,
                                         double tolerance,
                                         const double ray_point[3],
                                         const double unit_ray_dir[3],
                                         const double* ray_length = 0 );
    
    /**\brief Intersect ray with tree
     *
     * Return the tree nodes (as MBENTITYSET handles) for the leaf boxes
     * of the tree intersected by a ray.
     *\param boxes_out    The boxes intersected by the ray.
     *\param tolerance     The tolerance to use in intersection checks.
     *\param ray_point     The base point of the ray.
     *\param unit_ray_dir  The ray direction vector (must be unit length)
     *\param ray_length    Optional ray length (intersect segment instead of ray.)
     */
    MBErrorCode ray_intersect_boxes( MBRange& boxes_out,
                                     MBEntityHandle root_set,
                                     double tolerance,
                                     const double ray_point[3],
                                     const double unit_ray_dir[3],
                                     const double* ray_length = 0 );

    /**\brief Intersect ray with triangles contained in passed MBENTITYSETs */
    MBErrorCode ray_intersect_triangles( 
                          std::vector<double>& intersection_distances_out,
                          const MBRange& leaf_boxes_containing_tris,
                          double tolerance,
                          const double ray_point[3],
                          const double unit_ray_dir[3],
                          const double* ray_length = 0);
                          

    /**\brief Intersect a ray with the triangles contained within the tree
     *
     * Intersect a ray with the triangles contained in the tree and return
     * the distance at which the intersection occured.
     *\param distances_out The output list of intersection points on the ray.
     *\param sets_out      The contained set encountered during the tree traversal
     *                     (see 'set_build').  For the most common use, this is the
     *                     set corresponding to the geometric surface containing the
     *                     intersected triangle.
     *\param root_set      The MBENTITYSET representing the root of the tree.
     *\param min_tolerance_intersections This method returns all intersections
     *                     within 'tolerance' of the start of the ray and if 
     *                     the number of intersections within the 'tolerance' of the
     *                     ray start point is less than this number, the next closest
     *                     intersection.  If the desired result is only the closest
     *                     intersection, pass zero for this argument.
     *\param tolerance     The tolerance to use in intersection checks.
     *\param ray_point     The base point of the ray.
     *\param unit_ray_dir  The ray direction vector (must be unit length)
     *\param ray_length    Optional ray length (intersect segment instead of ray.)
     */
    MBErrorCode ray_intersect_sets( std::vector<double>& distances_out,
                                    std::vector<MBEntityHandle>& sets_out,
                                    MBEntityHandle root_set,
                                    double tolerance,
                                    unsigned min_tolerace_intersections,
                                    const double ray_point[3],
                                    const double unit_ray_dir[3],
                                    const double* ray_length = 0 );
    
    /**\brief Find closest surface, facet in surface, and location on facet
     *
     * Find the closest location in the tree to the specified location.
     *\param point Location to search from
     *\param point_out Closest location on closest facet
     *\param facet_out Closest 2D element to input position
     *\param set_out Set containing closest facet.  0 if tree was not 
     *               constructed using 'set_build'
     */
    MBErrorCode closest_to_location( const double* point,
                                     MBEntityHandle tree_root,
                                     double tolerance,
                                     double* point_out,
                                     MBEntityHandle& facet_out,
                                     MBEntityHandle* set_out = 0);
    
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
  
    /** \brief Implement this and pass instance to preorder_traverse
     * 
     * This interface may be implemented and an instance passed to
     * preorder_traverse to define some operation to do when traversing
     * the tree.
     */
    class Op {
      public:

        /**\brief Visit a node in the tree during a traversal.
         *
         * This method is called for each node in the tree visited
         * during a pre-order traversal.  
         *\param node The MBEntityHandle for the entity set for the tree node.
         *\param depth The current depth in the tree.
         *\param descend Output: if false, traversal will skip children
         *             of the current node, or if the current node is a
         *             leaf, the 'leaf' method will not be called.
         */
        virtual MBErrorCode visit( MBEntityHandle node,
                                   int depth,
                                   bool& descend ) = 0;
       
        /**\brief Process a leaf node during tree traversal */
        virtual MBErrorCode leaf( MBEntityHandle node ) = 0;

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
  
    struct SetData;
  private:
  
    MBErrorCode build_tree( const MBRange& entities, 
                            MBEntityHandle& set, 
                            int depth,
                            const Settings& settings );
  
    MBErrorCode build_sets( std::list<SetData>& sets,
                            MBEntityHandle& node_set,
                            int depth,
                            const Settings& settings );
  
    MBInterface* instance;
    MBTag tagHandle;
};

#endif
