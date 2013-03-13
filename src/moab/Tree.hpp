/**\file Tree.hpp
 * \class moab::Tree
 * \brief Parent class of various tree types in MOAB
 */

#ifndef MOAB_TREE_HPP
#define MOAB_TREE_HPP

#include "moab/Interface.hpp"
#include "moab/CartVect.hpp"
#include "moab/FileOptions.hpp"

#include <string>
#include <vector>
#include <math.h>
#include <assert.h>

namespace moab {

    class Interface;
    class Range;

    class Tree
    {
  public:
        /** \brief Constructor (bare)
         * \param iface MOAB instance 
         */
      Tree(Interface* iface);

        /** \brief Constructor (build the tree on construction)
         * Construct a tree object, and build the tree with entities input.  See comments
         * for build_tree() for detailed description of arguments.
         * \param iface MOAB instance 
         * \param entities Entities to build tree around
         * \param tree_root Root set for tree (see function description)
         * \param opts Options for tree (see function description)
         */
      Tree(Interface* iface, const Range &entities, 
           EntityHandle *tree_root_set = NULL, FileOptions *opts = NULL);

        /** \brief Destructor
         */
      virtual ~Tree();

        /** Build the tree
         * Build a tree with the entities input.  If a non-NULL tree_root_set pointer is input, 
         * use the pointed-to set as the root of this tree (*tree_root_set!=0) otherwise construct 
         * a new root set and pass its handle back in *tree_root_set.  Options vary by tree type, 
         * with a few common to all types of trees.  Common options:
         * MAX_PER_LEAF: max entities per leaf; default = 6
         * MAX_DEPTH: max depth of the tree; default = 30
         * MIN_WIDTH: minimum width of box, used like a tolerance; default = 1.0e-10
         * MESHSET_FLAGS: flags passed into meshset creation for tree nodes; should be a value from
         *          ENTITY_SET_PROPERTY (see Types.hpp); default = MESHSET_SET
         * CLEAN_UP: if false, do not delete tree sets upon tree class destruction; default = true
         * TAG_NAME: tag name to store tree information on tree nodes; default determined by tree type
         * \param entities Entities with which to build the tree
         * \param tree_root Root set for tree (see function description)
         * \param opts Options for tree (see function description)
         * \return Error is returned only on build failure
         */
      virtual ErrorCode build_tree(const Range& entities,
                                   EntityHandle *tree_root_set = NULL,
                                   FileOptions *options = NULL) = 0;

        /** \brief Destroy the tree maintained by this object, optionally checking we have the right root.
         * \param root If non-NULL, check that this is the root, return failure if not
         */
      virtual ErrorCode delete_tree(EntityHandle root = 0) = 0;

        /** \brief Get bounding box for entire tree
         * If no tree has been built yet, returns 3*0 for all dimensions.
         * \param box_min Minimum corner of box
         * \param box_max Maximum corner of box
         * \param tree_node If non-NULL, get bounding box for this node, otherwise for tree root
         * \return Only returns error on fatal condition
         */
      virtual ErrorCode get_bounding_box(double box_min[3], double box_max[3], EntityHandle *tree_node = NULL) const;
  
        /** \brief Get bounding box for entire tree
         * If no tree has been built yet, returns 3*0 for all dimensions.
         * \param box_min Minimum corner of box
         * \param box_max Maximum corner of box
         * \param tree_node If non-NULL, get bounding box for this node, otherwise for tree root
         * \return Only returns error on fatal condition
         */
      virtual ErrorCode get_bounding_box(CartVect &box_min, CartVect &box_max, EntityHandle *tree_node = NULL) const;
  
        /** \brief Return some basic information about the tree
         * Stats are returned for tree starting from input node or tree root (root = 0)
         * \param root If non-0, give stats below and including root
         * \param min Minimum corner of bounding box
         * \param max Maximum corner of bounding box
         * \param max_dep Maximum depth of tree below root
         */
      virtual ErrorCode get_info(EntityHandle root,
                                 double min[3], double max[3], 
                                 unsigned int &max_dep);
  
        /** \brief Find all trees, by bounding box tag
         */
      ErrorCode find_all_trees( Range& results );
      
        /** \brief Get leaf containing input position.
         *
         * Does not take into account global bounding box of tree.
         * - Therefore there is always one leaf containing the point.
         * - If caller wants to account for global bounding box, then
         * caller can test against that box and not call this method
         * at all if the point is outside the box, as there is no leaf
         * containing the point in that case.
         * \param point Point to be located in tree
         * \param leaf_out Leaf containing point
         * \param multiple_leaves Some tree types can have multiple leaves containing a point;
         *          if non-NULL, this parameter is returned true if multiple leaves contain
         *          the input point
         * \param start_node Start from this tree node (non-NULL) instead of tree root (NULL)
         * \return Non-success returned only in case of failure; not-found indicated by leaf_out=0
         */
      virtual ErrorCode point_search(const double *point,
                                     EntityHandle& leaf_out,
                                     bool *multiple_leaves = NULL,
                                     EntityHandle *start_node = NULL) = 0;

        /** \brief Find all leaves within a given distance from point
         * If dists_out input non-NULL, also returns distances from each leaf; if
         * point i is inside leaf, 0 is given as dists_out[i]
         * \param point Point to be located in tree
         * \param distance Distance within which to query
         * \param leaves Leaves within distance or containing point
         * \param dists If non-NULL, will contain distsances to leaves
         * \param start_node Start from this tree node (non-NULL) instead of tree root (NULL)
         */
      virtual ErrorCode distance_search(const double *point,
                                        const double distance,
                                        std::vector<EntityHandle>& leaves_out,
                                        std::vector<double> *dists_out = NULL,
                                        EntityHandle *start_node = NULL) = 0;

        /** \brief Compute bounding box of entities in elems
         * \param elems Entities for which bounding box is computed
         * \param box_min Minimum corner of box
         * \param box_max Maximum corner of box
         * \return This function returns error only under catastrophic; in the case of no entities, 
         *          it just returns a zero-extent box.
         */
      ErrorCode compute_bounding_box(const Range& elems, CartVect &box_min, CartVect &box_max) const;
      
        /** \brief Compute bounding box of entities in elems
         * \param elems Entities for which bounding box is computed
         * \param box_min Minimum corner of box
         * \param box_max Maximum corner of box
         * \return This function returns error only under catastrophic; in the case of no entities, 
         *          it just returns a zero-extent box.
         */
      ErrorCode compute_bounding_box(const Range& elems, double box_min[3], double box_max[3]) const;
      
        /** \brief Compute bounding box of an entity
         * \param ent Entity for which bounding box is computed
         * \param box_min Minimum corner of box
         * \param box_max Maximum corner of box
         * \return This function returns error only under catastrophic; in the case of no entities, 
         *          it just returns a zero-extent box.
         */
      ErrorCode compute_bounding_box(const EntityHandle ent, CartVect &box_min, CartVect &box_max) const;
      
        /** \brief Compute bounding box of an entity
         * \param ent Entity for which bounding box is computed
         * \param box_min Minimum corner of box
         * \param box_max Maximum corner of box
         * \return This function returns error only under catastrophic; in the case of no entities, 
         *          it just returns a zero-extent box.
         */
      ErrorCode compute_bounding_box(const EntityHandle ent, double box_min[3], double box_max[3]) const;
      
        /** \brief Return the MOAB interface associated with this tree
         */
      Interface* moab() { return mbImpl; }

        /** \brief Return the MOAB interface associated with this tree
         */
      const Interface* moab() const { return mbImpl; }

        /** \brief Get max depth set on tree */
      double get_max_depth() {return maxDepth;}
      
        /** \brief Get max entities per leaf set on tree */
      double get_max_per_leaf() {return maxPerLeaf;}
      
            
  protected:

        /** \brief Parse options common to all trees
         * \param options Options for representing tree; see Tree::build_tree() and subclass build_tree()
         *          functions for allowed options
         * \return Non-success returned from base class function only under catastrophic circumstances;
         *          derived classes also can recognize subclass-specific options
         */
      ErrorCode parse_common_options(FileOptions &options);

        /** \brief Get the box tag, possibly constructing it first
         * \param create_if_missing If true and it has not been made yet, make it
         */
      Tag get_box_tag(bool create_if_missing = true);

        // moab instance
      Interface *mbImpl;

        // bounding box corners for entire tree
      CartVect boxMin, boxMax;

        // max entities per leaf
      int maxPerLeaf;
      
        // max depth of tree
      int maxDepth;
      
        // tree depth, set by build_tree
      int treeDepth;
      
        // min width of box, handled like tolerance
      double minWidth;

        // meshset creation flags
      unsigned int meshsetFlags;

        // clean up flag
      bool cleanUp;

        // tree root
      EntityHandle myRoot;

        // tag used to mark bounding box of nodes
      Tag boxTag;

        // tag name used for boxTag
      std::string boxTagName;
      
  private:

    };

    inline Tree::Tree(Interface* iface) 
            : mbImpl(iface), maxPerLeaf(6), maxDepth(30), treeDepth(0), minWidth(1.0e-10),
              meshsetFlags(0), cleanUp(true), myRoot(0), boxTag(0)
    {}

    inline Tree::~Tree() 
    {
    }
    
    inline ErrorCode Tree::get_bounding_box(double box_min[3], 
                                            double box_max[3], EntityHandle *tree_node) const
    {
      return get_bounding_box(*reinterpret_cast<CartVect*>(box_min), *reinterpret_cast<CartVect*>(box_max), tree_node);
    }
  
    inline ErrorCode Tree::get_bounding_box(CartVect &box_min, CartVect &box_max, EntityHandle *tree_node) const
    {
      double tmp_tag[6];
      ErrorCode rval = moab()->tag_get_data(const_cast<Tree*>(this)->get_box_tag(), (tree_node ? tree_node : &myRoot), 1, tmp_tag);
      if (MB_SUCCESS != rval) return rval;
      box_min = tmp_tag;
      box_max = tmp_tag+3;
      return MB_SUCCESS;
    }
  
    inline ErrorCode Tree::compute_bounding_box(const Range &elems, double box_min[3], double box_max[3]) const
    {
      return compute_bounding_box(elems, *reinterpret_cast<CartVect*>(box_min), *reinterpret_cast<CartVect*>(box_max));
    }
  
    inline ErrorCode Tree::compute_bounding_box(const EntityHandle ent, double box_min[3], double box_max[3]) const
    {
      Range tmp_range(ent, ent);
      return compute_bounding_box(tmp_range, *reinterpret_cast<CartVect*>(box_min), *reinterpret_cast<CartVect*>(box_max));
    }
  
    inline ErrorCode Tree::compute_bounding_box(const EntityHandle ent, CartVect &box_min, CartVect &box_max) const
    {
      Range tmp_range(ent, ent);
      return compute_bounding_box(tmp_range, box_min, box_max);
    }
  
    inline ErrorCode Tree::get_info(EntityHandle /* root */,
                                    double * /*min[3]*/, double * /* max[3]*/, 
                                    unsigned int &/*dep*/) 
    {
      return MB_NOT_IMPLEMENTED;
    }

    inline Tag Tree::get_box_tag(bool create_if_missing) 
    {
      if (!boxTag && create_if_missing) {
        assert(boxTagName.length() > 0);
        ErrorCode rval = moab()->tag_get_handle(boxTagName.c_str(), 6, MB_TYPE_DOUBLE, boxTag, MB_TAG_CREAT | MB_TAG_SPARSE);
        if (MB_SUCCESS != rval) return 0;
      }
      
      return boxTag;
    }
    
} // namespace moab 

#endif
