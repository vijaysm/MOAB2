/**\file TreeStats.hpp
 * \class moab::TreeStats
 * \brief Traversal statistics accumulating and reporting
 *
 * Class to accumulate statistics on traversal performance. This structure contains the
 * count of nodes visited at each level in a tree, and the count of traversals that ended
 * at each level.  One TrvStats structure can be used with multiple trees or multiple
 * queries, or used on only a single tree or a single query.
 *
 * Note that these traversal statistics are not related to the stats() query below,
 * which calculates static information about a tree.  These statistics relate
 * to a tree's dynamic behavior on particular operations.
 */

#ifndef TREESTATS_HPP
#define TREESTATS_HPP

#include "moab/Interface.hpp"

#include <vector>
#include <iostream>
#include <string>

namespace moab 
{
    class TreeStats{
  public:
        //! constructor
      TreeStats() {reset();}

        /** \brief Given a root node, compute the stats for a tree 
         * \param impl MOAB instance pointer
         * \param root_node Root entity set for the tree
         */
      ErrorCode compute_stats(Interface *impl, EntityHandle root_node);
      
        //! reset traversal counters
      void reset_trav_stats();
              
        //! reset all counters
      void reset();
              
        //! print the contents of this structure
      void print() const ;

        //! output the contents of this structure on a single line
      void output() const ;

        // times
      double initTime;

        // tree stats that depend only on structure (not traversal)
      unsigned int maxDepth;
      unsigned int numNodes;
      unsigned int numLeaves;

        // traversal statistics
      unsigned int nodesVisited;
      unsigned int leavesVisited;
      unsigned int numTraversals;
      unsigned int leafObjectTests;
      unsigned int boxElemTests;

  private:
      ErrorCode traverse(Interface *impl, EntityHandle node, unsigned int &depth);
      
    };

    inline ErrorCode TreeStats::compute_stats(Interface *impl, EntityHandle root_node) 
    {
      maxDepth = 0;
      numNodes = 0;
      numLeaves = 0;
      return traverse(impl, root_node, maxDepth);
    }
      
    inline ErrorCode TreeStats::traverse(Interface *impl, EntityHandle node, unsigned int &depth) 
    {
      depth++;
      numNodes++;
      std::vector<EntityHandle> children;
      children.reserve(2);
      ErrorCode rval = impl->get_child_meshsets(node, children);
      if (MB_SUCCESS != rval) return rval;
      if (children.empty()) {
        numLeaves++;
        return MB_SUCCESS;
      }
      else {
        unsigned int right_depth = depth, left_depth = depth;
        rval = traverse(impl, children[0], left_depth);
        if (MB_SUCCESS != rval) return rval;
        rval = traverse(impl, children[1], right_depth);
        if (MB_SUCCESS != rval) return rval;
        depth = std::max(left_depth, right_depth);
        return MB_SUCCESS;
      }
    }

    inline void TreeStats::reset()
    {
      initTime = 0.0;
      
      maxDepth = 0;
      numNodes = 0;
      numLeaves = 0;
      
      reset_trav_stats();
    }
    
    inline void TreeStats::reset_trav_stats() 
    {
      nodesVisited = 0;
      leavesVisited = 0;
      numTraversals = 0;
      leafObjectTests = 0;
      boxElemTests = 0;
    }
    
    inline void TreeStats::print() const {
      std::cout << "Tree initialization time = " << initTime << " seconds" << std::endl;
      
      std::cout << "Num nodes         = " << numNodes << std::endl;
      std::cout << "Num leaves        = " << numLeaves << std::endl;
      std::cout << "Max depth         = " << maxDepth << std::endl << std::endl;

      std::cout << "NodesVisited      = " << nodesVisited << std::endl;
      std::cout << "LeavesVisited     = " << leavesVisited << std::endl;
      std::cout << "Num Traversals    = " << numTraversals << std::endl;
      std::cout << "Leaf Object Tests = " << leafObjectTests << std::endl;
      std::cout << "Box-Element Tests = " << boxElemTests << std::endl;
    }

    inline void TreeStats::output() const 
    {
      std::cout << initTime << " " << nodesVisited << " " << leavesVisited << " " << numTraversals << " " << leafObjectTests << " " << boxElemTests
                << " # initTime, nodesVisited, leavesVisited, numTraversals, leafObjectTests, boxElemTests" << std::endl;
    }
}

    
      


#endif
