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

#include <vector>
#include <iostream>
#include <string>

namespace moab 
{
    class TreeStats{
  public:
        //! constructor
      TreeStats() 
              : leafObjectTests(0) 
          {}      
      
        //! return counts of nodes visited, indexed by tree depth.  
        //! the counts include both leaves and interior nodes
      const unsigned int &nodes_visited() const
          {return nodesVisited;}

        //! return counts of tree leaves visited, indexed by tree depth
      const unsigned int &leaves_visited() const
          {return leavesVisited;}

        //! return counts of traversals ended, indexed by tree depth
      const unsigned int &num_traversals() const 
          {return numTraversals;}

        //! return counts of nodes visited, indexed by tree depth.  
        //! the counts include both leaves and interior nodes
      unsigned int &nodes_visited() 
          {return nodesVisited;}

        //! return counts of tree leaves visited, indexed by tree depth
      unsigned int &leaves_visited()
          {return leavesVisited;}

        //! return counts of traversals ended, indexed by tree depth
      unsigned int &num_traversals()
          {return numTraversals;}

        //! return total number of leaf-object tests (ray-triangle, point in elem, etc.) performed
      const unsigned int &leaf_object_tests() const
          {return leafObjectTests;}

      unsigned int &leaf_object_tests()
          {return leafObjectTests;}

        //! reset all counters on this structure
      void reset();
              
        //! print the contents of this structure to given stream
      void print() const ;

  private:
      unsigned int nodesVisited;
      unsigned int leavesVisited;
      unsigned int numTraversals;
      unsigned int leafObjectTests;
    };

    inline void TreeStats::reset()
    {
      nodesVisited = 0;
      leavesVisited = 0;
      numTraversals = 0;
      leafObjectTests = 0;
    }
    
    inline void TreeStats::print() const {
      std::cout << "NodesVisited      = " << nodesVisited << std::endl;
      std::cout << "LeavesVisited     = " << leavesVisited << std::endl;
      std::cout << "Num Traversals    = " << numTraversals << std::endl;
      std::cout << "Leaf Object Tests = " << leafObjectTests << std::endl;
    }
}

    
      


#endif
